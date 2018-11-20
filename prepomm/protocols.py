"""
Common protocols for simulation setup. If the protocol fails, you may just
have to do some manual setup. These aren't intended to work in all
situations, just in some common ones.
"""

from __future__ import print_function
import os
import tempfile
import collections
import copy
import numpy as np
import simtk.openmm as mm
import simtk.openmm.app as app
from simtk.unit import nanometer as nm
from simtk.unit import picosecond as ps
import simtk.unit as u
import mdtraj as md
import openmmtools


def topology_and_positions(input_data):
    """Extracts topology and positions from a wide range of inputs.

    Parameters
    ----------
    input_data :
        :class:`mdtraj.Trajectory` or :class:`openmm.app.Simulation` or name
        of PDB file or tuple or (topology, positions) or any class (e.g.,
        :class:`openmm.app.Modeller` with ``topology`` and ``positions``
        attributes

    Returns
    -------
    topology : :class:`openmm.Topology`
    positions : :class:`.unit.Quantity` (nanometers)
    """
    if isinstance(input_data, md.Trajectory):
        # using a temporary file seems to work best
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            input_data.save(tmp.name)
            pdb = app.PDBFile(tmp.name)
            return pdb.topology, pdb.positions
    elif isinstance(input_data, app.Simulation):
        state = input_data.context.getState(getPositions=True)
        return input_data.topology, state.getPositions()
    else:
        try:
            topology, positions = input_data.topology, input_data.positions
        except AttributeError:
            if isinstance(input_data, tuple) and len(input_data) == 2:
                return input_data
            elif os.path.isfile(input_data):
                pdb = app.PDBFile(input_data)
                return pdb.topology, pdb.positions
            else:
                return topology, positions
    # only reach here if nothing else worked
    raise RuntimeError("Unknown input %s. Not md.Trajecory or filename.",
                       input_data)

# TODO: clear out the namespace
_topology_and_positions = topology_and_positions
_modeller_args = _topology_and_positions

def make_simulation(pdbfile_or_mdtraj, ff_models, integrator=None,
                    platform=None, platform_properties=None):
    topology, positions = _topology_and_positions(pdbfile_or_mdtraj)
    forcefield = app.ForceField(*ff_models)
    if topology.getPeriodicBoxVectors():
        nonbonded=app.PME
    else:
        nonbonded=app.NoCutoff
    system = forcefield.createSystem(topology,
                                     nonbondedMethod=nonbonded)
    if integrator is None:
        integrator = mm.VerletIntegrator(0.001*ps)
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    return simulation

def make_modeller(pdbfile_or_mdtraj):
    return app.Modeller(*_topology_and_positions(pdbfile_or_mdtraj))

def minimize_vacuum(input_setup, ff_models, integrator=None):
    modeller = make_modeller(input_setup)
    forcefield = mm.app.ForceField(*ff_models)
    modeller.addHydrogens(forcefield)
    simulation = make_simulation((modeller.topology, modeller.positions),
                                 ff_models,
                                 integrator)
    simulation.minimizeEnergy()
    return simulation


def addH_and_solvate(input_setup, ff_models, box_vectors=None):
    # TODO: add more options for addHydrogens and addSolvent
    """
    This is largely stolen from the OpenMM docs, Example 5-2
    (http://docs.openmm.org/latest/userguide/application.html#saving-the-results)
    """
    print('Loading...')
    args = _modeller_args(input_setup)
    modeller = app.Modeller(*args)
    forcefield = app.ForceField(*ff_models)
    print('Adding hydrogens...')
    modeller.addHydrogens(forcefield)
    print('Adding solvent...')
    if box_vectors:
        print("... using box vectors provided:", box_vectors)
        solvent_kwargs = {'boxVectors': box_vectors}
    else:
        print("... using cubic box with 1 nm padding")
        solvent_kwargs = {'padding': 1*nm}
    modeller.addSolvent(forcefield, model=ff_models.water_model,
                        **solvent_kwargs)
    print(modeller.topology)
    print('Minimizing...')
    system = forcefield.createSystem(modeller.topology,
                                     nonbondedMethod=app.PME)
    integrator = mm.VerletIntegrator(0.001*ps)
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=100)
    print('Done')
    return simulation


# TODO: add position_constrained, equilibrate_nvt, equilibrate_npt
def run_position_constrained(simulation, constrained_atoms, duration):
    pass


def default_nvt_simulation(simulation, temperature=300.0*u.kelvin):
    integrator = openmmtools.integrators.VVVRIntegrator(
        temperature=temperature,
        collision_rate=1.0 / ps,
        timestep=0.002 * ps
    )
    topology, positions = _topology_and_positions(simulation)
    new_sim = app.Simulation(topology,
                             simulation.system,
                             integrator)
    new_sim.context.setPositions(positions)
    if isinstance(simulation, app.Simulation):
        state = simulation.context.getState(getVelocities=True)
        velocities = state.getVelocities()
        new_sim.context.setVelocities(velocities)
    return new_sim


def default_npt_simulation(simulation, temperature=300*u.kelvin,
                           pressure=1*u.atmosphere):
    integrator = openmmtools.integrators.VVVRIntegrator(
        temperature=temperature,
        collision_rate=1.0 / ps,
        timestep=0.002 * ps
    )
    system = copy.copy(simulation.system)
    system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))
    topology, positions = _topology_and_positions(simulation)
    new_sim = app.Simulation(topology, system, integrator)
    new_sim.context.setPositions(positions)
    if isinstance(simulation, app.Simulation):
        state = simulation.context.getState(getVelocities=True)
        velocities = state.getVelocities()
        new_sim.context.setVelocities(velocities)
    return new_sim
