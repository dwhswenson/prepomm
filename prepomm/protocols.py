"""
Common protocols for simulation setup. If the protocol fails, you may just
have to do some manual setup. These aren't intended to work in all
situations, just in some common ones.
"""

from __future__ import print_function
import sys
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
import parmed

from .tools import steps_for_duration

def default_integrator(temperature):
    if isinstance(temperature, mm.Integrator):
        integrator = temperature
        temperature = integrator.getTemperature()
    else:
        integrator = openmmtools.integrators.VVVRIntegrator(
            temperature=temperature,
            collision_rate=1.0 / ps,
            timestep=0.002 * ps
        )
    return integrator, temperature

def default_barostat(pressure, temperature):
    if isinstance(pressure, mm.Force):
        barostat = pressure
    else:
        barostat = mm.MonteCarloBarostat(pressure, temperature, 25)
    return barostat


def default_equilibration_reporters(simulation, base_name, infix="_equil"):
    interval = steps_for_duration(1.0*ps, simulation)
    traj_reporter = md.reporters.DCDReporter(base_name + infix + ".dcd",
                                             interval)
    state_data = app.StateDataReporter(
        base_name + infix + ".csv",
        interval,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True
    )
    progress = app.StateDataReporter(
        sys.stdout,
        interval,
        time=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
        speed=True,
        elapsedTime=True
    )
    return [traj_reporter, state_data, progress]


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


def addH_and_solvate(input_setup, ff_models, box_vectors=None, **kwargs):
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
    add_hydrogens_keywords = ['pH', 'variants', 'platform']
    add_hydrogens_kwargs = {k: v for k, v in kwargs.items()
                            if k in add_hydrogens_keywords}
    modeller.addHydrogens(forcefield, **add_hydrogens_kwargs)
    print('Adding solvent...')
    add_solvent_keywords = [
        'positiveIon', 'negativeIon', 'ionicStrength', 'neutralize'
    ]
    solvent_kwargs = {k: v for k, v in kwargs.items()
                      if k in add_solvent_keywords}
    if box_vectors:
        print("... using box vectors provided:", box_vectors)
        solvent_kwargs.update({'boxVectors': box_vectors})
    else:
        print("... using cubic box with 1 nm padding")
        solvent_kwargs.update({'padding': 1*nm})

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


# TODO: add position_constrained
def run_position_constrained(simulation, constrained_atoms, duration,
                             file_basename, temperature=300.0*u.kelvin,
                             reporters=None):
    pos_const_sim = default_nvt_simulation(simulation, temperature)
    system = pos_const_sim.system
    old_masses = {atom_idx: system.getParticleMass(atom_idx)
                  for atom_idx in constrained_atoms}
    for atom_idx in constrained_atoms:
        system.setParticleMass(atom_idx, 0.0)

    out_sim = run_equilibration(pos_const_sim, duraction, file_basename,
                                reporters)
    system = out_sim.system
    for (atom_idx, mass) in old_masses.items():
        system.setParticleMass(atom_idx, mass)
    return out_sim


def default_nvt_simulation(simulation, temperature=300.0*u.kelvin):
    # note that temperature can also be an integrator
    integrator, temperature = default_integrator(temperature)
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
    integrator, temperature = default_integrator(temperature)
    barostat = default_barostat(pressure, temperature)
    system = copy.copy(simulation.system)
    system.addForce(barostat)
    topology, positions = _topology_and_positions(simulation)
    new_sim = app.Simulation(topology, system, integrator)
    new_sim.context.setPositions(positions)
    if isinstance(simulation, app.Simulation):
        state = simulation.context.getState(getVelocities=True)
        velocities = state.getVelocities()
        new_sim.context.setVelocities(velocities)
    return new_sim


def run_equilibration(simulation, duration, file_basename, reporters=None):
    n_steps = steps_for_duration(duration, simulation)
    if reporters is None:
        reporters = default_equilibration_reporters(simulation,
                                                    file_basename)
    simulation.reporters.extend(reporters)
    simulation.step(n_steps)
    return simulation

