"""
Common protocols for simulation setup. If the protocol fails, you may just
have to do some manual setup.
"""

from __future__ import print_function
import simtk.openmm as mm
import simtk.openmm.app as app
from simtk.unit import nanometer as nm
from simtk.unit import picosecond as ps

def addH_and_solvate(pdb_infile, pdb_outfile, ff_models, box_vectors=None):
    # TODO: add more options for addHydrogens and addSolvent
    """
    This is largely stolen from the OpenMM docs, Example 5-2
    (http://docs.openmm.org/latest/userguide/application.html#saving-the-results)
    """
    print('Loading...')
    pdb = app.PDBFile(pdb_infile)
    forcefield = app.ForceField(*ff_models)
    modeller = app.Modeller(pdb.topology, pdb.positions)
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
    print('Minimizing...')
    system = forcefield.createSystem(modeller.topology,
                                     nonbondedMethod=app.PME)
    integrator = mm.VerletIntegrator(0.001*ps)
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=100)
    print('Saving...')
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(pdb_outfile, 'w') as pdb_out:
        app.PDBFile.writeFile(simulation.topology, positions, pdb_out)
    print('Done')
