"""
Tools used elsewhere in this package
"""
import os
import collections
import mdtraj as md
from simtk import openmm as mm
from simtk.openmm import app
from simtk import unit as u

def _traj_from_file_or_traj(file_or_traj):
    if isinstance(file_or_traj, md.Trajectory):
        traj = file_or_traj
    elif os.path.isfile(file_or_traj):
        traj = md.load(file_or_traj)
    else:
        raise TypeError("%s is neither a trajectory nor a filename",
                        file_or_traj)
    return traj

def steps_for_duration(duration, simulation):
    if isinstance(duration, u.Quantity):
        return int(duration / simulation.integrator.getStepSize())
    elif isinstance(duration, int):
        return duration
    else:
        raise RuntimeError("Unable to treat duration: %s", duration)

def simulation_write_pdb(simulation, pdb_outfile):
    """Write out the current state of the simulation as a PDB"""
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(pdb_outfile, 'w') as pdb_out:
        app.PDBFile.writeFile(simulation.topology, positions, pdb_out)


def simulation_to_mdtraj(simulation):
    topology, positions = _topology_and_positions(simulation)
    md_topology = md.Topology.from_openmm(topology)
    xyz = np.array([positions.value_in_unit(u.nanometer)])
    trajectory = md.Trajectory(xyz, topology)  # TODO unitcells
    return trajectory
    # with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
        # app.PDBFile.writeFile(topology, positions, tmp)
        # trajectory = md.load(tmp.name)
    return trajectory

def residue_type(res):
    if res.is_protein:
        return 'protein'
    elif res.is_nucleic:
        return 'nucleic'
    elif res.is_water:
        return 'water'
    else:
        return 'other'

def topology_describe(topology):
    total_str = ""
    for chain in topology.chains:
        restypes = collections.Counter([residue_type(res)
                                        for res in chain.residues])
        mystr = ", ".join([str(v) + " " + k + " residues"
                           for k, v in restypes.items()])
        total_str += mystr + "\n"
    return total_str
