"""
Miscellaneous analysis functions
"""
import itertools
import os.path

import scipy.spatial
import mdtraj as md
from simtk import unit as u

from .tools import _traj_from_file_or_traj

def max_atom_distance(file_or_trajectory):
    trajectory = _traj_from_file_or_traj(file_or_trajectory)
    positions = trajectory.xyz[0]
    hull = scipy.spatial.ConvexHull(positions)
    atom_pairs = list(itertools.combinations(hull.vertices, 2))
    distances = md.compute_distances(trajectory, atom_pairs)[0]
    return distances.max()

def concentration(file_or_trajectory, ion_name):
    trajectory = _traj_from_file_or_traj(file_or_trajectory)
    n_ions = len(trajectory.topology.select("name == " + ion_name))
    volumes = trajectory.unitcell_volumes * u.nanometer**3
    conc = (n_ions / volumes / u.AVOGADRO_CONSTANT_NA).in_units_of(u.molar)
    return conc

