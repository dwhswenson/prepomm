"""
Miscellaneous analysis functions
"""
import itertools
import scipy.spatial
import mdtraj as md

def max_atom_distance(trajectory):
    positions = trajectory.xyz[0]
    hull = scipy.spatial.ConvexHull(positions)
    atom_pairs = list(itertools.combinations(hull.vertices, 2))
    distances = md.compute_distances(trajectory, atom_pairs)[0]
    return distances.max()
