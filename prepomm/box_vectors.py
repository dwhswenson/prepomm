"""
Tools to set up periodic box vectors.
"""
from __future__ import print_function
from math import sqrt

import mdtraj as md
from simtk.openmm.vec3 import Vec3
from simtk.unit import nanometer as nm

from .analysis import max_atom_distance

def dodecahedral_box_vectors(base_size, padding=0.0):
    vectors = Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0.5, 0.5, sqrt(2)/2)
    return (float(base_size) + padding) * (vectors * nm)

def pdb_to_padded_dodecahedral_box(filename, padding, verbose=True):
    traj = md.load(filename)
    base_size = max_atom_distance(traj)
    if verbose:
        print(base_size)
    box_vectors = dodecahedral_box_vectors(base_size, padding)
    return box_vectors
