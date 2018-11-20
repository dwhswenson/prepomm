"""
Tools for troubleshooting setup problems in OpenMM
"""

def residue_bonds_graph(topology, residue_numbers):
    """
    Requires networkx.

    Parameters
    ----------
    topology : mdtraj.Topology
    residue_numbers : int or list of int
        count from 0
    """
    try:
        residue_numbers = list(residue_numbers)
    except TypeError:
        residue_numbers = [residue_numbers]
    graph = topology.to_bondgraph()
    res_atoms = sum([list(topology.residue(res).atoms)
                     for res in residue_numbers], [])
    subgraph = graph.subgraph(res_atoms)
    return subgraph


