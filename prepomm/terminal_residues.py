def get_5prime_phosphate(topology, chainid):
    fiveprime_res = list(topology.chain(chainid).residues)[0]
    fiveprime_res_sel = "resid " + str(fiveprime_res.index)
    phosphate_atom_sel = "((name =~ 'OP[1-2]') or name == 'P')"
    selection_str = fiveprime_res_sel + " and " + phosphate_atom_sel
    phosphate_5prime = list(topology.select(selection_str))
    return phosphate_5prime


def remove_5prime_phosphate(trajectory, chains):
    topology = trajectory.topology
    fiveprime_phosphates = set(sum([get_5prime_phosphate(topology, chain)
                                    for chain in chains], []))
    atoms_to_keep = set(range(trajectory.n_atoms)) - fiveprime_phosphates
    return trajectory.atom_slice(list(atoms_to_keep))
