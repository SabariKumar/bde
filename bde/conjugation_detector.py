#Yeonjoon's dirad conjugation code.


import networkx as nx
from rdkit import Chem
from itertools import islice, chain

# Just for drawing a molecule with atom indices
def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol


# Method #2
# 1. Find bonds consisting of non-SP3 atoms, allylic radical centers
# 2. Construct a graph consisting of the bonds identified from #1
# 3. Yen's Kth shortest path algorithm
# 4. Keep '1.5-1.5-1.5....' or 'multiple-single-multiple-single' paths consisting of at least 4 atoms
# (Rule out double-double-double-double ....)
def conj_paths(smi):
    BO_to_str = {1.0: 'S', 1.5: 'R', 2.0: 'D', 3.0: 'T'}

    # Steps 1-2
    G, G_no_W = nxGraph_from_smiles(smi)

    # Step 3
    unweighted_distances = nx.all_pairs_shortest_path_length(G_no_W)
    pairwise_d = {}

    for unweighted_distance in unweighted_distances:
        source_node, distance_to_each_node = unweighted_distance

        for target_node in distance_to_each_node.keys():
            if source_node != target_node:
                pairwise_d[(source_node, target_node)] = distance_to_each_node[target_node]

    # from longest to shoretst
    pairwise_d = {k: v for k, v in sorted(pairwise_d.items(), key=lambda item: item[1], reverse=True)}

    conj_paths_candidates = []
    for atom_pair in pairwise_d.keys():
        source_node, target_node = atom_pair

        if source_node > target_node:
            continue

        longest_path = [str(x) for x in k_shortest_paths(G_no_W, source_node, target_node, 10)[-1]]

        if len(longest_path) < 4:  # at least four atoms to define conjugation
            continue

        bond_orders = [G[int(longest_path[i])][int(longest_path[i + 1])]['weight'] \
                       for i in range(len(longest_path) - 1)]

        bond_orders_str = ''.join([BO_to_str[bo] for bo in bond_orders])

        if ('DD' not in bond_orders_str) and ('SS' not in bond_orders_str):
            conj_paths_candidates.append(longest_path)
    conj_paths_candidates.sort(key=len, reverse=True)

    conj_paths_final = []
    for conj_path in conj_paths_candidates:
        is_new = True
        for conj_path_fin in conj_paths_final:
            if len(set(conj_path) - set(conj_path_fin)) == 0:
                is_new = False
                break
        if is_new:
            conj_paths_final.append(conj_path)

    # print(conj_paths_final)
    conj_atoms = [int(x) for x in sorted(list(set(list(chain(*conj_paths_final)))))]
    return conj_paths_final, conj_atoms


def k_shortest_paths(G, source, target, k):
    return list(islice(nx.shortest_simple_paths(G, source, target, weight='weight'), k))


def nxGraph_from_smiles(smi):
    mol = Chem.MolFromSmiles(smi)
    edges_with_bond_orders = []
    edges = []

    hyb = [x.GetHybridization() for x in mol.GetAtoms()]
    # Chem.Kekulize(mol, clearAromaticFlags=True)

    # Conjugation: typically, a series of non-SP3 atoms
    # One exception: Cumulated systems (e.g., allene) - they will be ruled out later
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()

        if hyb[a1] != Chem.HybridizationType.SP3 and \
                hyb[a2] != Chem.HybridizationType.SP3:
            edges_with_bond_orders.append((a1, a2, bond.GetBondTypeAsDouble()))
            edges.append((a1, a2))

    # Check allylic radicals
    for allyl_center in mol.GetSubstructMatches(Chem.MolFromSmarts('[#6;X3v3+0]-[#6]=[#6X3]')):
        # Assign 1.5 bond order for the allylic radical C-C bond
        edges_with_bond_orders.append((allyl_center[0], allyl_center[1], 1.5))
        edges.append((allyl_center[0], allyl_center[1]))

    # just in case
    edges_with_bond_orders = list(set(edges_with_bond_orders))
    edges = list(set(edges))

    G = nx.Graph()
    G.add_weighted_edges_from(edges_with_bond_orders)

    G_no_W = nx.Graph()
    G_no_W.add_edges_from(edges)

    return G, G_no_W