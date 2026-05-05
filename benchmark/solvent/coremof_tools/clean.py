import os, functools, itertools, collections
from ase.io import read
from scipy.sparse.csgraph import connected_components
from ase.neighborlist import NeighborList
import numpy as np
from utils import ALLIONS, COVALENTRADII, METAL


ions_list = set(ALLIONS)
metal_list = [element for element, is_metal in METAL.items() if is_metal]

def build_ASE_neighborlist(cif, skin):
    radii = [COVALENTRADII[i] for i in cif.get_chemical_symbols()]
    ASE_neighborlist = NeighborList(radii, self_interaction=False, bothways=True, skin=skin)
    ASE_neighborlist.update(cif)
    return ASE_neighborlist

def find_clusters(adjacency_matrix, atom_count):
    clusters = []
    cluster_count, clusterIDs = connected_components(adjacency_matrix, directed=True)
    for n in range(cluster_count):
        clusters.append([i for i in range(atom_count) if clusterIDs[i] == n])
    return clusters

def find_metal_connected_atoms(structure, neighborlist):
    metal_connected_atoms = []
    metal_atoms = []
    for i, elem in enumerate(structure.get_chemical_symbols()):
        if elem in metal_list:
            neighbors, _ = neighborlist.get_neighbors(i)
            metal_connected_atoms.append(neighbors)
            metal_atoms.append(i)
    return metal_connected_atoms, metal_atoms, structure

def CustomMatrix(neighborlist, atom_count):
    matrix = np.zeros((atom_count, atom_count), dtype=int)
    for i in range(atom_count):
        neighbors, _ = neighborlist.get_neighbors(i)
        for j in neighbors:
            matrix[i][j] = 1
    return matrix

def mod_adjacency_matrix(adj_matrix, MetalConAtoms, MetalAtoms, atom_count, struct):
    clusters = find_clusters(adj_matrix, atom_count)
    for i, element_1 in enumerate(MetalAtoms):
        for _, element_2 in enumerate(MetalConAtoms[i]):
            if struct[element_2].symbol == "O":
                tmp = len(find_clusters(adj_matrix, atom_count))
                adj_matrix[element_2][element_1] = 0
                adj_matrix[element_1][element_2] = 0
                new_clusters = find_clusters(adj_matrix, atom_count)
                if tmp == len(new_clusters):
                    adj_matrix[element_2][element_1] = 1
                    adj_matrix[element_1][element_2] = 1
                for ligand in new_clusters:
                    if ligand not in clusters:
                        tmp3 = struct[ligand].get_chemical_symbols()
                        if "O" and "H" in tmp3 and len(tmp3) == 2:
                            adj_matrix[element_2][element_1] = 1
                            adj_matrix[element_1][element_2] = 1
    return adj_matrix

def cmp(x, y):
    return (x > y) - (x < y)

def cluster_to_formula(cluster, cif):
    symbols = [cif[i].symbol for i in cluster]
    count = collections.Counter(symbols)
    formula = ''.join([atom + (str(count[atom]) if count[atom] > 1 else '') for atom in sorted(count)])
    return formula

def free_clean(cif_path, skin):
    cif = read(cif_path)
    atom_count = len(cif.get_chemical_symbols())
    ASE_neighborlist = build_ASE_neighborlist(cif, skin)
    a = CustomMatrix(ASE_neighborlist,atom_count)
    b = find_clusters(a,atom_count)
    b.sort(key=functools.cmp_to_key(lambda x,y: cmp(len(x), len(y))))
    b.reverse()
    cluster_length=[]
    solvated_cluster = []
    ions_cluster = []
    printed_formulas = []
    for index, _ in enumerate(b):
        cluster_formula = cluster_to_formula(b[index], cif) 
        if cluster_formula in ions_list:
            ions_cluster.append(b[index])
        else:
            tmp = len(b[index])
            if len(cluster_length) > 0:
                if tmp > max(cluster_length):
                    cluster_length = []
                    solvated_cluster = []
                    solvated_cluster.append(b[index])
                    cluster_length.append(tmp)
                elif tmp > 0.5 * max(cluster_length):
                    solvated_cluster.append(b[index])
                    cluster_length.append(tmp)
                else:
                    formula = cluster_to_formula(b[index], cif)
                    if formula not in printed_formulas:
                        printed_formulas.append(formula)
            else:
                solvated_cluster.append(b[index])
                cluster_length.append(tmp)
    solvated_cluster = solvated_cluster + ions_cluster
    solvated_merged = list(itertools.chain.from_iterable(solvated_cluster))
    atom_count = len(cif[solvated_merged].get_chemical_symbols())
    return printed_formulas, cif[solvated_merged]

def all_clean(cif_path, skin):
    cif = read(cif_path)
    atom_count = len(cif.get_chemical_symbols())
    ASE_neighborlist = build_ASE_neighborlist(cif,skin)
    a = CustomMatrix(ASE_neighborlist,atom_count)
    b = find_clusters(a,atom_count)
    b.sort(key=functools.cmp_to_key(lambda x,y: cmp(len(x), len(y))))
    b.reverse()
    cluster_length=[]
    solvated_cluster = []
    ions_cluster = []
    printed_formulas = []
    for index, _ in enumerate(b):
        cluster_formula = cluster_to_formula(b[index], cif) 
        if cluster_formula in ions_list:
            ions_cluster.append(b[index])
            solvated_cluster.append(b[index])
        else:
            tmp = len(b[index])
            if len(cluster_length) > 0:
                if tmp > max(cluster_length):
                    cluster_length = []
                    solvated_cluster = []
                    solvated_cluster.append(b[index])
                    cluster_length.append(tmp)
                if tmp > 0.5 * max(cluster_length):
                    solvated_cluster.append(b[index])
                    cluster_length.append(tmp)
                else:
                    formula = cluster_to_formula(b[index], cif)
                    if formula not in printed_formulas:
                        printed_formulas.append(formula)
            else:
                solvated_cluster.append(b[index])
                cluster_length.append(tmp)
    solvated_merged = list(itertools.chain.from_iterable(solvated_cluster))
    atom_count = len(cif[solvated_merged].get_chemical_symbols())
    newASE_neighborlist = build_ASE_neighborlist(cif[solvated_merged],skin)
    MetalCon, MetalAtoms, struct = find_metal_connected_atoms(cif[solvated_merged], newASE_neighborlist)
    c = CustomMatrix(newASE_neighborlist,atom_count)
    d = mod_adjacency_matrix(c, MetalCon, MetalAtoms,atom_count,struct)
    solvated_clusters2 = find_clusters(d,atom_count)
    solvated_clusters2.sort(key=functools.cmp_to_key(lambda x,y: cmp(len(x), len(y))))
    solvated_clusters2.reverse()
    cluster_length=[]
    final_clusters = []
    for index, _ in enumerate(solvated_clusters2):
        cluster_formula2 = cluster_to_formula(solvated_clusters2[index], struct) 
        if cluster_formula2 in ions_list:
            final_clusters.append(solvated_clusters2[index])
        else:
            tmp = len(solvated_clusters2[index])
            if len(cluster_length) > 0:
                if tmp > max(cluster_length):
                    cluster_length = []
                    final_clusters = []
                    final_clusters.append(solvated_clusters2[index])
                    cluster_length.append(tmp)
                if tmp > 0.5 * max(cluster_length):
                    final_clusters.append(solvated_clusters2[index])
                    cluster_length.append(tmp)
                else:
                    formula = cluster_to_formula(solvated_clusters2[index], struct)
                    if formula not in printed_formulas:
                        printed_formulas.append(formula)
            else:
                final_clusters.append(solvated_clusters2[index])
                cluster_length.append(tmp)
    final_clusters = final_clusters
    final_merged = list(itertools.chain.from_iterable(final_clusters))
    tmp = struct[final_merged].get_chemical_symbols()
    tmp.sort()
    
    return printed_formulas, struct[final_merged]