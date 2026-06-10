from gemmi import cif
import numpy as np
from pymatgen.core import Structure


def get_label(cif_file):
    doc = cif.read_file(cif_file)
    block = doc.sole_block()
    labels = block.find_loop("_atom_site_type_symbol")
    return labels

def get_pos(cif_file):
    doc = cif.read_file(cif_file)
    block = doc.sole_block()
    x = block.find_loop("_atom_site_fract_x")
    y = block.find_loop("_atom_site_fract_y")
    z = block.find_loop("_atom_site_fract_z")
    pos =[]
    for i in range(len(x)):
        pos.append([float(x[i]),float(y[i]),float(z[i])])
    return pos

def get_charge(cif_file):
    doc = cif.read_file(cif_file)
    block = doc.sole_block()
    charges = block.find_loop("_atom_site_charge")
    return charges

def unit_cell(cif_file, cutoff):
    struct = Structure.from_file(cif_file)
    lattice = struct.lattice
    A = lattice.matrix[0]
    B = lattice.matrix[1]
    C = lattice.matrix[2]
    Wa = np.divide(np.linalg.norm(np.dot(np.cross(B, C), A)), np.linalg.norm(np.cross(B, C)))
    Wb = np.divide(np.linalg.norm(np.dot(np.cross(C, A), B)), np.linalg.norm(np.cross(C, A)))
    Wc = np.divide(np.linalg.norm(np.dot(np.cross(A, B), C)), np.linalg.norm(np.cross(A, B)))
    uc_x = int(np.ceil(cutoff / (0.5 * Wa)))
    uc_y = int(np.ceil(cutoff / (0.5 * Wb)))
    uc_z = int(np.ceil(cutoff / (0.5 * Wc)))
    return [uc_x, uc_y, uc_z]