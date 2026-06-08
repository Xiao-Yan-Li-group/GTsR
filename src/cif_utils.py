from __future__ import annotations

import os
from pathlib import Path
from pymatgen.core import Structure
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from ase.data import atomic_numbers
import numpy as np
from pymatgen.io.cif import CifWriter


def cif2graph(
                cif_path: str | Path,
                radius: float = 8.0,
                max_neighbors: int = 12,
            ) -> dict:
    
    structure = read(cif_path)
    struct = AseAtomsAdaptor.get_structure(structure)

    ele_list = [sp.symbol for sp in struct.species]
    numbers = [atomic_numbers[sym] for sym in ele_list]
    c_index, n_index, _, dists = struct.get_neighbor_list(
        r=radius, numerical_tol=0, exclude_self=True)
    n_sites = len(struct)
    selected_idx_chunks = []
    for i in range(n_sites):
        mask = (c_index == i)
        if np.any(mask):
            idx = np.where(mask)[0]
            order = np.argsort(dists[idx])[:max_neighbors]
            selected_idx_chunks.append(idx[order])
    if selected_idx_chunks:
        selected_idx = np.concatenate(selected_idx_chunks)
        index1 = c_index[selected_idx]
        index2 = n_index[selected_idx]
        dij    = dists[selected_idx]
    else:
        index1 = np.array([], dtype=int)
        index2 = np.array([], dtype=int)
        dij    = np.array([], dtype=float)

    return {
            "rcut": float(radius),
            "max_neighbors": int(max_neighbors),
            "numbers": numbers,
            "index1": index1.tolist(),
            "index2": index2.tolist(),
            "dij": dij.tolist()
            }
        

def cif2pos(cif_path: str | Path) -> np.ndarray:
    structure = read(cif_path)
    return structure.get_positions()


def make_label(cif_path1: str | Path, cif_path2: str | Path, tol=1e-3) -> np.ndarray:

    struct1 = Structure.from_file(cif_path1)
    struct2 = Structure.from_file(cif_path2)

    dists = struct1.lattice.get_all_distances(struct1.frac_coords, struct2.frac_coords)
    min_dists = np.min(dists, axis=1)
    labels = np.where(min_dists < tol, 0, 1)
        
    return labels
    

def label2cif(cif_path1: str | Path, label: np.ndarray, output_dir: str) -> Structure:
    struct = Structure.from_file(cif_path1)
    
    solvent = label == 1
    framework = ~solvent

    def extract(mask):
        indices = np.where(mask)[0]
        sites   = [struct[i] for i in indices]
        return Structure.from_sites(sites)
    
    os.makedirs(output_dir, exist_ok=True)

    if sum(label) == 0:
        framework_atoms = extract(framework)
        output_path = os.path.join(output_dir, Path(cif_path1).stem+"_framework.cif")
        CifWriter(framework_atoms).write_file(output_path)
    else:
        solvent_atoms = extract(solvent)
        framework_atoms = extract(framework)
        solvent_path = os.path.join(output_dir, Path(cif_path1).stem+"_solvent.cif")
        framework_path = os.path.join(output_dir, Path(cif_path1).stem+"_framework.cif")
        CifWriter(solvent_atoms).write_file(solvent_path)
        CifWriter(framework_atoms).write_file(framework_path)
