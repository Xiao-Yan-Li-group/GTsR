import pymatgen.core as mg
from molSimplify.Informatics.MOF.MOF_descriptors import get_MOF_descriptors
import os, stat, shutil


def remove_dir_with_permissions(dir_path):
    def handle_permission_error(func, path, exc_info):
        os.chmod(path, stat.S_IWUSR)
        func(path)
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path, onerror=handle_permission_error)


def RACs(structure):
    os.makedirs("tmp_rac", exist_ok=True)
    name = os.path.basename(structure).replace(".cif", "")
    full_names, full_descriptors = get_MOF_descriptors(
        structure, 3,
        path='tmp_rac',
        xyz_path=f'tmp_rac/{name}.xyz',
        max_num_atoms=6000
    )
    descriptor_data = dict(zip(full_names, full_descriptors))
    remove_dir_with_permissions("tmp_rac")
    return descriptor_data


def get_cell(cif_path):
    structure = mg.Structure.from_file(cif_path)
    return structure.lattice.matrix


def flatten_rac(descriptor_data, suffix):
    return {f"{k}_{suffix}": v for k, v in descriptor_data.items()}


def flatten_cell(cell, suffix):
    flat = {}
    for r in range(3):
        for c in range(3):
            flat[f"cell_{r}{c}_{suffix}"] = cell[r, c]
    return flat
