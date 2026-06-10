import pandas as pd
from tqdm import tqdm
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

PORE_COLS = ["LCD", "PLD", "LFPD", "volume", "density", "PV"]

def load_pore(csv_path):
    df = pd.read_csv(csv_path)
    df["suffix"] = df["name"].str.extract(r"_([01])$")
    df["mof"]    = df["name"].str.replace(r"_[01]$", "", regex=True)

    pore_list = []
    for suffix in ["0", "1"]:
        sub = df[df["suffix"] == suffix].set_index("mof")[PORE_COLS]
        sub.columns = [f"{c}_{suffix}" for c in PORE_COLS]
        pore_list.append(sub)

    return pd.concat(pore_list, axis=1)


pore     = load_pore("pores.csv")
data     = pd.read_csv("/mnt/d/Project/GTSR/benchmark/stability/ml/data/dataset_clean.csv")
cif_path = "/mnt/d/Project/GTSR/benchmark/stability/md/clean/npt_cifs"

rows = []

for i, name in tqdm(enumerate(data["filename"][:2]), total=len(data)):
    target = data["target"][i]
    folder = "Y" if target == 1 else "N"

    cif_0 = os.path.join(cif_path, folder, name + "_0.cif")
    cif_1 = os.path.join(cif_path, folder, name + "_1.cif")
    cif_ori = os.path.join("/scratch/guobinzhao/calculation/ASMR/md/structures/", folder, name + ".cif")

    try:
        rac  = RACs(structure=cif_ori)
        cell_0 = get_cell(cif_0)
        cell_1 = get_cell(cif_1)
    except Exception as e:
        print(f"{name} failed")
        continue

    row = {"filename": name}
    row.update(flatten_rac(rac,  suffix="0"))
    row.update(flatten_cell(cell_0, suffix="0"))
    row.update(flatten_cell(cell_1, suffix="1"))

    pore_key = name
    for col in pore.columns:
        row[col] = pore.loc[pore_key, col]

    row["target"] = target
    rows.append(row)

df_out = pd.DataFrame(rows)
df_out.to_csv("features.csv", index=False)