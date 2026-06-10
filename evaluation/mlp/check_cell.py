from ase.io import read
import os
import pandas as pd
from tqdm import tqdm


def cell(cif_ref, cif_opt):
    atoms0 = read(cif_ref)
    atoms1 = read(cif_opt)

    m0 = atoms0.get_cell()[:]
    m1 = atoms1.get_cell()[:]

    vol0 = atoms0.get_volume()
    vol1 = atoms1.get_volume()

    row = {}
    for i in range(3):
        for j in range(3):
            row[f"cell_{i}{j}_0"] = m0[i, j]
    for i in range(3):
        for j in range(3):
            row[f"cell_{i}{j}_1"] = m1[i, j]

    row["volume_0"] = vol0
    row["volume_1"] = vol1

    return row


types = ["n", "y"]

for t in types:
    folder1 = f"./clean_{t}"
    folder2 = f"./results_{t}/"
    n_data = pd.read_csv(f"./results_{t}/mlp_results.csv")

    data = []
    for cifname in tqdm(n_data["name"][:]):
        cif_ref = os.path.join(folder1, cifname + ".cif")
        cif_opt = os.path.join(folder2, cifname, cifname + "_opt.cif")
        try:
            row = cell(cif_ref, cif_opt)
            row["name"] = cifname
        except Exception as e:
            print(f"{cifname} failed: {e}")
            row = {"name": cifname}
        data.append(row)

    cols = ["name"]
    for suffix in ["_0", "_1"]:
        for i in range(3):
            for j in range(3):
                cols.append(f"cell_{i}{j}{suffix}")
    cols += ["volume_0", "volume_1"]

    pd.DataFrame(data, columns=cols).to_csv(f"./results_{t}/cell_diff.csv", index=False)
