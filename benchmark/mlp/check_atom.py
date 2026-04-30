import numpy as np
from gemmi import cif
import os
import pandas as pd
from tqdm import tqdm


def get_lattice_matrix(block):
    a = float(block.find_value("_cell_length_a").split("(")[0])
    b = float(block.find_value("_cell_length_b").split("(")[0])
    c = float(block.find_value("_cell_length_c").split("(")[0])
    alpha = np.radians(float(block.find_value("_cell_angle_alpha").split("(")[0]))
    beta  = np.radians(float(block.find_value("_cell_angle_beta").split("(")[0]))
    gamma = np.radians(float(block.find_value("_cell_angle_gamma").split("(")[0]))

    cx = np.cos(beta)
    cy = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(max(1.0 - cx**2 - cy**2, 0.0))

    matrix = np.array([
        [a,                b * np.cos(gamma), c * cx],
        [0,                b * np.sin(gamma), c * cy],
        [0,                0,                 c * cz],
    ])
    return matrix


def read_structure(cif_path):
    doc = cif.read_file(cif_path)
    block = doc.sole_block()

    lattice = get_lattice_matrix(block)

    labels  = list(block.find_loop("_atom_site_label"))
    x = [float(v.split("(")[0]) for v in block.find_loop("_atom_site_fract_x")]
    y = [float(v.split("(")[0]) for v in block.find_loop("_atom_site_fract_y")]
    z = [float(v.split("(")[0]) for v in block.find_loop("_atom_site_fract_z")]
    frac_coords = np.array(list(zip(x, y, z)))  # (N, 3)

    charge_loop = block.find_loop("_atom_site_charge")
    if len(charge_loop) > 0:
        charges = np.array([float(v) for v in charge_loop])
    else:
        charges = None

    return labels, frac_coords, charges, lattice


def frac_displacement_cart(frac1, frac2, lattice):
    delta = frac2 - frac1
    delta -= np.round(delta)
    delta_cart = delta @ lattice
    distances = np.linalg.norm(delta_cart, axis=1)
    return distances


def compare(cif_ref, cif_opt):
    labels1, frac1, charges1, lat1 = read_structure(cif_ref)
    labels2, frac2, charges2, _    = read_structure(cif_opt)

    if len(labels1) != len(labels2):
        print(f"[WARNING] Atom count mismatch: {len(labels1)} vs {len(labels2)}")
        return None

    displacements = frac_displacement_cart(frac1, frac2, lat1)

    result = {
        "labels": labels1,
        "displacement_A": displacements,
        "mean_displacement": float(np.mean(displacements)),
        "max_displacement": float(np.max(displacements)),
    }

    if charges1 is not None and charges2 is not None:
        charge_diff = charges2 - charges1
        result["charge_diff"] = charge_diff
        result["mean_abs_charge_diff"] = float(np.mean(np.abs(charge_diff)))
        result["max_abs_charge_diff"] = float(np.max(np.abs(charge_diff)))
    else:
        result["charge_diff"] = None
        result["mean_abs_charge_diff"] = None
        result["max_abs_charge_diff"] = None

    return result



folder = "./pacman_data"
types = ["n", "y"]


for t in types:
    n_data = pd.read_csv("./results_"+t+"/mlp_results.csv")

    data = []
    for cifname in tqdm(n_data["name"][:]):
        ref_cif = os.path.join(folder, cifname + "_pacman.cif")
        opt_cif = os.path.join(folder, cifname + "_opt_pacman.cif")

        result = compare(ref_cif, opt_cif)

        if result is not None:
            row = [
                cifname,
                result["mean_displacement"],
                result["max_displacement"],
                result["mean_abs_charge_diff"],
                result["max_abs_charge_diff"],
            ]
            print(
                f"{cifname}  mean_disp={result['mean_displacement']:.4f} Å  "
                f"max_disp={result['max_displacement']:.4f} Å  "
                f"mean_Δq={result['mean_abs_charge_diff']}  "
                f"max_Δq={result['max_abs_charge_diff']}"
            )
        else:
            row = [cifname, None, None, None, None]
            print(f"{cifname}  compare failed")

        data.append(row)

    pd.DataFrame(
        data,
        columns=["name", "mean_displacement_A", "max_displacement_A",
                "mean_abs_charge_diff", "max_abs_charge_diff"],
    ).to_csv("./results_"+t+"/atom_diff.csv", index=False)
