ATR = {
        'H':0.38,
        'Li':0.86,
        'Be':0.53,
        'B':1.01,
        'C':0.88,
        'N':0.86,
        'O':0.89,
        'F':0.82,
        'Na':1.15,
        'Mg':1.28,
        'Al':1.53,
        'Si':1.38,
        'P':1.28,
        'S':1.20,
        'Cl':1.17,
        'K':1.44,
        'Ca':1.17,
        'Sc':1.62,
        'Ti':1.65,
        'V':1.51,
        'Cr':1.53,
        'Mn':1.53,
        'Fe':1.43,
        'Co':1.31,
        'Ni':1.33,
        'Cu':1.31,
        'Zn':1.41,
        'Ga':1.40,
        'Ge':1.35,
        'As':1.39,
        'Se':1.40,
        'Br':1.39,
        'Rb':1.65,
        'Sr':1.30,
        'Y':1.84,
        'Zr':1.73,
        'Nb':1.66,
        'Mo':1.57,
        'Ru':1.58,
        'Rh':1.63,
        'Pd':1.68,
        'Ag':1.56,
        'Cd':1.56,
        'In':1.53,
        'Sn':1.64,
        'Sb':1.64,
        'Te':1.65,
        'I':1.58,
        'Cs':1.85,
        'Ba':1.52,
        'La':1.91,
        'Ce':1.98,
        'Pr':1.75,
        'Nd':1.92,
        'Sm':1.89,
        'Eu':1.83,
        'Gd':1.79,
        'Tb':1.82,
        'Dy':1.79,
        'Ho':1.63,
        'Er':1.80,
        'Tm':1.84,
        'Yb':1.80,
        'Lu':1.86,
        'Hf':1.73,
        'W':1.33,
        'Re':1.29,
        'Ir':1.50,
        'Pt':1.66,
        'Au':1.68,
        'Hg':1.88,
        'Pb':1.72,
        'Bi':1.72,
        'Th':1.97,
        'U':1.76,
        'Np':1.73,
        'Pu':1.71
        }

from ase.io import read
import sys
import os
from tqdm import tqdm
from glob import glob
from pathlib import Path
from multiprocessing import Pool


def check(cif_path):
    atoms = read(cif_path)
    sym = atoms.get_chemical_symbols()
    for a in range(len(atoms)):
        nl = []
        for b in range(len(atoms)):
            if a == b:
                continue
            d = atoms.get_distance(a, b, mic=True)
            if d <= (ATR[sym[a]] + ATR[sym[b]]):
                nl.append(b)
        if len(nl) == 0:
            return True
    return False


def worker(cif_path):
    try:
        result = check(cif_path)
        return ("free" if result else "ok", Path(cif_path).stem)
    except Exception:
        return ("error", Path(cif_path).stem)


if __name__ == "__main__":
    cif_folder = sys.argv[1]
    n_processes = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()

    cif_files = glob(cif_folder + "/*.cif")

    with Pool(processes=n_processes) as pool:
        results = list(tqdm(pool.imap_unordered(worker, cif_files), total=len(cif_files)))

    has_free = [name for status, name in results if status == "free"]
    errors   = [name for status, name in results if status == "error"]

    with open("has_free_list.txt", "w") as f:
        for name in has_free:
            f.write(f"{name}\n")

    with open("error_list.txt", "w") as f:
        for name in errors:
            f.write(f"{name}\n")
