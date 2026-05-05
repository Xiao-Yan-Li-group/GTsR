"pip install git+https://github.com/sxm13/mofchecker_2.0.git@main"

from mofchecker import MOFChecker


def check(cif_path):
    checker = MOFChecker.from_cif(cif_path)
    check_result = checker.get_mof_descriptors()

    error = "has_lone_molecule"

    if check_result.get(error, False):
        return True
    else:
        return False


import sys
from tqdm import tqdm
from glob import glob
from pathlib import Path
from multiprocessing import Pool
import os


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
