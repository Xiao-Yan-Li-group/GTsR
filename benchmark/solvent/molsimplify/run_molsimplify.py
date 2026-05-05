from molSimplify.Informatics.MOF.PBC_functions import solvent_removal
from glob import glob
import os
import sys
import time
import multiprocessing
from tqdm import tqdm
from pathlib import Path


# memory killed


TIMEOUT = 300


def _run(cif_path, save_path, wiggle_room, error_q):
    try:
        solvent_removal(cif_path=cif_path, new_cif_path=save_path, wiggle_room=wiggle_room)
    except Exception as e:
        error_q.put(str(e))


def process_one(cif_path, save_path, wiggle_room):
    """Run solvent_removal in a subprocess; return error string or None on success."""
    error_q = multiprocessing.Queue()
    p = multiprocessing.Process(target=_run, args=(cif_path, save_path, wiggle_room, error_q))
    p.start()
    p.join(TIMEOUT)

    if p.is_alive():
        p.kill()
        p.join()
        return "killed: timeout/OOM"

    if p.exitcode != 0:
        try:
            return error_q.get_nowait()
        except Exception:
            return f"killed: exitcode={p.exitcode}"

    return None


if __name__ == "__main__":
    cif_folder = sys.argv[1]
    save_folder = sys.argv[2]
    wiggle_room = 1

    os.makedirs(save_folder, exist_ok=True)

    fails = []
    t_start = time.time()

    for cif_path in tqdm(sorted(glob(os.path.join(cif_folder, "*.cif")))):
        name = Path(cif_path).stem
        new_path = os.path.join(save_folder, f"{name}_molsimplify_{wiggle_room}_F.cif")

        if os.path.exists(new_path):
            continue

        err = process_one(cif_path, new_path, wiggle_room)
        if err:
            fails.append((name, err))

    with open("failed_list.txt", "w", encoding="utf-8") as f:
        f.write(f"time: {time.time() - t_start}\n")
        for name, err in fails:
            f.write(f"{name}\t{err}\n")
