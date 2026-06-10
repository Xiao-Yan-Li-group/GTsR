# -*- coding: utf-8 -*-
import os
import sys
import csv
import time
import logging
import numpy as np
from pathlib import Path
from glob import glob

from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.optimize import LBFGS
from ase.filters import FrechetCellFilter
from ase.md.langevin import Langevin
from ase.md.nptberendsen import NPTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units

from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
from pymatgen.io.ase import AseAtomsAdaptor

from huggingface_hub import login
from fairchem.core import pretrained_mlip, FAIRChemCalculator

_bar = getattr(units, "bar", 1e-4 * units.GPa)

PARAMS = {
    "opt": {
        "fmax": 1e-3,
        "steps": 5000,
        "optimizer": "LBFGS",
        "filter": FrechetCellFilter,
    },
    "nvt": {
        "total_steps": 1000,
        "trajectory_interval": 100,
        "log_interval": 500,
        "timestep_fs": 1,
        "temperature_K": 300,
        "friction_per_fs": 0.01,
    },
    "npt": {
        "total_steps": 50000,
        "trajectory_interval": 100,
        "log_interval": 500,
        "timestep_fs": 1,
        "temperature_K": 300,
        "externalstress_bar": 1,
        "ttime_fs": 100,
        "ptime_fs": 1000,
        "bulk_modulus_GPa": 5,
    },
}


# ─────────────────────────── utilities ────────────────────────────────────────

def setup_logger(log_path: str) -> logging.Logger:
    logger = logging.getLogger(log_path)
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(log_path, mode="w")
    fh.setFormatter(logging.Formatter("%(asctime)s  %(message)s", "%H:%M:%S"))
    logger.addHandler(fh)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(sh)
    return logger


_sm = StructureMatcher(
    comparator=ElementComparator(),
    ltol=0.2, stol=0.3, angle_tol=5,
    primitive_cell=True, scale=True,
    attempt_supercell=False,
)

def compute_rmsd(atoms_ref, atoms_fin) -> float:
    try:
        s1 = AseAtomsAdaptor.get_structure(atoms_ref)
        s2 = AseAtomsAdaptor.get_structure(atoms_fin)
        d = _sm.get_rms_dist(s1, s2)
        return float(d[0]) if d is not None else float("nan")
    except Exception:
        return float("nan")


def _fmax(forces: np.ndarray) -> float:
    return float(np.sqrt((forces ** 2).sum(axis=1).max()))


def _write_csv(path: str, rows: list, fieldnames: list) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


# ─────────────────────────── OPT ──────────────────────────────────────────────

def run_opt(atoms, out_dir: Path, params: dict) -> None:

    p = params["opt"]
    traj_path = str(out_dir / "opt.traj")
    csv_path  = str(out_dir / "opt.csv")

    filtered = p["filter"](atoms)
    opt = LBFGS(filtered, trajectory=traj_path, logfile=os.devnull)

    rows: list = []
    t0 = time.time()

    def _cb():
        step = opt.get_number_of_steps()
        e    = atoms.get_potential_energy()
        fm   = _fmax(filtered.get_forces())
        rows.append({
            "step":      step,
            "time_s":    round(time.time() - t0, 3),
            "energy_eV": round(e, 6),
            "fmax_eV_A": round(fm, 6),
        })

    opt.attach(_cb, interval=1)

    try:
        opt.run(fmax=p["fmax"], steps=p["steps"])
    except Exception:
        pass

    _write_csv(csv_path, rows, ["step", "time_s", "energy_eV", "fmax_eV_A"])


# ─────────────────────────── NVT ──────────────────────────────────────────────

def run_nvt(atoms, out_dir: Path, params: dict) -> None:
    p = params["nvt"]
    traj_path = str(out_dir / "nvt.traj")
    csv_path  = str(out_dir / "nvt.csv")

    MaxwellBoltzmannDistribution(atoms, temperature_K=p["temperature_K"])

    dyn = Langevin(
        atoms,
        timestep=p["timestep_fs"] * units.fs,
        temperature_K=p["temperature_K"],
        friction=p["friction_per_fs"] / units.fs,
        fixcm=False,
    )

    traj = Trajectory(traj_path, "w", atoms)
    dyn.attach(traj.write, interval=p["trajectory_interval"])

    rows: list = []

    def _log():
        step = dyn.get_number_of_steps()
        e    = atoms.get_potential_energy()
        fm   = _fmax(atoms.get_forces())
        rows.append({
            "step":      step,
            "time_fs":   step * p["timestep_fs"],
            "energy_eV": round(e, 6),
            "fmax_eV_A": round(fm, 6),
        })

    dyn.attach(_log, interval=p["log_interval"])
    dyn.run(p["total_steps"])

    final = dyn.get_number_of_steps()
    if not rows or rows[-1]["step"] != final:
        _log()

    traj.close()
    _write_csv(csv_path, rows, ["step", "time_fs", "energy_eV", "fmax_eV_A"])


# ─────────────────────────── NPT ──────────────────────────────────────────────

def run_npt(atoms, out_dir: Path, params: dict) -> None:
    p = params["npt"]
    traj_path = str(out_dir / "npt.traj")
    csv_path  = str(out_dir / "npt.csv")

    dyn = NPTBerendsen(
        atoms,
        timestep=p["timestep_fs"] * units.fs,
        temperature_K=p["temperature_K"],
        pressure_au=p["externalstress_bar"] * _bar,
        taut=p["ttime_fs"] * units.fs,
        taup=p["ptime_fs"] * units.fs,
        compressibility_au=1.0 / (p["bulk_modulus_GPa"] * units.GPa),
    )

    traj = Trajectory(traj_path, "w", atoms)
    dyn.attach(traj.write, interval=p["trajectory_interval"])

    rows: list = []

    def _log():
        step = dyn.get_number_of_steps()
        e    = atoms.get_potential_energy()
        fm   = _fmax(atoms.get_forces())
        rows.append({
            "step":      step,
            "time_fs":   step * p["timestep_fs"],
            "energy_eV": round(e, 6),
            "fmax_eV_A": round(fm, 6),
        })

    dyn.attach(_log, interval=p["log_interval"])
    dyn.run(p["total_steps"])

    final = dyn.get_number_of_steps()
    if not rows or rows[-1]["step"] != final:
        _log()

    traj.close()
    _write_csv(csv_path, rows, ["step", "time_fs", "energy_eV", "fmax_eV_A"])



def simulate_cif(cif_path: str, out_dir: str, calc, params: dict) -> dict:
    nan = float("nan")
    name    = Path(cif_path).stem
    out_dir = Path(out_dir) / name
    out_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logger(str(out_dir / f"{name}.log"))
    logger.info(f"Structure : {cif_path}")

    atoms = read(cif_path)
    atoms.calc = calc
    n_atoms = len(atoms)
    logger.info(f"N atoms   : {n_atoms}")

    result = {
        "name":       name,
        "init_E":     nan,
        "opt_E":      nan,
        "opt_RMSD":   nan,
        "init_NPT_E": nan,
        "NPT_E":      nan,
        "NPT_RMSD":   nan,
        "N_atoms":    n_atoms,
    }

    atoms_init = atoms.copy()
    result["init_E"] = atoms.get_potential_energy()
    logger.info(f"init E    : {result['init_E']:.6f} eV")

    logger.info("=== OPT ===")
    try:
        run_opt(atoms, out_dir, params)
        result["opt_E"]    = atoms.get_potential_energy()
        result["opt_RMSD"] = compute_rmsd(atoms_init, atoms)
        write(str(out_dir / f"{name}_opt.cif"), atoms)
        logger.info(f"opt E     : {result['opt_E']:.6f} eV  |  RMSD: {result['opt_RMSD']}")
    except Exception as e:
        logger.error(f"OPT failed: {e}")

    logger.info("=== NVT ===")
    try:
        run_nvt(atoms, out_dir, params)
        logger.info("NVT done")
    except Exception as e:
        logger.error(f"NVT failed: {e}")

    logger.info("=== NPT ===")
    try:
        atoms_pre_npt        = atoms.copy()
        result["init_NPT_E"] = atoms.get_potential_energy()
        run_npt(atoms, out_dir, params)
        result["NPT_E"]    = atoms.get_potential_energy()
        result["NPT_RMSD"] = compute_rmsd(atoms_pre_npt, atoms)
        write(str(out_dir / f"{name}_npt.cif"), atoms)
        logger.info(f"NPT E     : {result['NPT_E']:.6f} eV  |  RMSD: {result['NPT_RMSD']}")
    except Exception as e:
        logger.error(f"NPT failed: {e}")

    return result


def main():
    import pandas as pd
    from tqdm import tqdm

    cif_dir = sys.argv[1]
    out_dir = sys.argv[2]

    cif_files = sorted(glob(os.path.join(cif_dir, "*.cif")))
    login(token="")
    print("Loading uma-s-1.1 ...")
    predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
    calc = FAIRChemCalculator(predictor, task_name="odac")

    all_results: list = []
    failed:      list = []

    for cif_path in tqdm(cif_files):
        name = Path(cif_path).stem
        try:
            result = simulate_cif(
                cif_path=cif_path,
                out_dir=out_dir,
                calc=calc,
                params=PARAMS,
            )
            all_results.append(result)
        except Exception as e:
            print(f"  FAILED: {name}  →  {e}")
            failed.append({"name": name, "error": str(e)})

    if all_results:
        cols = ["name", "init_E", "opt_E", "opt_RMSD",
                "init_NPT_E", "NPT_E", "NPT_RMSD", "N_atoms"]
        pd.DataFrame(all_results)[cols].to_csv(
            os.path.join(out_dir, "mlp_results.csv"), index=False)

    if failed:
        pd.DataFrame(failed).to_csv(
            os.path.join(out_dir, "failed.csv"), index=False)
        print(f"Failed list -> {out_dir}/failed.csv")


if __name__ == "__main__":
    main()
