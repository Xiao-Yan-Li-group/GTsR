# MLP Benchmark — Structure Relaxation & MD with UMA

This module benchmarks porous crystal structures using the **UMA-S-1.1** machine learning interatomic potential (MLIP) from [FairChem](https://github.com/FAIR-Chem/fairchem). For each input CIF file, it runs a three-stage simulation pipeline: geometry optimization → NVT MD → NPT MD, and outputs energy/force trajectories plus a summary CSV.

---

## Pipeline Overview

```
CIF file
   │
   ├─ [OPT]  LBFGS geometry optimization  →  *_opt.cif  +  opt.csv / opt.traj
   │
   ├─ [NVT]  Langevin MD (300 K)          →  nvt.csv  /  nvt.traj
   │
   └─ [NPT]  Berendsen NPT (300 K, 1 bar) →  *_npt.cif  +  npt.csv / npt.traj
```

Results for all structures are aggregated into `mlp_results.csv`.

---

## Quick Start

```bash
# 1. Create and activate the environment
conda env create -f environment.yml
conda activate uma

# 2. Run the benchmark
python run_mlp.py <cif_dir> <out_dir>
```

| Argument | Description |
|---|---|
| `cif_dir` | Directory containing input `.cif` files |
| `out_dir` | Root output directory (subdirectory per structure) |

**Example**

```bash
python run_mlp.py data/cifs/ results/mlp/
```

---

## Output Structure

```
out_dir/
├── mlp_results.csv          # summary table for all structures
├── failed.csv               # structures that raised exceptions (if any)
└── <structure_name>/
    ├── <name>.log           # per-structure log
    ├── opt.csv              # step, time, energy, fmax  (OPT)
    ├── opt.traj             # ASE trajectory (OPT)
    ├── <name>_opt.cif       # optimized structure
    ├── nvt.csv              # step, time_fs, energy, fmax  (NVT)
    ├── nvt.traj             # ASE trajectory (NVT)
    ├── npt.csv              # step, time_fs, energy, fmax  (NPT)
    ├── npt.traj             # ASE trajectory (NPT)
    └── <name>_npt.cif       # final NPT structure
```

### `mlp_results.csv` columns

| Column | Description |
|---|---|
| `name` | Structure name (CIF stem) |
| `N_atoms` | Number of atoms in the unit cell |
| `init_E` | Potential energy before optimization (eV) |
| `opt_E` | Potential energy after OPT (eV) |
| `opt_RMSD` | RMS displacement OPT vs. initial structure (fractional) |
| `init_NPT_E` | Energy at start of NPT (after NVT) (eV) |
| `NPT_E` | Energy at end of NPT (eV) |
| `NPT_RMSD` | RMS displacement NPT final vs. NPT initial (fractional) |

---

## Simulation Parameters

All parameters are defined in `PARAMS` at the top of [run_mlp.py](run_mlp.py) and can be edited directly.

### Geometry Optimization (`opt`)

| Parameter | Default | Description |
|---|---|---|
| `fmax` | `1e-3 eV/Å` | Convergence threshold on max force |
| `steps` | `5000` | Maximum number of LBFGS steps |
| `optimizer` | `LBFGS` | ASE optimizer class |
| `filter` | `FrechetCellFilter` | Allows cell shape/volume to relax |

### NVT Molecular Dynamics (`nvt`)

| Parameter | Default | Description |
|---|---|---|
| `total_steps` | `1000` | Total MD steps |
| `timestep_fs` | `1 fs` | Integration timestep |
| `temperature_K` | `300 K` | Target temperature |
| `friction_per_fs` | `0.01 /fs` | Langevin friction coefficient |
| `trajectory_interval` | `100` | Steps between trajectory writes |
| `log_interval` | `500` | Steps between CSV log entries |

### NPT Molecular Dynamics (`npt`)

| Parameter | Default | Description |
|---|---|---|
| `total_steps` | `50000` | Total MD steps |
| `timestep_fs` | `1 fs` | Integration timestep |
| `temperature_K` | `300 K` | Target temperature |
| `externalstress_bar` | `1 bar` | External pressure |
| `ttime_fs` | `100 fs` | Thermostat time constant |
| `ptime_fs` | `1000 fs` | Barostat time constant |
| `bulk_modulus_GPa` | `5 GPa` | Estimated bulk modulus (for Berendsen) |
| `trajectory_interval` | `100` | Steps between trajectory writes |
| `log_interval` | `500` | Steps between CSV log entries |

---

## Model

| Item | Value |
|---|---|
| Model | `uma-s-1p1` (UMA-S 1.1) |
| Task | `odac` (Open Direct Air Capture) |
| Backend | FairChem `FAIRChemCalculator` |
| Device | CUDA (GPU required) |
| Source | HuggingFace Hub (auto-downloaded on first run) |

---

## Dependencies

Key packages (see [environment.yml](environment.yml) for the full pinned list):

| Package | Role |
|---|---|
| `fairchem-core` | MLIP model and calculator |
| `ase` | Atoms object, optimizers, MD drivers |
| `pymatgen` | RMSD calculation via `StructureMatcher` |
| `torch` + CUDA 12 | GPU inference |
| `pandas` | CSV aggregation |
| `tqdm` | Progress bar |

---

## Helper Scripts

| Script | Description |
|---|---|
| [check_dist.py](check_dist.py) | Check interatomic distances in structures |
| [check_atom.py](check_atom.py) | Inspect atom types and counts |
| [check_cell.py](check_cell.py) | Inspect unit cell parameters |
| [paras.py](paras.py) | Elemental property tables (electronegativity, ionization energy) |
