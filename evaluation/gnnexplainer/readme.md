# GNNExplainer for GTsR

`explain_all.py` runs [GNNExplainer](https://arxiv.org/abs/1903.03894) on **every atom** of a MOF and writes one CSV row per atom. Each row shows how strongly the atom is predicted as solvent and how much the GTsR model relies on that atom to make its predictions.

For each target atom, a soft mask on the graph edges is optimised so that the model keeps its own prediction for that atom while using as few edges as possible. Several atoms are explained in one batch (one copy of the crystal per atom), so a GPU is recommended.

---

## Quick Start

```bash
conda activate gtsr                      # environment.yml in the repository root
cd evaluation/gnnexplainer

python explain_all.py --cif ADABEO.cif --model free
```

The results are written to `output/free/ADABEO/atoms.csv`.

---

## Usage

```bash
# one structure, free-solvent model (ckpt/free_best.pth)
python explain_all.py --cif ADABEO.cif --model free

# one structure, all-solvent model (ckpt/all_best.pth)
python explain_all.py --cif ADABEO.cif --model all

# several structures in one run (the model is loaded once)
python explain_all.py --cif ADABEO.cif MOF2.cif MOF3.cif --model free

# every CIF in a folder
python explain_all.py --cif path/to/cifs/*.cif --model all

# your own checkpoint instead of ckpt/<model>_best.pth
python explain_all.py --cif ADABEO.cif --model free --checkpoint path/to/best.pth

# write the results somewhere else (default: ./output)
python explain_all.py --cif ADABEO.cif --model free --out results/gnnexplainer

# choose the device (default: GPU if available, otherwise CPU)
python explain_all.py --cif ADABEO.cif --model free --device cpu
python explain_all.py --cif ADABEO.cif --model free --device cuda:1

# quicker, rougher explanations / smaller batches for a small GPU
python explain_all.py --cif ADABEO.cif --model free --epochs 100 --batch 8
```

Both models can be run on the same structure. Their results go to separate folders (`output/free/` and `output/all/`).

Structures that already have an `atoms.csv` are skipped, so an interrupted run can simply be started again. Delete a structure's folder to recompute it.

---

## Options

| Option | Default | Description |
|---|---|---|
| `--cif` | required | One or more CIF files |
| `--model` | required | `free` (remove free solvent) or `all` (remove all solvent); also the output subfolder |
| `--checkpoint` | `ckpt/<model>_best.pth` | Model file to explain |
| `--out` | `./output` | Output root; results go to `<out>/<model>/<cif name>/` |
| `--batch` | `32` | Atoms explained together; GPU memory grows with it |
| `--epochs` | `300` | Optimisation steps per batch |
| `--lr` | `0.01` | Learning rate of the edge masks |
| `--size-coef` | `0.005` | Penalty on the total mask size (sparser explanations when larger) |
| `--ent-coef` | `1.0` | Penalty on mask entropy (pushes masks towards 0 or 1) |
| `--device` | `auto` | `auto`, `cpu`, `cuda`, `cuda:1`, ... |

---

## Output

```
output/<model>/<cif name>/
├── atoms.csv           one row per atom
├── masks.npz           edge masks (atom × edge) and the graph
└── <cif name>.extxyz   structure with the atoms.csv columns as per-atom properties (e.g. for OVITO)
```

| Column | Description |
|---|---|
| `atom` | Atom index (0-based, CIF order) |
| `element` | Element symbol |
| `probability` | Solvent probability from the model |
| `pred` | 1 = solvent, 0 = framework (threshold stored in the checkpoint) |
| `p_no_edges` | Probability when all edges are removed |
| `edge_dependence` | `probability - p_no_edges`: how much the prediction depends on the neighbours |
| `p_soft_mask` | Probability with the learned soft mask |
| `p_kept_edges` | Probability with only the edges whose mask is > 0.5 |
| `n_edges_kept` | Number of edges with mask > 0.5 |
| `explanation_complete` | 1 if the kept edges alone reproduce the prediction |
| `influence` | How much this atom is used to explain **all** atoms |
| `influence_solvent` | How much this atom is used to explain the **predicted-solvent** atoms |

The influence of atom *j* is the sum, over all explained atoms *t* and over all edges *e* sent by *j*, of max(0, mask<sub>t</sub>(e) − median(mask<sub>t</sub>)). It has no unit and can only be compared between atoms of the same structure.

---

## Runtime

A cell with about 300 atoms takes about 1 min on one GPU with the default settings. On a single CPU core the same cell takes several hours.
