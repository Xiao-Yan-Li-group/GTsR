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


