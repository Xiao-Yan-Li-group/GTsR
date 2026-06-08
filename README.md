# GTsR3

GTsR3 is the ASE/pymatgen-reading variant of GTsR2. It prepares two atom-level solvent-removal datasets from
`clean_dataset_test` using `dataset_labels.csv`.

## CIF Reading

`GTsR3/model/cif_utils.py` reads CIF structures with this order:

1. `ase.io.read(...)`
2. `pymatgen.io.ase.AseAtomsAdaptor.get_structure(...)`
3. fallback to `pymatgen.core.Structure.from_file(...)`
4. fallback to `pymatgen.io.cif.CifParser(..., occupancy_tolerance=10)`

The same read structure is used for labels, fractional positions, atomic numbers, and pymatgen neighbor lists.

Install the reader dependencies before running:

```bash
pip install -r GTsR3/requirements.txt
```

## Tasks

- `free`: if `r2f=1`, labels atoms removed from `R -> F`; if `r2f=0`, uses `R` with all-zero labels.
- `coordinated`: uses `F` when available, otherwise `R`; if `f2a=1`, labels atoms removed from source `-> A`; if `f2a=0`, uses all-zero labels.

Both tasks build crystal graphs with cutoff `8.0` and `12` nearest neighbors.

## Prepare

```bash
python GTsR3/prepare_dataset.py --num-workers 8
```

Outputs:

- `GTsR3/data/free/...`
- `GTsR3/data/coordinated/...`
- `GTsR3/data/splits/free_*.csv`
- `GTsR3/data/splits/coordinated_*.csv`
- `GTsR3/data/label_summary.csv`
- `GTsR3/data/prepare_errors.csv` if any samples fail

Debug run:

```bash
python GTsR3/prepare_dataset.py --limit 20 --write-marked-cif
```

By default, compared structures are kept when at least 50% of target atoms align. Samples below that threshold are
written to `prepare_errors.csv`. Adjust with `--min-alignment-ratio`.

## Train

```bash
python GTsR3/train_solvent.py --task free --epochs 100 --batch-size 16
python GTsR3/train_solvent.py --task coordinated --epochs 100 --batch-size 16
```

Best checkpoints are written under `GTsR3/pth/<task>/best.pth`.

## Label Check

Split prepared labels into solvent and MOF CIFs:

```bash
python GTsR3/label_check.py --task both --cif-id 1101175
```

Outputs go to `GTsR3/label_check/<task>/`.
