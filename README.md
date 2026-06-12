# GTsR

<div align="center">
        <img src="https://raw.githubusercontent.com/Xiao-Yan-Li-group/GTsR/main/webapp/imgs/gtsr_logo.png" alt="GTsR logo" width="500"/>
</div> 

**GTsR (GNN Tool for Solvent Removal)** is a tool for solvent identification, solvent removal, and activation-stability prediction in metal-organic frameworks (MOFs).

GTsR uses graph neural networks to classify atoms in CIF structures and generate solvent-free framework CIF files. It also provides a random forest model that predicts the activation stability of cleaned MOFs using structural, pore, and RAC descriptors.

## Bundled Models

| `checkpoint` | Model file | Purpose |
| --- | --- | --- |
| `free` (default) | `ckpt/free_best.pth` | Remove free solvent |
| `all` | `ckpt/all_best.pth` | Remove all solvent |
| `stability` | `ckpt/stability_best.pkl` | Predict activation stability |

The `free` and `all` checkpoints are atom-level GNN classifiers. The `stability` checkpoint is a random forest model bundled with its missing-value imputer.

## Requirements

- Python 3.9 or later
- PyTorch
- ASE
- pymatgen
- NumPy
- scikit-learn
- molSimplify 1.8.0
- RDKit
- NetworkX
- Zeo++, required only for stability prediction
- Streamlit and stmol, required only for the web interface

Stability prediction calls the Zeo++ `network` executable. GTsR searches for it in the system `PATH` and next to the active Python interpreter.

## Installation

Creating an isolated Conda environment is recommended. Install scientific dependencies such as Zeo++ and RDKit, then install GTsR from source:

```bash
git clone https://github.com/coollkr/GTsR.git
cd GTsR
pip install -e .
pip install -r requirements.txt
```

Confirm that the Zeo++ executable required for stability prediction is available:

```bash
network
```

## Python API

### Remove Free Solvent

```python
from gtsr import GTsRunner

runner = GTsRunner(checkpoint="free")
result = runner.clean(
    cif="input.cif",
    output="prediction",
    threshold=0.5,
)

print("Framework:", result["framework"])
print("Solvent:", result["solvent"])
print("Solvent SMILES:", result["solvent_smiles"])
```

Omitting `checkpoint` or passing an empty string also selects the `free` model:

```python
runner = GTsRunner()
```

### Remove All Solvent

```python
from gtsr import GTsRunner

runner = GTsRunner(checkpoint="all")
result = runner.clean(
    cif="input.cif",
    output="prediction",
    threshold=0.5,
)
```

`threshold` must be between `0.0` and `1.0`. If `output` is omitted, results are written to a `<filename>_gtsr` directory beside the input CIF.

### Use a Custom GNN Checkpoint

```python
from gtsr import GTsRunner

runner = GTsRunner(checkpoint="path/to/custom_best.pth", device="cpu")
result = runner.clean(cif="input.cif")
```

`device` accepts `cpu`, `cuda`, `auto`, or a PyTorch `torch.device`. When omitted, GTsR automatically uses a GPU if CUDA is available.

### Predict Activation Stability

The stability model uses 191 input features derived from the atom count, unit cell, Zeo++ pore descriptors, and molSimplify RAC descriptors.

```python
from gtsr import GTsRunner

runner = GTsRunner(checkpoint="stability")
score = runner.stability(cif="cleaned_framework.cif")

if score == 1:
    print("The cleaned structure is stable.")
else:
    print("The cleaned structure is not stable.")
```

Stability prediction returns a single class:

- `1`: predicted to be stable.
- `0`: predicted to be unstable.

Stability prediction should generally be performed on the framework CIF generated after solvent removal.

## `clean()` Result

`clean()` returns a dictionary containing the following fields:

| Field | Description |
| --- | --- |
| `input` | Absolute path to the input CIF |
| `output` | Output directory |
| `framework` | Path to the cleaned framework CIF |
| `solvent` | Path to the solvent CIF, or `None` if no file was generated |
| `checkpoint` | Path to the checkpoint used for prediction |
| `task` | Task name stored in the checkpoint |
| `threshold` | Atom-classification threshold |
| `num_atoms` | Total number of atoms |
| `num_framework_atoms` | Number of framework atoms |
| `num_solvent_atoms` | Number of solvent atoms |
| `probabilities` | Solvent probability for each atom |
| `labels` | Predicted class label for each atom |
| `solvent_smiles` | SMILES strings of identified solvents |

## Web Interface

Start the Streamlit web interface with:

```bash
streamlit run webapp/Home.py
```

The web interface supports:

- Uploading and visualizing CIF files.
- Selecting the `Free Only` or `All` solvent model.
- Adjusting the classification threshold.
- Viewing framework and solvent structures.
- Downloading cleaned MOFs.
- Predicting the activation stability of cleaned structures.

## Project Structure

```text
GTsR/
├── ckpt/                    # Bundled GNN and random forest models
├── gtsr/
│   ├── __init__.py          # Public Python API
│   └── runner.py            # GTsRunner inference interface
├── src/
│   ├── GCN.py               # GNN model
│   ├── cif_utils.py         # CIF, pore, and RAC descriptor utilities
│   ├── data.py              # Graph data and distance-expansion utilities
│   └── utils.py             # Checkpoint utilities
├── webapp/                  # Streamlit web interface
├── experiments/             # Training and experiment scripts
├── evaluation/              # Model evaluation and data processing
├── requirements.txt
└── setup.py
```

## Build the Package

```bash
python setup.py sdist bdist_wheel
```

The generated package includes:

- `free_best.pth`
- `all_best.pth`
- `stability_best.pkl`

## Troubleshooting

### `Zeo++ executable 'network' was not found`

Stability prediction requires Zeo++. Install Zeo++ and ensure that `network` is available in the system `PATH` or in the active Python environment's `bin` directory.

### `CUDA was requested but is not available`

CUDA is unavailable in the current environment. Initialize the GNN runner with the CPU device:

```python
runner = GTsRunner(checkpoint="free", device="cpu")
```

### `Checkpoint not found`

Confirm that the custom checkpoint path exists. To use a bundled model, pass `free`, `all`, or `stability`.

## Citation

Update the following entry when the associated publication becomes available:

```bibtex
@article{gtsr-xyl-group,
  title   = {GTSR: A GNN Based Tool for Solvent Removal from MOF with Stability Check},
  author  = {Liang, Kairui and Zhao, Guobin and Li, Xiao-Yan},
  year    = {2026}
}
```

## License

The repository's [`LICENSE`](LICENSE) file currently uses the MIT License.
