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

## Installation

```bash
git clone https://github.com/coollkr/GTsR.git
cd GTsR
pip install -e .
```

## Usage

### Solvent Removal

```python
from gtsr import GTsRunner

runner = GTsRunner(checkpoint="free") ### for free solvent removal
runner = GTsRunner(checkpoint="all") ### for all solvent removal
runner = GTsRunner(checkpoint="path/to/ckpt.pth", device="cpu") #### use your model
result = runner.clean(
    cif="input.cif",
    output="prediction",
    threshold=0.5,
)
```

### Predict Activation Stability

```python
from gtsr import GTsRunner

runner = GTsRunner(checkpoint="stability")
score = runner.stability(cif="cleaned_framework.cif")

if score == 1:
    print("The cleaned structure is stable.")
else:
    print("The cleaned structure is not stable.")
```

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

[Host on Streamlit](https://xiao-yan-li-group.streamlit.app/GTsR)
or in your location
```bash
streamlit run webapp/Home.py
```

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
