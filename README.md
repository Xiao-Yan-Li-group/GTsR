# GTsR3

## Python API

Install the package and use the bundled free-solvent checkpoint:

```python
from gtsr import GTsRunner

runner = GTsRunner(checkpoint="")
result = runner.predict(
    cif="input.cif",
    output="prediction",
    threshold=0.5,
)

print(result["framework"])
print(result["solvent"])
```

Use the coordinated-solvent checkpoint when needed:

```python
free_runner = GTsRunner(checkpoint="free")
all_runner = GTsRunner(checkpoint="all")
```

You can still pass a custom checkpoint path:

```python
runner = GTsRunner(checkpoint="path/to/custom_best.pth")
```

Build the package locally:

```bash
python setup.py sdist bdist_wheel
```
