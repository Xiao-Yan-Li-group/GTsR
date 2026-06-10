from pathlib import Path

from setuptools import setup


ROOT = Path(__file__).resolve().parent


setup(
    name="gtsr-mof",
    version="0.0.1",
    description="Graph neural network tool for solvent removal from MOF structures",
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    author="Xiao-Yan Li Group",
    license="CC-BY-4.0",
    python_requires=">=3.9",
    packages=["gtsr", "gtsr.src", "gtsr.ckpt"],
    py_modules=["GTsRunner"],
    package_dir={
        "gtsr": "gtsr",
        "gtsr.src": "src",
        "gtsr.ckpt": "ckpt",
    },
    package_data={
        "gtsr.ckpt": [
            "free_best.pth",
            "all_best.pth",
            "stability_best.pkl",
        ],
    },
    include_package_data=True,
    install_requires=[
        "ase>=3.19",
        "numpy>=1.21",
        "pymatgen>=2018.6.11",
        "scikit-learn>=1.0",
        "torch>=1.12",
        "molSimplify==1.8.0"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
