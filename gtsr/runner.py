from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import torch

try:
    from .src.GCN import SolventAtomClassifier
    from .src.cif_utils import cif2graph, cif2pos, label2cif, get_sol_smi
    from .src.data import GaussianDistance
    from .src.utils import load_checkpoint
except ImportError:
    from src.GCN import SolventAtomClassifier
    from src.cif_utils import cif2graph, cif2pos, label2cif, get_sol_smi
    from src.data import GaussianDistance
    from src.utils import load_checkpoint


PACKAGE_DIR = Path(__file__).resolve().parent


def _bundled_checkpoint(name: str) -> Path:
    filename = f"{name}_best.pth"
    candidates = (
        PACKAGE_DIR / "ckpt" / filename,
        PACKAGE_DIR.parent / "ckpt" / filename,
    )
    return next((path for path in candidates if path.is_file()), candidates[0])


CHECKPOINTS = {
    "free": _bundled_checkpoint("free"),
    "all": _bundled_checkpoint("all"),
}
DEFAULT_CHECKPOINT = CHECKPOINTS["free"]


class GTsRunner:

    def __init__(
        self,
        checkpoint: str | Path = "",
        device: str | torch.device | None = None,
    ) -> None:
        self.device = self._resolve_device(device)
        self.checkpoint_path = self._resolve_checkpoint(checkpoint)
        self.checkpoint = load_checkpoint(self.checkpoint_path, device=self.device)

        model_config = self.checkpoint.get("model_config")
        if not isinstance(model_config, dict):
            raise ValueError(
                f"Checkpoint does not contain a valid model_config: {self.checkpoint_path}"
            )

        self.model = SolventAtomClassifier(**model_config).to(self.device)
        self.model.load_state_dict(self.checkpoint["state_dict"])
        self.model.eval()

        self.radius = float(self.checkpoint.get("radius", 8.0))
        self.dmin = float(self.checkpoint.get("dmin", 0.0))
        self.step = float(self.checkpoint.get("step", 0.2))
        self.default_threshold = float(self.checkpoint.get("threshold", 0.5))
        self.task = str(self.checkpoint.get("task", "unknown"))
        self.max_atomic_number = 118
        self.gdf = GaussianDistance(
            dmin=self.dmin,
            dmax=self.radius,
            step=self.step,
        )

    def predict(
        self,
        cif: str | Path = "",
        output: str | Path = "",
        threshold: float | None = None,
    ) -> dict[str, Any]:

        cif_path = self._resolve_cif(cif)
        output_dir = self._resolve_output(cif_path, output)
        cutoff = self.default_threshold if threshold is None else float(threshold)
        if not 0.0 <= cutoff <= 1.0:
            raise ValueError(f"threshold must be between 0 and 1, got {cutoff}")

        tensors = self._build_tensors(cif_path)
        with torch.inference_mode():
            probabilities = torch.sigmoid(self.model(*tensors)).cpu().numpy()

        labels = (probabilities >= cutoff).astype(np.int64)
        label2cif(cif_path, labels, str(output_dir))

        stem = cif_path.stem
        framework_path = output_dir / f"{stem}_framework.cif"
        solvent_path = output_dir / f"{stem}_solvent.cif"

        sol_smis = get_sol_smi(solvent_path)

        return {
            "input": str(cif_path),
            "output": str(output_dir),
            "framework": str(framework_path),
            "solvent": str(solvent_path) if solvent_path.exists() else None,
            "checkpoint": str(self.checkpoint_path),
            "task": self.task,
            "threshold": cutoff,
            "num_atoms": int(labels.size),
            "num_framework_atoms": int((labels == 0).sum()),
            "num_solvent_atoms": int((labels == 1).sum()),
            "probabilities": probabilities.tolist(),
            "labels": labels.tolist(),
            "solvent_smiles": sol_smis
        }

    def _build_tensors(self, cif_path: Path) -> tuple[torch.Tensor, ...]:
        graph = cif2graph(cif_path, radius=self.radius)
        positions = np.asarray(cif2pos(cif_path), dtype=np.float32)
        numbers = np.asarray(graph["numbers"], dtype=np.int64)

        if numbers.size == 0:
            raise ValueError(f"CIF contains no atoms: {cif_path}")
        if numbers.min() < 1 or numbers.max() > self.max_atomic_number:
            raise ValueError(
                f"CIF contains an unsupported atomic number; supported range is "
                f"1-{self.max_atomic_number}"
            )
        if len(positions) != len(numbers):
            raise ValueError(
                f"Position/atom mismatch in {cif_path}: "
                f"{len(positions)} positions for {len(numbers)} atoms"
            )

        atom_features = np.eye(
            self.max_atomic_number + 1,
            dtype=np.float32,
        )[numbers]
        atom_features = np.concatenate([atom_features, positions], axis=1)

        distances = np.asarray(graph["dij"], dtype=np.float32)
        neighbor_features = self.gdf.expand(distances).astype(np.float32)
        index1 = np.asarray(graph["index1"], dtype=np.int64)
        index2 = np.asarray(graph["index2"], dtype=np.int64)
        atom_index = np.zeros(len(numbers), dtype=np.int64)

        tensors = (
            torch.from_numpy(atom_features),
            torch.from_numpy(neighbor_features),
            torch.from_numpy(index1),
            torch.from_numpy(index2),
            torch.from_numpy(atom_index),
        )
        return tuple(tensor.to(self.device) for tensor in tensors)

    @staticmethod
    def _resolve_device(device: str | torch.device | None) -> torch.device:
        if device is None or str(device).lower() == "auto":
            return torch.device("cuda" if torch.cuda.is_available() else "cpu")
        resolved = torch.device(device)
        if resolved.type == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("CUDA was requested but is not available")
        return resolved

    @staticmethod
    def _resolve_checkpoint(checkpoint: str | Path) -> Path:
        checkpoint_name = str(checkpoint).strip().lower()
        if not checkpoint_name:
            path = DEFAULT_CHECKPOINT
        elif checkpoint_name in CHECKPOINTS:
            path = CHECKPOINTS[checkpoint_name]
        else:
            path = Path(checkpoint).expanduser()
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Checkpoint not found: {path}")
        return path

    @staticmethod
    def _resolve_cif(cif: str | Path) -> Path:
        if not cif:
            raise ValueError("cif must be a path to an input CIF file")
        path = Path(cif).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"CIF not found: {path}")
        return path

    @staticmethod
    def _resolve_output(cif_path: Path, output: str | Path) -> Path:
        path = (
            Path(output).expanduser()
            if output
            else cif_path.parent / f"{cif_path.stem}_gtsr"
        )
        path = path.resolve()
        path.mkdir(parents=True, exist_ok=True)
        return path
