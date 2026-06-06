from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import torch

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from model.GCN_solvent import SolventAtomClassifier
from model.cif_utils import build_crystal_graph, read_cif_structure, write_labeled_cif
from model.data_solvent import GaussianDistance


def choose_device(device_arg: str) -> torch.device:
    if device_arg == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(device_arg)


def load_checkpoint(path: Path, device: torch.device) -> dict:
    try:
        return torch.load(path, map_location=device, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=device)


def one_hot(atomic_number: int, max_atomic_number: int = 118) -> np.ndarray:
    fea = np.zeros((max_atomic_number + 1,), dtype=np.float32)
    if atomic_number < 1 or atomic_number > max_atomic_number:
        raise ValueError(f"Unsupported atomic number: {atomic_number}")
    fea[atomic_number] = 1.0
    return fea


def graph_to_tensors(graph: dict, pos: np.ndarray, radius: float, dmin: float, step: float):
    gdf = GaussianDistance(dmin=dmin, dmax=radius, step=step)
    atom_fea = np.vstack([one_hot(int(number)) for number in graph["numbers"]])
    index1 = np.array(graph["index1"], dtype=np.int64)
    index2 = np.array(graph["index2"], dtype=np.int64)
    dij = np.array(graph["dij"], dtype=np.float32)
    nbr_fea = gdf.expand(dij).astype(np.float32)

    atom_fea_tensor = torch.Tensor(np.concatenate([atom_fea, pos.astype(np.float32)], axis=1))
    return (
        atom_fea_tensor,
        torch.Tensor(nbr_fea),
        torch.LongTensor(index1),
        torch.LongTensor(index2),
        torch.LongTensor([0] * len(atom_fea)),
    )


def predict_solvent_atoms(
    cif_path: Path,
    checkpoint_path: Path,
    threshold: float | None = None,
    radius: float | None = None,
    max_neighbors: int = 12,
    device_name: str = "auto",
) -> dict:
    device = choose_device(device_name)
    checkpoint = load_checkpoint(checkpoint_path, device)
    model_config = checkpoint["model_config"]
    threshold_value = threshold if threshold is not None else float(checkpoint.get("threshold", 0.5))
    radius_value = radius if radius is not None else float(checkpoint.get("radius", 8.0))
    dmin = float(checkpoint.get("dmin", 0.0))
    step = float(checkpoint.get("step", 0.2))

    graph = build_crystal_graph(cif_path, radius=radius_value, max_neighbors=max_neighbors)
    structure = read_cif_structure(cif_path)
    pos = structure.frac_coords.astype(np.float32)
    tensors = tuple(tensor.to(device) for tensor in graph_to_tensors(graph, pos, radius_value, dmin, step))

    model = SolventAtomClassifier(**model_config).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    with torch.no_grad():
        probabilities = torch.sigmoid(model(*tensors)).cpu().numpy()
    labels = (probabilities >= threshold_value).astype(int)

    return {
        "structure": structure,
        "probabilities": probabilities,
        "labels": labels,
        "threshold": float(threshold_value),
        "radius": float(radius_value),
        "dmin": dmin,
        "step": step,
        "device": str(device),
    }


def write_prediction_csv(
    output_csv: Path,
    structure,
    probabilities: np.ndarray,
    labels: np.ndarray,
    include_kept: bool = False,
) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        header = ["atom_index", "atom_label", "element", "fract_x", "fract_y", "fract_z", "probability", "pred_label"]
        if include_kept:
            header.append("kept")
        writer.writerow(header)
        for i, atom in enumerate(structure.atoms):
            row = [
                i,
                atom.label,
                atom.symbol,
                f"{atom.frac[0]:.8f}",
                f"{atom.frac[1]:.8f}",
                f"{atom.frac[2]:.8f}",
                f"{float(probabilities[i]):.6f}",
                int(labels[i]),
            ]
            if include_kept:
                row.append(int(labels[i] == 0))
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict GTsR3 solvent atoms in a CIF file.")
    parser.add_argument("--cif", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--marked-cif", type=Path, default=None)
    parser.add_argument("--threshold", type=float, default=None)
    parser.add_argument("--radius", type=float, default=None)
    parser.add_argument("--max-neighbors", type=int, default=12)
    parser.add_argument("--device", default="auto")
    args = parser.parse_args()

    result = predict_solvent_atoms(
        cif_path=args.cif,
        checkpoint_path=args.checkpoint,
        threshold=args.threshold,
        radius=args.radius,
        max_neighbors=args.max_neighbors,
        device_name=args.device,
    )
    structure = result["structure"]
    probabilities = result["probabilities"]
    labels = result["labels"]
    write_prediction_csv(args.output_csv, structure, probabilities, labels)

    if args.marked_cif is not None:
        write_labeled_cif(
            args.cif,
            args.marked_cif,
            [
                ("_atom_site_solvent_probability", probabilities),
                ("_atom_site_solvent_pred_label", labels),
            ],
        )

    print(f"predicted atoms: {len(labels)}")
    print(f"predicted solvent atoms: {int(labels.sum())}")
    print(f"threshold: {result['threshold']}")
    print(f"csv: {args.output_csv}")
    if args.marked_cif is not None:
        print(f"marked cif: {args.marked_cif}")


if __name__ == "__main__":
    main()
