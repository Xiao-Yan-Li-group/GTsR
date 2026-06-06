from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from model.cif_utils import read_cif_structure, write_structure_subset_cif


def resolve_cif_dir(dataset_dir: Path) -> Path:
    nested = dataset_dir / "clean_dataset"
    return nested if nested.is_dir() else dataset_dir


def source_cif_from_json(cif_dir: Path, graph_json: Path, cif_id: str) -> Path:
    if graph_json.exists():
        with graph_json.open(encoding="utf-8") as handle:
            graph = json.load(handle)
        source_name = graph.get("metadata", {}).get("source_file")
        if source_name:
            source_file = cif_dir / source_name
            if source_file.exists():
                return source_file
    fallback = cif_dir / f"{cif_id}_R.cif"
    if fallback.exists():
        return fallback
    raise FileNotFoundError(f"Cannot find source CIF for {cif_id}")


def read_ids(label_dir: Path, requested_ids: list[str] | None, limit: int | None) -> list[str]:
    if requested_ids:
        ids = requested_ids
    else:
        ids = sorted(path.stem for path in label_dir.glob("*.npy"))
    if limit is not None:
        ids = ids[:limit]
    return ids


def write_split_cif(
    template_cif: Path,
    output_cif: Path,
    labels: np.ndarray,
    keep_solvent: bool,
    label_column: str,
) -> int:
    structure = read_cif_structure(template_cif)
    if labels.shape[0] != len(structure.atoms):
        raise ValueError(
            f"Label length mismatch for {template_cif}: {labels.shape[0]} labels, "
            f"{len(structure.atoms)} atoms"
        )

    target_value = 1 if keep_solvent else 0
    keep_mask = np.array([int(value >= 0.5) == target_value for value in labels], dtype=bool)
    return write_structure_subset_cif(
        template_cif,
        output_cif,
        keep_mask,
        [(label_column, labels.astype(int))],
    )


def split_one(
    cif_id: str,
    task: str,
    dataset_dir: Path,
    data_dir: Path,
    output_dir: Path,
) -> tuple[int, int]:
    cif_dir = resolve_cif_dir(dataset_dir)
    task_data_dir = data_dir / task
    source_file = source_cif_from_json(cif_dir, task_data_dir / "json" / f"{cif_id}.json", cif_id)
    labels = np.load(task_data_dir / "npy" / "label" / f"{cif_id}.npy").astype(np.float32)

    task_output = output_dir / task
    label_column = f"_atom_site_{task}_solvent_label"
    solvent_count = write_split_cif(
        source_file,
        task_output / f"{cif_id}_{task}_solvent.cif",
        labels,
        keep_solvent=True,
        label_column=label_column,
    )
    mof_count = write_split_cif(
        source_file,
        task_output / f"{cif_id}_{task}_mof.cif",
        labels,
        keep_solvent=False,
        label_column=label_column,
    )
    return solvent_count, mof_count


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Split GTsR3 source CIFs into solvent and MOF CIFs with prepared atom labels."
    )
    parser.add_argument("--dataset-dir", type=Path, default=PROJECT_ROOT / "clean_dataset_test")
    parser.add_argument("--data-dir", type=Path, default=SCRIPT_DIR / "data")
    parser.add_argument("--output-dir", type=Path, default=SCRIPT_DIR / "label_check")
    parser.add_argument("--task", choices=("free", "coordinated", "both"), default="both")
    parser.add_argument("--cif-id", action="append", default=None, help="CIF id to split. Repeat for many ids.")
    parser.add_argument("--limit", type=int, default=None, help="Process only the first N ids.")
    args = parser.parse_args()

    tasks = ("free", "coordinated") if args.task == "both" else (args.task,)
    total = 0
    errors: list[tuple[str, str, str]] = []
    for task in tasks:
        label_dir = args.data_dir / task / "npy" / "label"
        ids = read_ids(label_dir, args.cif_id, args.limit)
        for cif_id in ids:
            try:
                solvent_count, mof_count = split_one(
                    cif_id=cif_id,
                    task=task,
                    dataset_dir=args.dataset_dir,
                    data_dir=args.data_dir,
                    output_dir=args.output_dir,
                )
                total += 1
                print(f"{task} {cif_id}: solvent={solvent_count}, mof={mof_count}")
            except Exception as exc:
                errors.append((task, cif_id, str(exc)))

    if errors:
        print("errors:")
        for task, cif_id, error in errors:
            print(f"{task} {cif_id}: {error}")
    print(f"written pairs: {total}")
    print(f"output: {args.output_dir}")


if __name__ == "__main__":
    main()
