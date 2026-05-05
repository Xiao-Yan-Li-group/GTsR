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

from model.cif_utils import read_cif_structure


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


def data_block_name(path: Path) -> str:
    return "".join(char if char.isalnum() or char == "_" else "_" for char in path.stem) or "structure"


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
    output_lines = [
        f"data_{data_block_name(template_cif)}",
        f"_cell_length_a    {structure.cell['_cell_length_a']:.8f}",
        f"_cell_length_b    {structure.cell['_cell_length_b']:.8f}",
        f"_cell_length_c    {structure.cell['_cell_length_c']:.8f}",
        f"_cell_angle_alpha {structure.cell['_cell_angle_alpha']:.8f}",
        f"_cell_angle_beta  {structure.cell['_cell_angle_beta']:.8f}",
        f"_cell_angle_gamma {structure.cell['_cell_angle_gamma']:.8f}",
        "_symmetry_space_group_name_H-M 'P 1'",
        "_symmetry_Int_Tables_number 1",
        "loop_",
        " _atom_site_label",
        " _atom_site_type_symbol",
        " _atom_site_fract_x",
        " _atom_site_fract_y",
        " _atom_site_fract_z",
        f" {label_column}",
    ]

    kept_atoms = 0
    for atom_index, atom in enumerate(structure.atoms):
        label = int(labels[atom_index] >= 0.5)
        if label == target_value:
            output_lines.append(
                f" {atom.label} {atom.symbol} "
                f"{atom.frac[0]:.8f} {atom.frac[1]:.8f} {atom.frac[2]:.8f} {label}"
            )
            kept_atoms += 1

    output_cif.parent.mkdir(parents=True, exist_ok=True)
    output_cif.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
    return kept_atoms


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
