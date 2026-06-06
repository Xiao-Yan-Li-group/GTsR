from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from model.cif_utils import (  # noqa: E402
    label_removed_atom_lists,
    pymatgen_structure_to_cif_atoms,
    read_cif_structure,
    read_pymatgen_structure,
    subset_pymatgen_structure,
    write_structure_subset_cif,
)

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    def tqdm(items, **_kwargs):
        return items


SUMMARY_FIELDS = [
    "refcode",
    "task",
    "source_file",
    "target_file",
    "source_atoms",
    "target_atoms",
    "label_kept",
    "label_removed",
    "pos_count_ok",
    "json_count_ok",
    "passed",
    "match_digits_used",
    "reason",
]


def resolve_cif_dir(dataset_dir: Path) -> Path:
    nested = dataset_dir / "clean_dataset"
    return nested if nested.is_dir() else dataset_dir


def read_summary_rows(path: Path, task: str, limit: int | None) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if task != "both":
        rows = [row for row in rows if row["task"] == task]
    return rows[:limit] if limit is not None else rows


def count_json_atoms(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open(encoding="utf-8") as handle:
        return len(json.load(handle).get("numbers", []))


def check_kept_against_target(kept_structure, target_file: Path, match_digits: list[int]) -> dict:
    if kept_structure is None or len(kept_structure) == 0:
        return {"passed": False, "match_digits_used": "", "reason": "label removed all source atoms", "attempts": []}

    reference_atoms = pymatgen_structure_to_cif_atoms(kept_structure)
    target_atoms = read_cif_structure(target_file).atoms
    attempts = []
    for digits in match_digits:
        result = label_removed_atom_lists(reference_atoms, target_atoms, match_digits=digits)
        passed = result.kept_count == result.target_count and result.removed_count == 0
        attempt = {
            "match_digits": int(digits),
            "kept_atoms": int(result.kept_count),
            "target_atoms": int(result.target_count),
            "extra_atoms_after_label": int(result.removed_count),
            "passed": bool(passed),
        }
        attempts.append(attempt)
        if passed:
            return {"passed": True, "match_digits_used": int(digits), "reason": "", "attempts": attempts}

    best = max(attempts, key=lambda item: (item["kept_atoms"], -item["extra_atoms_after_label"]))
    return {
        "passed": False,
        "match_digits_used": best["match_digits"],
        "reason": (
            f"kept_from_label does not match target; best kept {best['kept_atoms']}/"
            f"{best['target_atoms']}, extra_atoms_after_label={best['extra_atoms_after_label']}"
        ),
        "attempts": attempts,
    }


def verify_row(row: dict, cif_dir: Path, data_dir: Path, match_digits: list[int]) -> dict:
    refcode = row["refcode"]
    task = row["task"]
    source_file = cif_dir / row["source_file"]
    target_file = cif_dir / row["target_file"] if row.get("target_file") else None
    label_file = data_dir / task / "npy" / "label" / f"{refcode}.npy"
    pos_file = data_dir / task / "npy" / "pos" / f"{refcode}.npy"
    json_file = data_dir / task / "json" / f"{refcode}.json"

    labels = np.load(label_file).astype(np.float32)
    source = read_cif_structure(source_file)
    if len(labels) != len(source.atoms):
        raise ValueError(f"Label length mismatch: {len(labels)} labels, {len(source.atoms)} source atoms")

    structure = read_pymatgen_structure(source_file)
    solvent_mask = labels >= 0.5
    kept_structure = subset_pymatgen_structure(structure, ~solvent_mask)
    label_removed = int(solvent_mask.sum())
    label_kept = int((~solvent_mask).sum())

    if target_file is None:
        check = {
            "passed": label_removed == 0,
            "match_digits_used": "",
            "reason": "" if label_removed == 0 else f"no target file but label removes {label_removed} atoms",
            "attempts": [],
        }
        target_atoms = ""
    else:
        check = check_kept_against_target(kept_structure, target_file, match_digits)
        target_atoms = len(read_cif_structure(target_file).atoms)

    pos_atoms = np.load(pos_file).shape[0] if pos_file.exists() else None
    json_atoms = count_json_atoms(json_file)
    return {
        "refcode": refcode,
        "task": task,
        "source_file": source_file.name,
        "target_file": target_file.name if target_file is not None else "",
        "source_atoms": len(source.atoms),
        "target_atoms": target_atoms,
        "label_kept": label_kept,
        "label_removed": label_removed,
        "pos_count_ok": "" if pos_atoms is None else int(pos_atoms == len(labels)),
        "json_count_ok": "" if json_atoms is None else int(json_atoms == len(labels)),
        "passed": int(check["passed"]),
        "match_digits_used": check["match_digits_used"],
        "reason": check["reason"],
        "_attempts": check["attempts"],
        "_source_path": source_file,
        "_target_path": target_file,
        "_label_path": label_file,
    }


def write_failure(output_dir: Path, record: dict) -> None:
    case_dir = output_dir / "failed" / record["task"] / record["refcode"]
    case_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(record["_source_path"], case_dir / "source.cif")
    if record["_target_path"] is not None:
        shutil.copy2(record["_target_path"], case_dir / "target.cif")

    labels = np.load(record["_label_path"]).astype(np.float32)
    label_column = f"_atom_site_{record['task']}_ver_label"
    write_structure_subset_cif(
        record["_source_path"],
        case_dir / "kept_from_label.cif",
        labels < 0.5,
        [(label_column, labels.astype(int))],
    )
    write_structure_subset_cif(
        record["_source_path"],
        case_dir / "solvent_from_label.cif",
        labels >= 0.5,
        [(label_column, labels.astype(int))],
    )
    (case_dir / "reason.json").write_text(
        json.dumps(
            {
                "refcode": record["refcode"],
                "task": record["task"],
                "reason": record["reason"],
                "attempts": record["_attempts"],
                "label_kept": record["label_kept"],
                "label_removed": record["label_removed"],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def write_summary(path: Path, records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for record in records:
            writer.writerow({field: record.get(field, "") for field in SUMMARY_FIELDS})


def error_record(row: dict, error: Exception) -> dict:
    record = {field: row.get(field, "") for field in SUMMARY_FIELDS}
    record.update({"passed": 0, "reason": str(error)})
    return record


def process_job(job: dict) -> dict:
    try:
        return verify_row(
            row=job["row"],
            cif_dir=job["cif_dir"],
            data_dir=job["data_dir"],
            match_digits=job["match_digits"],
        )
    except Exception as exc:
        return error_record(job["row"], exc)


def collect_records(rows: list[dict], cif_dir: Path, data_dir: Path, match_digits: list[int], num_workers: int) -> list[dict]:
    jobs = [
        {
            "row": row,
            "cif_dir": cif_dir,
            "data_dir": data_dir,
            "match_digits": match_digits,
        }
        for row in rows
    ]
    records: list[dict] = []

    if num_workers <= 1:
        for job in tqdm(jobs, total=len(jobs), desc="Verifying labels", unit="sample", dynamic_ncols=True):
            records.append(process_job(job))
        return records

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        future_to_job = {executor.submit(process_job, job): job for job in jobs}
        for future in tqdm(
            as_completed(future_to_job),
            total=len(future_to_job),
            desc=f"Verifying labels ({num_workers} workers)",
            unit="sample",
            dynamic_ncols=True,
        ):
            records.append(future.result())
    records.sort(key=lambda record: (record.get("task", ""), record.get("refcode", "")))
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify prepared GTsR labels against target CIFs.")
    parser.add_argument("--dataset-dir", type=Path, default=PROJECT_ROOT / "clean_dataset_test")
    parser.add_argument("--data-dir", type=Path, default=SCRIPT_DIR / "data")
    parser.add_argument("--summary-csv", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "ver_label")
    parser.add_argument("--task", choices=("free", "coordinated", "both"), default="both")
    parser.add_argument("--match-digits", type=int, nargs="+", default=[3, 2])
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument(
        "--num-workers",
        type=int,
        default=min(8, max(1, os.cpu_count() or 1)),
        help="Parallel worker processes. Use 1 for serial mode.",
    )
    args = parser.parse_args()

    summary_csv = args.summary_csv or args.data_dir / "label_summary.csv"
    rows = read_summary_rows(summary_csv, args.task, args.limit)
    cif_dir = resolve_cif_dir(args.dataset_dir)

    records = collect_records(rows, cif_dir, args.data_dir, args.match_digits, max(1, args.num_workers))
    for record in tqdm(records, total=len(records), desc="Writing failed cases", unit="sample", dynamic_ncols=True):
        if not int(record["passed"]) and "_source_path" in record:
            write_failure(args.output_dir, record)

    write_summary(args.output_dir / "summary.csv", records)
    failed = sum(1 for record in records if not int(record["passed"]))
    print(f"checked samples: {len(records)}")
    print(f"failed samples: {failed}")
    print(f"summary: {args.output_dir / 'summary.csv'}")
    print(f"failed output: {args.output_dir / 'failed'}")


if __name__ == "__main__":
    main()
