from __future__ import annotations

import argparse
import csv
import json
import random
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from model.cif_utils import build_crystal_graph, label_removed_atoms, read_cif_structure, write_labeled_cif

try:
    from tqdm import tqdm
except Exception:  # pragma: no cover
    def tqdm(items, **_kwargs):
        return items


TASKS = ("free", "coordinated")


def resolve_cif_dir(dataset_dir: Path) -> Path:
    nested = dataset_dir / "clean_dataset"
    return nested if nested.is_dir() else dataset_dir


def pick_cif(cif_dir: Path, refcode: str, suffixes: tuple[str, ...], required: bool = True) -> Path | None:
    matched = [cif_dir / f"{refcode}{suffix}" for suffix in suffixes if (cif_dir / f"{refcode}{suffix}").exists()]
    if len(matched) == 1:
        return matched[0]
    if not matched and not required:
        return None
    raise FileNotFoundError(f"Cannot uniquely pick {suffixes} for {refcode} in {cif_dir}")


def read_label_rows(labels_csv: Path) -> list[dict]:
    with labels_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        if "refcode" not in row:
            raise ValueError(f"dataset label row is missing refcode: {row}")
        row["refcode"] = row["refcode"].strip()
        row["r2f"] = str(row.get("r2f", "0")).strip()
        row["f2a"] = str(row.get("f2a", "0")).strip()
    return [row for row in rows if row["refcode"]]


def flag(row: dict, name: str) -> int:
    value = str(row.get(name, "0")).strip()
    return 1 if value in {"1", "true", "True", "yes", "Y"} else 0


def write_id_csv(path: Path, ids: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        for cif_id in ids:
            writer.writerow([cif_id])


def format_transform(axis_order: tuple[int, int, int], axis_signs: tuple[int, int, int]) -> str:
    axis_names = ("x", "y", "z")
    return ",".join(("" if sign > 0 else "-") + axis_names[axis] for axis, sign in zip(axis_order, axis_signs))


def zero_alignment(atom_count: int) -> dict:
    return {
        "compared": False,
        "kept_atoms": int(atom_count),
        "target_atoms": int(atom_count),
        "alignment_ratio": 1.0,
        "removed_atoms": 0,
        "shift": [0.0, 0.0, 0.0],
        "axis_order": [0, 1, 2],
        "axis_signs": [1, 1, 1],
        "transform": "x,y,z",
        "complete": True,
        "match_digits_used": "",
    }


def alignment_summary(label_result, match_digits: int) -> dict:
    alignment_ratio = 1.0
    if label_result.target_count > 0:
        alignment_ratio = label_result.kept_count / label_result.target_count
    return {
        "compared": True,
        "kept_atoms": int(label_result.kept_count),
        "target_atoms": int(label_result.target_count),
        "alignment_ratio": float(alignment_ratio),
        "removed_atoms": int(label_result.removed_count),
        "shift": [float(value) for value in label_result.shift],
        "axis_order": [int(value) for value in label_result.axis_order],
        "axis_signs": [int(value) for value in label_result.axis_signs],
        "transform": format_transform(label_result.axis_order, label_result.axis_signs),
        "complete": bool(label_result.kept_count == label_result.target_count),
        "match_digits_used": int(match_digits),
    }


def split_ids(records: list[dict], val_ratio: float, test_ratio: float, seed: int) -> dict[str, list[str]]:
    rng = random.Random(seed)
    pos = [record["refcode"] for record in records if int(record["removed"]) > 0]
    neg = [record["refcode"] for record in records if int(record["removed"]) == 0]
    rng.shuffle(pos)
    rng.shuffle(neg)

    def split_group(ids: list[str]) -> tuple[list[str], list[str], list[str]]:
        n_total = len(ids)
        n_test = int(round(n_total * test_ratio))
        n_val = int(round(n_total * val_ratio))
        return ids[n_test + n_val:], ids[n_test:n_test + n_val], ids[:n_test]

    train_pos, val_pos, test_pos = split_group(pos)
    train_neg, val_neg, test_neg = split_group(neg)
    splits = {
        "train": train_pos + train_neg,
        "val": val_pos + val_neg,
        "test": test_pos + test_neg,
    }
    for ids in splits.values():
        rng.shuffle(ids)
    return splits


def task_files(cif_dir: Path, row: dict, task: str) -> tuple[Path, str, Path | None, str | None]:
    refcode = row["refcode"]
    r_file = pick_cif(cif_dir, refcode, ("_R.cif",), required=True)
    f_file = pick_cif(cif_dir, refcode, ("_F.cif", "_IF.cif"), required=False)
    a_file = pick_cif(cif_dir, refcode, ("_A.cif", "_IA.cif"), required=False)

    if task == "free":
        if flag(row, "r2f"):
            if f_file is None:
                raise FileNotFoundError(f"{refcode} has r2f=1 but no F/IF CIF")
            return r_file, "R", f_file, "F"
        return r_file, "R", None, None

    if task == "coordinated":
        source_file = f_file if f_file is not None else r_file
        source_role = "F" if f_file is not None else "R"
        if flag(row, "f2a"):
            if a_file is None:
                raise FileNotFoundError(f"{refcode} has f2a=1 but no A/IA CIF")
            return source_file, source_role, a_file, "A"
        return source_file, source_role, None, None

    raise ValueError(f"Unknown task: {task}")


def make_labels(
    source_file: Path,
    target_file: Path | None,
    match_digits_options: list[int],
) -> tuple[np.ndarray, dict]:
    if isinstance(match_digits_options, int):
        match_digits_options = [match_digits_options]
    structure = read_cif_structure(source_file)
    if target_file is None:
        labels = np.zeros((len(structure.atoms),), dtype=np.float32)
        return labels, zero_alignment(len(structure.atoms))

    attempts: list[tuple[np.ndarray, dict]] = []
    for match_digits in match_digits_options:
        result = label_removed_atoms(source_file, target_file, match_digits=match_digits)
        summary = alignment_summary(result, match_digits)
        attempts.append((result.labels, summary))
        if summary["complete"]:
            return result.labels, summary

    best_labels, best_summary = max(
        attempts,
        key=lambda item: (int(item[1]["kept_atoms"]), float(item[1]["alignment_ratio"])),
    )
    best_summary = dict(best_summary)
    best_summary["attempts"] = [
        {
            "match_digits": int(summary["match_digits_used"]),
            "kept_atoms": int(summary["kept_atoms"]),
            "target_atoms": int(summary["target_atoms"]),
            "alignment_ratio": float(summary["alignment_ratio"]),
            "transform": summary["transform"],
            "shift": summary["shift"],
        }
        for _, summary in attempts
    ]
    return best_labels, best_summary


def process_one(
    row: dict,
    task: str,
    cif_dir: Path,
    output_dir: Path,
    radius: float,
    max_neighbors: int,
    match_digits: list[int],
    write_marked_cif: bool,
    allow_partial_alignment: bool,
    min_alignment_ratio: float,
) -> dict:
    refcode = row["refcode"]
    source_file, source_role, target_file, target_role = task_files(cif_dir, row, task)
    structure = read_cif_structure(source_file)
    labels, alignment = make_labels(source_file, target_file, match_digits)
    if labels.shape[0] != len(structure.atoms):
        raise ValueError(f"Label length mismatch for {refcode}/{task}")
    if alignment["compared"] and not alignment["complete"] and not allow_partial_alignment:
        attempts = alignment.get("attempts", [])
        raise ValueError(
            f"{task} target alignment incomplete after match_digits={match_digits}: "
            f"best kept {alignment['kept_atoms']}/{alignment['target_atoms']} "
            f"({float(alignment['alignment_ratio']):.3f}), "
            f"best_digits={alignment['match_digits_used']}, attempts={attempts}"
        )
    if alignment["compared"] and allow_partial_alignment and float(alignment["alignment_ratio"]) < min_alignment_ratio:
        raise ValueError(
            f"{task} target alignment below threshold: kept {alignment['kept_atoms']}/"
            f"{alignment['target_atoms']} ({float(alignment['alignment_ratio']):.3f}), "
            f"threshold={min_alignment_ratio:.3f}, transform={alignment['transform']}, "
            f"shift={alignment['shift']}"
        )

    graph = build_crystal_graph(source_file, radius=radius, max_neighbors=max_neighbors)
    graph["metadata"] = {
        "refcode": refcode,
        "task": task,
        "source_role": source_role,
        "source_file": source_file.name,
        "target_role": target_role,
        "target_file": target_file.name if target_file is not None else None,
        "r2f": flag(row, "r2f"),
        "f2a": flag(row, "f2a"),
        "target_alignment": alignment,
    }

    task_dir = output_dir / task
    (task_dir / "json").mkdir(parents=True, exist_ok=True)
    (task_dir / "npy" / "pos").mkdir(parents=True, exist_ok=True)
    (task_dir / "npy" / "label").mkdir(parents=True, exist_ok=True)

    with (task_dir / "json" / f"{refcode}.json").open("w", encoding="utf-8") as handle:
        json.dump(graph, handle)
    np.save(task_dir / "npy" / "pos" / f"{refcode}.npy", structure.frac_coords.astype(np.float32))
    np.save(task_dir / "npy" / "label" / f"{refcode}.npy", labels.astype(np.float32))

    if write_marked_cif:
        write_labeled_cif(
            source_file,
            output_dir / "labeled_cif" / task / f"{refcode}_{source_role}_labeled.cif",
            [(f"_atom_site_{task}_solvent_label", labels)],
        )

    return {
        "refcode": refcode,
        "task": task,
        "source_role": source_role,
        "source_file": source_file.name,
        "target_role": target_role or "",
        "target_file": target_file.name if target_file is not None else "",
        "atom_count": len(structure.atoms),
        "edge_count": len(graph["index1"]),
        "removed": int(labels.sum()),
        "kept_atoms": int(alignment["kept_atoms"]),
        "target_atoms": int(alignment["target_atoms"]),
        "alignment_ratio": f"{float(alignment['alignment_ratio']):.6f}",
        "alignment_complete": int(alignment["complete"]),
        "match_digits_used": alignment["match_digits_used"],
        "shift": ";".join(f"{float(value):.4f}" for value in alignment["shift"]),
        "transform": alignment["transform"],
        "r2f": flag(row, "r2f"),
        "f2a": flag(row, "f2a"),
    }


def process_job(job: dict) -> tuple[bool, dict]:
    row = job["row"]
    task = job["task"]
    try:
        record = process_one(
            row=row,
            task=task,
            cif_dir=job["cif_dir"],
            output_dir=job["output_dir"],
            radius=job["radius"],
            max_neighbors=job["max_neighbors"],
            match_digits=job["match_digits"],
            write_marked_cif=job["write_marked_cif"],
            allow_partial_alignment=job["allow_partial_alignment"],
            min_alignment_ratio=job["min_alignment_ratio"],
        )
        return True, record
    except Exception as exc:
        payload = {"refcode": row.get("refcode", ""), "task": task, "error": str(exc)}
        try:
            source_file, source_role, target_file, target_role = task_files(job["cif_dir"], row, task)
            payload.update(
                {
                    "source_role": source_role,
                    "source_file": str(source_file),
                    "target_role": target_role or "",
                    "target_file": str(target_file) if target_file is not None else "",
                }
            )
        except Exception:
            payload.update({"source_role": "", "source_file": "", "target_role": "", "target_file": ""})
        return False, payload


def build_jobs(rows: list[dict], args, cif_dir: Path) -> list[dict]:
    jobs: list[dict] = []
    for row in rows:
        for task in TASKS:
            jobs.append(
                {
                    "row": row,
                    "task": task,
                    "cif_dir": cif_dir,
                    "output_dir": args.output_dir,
                    "radius": args.radius,
                    "max_neighbors": args.max_neighbors,
                    "match_digits": args.match_digits,
                    "write_marked_cif": args.write_marked_cif,
                    "allow_partial_alignment": args.allow_partial_alignment,
                    "min_alignment_ratio": args.min_alignment_ratio,
                }
            )
    return jobs


def collect_jobs(jobs: list[dict], num_workers: int) -> tuple[list[dict], list[dict]]:
    all_records: list[dict] = []
    errors: list[dict] = []

    if num_workers <= 1:
        for job in tqdm(jobs, total=len(jobs), desc="Preparing GTsR3 data", unit="job", dynamic_ncols=True):
            success, payload = process_job(job)
            if success:
                all_records.append(payload)
            else:
                errors.append(payload)
    else:
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            future_to_job = {executor.submit(process_job, job): job for job in jobs}
            futures = as_completed(future_to_job)
            for future in tqdm(
                futures,
                total=len(future_to_job),
                desc=f"Preparing GTsR3 data ({num_workers} workers)",
                unit="job",
                dynamic_ncols=True,
            ):
                job = future_to_job[future]
                try:
                    success, payload = future.result()
                except Exception as exc:
                    payload = {"refcode": job["row"].get("refcode", ""), "task": job["task"], "error": str(exc)}
                    success = False
                if success:
                    all_records.append(payload)
                else:
                    errors.append(payload)

    all_records.sort(key=lambda record: (record["task"], record["refcode"]))
    errors.sort(key=lambda error: (error["task"], error["refcode"]))
    return all_records, errors


def write_summary(path: Path, records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "refcode",
        "task",
        "source_role",
        "source_file",
        "target_role",
        "target_file",
        "atom_count",
        "edge_count",
        "removed",
        "kept_atoms",
        "target_atoms",
        "alignment_ratio",
        "alignment_complete",
        "match_digits_used",
        "shift",
        "transform",
        "r2f",
        "f2a",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def write_error_report(output_dir: Path, errors: list[dict]) -> Path:
    error_dir = output_dir / "error"
    error_dir.mkdir(parents=True, exist_ok=True)
    error_csv = error_dir / "alignment_errors.csv"
    fieldnames = ["refcode", "task", "source_role", "source_file", "target_role", "target_file", "error"]
    with error_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for error in errors:
            writer.writerow({name: error.get(name, "") for name in fieldnames})

    for error in errors:
        refcode = error.get("refcode") or "unknown"
        task = error.get("task") or "unknown"
        case_dir = error_dir / refcode / task
        case_dir.mkdir(parents=True, exist_ok=True)
        source_file = Path(error["source_file"]) if error.get("source_file") else None
        target_file = Path(error["target_file"]) if error.get("target_file") else None
        if source_file is not None and source_file.exists():
            shutil.copy2(source_file, case_dir / "source.cif")
        if target_file is not None and target_file.exists():
            shutil.copy2(target_file, case_dir / "target.cif")
        (case_dir / "reason.txt").write_text(str(error.get("error", "")) + "\n", encoding="utf-8")
    return error_csv


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare GTsR3 datasets from clean_dataset_test.")
    parser.add_argument("--dataset-dir", type=Path, default=PROJECT_ROOT / "clean_dataset_test")
    parser.add_argument("--labels-csv", type=Path, default=PROJECT_ROOT / "dataset_labels.csv")
    parser.add_argument("--output-dir", type=Path, default=SCRIPT_DIR / "data")
    parser.add_argument("--radius", type=float, default=8.0)
    parser.add_argument("--max-neighbors", type=int, default=12)
    parser.add_argument(
        "--match-digits",
        type=int,
        nargs="+",
        default=[3, 2],
        help="Coordinate matching precision attempts. Default tries 3 digits, then 2 digits.",
    )
    parser.add_argument("--val-ratio", type=float, default=0.1)
    parser.add_argument("--test-ratio", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=1126)
    parser.add_argument("--write-marked-cif", action="store_true")
    parser.add_argument("--allow-partial-alignment", action="store_true")
    parser.add_argument(
        "--min-alignment-ratio",
        type=float,
        default=0.5,
        help="Only fail compared samples when kept_atoms / target_atoms is below this value.",
    )
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--num-workers", type=int, default=1, help="Parallel worker processes. Use 1 for serial mode.")
    args = parser.parse_args()

    cif_dir = resolve_cif_dir(args.dataset_dir)
    rows = read_label_rows(args.labels_csv)
    if args.limit is not None:
        rows = rows[:args.limit]

    jobs = build_jobs(rows, args, cif_dir)
    all_records, errors = collect_jobs(jobs, max(1, args.num_workers))
    records_by_task: dict[str, list[dict]] = {
        task: [record for record in all_records if record["task"] == task]
        for task in TASKS
    }

    write_summary(args.output_dir / "label_summary.csv", all_records)
    errors_path = args.output_dir / "prepare_errors.csv"
    error_report_path = None
    if errors:
        error_report_path = write_error_report(args.output_dir, errors)
        with errors_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["refcode", "task", "source_file", "target_file", "error"])
            writer.writeheader()
            writer.writerows({name: error.get(name, "") for name in writer.fieldnames} for error in errors)
    elif errors_path.exists():
        errors_path.unlink()

    split_outputs = []
    for task, records in records_by_task.items():
        splits = split_ids(records, args.val_ratio, args.test_ratio, args.seed)
        for split_name, split_ids_ in splits.items():
            split_outputs.append((task, split_name, split_ids_))

    for task, split_name, split_ids_ in tqdm(
        split_outputs,
        total=len(split_outputs),
        desc="Writing split CSVs",
        unit="split",
        dynamic_ncols=True,
    ):
        write_id_csv(args.output_dir / "splits" / f"{task}_{split_name}.csv", split_ids_)

    print(f"prepared samples: {len(all_records)}")
    print(f"errors: {len(errors)}")
    print(f"summary: {args.output_dir / 'label_summary.csv'}")
    if error_report_path is not None:
        print(f"alignment errors: {error_report_path}")
    print(f"splits: {args.output_dir / 'splits'}")


if __name__ == "__main__":
    main()
