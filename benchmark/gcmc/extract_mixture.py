import re
import sys
import argparse
import pandas as pd
from pathlib import Path
from tqdm import tqdm


UNIT_ALIASES = {
    "mol/kg":            "mol/kg-framework",
    "mol/kg-framework":  "mol/kg-framework",
    "molecules/cell":    "molecules/cell",
    "molecules/uc":      "molecules/uc",
    "mg/g":              "mg/g-framework",
    "mg/g-framework":    "mg/g-framework",
}

RE_COMPONENT = re.compile(r"^Component\s+(\d+)\s+\((\w+)\)", re.MULTILINE)

RE_LOADING = re.compile(
    r"(Abs\.|Excess)\s+loading average\s+([\d.e+\-]+)\s+\+/-\s+([\d.e+\-]+)\s+\[([^\]]+)\]"
)


def parse_components(text: str) -> list[str]:
    matches = RE_COMPONENT.findall(text)
    by_idx: dict[int, str] = {}
    for idx, name in matches:
        by_idx[int(idx)] = name.lower()
    return [by_idx[i] for i in sorted(by_idx)]


def parse_loadings(text: str, unit: str, loading_type: str) -> dict[str, tuple[float, float]]:
    sections: list[tuple[str, str]] = []
    boundaries = list(RE_COMPONENT.finditer(text))
    for i, m in enumerate(boundaries):
        name = m.group(2).lower()
        start = m.start()
        end = boundaries[i + 1].start() if i + 1 < len(boundaries) else len(text)
        sections.append((name, text[start:end]))

    prefix = "Abs." if loading_type.lower() == "abs" else "Excess"
    result: dict[str, tuple[float, float]] = {}
    for name, sec in sections:
        for m in RE_LOADING.finditer(sec):
            ltype, val, err, u = m.group(1), m.group(2), m.group(3), m.group(4)
            if ltype == prefix and u.strip() == unit:
                result[name] = (float(val), float(err))
    return result


def extract_file(txt_path: Path, unit: str, loading_type: str) -> dict[str, tuple[float, float]] | None:
    text = txt_path.read_text(errors="ignore")
    loadings = parse_loadings(text, unit, loading_type)
    return loadings if loadings else None


def build_columns(components: list[str]) -> list[str]:
    cols = ["name"]
    for c in components:
        cols += [f"{c} uptake", f"{c} error"]
    return cols


def main():
    parser = argparse.ArgumentParser(
        description="Extract mixture GCMC loadings from RASPA output files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("jobs_dir", help="Root folder containing per-MOF sub-folders")
    parser.add_argument("--unit", default="mol/kg-framework",
                        help="Loading unit (default: mol/kg-framework)")
    parser.add_argument("--loading", default="abs", choices=["abs", "excess"],
                        help="Abs or excess loading (default: abs)")
    parser.add_argument("--out", default=None,
                        help="Output CSV path (default: <jobs_dir>/mixture_results.csv)")
    args = parser.parse_args()

    unit = UNIT_ALIASES.get(args.unit.lower().strip(), args.unit.strip())
    root = Path(args.jobs_dir).resolve()
    out_path = Path(args.out) if args.out else root / "mixture_results.csv"

    records: list[dict] = []
    all_components: list[str] | None = None

    subdirs = sorted(d for d in root.iterdir() if d.is_dir())
    for sub in tqdm(subdirs, desc="Parsing"):
        txts = list(sub.glob("*.txt")) or list(sub.glob("output/*.txt"))
        if not txts:
            continue

        txt = txts[0]
        text = txt.read_text(errors="ignore")

        components = parse_components(text)
        if not components:
            tqdm.write(f"[SKIP]: {sub.name}")
            continue
        if all_components is None:
            all_components = components
            print(f"[Detected]: {components}")

        loadings = parse_loadings(text, unit, args.loading)
        if not loadings:
            tqdm.write(f"[SKIP]: {sub.name}")
            continue

        row: dict = {"name": sub.name}
        for comp in components:
            if comp in loadings:
                val, err = loadings[comp]
                row[f"{comp} uptake"] = val
                row[f"{comp} error"] = err
            else:
                row[f"{comp} uptake"] = None
                row[f"{comp} error"] = None
        records.append(row)

    cols = build_columns(all_components or [])
    df = pd.DataFrame(records, columns=[c for c in cols if c in records[0] or c == "name"])
    df.to_csv(out_path, index=False)


if __name__ == "__main__":
    main()
