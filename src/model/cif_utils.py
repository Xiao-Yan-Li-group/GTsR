from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import permutations, product
from pathlib import Path
from typing import Iterable

import numpy as np


PERIODIC_TABLE_SYMBOLS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs",
    "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]
SYMBOL_TO_Z = {symbol: i + 1 for i, symbol in enumerate(PERIODIC_TABLE_SYMBOLS)}


@dataclass(frozen=True)
class CifAtom:
    symbol: str
    label: str
    frac: tuple[float, float, float]


@dataclass(frozen=True)
class CifStructure:
    path: Path
    cell: dict[str, float]
    atoms: list[CifAtom]
    pmg_structure: object | None = None

    @property
    def frac_coords(self) -> np.ndarray:
        return np.array([atom.frac for atom in self.atoms], dtype=float)

    @property
    def symbols(self) -> list[str]:
        return [atom.symbol for atom in self.atoms]

    @property
    def numbers(self) -> list[int]:
        return [atomic_number(atom.symbol) for atom in self.atoms]


@dataclass(frozen=True)
class LabelResult:
    labels: np.ndarray
    shift: tuple[float, float, float]
    axis_order: tuple[int, int, int]
    axis_signs: tuple[int, int, int]
    kept_count: int
    target_count: int

    @property
    def removed_count(self) -> int:
        return int(self.labels.sum())


@dataclass(frozen=True)
class Alignment:
    shift: tuple[float, float, float]
    axis_order: tuple[int, int, int]
    axis_signs: tuple[int, int, int]
    score: int


def atomic_number(symbol: str) -> int:
    try:
        return SYMBOL_TO_Z[symbol]
    except KeyError as exc:
        raise ValueError(f"Unsupported element symbol: {symbol}") from exc


def clean_symbol(raw: str) -> str:
    raw = raw.strip().strip("'\"")
    match = re.search(r"[A-Za-z]+", raw)
    if not match:
        raise ValueError(f"Cannot infer element symbol from {raw!r}")
    symbol = match.group(0)
    if len(symbol) == 1:
        return symbol.upper()
    return symbol[0].upper() + symbol[1:].lower()


def _structure_from_ase_pymatgen(path: Path):
    errors: list[str] = []

    try:
        from ase.io import read  # type: ignore
        from pymatgen.io.ase import AseAtomsAdaptor  # type: ignore

        atoms = read(str(path))
        return AseAtomsAdaptor.get_structure(atoms)
    except Exception as exc:
        errors.append(f"ASE read + pymatgen conversion failed: {type(exc).__name__}: {exc}")

    try:
        from pymatgen.core import Structure  # type: ignore

        return Structure.from_file(str(path))
    except Exception as exc:
        errors.append(f"Structure.from_file failed: {type(exc).__name__}: {exc}")

    raise ValueError(f"Could not read {path}. Attempts: {' | '.join(errors)}")


def _cell_from_pymatgen_structure(structure) -> dict[str, float]:
    lengths = structure.lattice.abc
    angles = structure.lattice.angles
    return {
        "_cell_length_a": float(lengths[0]),
        "_cell_length_b": float(lengths[1]),
        "_cell_length_c": float(lengths[2]),
        "_cell_angle_alpha": float(angles[0]),
        "_cell_angle_beta": float(angles[1]),
        "_cell_angle_gamma": float(angles[2]),
    }


def _site_symbol(site) -> str:
    try:
        raw = str(site.specie)
    except Exception:
        raw = str(site.species)
    return clean_symbol(raw)


def read_cif_structure(cif_path: str | Path) -> CifStructure:
    path = Path(cif_path)
    structure = _structure_from_ase_pymatgen(path)
    atoms = pymatgen_structure_to_cif_atoms(structure)
    if not atoms:
        raise ValueError(f"No atoms parsed from {path}")
    return CifStructure(
        path=path,
        cell=_cell_from_pymatgen_structure(structure),
        atoms=atoms,
        pmg_structure=structure,
    )


def pymatgen_structure_to_cif_atoms(structure) -> list[CifAtom]:
    atoms: list[CifAtom] = []
    for atom_index, site in enumerate(structure.sites):
        symbol = _site_symbol(site)
        label = getattr(site, "label", None) or f"{symbol}{atom_index + 1}"
        frac = (
            float(site.frac_coords[0]) % 1.0,
            float(site.frac_coords[1]) % 1.0,
            float(site.frac_coords[2]) % 1.0,
        )
        atoms.append(CifAtom(symbol=symbol, label=str(label), frac=frac))
    return atoms


def read_pymatgen_structure(cif_path: str | Path):
    return _structure_from_ase_pymatgen(Path(cif_path))


def _site_property_key(name: str) -> str:
    return name[len("_atom_site_"):] if name.startswith("_atom_site_") else name


def _copy_with_site_properties(structure, label_columns: list[tuple[str, Iterable[float | int]]]):
    copied = structure.copy()
    atom_count = len(copied)
    for name, values in label_columns:
        values = list(values)
        if len(values) != atom_count:
            raise ValueError(f"Column {name} has {len(values)} values, but structure has {atom_count} atoms")
        copied.add_site_property(_site_property_key(name), values)
    return copied


def subset_pymatgen_structure(structure, keep_mask: Iterable[bool]):
    mask = np.array(list(keep_mask), dtype=bool)
    if len(mask) != len(structure):
        raise ValueError(f"keep_mask has {len(mask)} values, but structure has {len(structure)} atoms")
    if not mask.any():
        return None
    from pymatgen.core import Structure  # type: ignore

    return Structure.from_sites([structure[int(i)] for i in np.where(mask)[0]])


def write_pymatgen_cif(structure, output_cif: str | Path, write_site_properties: bool = False) -> None:
    if len(structure) == 0:
        raise ValueError("pymatgen CifWriter cannot write an empty structure")
    from pymatgen.io.cif import CifWriter  # type: ignore

    out_path = Path(output_cif)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = CifWriter(
        structure,
        refine_struct=False,
        significant_figures=8,
        write_site_properties=write_site_properties,
    )
    out_path.write_text(str(writer), encoding="utf-8")


def write_structure_subset_cif(
    template_cif: str | Path,
    output_cif: str | Path,
    keep_mask: Iterable[bool],
    label_columns: list[tuple[str, Iterable[float | int]]] | None = None,
) -> int:
    structure = read_pymatgen_structure(template_cif)
    if label_columns:
        structure = _copy_with_site_properties(structure, label_columns)
    subset = subset_pymatgen_structure(structure, keep_mask)
    if subset is None:
        out_path = Path(output_cif)
        if out_path.exists():
            out_path.unlink()
        return 0
    write_pymatgen_cif(subset, output_cif, write_site_properties=bool(label_columns))
    return len(subset)


def cell_to_lattice(cell: dict[str, float]) -> np.ndarray:
    a = cell["_cell_length_a"]
    b = cell["_cell_length_b"]
    c = cell["_cell_length_c"]
    alpha = math.radians(cell["_cell_angle_alpha"])
    beta = math.radians(cell["_cell_angle_beta"])
    gamma = math.radians(cell["_cell_angle_gamma"])

    cos_alpha = math.cos(alpha)
    cos_beta = math.cos(beta)
    cos_gamma = math.cos(gamma)
    sin_gamma = math.sin(gamma)
    if abs(sin_gamma) < 1e-8:
        raise ValueError("Invalid lattice: sin(gamma) is too small")

    ax, ay, az = a, 0.0, 0.0
    bx, by, bz = b * cos_gamma, b * sin_gamma, 0.0
    cx = c * cos_beta
    cy = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cz_sq = c * c - cx * cx - cy * cy
    cz = math.sqrt(max(cz_sq, 0.0))
    return np.array([[ax, ay, az], [bx, by, bz], [cx, cy, cz]], dtype=float)


def pbc_fractional_delta(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    delta = a - b
    return delta - np.round(delta)


def _graph_with_pymatgen(
    structure,
    atom_count: int,
    radius: float,
    max_neighbors: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    try:
        if len(structure) != atom_count:
            return None
        center_idx, nbr_idx, _, distances = structure.get_neighbor_list(
            r=radius,
            numerical_tol=0,
            exclude_self=True,
        )
    except Exception:
        return None

    selected: list[np.ndarray] = []
    for atom_i in range(atom_count):
        edge_idx = np.nonzero(center_idx == atom_i)[0]
        if len(edge_idx) == 0:
            continue
        selected.append(edge_idx[np.argsort(distances[edge_idx])[:max_neighbors]])
    if not selected:
        return (
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            np.array([], dtype=float),
        )
    keep = np.concatenate(selected)
    return center_idx[keep], nbr_idx[keep], distances[keep]


def _graph_with_minimum_image(
    structure: CifStructure,
    radius: float,
    max_neighbors: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    frac = structure.frac_coords
    lattice = cell_to_lattice(structure.cell)
    index1: list[int] = []
    index2: list[int] = []
    dij: list[float] = []

    for i in range(len(frac)):
        candidates: list[tuple[float, int]] = []
        for j in range(len(frac)):
            if i == j:
                continue
            delta = pbc_fractional_delta(frac[j], frac[i])
            distance = float(np.linalg.norm(delta @ lattice))
            if 1e-8 < distance <= radius:
                candidates.append((distance, j))
        candidates.sort(key=lambda item: item[0])
        for distance, j in candidates[:max_neighbors]:
            index1.append(i)
            index2.append(j)
            dij.append(distance)

    return (
        np.array(index1, dtype=np.int64),
        np.array(index2, dtype=np.int64),
        np.array(dij, dtype=float),
    )


def build_crystal_graph(
    cif_path: str | Path,
    radius: float = 8.0,
    max_neighbors: int = 12,
    prefer_pymatgen: bool = True,
) -> dict:
    structure = read_cif_structure(cif_path)
    atom_count = len(structure.atoms)

    graph_edges = None
    if prefer_pymatgen and structure.pmg_structure is not None:
        graph_edges = _graph_with_pymatgen(structure.pmg_structure, atom_count, radius, max_neighbors)
    if graph_edges is None:
        graph_edges = _graph_with_minimum_image(structure, radius, max_neighbors)
    index1, index2, dij = graph_edges

    return {
        "rcut": float(radius),
        "max_neighbors": int(max_neighbors),
        "numbers": structure.numbers,
        "index1": index1.astype(int).tolist(),
        "index2": index2.astype(int).tolist(),
        "dij": dij.astype(float).tolist(),
    }


def _coord_key(frac: Iterable[float], digits: int) -> tuple[int, int, int]:
    scale = 10**digits
    return tuple(
        int(math.floor((float(value) % 1.0) * scale + 0.5 + 1e-8)) % scale
        for value in frac
    )


def _shift_key(shift: Iterable[float], digits: int) -> tuple[int, int, int]:
    return _coord_key(shift, digits)


def _key_to_shift(key: tuple[int, int, int], digits: int) -> tuple[float, float, float]:
    scale = float(10**digits)
    return tuple(value / scale for value in key)


def _atom_counter(
    atoms: list[CifAtom],
    shift: tuple[float, float, float],
    digits: int,
) -> Counter:
    shift_array = np.array(shift, dtype=float)
    counter: Counter = Counter()
    for atom in atoms:
        shifted = (np.array(atom.frac, dtype=float) + shift_array) % 1.0
        counter[(atom.symbol, _coord_key(shifted, digits))] += 1
    return counter


def _transform_atoms(
    atoms: list[CifAtom],
    axis_order: tuple[int, int, int],
    axis_signs: tuple[int, int, int],
) -> list[CifAtom]:
    order = np.array(axis_order, dtype=int)
    signs = np.array(axis_signs, dtype=float)
    transformed: list[CifAtom] = []
    for atom in atoms:
        frac = (np.array(atom.frac, dtype=float)[order] * signs) % 1.0
        transformed.append(CifAtom(symbol=atom.symbol, label=atom.label, frac=tuple(frac.tolist())))
    return transformed


def _score_shift(
    reference_atoms: list[CifAtom],
    target_atoms: list[CifAtom],
    shift: tuple[float, float, float],
    digits: int,
) -> int:
    target_counter = _atom_counter(target_atoms, shift, digits)
    score = 0
    for atom in reference_atoms:
        key = (atom.symbol, _coord_key(atom.frac, digits))
        if target_counter[key] > 0:
            target_counter[key] -= 1
            score += 1
    return score


def infer_alignment_shift(
    reference_atoms: list[CifAtom],
    target_atoms: list[CifAtom],
    digits: int = 4,
    max_candidates: int = 64,
    pair_budget: int = 100_000,
) -> tuple[float, float, float]:
    if not reference_atoms or not target_atoms:
        return (0.0, 0.0, 0.0)

    ref_by_symbol: dict[str, list[np.ndarray]] = defaultdict(list)
    tgt_by_symbol: dict[str, list[np.ndarray]] = defaultdict(list)
    for atom in reference_atoms:
        ref_by_symbol[atom.symbol].append(np.array(atom.frac, dtype=float))
    for atom in target_atoms:
        tgt_by_symbol[atom.symbol].append(np.array(atom.frac, dtype=float))

    common = [
        (symbol, len(ref_by_symbol[symbol]) * len(tgt_by_symbol[symbol]))
        for symbol in ref_by_symbol.keys() & tgt_by_symbol.keys()
    ]
    common.sort(key=lambda item: item[1])

    candidate_counts: Counter = Counter()
    candidate_counts[_shift_key((0.0, 0.0, 0.0), digits)] += 1
    pairs_seen = 0
    for symbol, product in common:
        refs = ref_by_symbol[symbol]
        tgts = tgt_by_symbol[symbol]
        if product == 0:
            continue
        if pairs_seen >= pair_budget:
            break

        if product > pair_budget - pairs_seen:
            limit = max(1, int(math.sqrt(pair_budget - pairs_seen)))
            refs = refs[:limit]
            tgts = tgts[:limit]

        for ref_frac in refs:
            for tgt_frac in tgts:
                shift = (ref_frac - tgt_frac) % 1.0
                candidate_counts[_shift_key(shift, digits)] += 1
                pairs_seen += 1
                if pairs_seen >= pair_budget:
                    break
            if pairs_seen >= pair_budget:
                break

    candidates = [
        _key_to_shift(key, digits)
        for key, _ in candidate_counts.most_common(max_candidates)
    ]
    best_shift = candidates[0]
    best_score = -1
    for shift in candidates:
        score = _score_shift(reference_atoms, target_atoms, shift, digits)
        if score > best_score:
            best_score = score
            best_shift = shift
    return best_shift


def infer_coordinate_alignment(
    reference_atoms: list[CifAtom],
    target_atoms: list[CifAtom],
    digits: int = 4,
) -> Alignment:
    best = Alignment(
        shift=(0.0, 0.0, 0.0),
        axis_order=(0, 1, 2),
        axis_signs=(1, 1, 1),
        score=-1,
    )
    target_possible = min(len(reference_atoms), len(target_atoms))

    def try_alignment(axis_order: tuple[int, int, int], axis_signs: tuple[int, int, int]) -> Alignment:
        nonlocal best
        transformed = _transform_atoms(target_atoms, axis_order, axis_signs)
        shift = infer_alignment_shift(reference_atoms, transformed, digits=digits)
        score = _score_shift(reference_atoms, transformed, shift, digits)
        alignment = Alignment(
            shift=shift,
            axis_order=axis_order,
            axis_signs=axis_signs,
            score=score,
        )
        if score > best.score:
            best = alignment
        return alignment

    first = try_alignment((0, 1, 2), (1, 1, 1))
    if first.score == target_possible:
        return first

    for axis_order in permutations((0, 1, 2)):
        axis_order = tuple(axis_order)
        if axis_order == (0, 1, 2):
            continue
        alignment = try_alignment(axis_order, (1, 1, 1))
        if alignment.score == target_possible:
            return alignment

    if best.score < int(0.8 * target_possible):
        for axis_order in permutations((0, 1, 2)):
            for axis_signs in product((1, -1), repeat=3):
                if axis_signs == (1, 1, 1):
                    continue
                alignment = try_alignment(tuple(axis_order), tuple(axis_signs))
                if alignment.score == target_possible:
                    return alignment

    return best


def label_removed_atoms(
    reference_cif: str | Path,
    target_cif: str | Path,
    match_digits: int = 4,
) -> LabelResult:
    reference = read_cif_structure(reference_cif)
    target = read_cif_structure(target_cif)
    return label_removed_atom_lists(reference.atoms, target.atoms, match_digits=match_digits)


def label_removed_atom_lists(
    reference_atoms: list[CifAtom],
    target_atoms: list[CifAtom],
    match_digits: int = 4,
) -> LabelResult:
    alignment = infer_coordinate_alignment(reference_atoms, target_atoms, digits=match_digits)
    transformed_target_atoms = _transform_atoms(
        target_atoms,
        alignment.axis_order,
        alignment.axis_signs,
    )
    target_counter = _atom_counter(transformed_target_atoms, alignment.shift, match_digits)

    labels: list[int] = []
    kept_count = 0
    for atom in reference_atoms:
        key = (atom.symbol, _coord_key(atom.frac, match_digits))
        if target_counter[key] > 0:
            labels.append(0)
            kept_count += 1
            target_counter[key] -= 1
        else:
            labels.append(1)

    return LabelResult(
        labels=np.array(labels, dtype=np.float32),
        shift=alignment.shift,
        axis_order=alignment.axis_order,
        axis_signs=alignment.axis_signs,
        kept_count=kept_count,
        target_count=len(target_atoms),
    )


def write_labeled_cif(
    template_cif: str | Path,
    output_cif: str | Path,
    label_columns: list[tuple[str, Iterable[float | int]]],
) -> None:
    structure = read_pymatgen_structure(template_cif)
    labeled = _copy_with_site_properties(structure, label_columns)
    write_pymatgen_cif(labeled, output_cif, write_site_properties=True)
