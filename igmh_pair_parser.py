"""Parsers for Multiwfn atom-pair delta-g and IBSIW text outputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import math
import re

BOHR_TO_ANGSTROM = 0.529177210903

ATOMIC_NUMBERS = {
    1: "H", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 14: "Si", 15: "P",
    16: "S", 17: "Cl", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 35: "Br", 53: "I",
}
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?")
ATOM_TOKEN_RE = re.compile(r"(?:^|\s)(\d+)\s*(?:\(([A-Z][a-z]?)\)|([A-Z][a-z]?))?")


@dataclass(frozen=True)
class AtomRecord:
    index: int
    element: str
    xyz_A: tuple[float, float, float]


@dataclass(frozen=True)
class PairRecord:
    atom1: int
    atom2: int
    element1: str | None = None
    element2: str | None = None
    distance_A: float | None = None
    deltaGpair: float | None = None
    IBSIW: float | None = None
    interaction_character: str | None = None

    @property
    def key(self) -> tuple[int, int]:
        return tuple(sorted((self.atom1, self.atom2)))


def _read_array_after_header(path: Path, header: str, expected_count: int, cast):
    values = []
    found = False
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not found:
                if line.lstrip().startswith(header):
                    found = True
                continue
            for token in line.split():
                try:
                    values.append(cast(token))
                except ValueError:
                    pass
                if len(values) == expected_count:
                    return values
    if len(values) != expected_count:
        raise ValueError(f"Could not read {expected_count} values after {header!r} in {path}")
    return values


def read_fchk_geometry(path: Path) -> dict[int, AtomRecord]:
    atom_count = None
    count_re = re.compile(r"^\s*Number of atoms\s+I\s+(?:N=\s*)?(\d+)")
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = count_re.search(line)
            if match:
                atom_count = int(match.group(1))
                break
    if not atom_count:
        raise ValueError(f"Could not read atom count from {path}")

    atomic_numbers = _read_array_after_header(path, "Atomic numbers", atom_count, int)
    coords = _read_array_after_header(path, "Current cartesian coordinates", atom_count * 3, float)
    atoms: dict[int, AtomRecord] = {}
    for idx, atomic_number in enumerate(atomic_numbers, start=1):
        xyz_values = coords[(idx - 1) * 3: idx * 3]
        xyz = (xyz_values[0] * BOHR_TO_ANGSTROM, xyz_values[1] * BOHR_TO_ANGSTROM, xyz_values[2] * BOHR_TO_ANGSTROM)
        atoms[idx] = AtomRecord(idx, ATOMIC_NUMBERS.get(atomic_number, f"Z{atomic_number}"), xyz)
    return atoms


def distance_A(atoms: dict[int, AtomRecord], atom1: int, atom2: int) -> float:
    xyz1 = atoms[atom1].xyz_A
    xyz2 = atoms[atom2].xyz_A
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(xyz1, xyz2)))


def _pair_and_value_from_line(line: str) -> tuple[int, int, str | None, str | None, float] | None:
    atom_matches = list(ATOM_TOKEN_RE.finditer(line))
    if len(atom_matches) < 2:
        return None
    try:
        atom1 = int(atom_matches[0].group(1))
        atom2 = int(atom_matches[1].group(1))
    except ValueError:
        return None
    if atom1 == atom2:
        return None
    element1 = atom_matches[0].group(2) or atom_matches[0].group(3)
    element2 = atom_matches[1].group(2) or atom_matches[1].group(3)

    tail_start = atom_matches[1].end()
    values = [float(match.group(0)) for match in FLOAT_RE.finditer(line[tail_start:])]
    if not values:
        return None
    return atom1, atom2, element1, element2, values[-1]


def parse_pair_value_file(path: Path) -> dict[tuple[int, int], tuple[str | None, str | None, float]]:
    """Parse the atom-pair section of a Multiwfn atmdg.txt or IBSIW.txt file."""
    pairs: dict[tuple[int, int], tuple[str | None, str | None, float]] = {}
    in_pair_section = False
    pair_line_re = re.compile(r"^\s*(\d+)\s+(\d+)\s*:\s*(" + FLOAT_RE.pattern + r")")
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if "Atomic pair delta-g indices" in line or "IBSIW index" in line:
                in_pair_section = True
                continue
            if not in_pair_section:
                parsed = _pair_and_value_from_line(line)
                if parsed and "(" in line:
                    atom1, atom2, element1, element2, value = parsed
                    pairs[tuple(sorted((atom1, atom2)))] = (element1, element2, value)
                continue
            if line.lstrip().startswith("Sum of"):
                break
            match = pair_line_re.search(line)
            if not match:
                parsed = _pair_and_value_from_line(line)
                if not parsed:
                    continue
                atom1, atom2, element1, element2, value = parsed
            else:
                atom1 = int(match.group(1))
                atom2 = int(match.group(2))
                element1 = element2 = None
                value = float(match.group(3))
            if atom1 != atom2:
                pairs[tuple(sorted((atom1, atom2)))] = (element1, element2, value)
    return pairs


def read_pair_outputs(directory: Path, fchk: Path | None = None) -> list[PairRecord]:
    atmdg_path = directory / "atmdg.txt"
    ibsiw_path = directory / "IBSIW.txt"
    delta_pairs = parse_pair_value_file(atmdg_path) if atmdg_path.exists() else {}
    ibsiw_pairs = parse_pair_value_file(ibsiw_path) if ibsiw_path.exists() else {}
    geometry = read_fchk_geometry(fchk) if fchk else {}

    records: list[PairRecord] = []
    for key in sorted(set(delta_pairs) | set(ibsiw_pairs)):
        atom1, atom2 = key
        element1 = element2 = None
        delta = None
        ibsiw = None
        if key in delta_pairs:
            element1, element2, delta = delta_pairs[key]
        if key in ibsiw_pairs:
            ie1, ie2, ibsiw = ibsiw_pairs[key]
            element1 = element1 or ie1
            element2 = element2 or ie2
        dist = None
        if geometry and atom1 in geometry and atom2 in geometry:
            element1 = geometry[atom1].element
            element2 = geometry[atom2].element
            dist = distance_A(geometry, atom1, atom2)
        records.append(PairRecord(atom1, atom2, element1, element2, dist, delta, ibsiw))
    return records


def rank_pairs(records: list[PairRecord], top_n: int = 10, exclude_pairs: set[tuple[int, int]] | None = None) -> list[PairRecord]:
    excluded = {tuple(sorted(pair)) for pair in (exclude_pairs or set())}
    candidates = [record for record in records if record.deltaGpair is not None and record.key not in excluded]
    return sorted(candidates, key=lambda record: abs(record.deltaGpair or 0.0), reverse=True)[:top_n]


def write_pair_csv(path: Path, system: str, records: list[PairRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "system", "atom1", "atom2", "element1", "element2", "fragment1_or_2",
            "distance_A", "deltaGpair", "IBSIW", "interaction_character",
        ])
        for record in records:
            writer.writerow([
                system,
                record.atom1,
                record.atom2,
                record.element1 or "",
                record.element2 or "",
                "interfragment",
                "" if record.distance_A is None else f"{record.distance_A:.6f}",
                "" if record.deltaGpair is None else f"{record.deltaGpair:.10g}",
                "" if record.IBSIW is None else f"{record.IBSIW:.10g}",
                record.interaction_character or "",
            ])
