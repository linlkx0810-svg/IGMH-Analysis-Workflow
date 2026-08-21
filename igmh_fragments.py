"""Fragment parsing and validation helpers for Multiwfn/VMD workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import csv
import re

DASH_TRANSLATION = str.maketrans({0x2013: "-", 0x2014: "-", 0x2212: "-"})


@dataclass(frozen=True)
class FragmentConfig:
    fragment1: str
    fragment2: str
    center_selection: str | None = None


@dataclass(frozen=True)
class FragmentValidationResult:
    atom_count: int
    fragment1_atoms: tuple[int, ...]
    fragment2_atoms: tuple[int, ...]
    uncovered_atoms: tuple[int, ...] = field(default_factory=tuple)
    warnings: tuple[str, ...] = field(default_factory=tuple)


class FragmentValidationError(ValueError):
    """Raised when fragment definitions are not scientifically usable."""


def _tokens(selection: str) -> list[str]:
    normalized = selection.translate(DASH_TRANSLATION)
    return [token for token in normalized.replace(",", " ").split() if token]


def parse_atom_selection(selection: str) -> list[int]:
    """Parse a one-based atom selection such as '1-4,8,10-12'."""
    if not selection or not selection.strip():
        raise FragmentValidationError("Atom selection is empty")

    atoms: list[int] = []
    seen: set[int] = set()
    duplicates: list[int] = []
    for token in _tokens(selection):
        if "-" in token:
            start_text, end_text = token.split("-", 1)
            if not start_text.isdigit() or not end_text.isdigit():
                raise FragmentValidationError(f"Invalid atom range: {token}")
            start = int(start_text)
            end = int(end_text)
            if start <= 0 or end <= 0 or end < start:
                raise FragmentValidationError(f"Invalid atom range: {token}")
            expanded = range(start, end + 1)
        else:
            if not token.isdigit():
                raise FragmentValidationError(f"Invalid atom index: {token}")
            atom = int(token)
            if atom <= 0:
                raise FragmentValidationError(f"Invalid atom index: {token}")
            expanded = (atom,)

        for atom in expanded:
            if atom in seen:
                duplicates.append(atom)
            seen.add(atom)
            atoms.append(atom)

    if duplicates:
        dup_text = ", ".join(str(atom) for atom in sorted(set(duplicates)))
        raise FragmentValidationError(f"Duplicated atom indices: {dup_text}")
    if not atoms:
        raise FragmentValidationError("Atom selection is empty")
    return atoms


def validate_fragments(fragment1: str, fragment2: str, atom_count: int) -> FragmentValidationResult:
    """Validate two one-based Multiwfn fragment selections against atom_count."""
    if atom_count <= 0:
        raise FragmentValidationError(f"Invalid atom count: {atom_count}")
    frag1 = parse_atom_selection(fragment1)
    frag2 = parse_atom_selection(fragment2)
    set1 = set(frag1)
    set2 = set(frag2)

    overlap = sorted(set1 & set2)
    if overlap:
        raise FragmentValidationError("Fragments overlap at atom(s): " + ", ".join(map(str, overlap)))

    out_of_range = sorted(atom for atom in set1 | set2 if atom > atom_count)
    if out_of_range:
        raise FragmentValidationError(
            "Atom index exceeds .fchk atom count "
            f"({atom_count}): " + ", ".join(map(str, out_of_range))
        )

    covered = set1 | set2
    uncovered = tuple(atom for atom in range(1, atom_count + 1) if atom not in covered)
    warnings: list[str] = []
    if uncovered:
        warnings.append(f"{len(uncovered)} atom(s) are not assigned to either fragment")

    return FragmentValidationResult(
        atom_count=atom_count,
        fragment1_atoms=tuple(frag1),
        fragment2_atoms=tuple(frag2),
        uncovered_atoms=uncovered,
        warnings=tuple(warnings),
    )


def fragment_to_vmd_indices(selection: str) -> str:
    atoms = parse_atom_selection(selection)
    parts: list[str] = []
    ranges: list[tuple[int, int]] = []
    start = prev = atoms[0]
    for atom in atoms[1:]:
        if atom == prev + 1:
            prev = atom
        else:
            ranges.append((start, prev))
            start = prev = atom
    ranges.append((start, prev))

    for start, end in ranges:
        if start == end:
            parts.append(str(start - 1))
        else:
            parts.append(f"{start - 1} to {end - 1}")
    return "index " + " ".join(parts)


def read_fragment_file(path: Path) -> dict[str, FragmentConfig]:
    configs: dict[str, FragmentConfig] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"system", "fragment1", "fragment2"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: system, fragment1, fragment2")
        for row in reader:
            system = (row.get("system") or "").strip()
            if not system:
                continue
            center_selection = (row.get("center_selection") or "").strip() or None
            configs[system] = FragmentConfig(
                fragment1=(row["fragment1"] or "").strip(),
                fragment2=(row["fragment2"] or "").strip(),
                center_selection=center_selection,
            )
    return configs


def read_fchk_atom_count(path: Path) -> int:
    pattern = re.compile(r"^\s*Number of atoms\s+I\s+(?:N=\s*)?(\d+)")
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = pattern.search(line)
            if match:
                return int(match.group(1))
    raise ValueError(f"Could not read atom count from {path}")
