"""Comparison configuration and shared-view transform helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import json

import numpy as np

from igmh_alignment import kabsch_align

BOHR_TO_ANGSTROM = 0.529177210903


class ComparisonConfigError(ValueError):
    """Raised when a comparison-rendering configuration is incomplete."""


@dataclass(frozen=True)
class CubeAtoms:
    coords: dict[int, np.ndarray]
    unit: str


@dataclass(frozen=True)
class ComparisonPlan:
    panel_order: list[str]
    transform_tcl: Path
    alignment_rmsd: dict[str, float]
    alignment_rmsd_A: dict[str, float]
    alignment_atom_mapping: dict[str, dict[str, list[int]]]
    camera_matrix: list[list[float]]
    scale: float
    projection: str
    screen_translate: list[float]
    per_panel_screen_translate: dict[str, list[float]]
    volume_registration: str


def load_comparison_config(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    required = ["reference_system", "panel_order"]
    missing = [key for key in required if key not in data]
    if missing:
        raise ComparisonConfigError("Comparison config missing key(s): " + ", ".join(missing))
    if "alignment_mappings" not in data and "alignment" not in data:
        raise ComparisonConfigError("Comparison config requires alignment_mappings")
    return data


def read_cube_atoms(path: Path) -> CubeAtoms:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        handle.readline()
        handle.readline()
        parts = handle.readline().split()
        natoms_signed = int(float(parts[0]))
        natoms = abs(natoms_signed)
        unit = "angstrom" if natoms_signed < 0 else "bohr"
        handle.readline()
        handle.readline()
        handle.readline()
        coords: dict[int, np.ndarray] = {}
        for idx in range(1, natoms + 1):
            atom_parts = handle.readline().split()
            coords[idx] = np.array([float(atom_parts[2]), float(atom_parts[3]), float(atom_parts[4])], dtype=float)
    return CubeAtoms(coords=coords, unit=unit)


def _coords(cube: CubeAtoms, atoms: list[int]) -> np.ndarray:
    missing = [atom for atom in atoms if atom not in cube.coords]
    if missing:
        raise ComparisonConfigError("Alignment atom index not present in cube: " + ", ".join(map(str, missing)))
    return np.vstack([cube.coords[atom] for atom in atoms])


def _identity4() -> np.ndarray:
    return np.eye(4, dtype=float)


def _vmd_matrix_from_row_transform(rotation_row: list[list[float]], translation: list[float]) -> np.ndarray:
    matrix = np.eye(4, dtype=float)
    matrix[:3, :3] = np.asarray(rotation_row, dtype=float).T
    matrix[:3, 3] = np.asarray(translation, dtype=float)
    return matrix


def _row_transform_coords(coords: np.ndarray, matrix_vmd: np.ndarray) -> np.ndarray:
    return (matrix_vmd[:3, :3] @ coords.T).T + matrix_vmd[:3, 3]


def _matrix4_from_config(value: Any) -> np.ndarray:
    matrix = np.asarray(value, dtype=float)
    if matrix.shape == (3, 3):
        out = np.eye(4, dtype=float)
        out[:3, :3] = matrix
        return out
    if matrix.shape == (4, 4):
        return matrix
    raise ComparisonConfigError("view_matrix must be 3x3 or 4x4")


def _matrix_tcl(matrix: np.ndarray) -> str:
    rows = []
    for row in matrix:
        rows.append("{" + " ".join(f"{float(value):.10f}" for value in row) + "}")
    return "{" + " ".join(rows) + "}"


def _center_matrix_tcl(center: np.ndarray) -> str:
    matrix = np.eye(4, dtype=float)
    matrix[:3, 3] = -center
    return _matrix_tcl(matrix)


def build_comparison_plan(config_path: Path, cube_dirs: dict[str, Path], work_dir: Path) -> ComparisonPlan:
    config = load_comparison_config(config_path)
    reference = config["reference_system"]
    panel_order = list(config["panel_order"])
    if reference not in cube_dirs:
        raise ComparisonConfigError(f"Reference system {reference!r} is not among configured systems")
    missing_panels = [label for label in panel_order if label not in cube_dirs]
    if missing_panels:
        raise ComparisonConfigError("Panel order references unknown system(s): " + ", ".join(missing_panels))

    mappings = config.get("alignment_mappings") or config.get("alignment")
    if not isinstance(mappings, dict):
        raise ComparisonConfigError("alignment_mappings must be an object keyed by target system")

    cubes = {label: read_cube_atoms(cube_dirs[label] / "sl2r.cub") for label in cube_dirs}
    reference_cube = cubes[reference]
    globals_by_system: dict[str, np.ndarray] = {reference: _identity4()}
    rmsd: dict[str, float] = {}
    rmsd_A: dict[str, float] = {}
    atom_mapping: dict[str, dict[str, list[int]]] = {}

    unit_to_A = BOHR_TO_ANGSTROM if reference_cube.unit == "bohr" else 1.0
    first_reference_atoms: list[int] | None = None
    for target, mapping in mappings.items():
        if target not in cube_dirs:
            raise ComparisonConfigError(f"Alignment target {target!r} is not among configured systems")
        reference_atoms = list(mapping.get("reference_atoms", []))
        target_atoms = list(mapping.get("target_atoms", []))
        if len(reference_atoms) != len(target_atoms):
            raise ComparisonConfigError(f"{target}: reference_atoms and target_atoms must have the same length")
        if len(reference_atoms) < 3:
            raise ComparisonConfigError(f"{target}: at least three corresponding atoms are required")
        first_reference_atoms = first_reference_atoms or reference_atoms
        result = kabsch_align(_coords(reference_cube, reference_atoms), _coords(cubes[target], target_atoms))
        globals_by_system[target] = _vmd_matrix_from_row_transform(result.rotation, result.translation)
        rmsd[target] = result.rmsd
        rmsd_A[target] = result.rmsd * unit_to_A
        atom_mapping[target] = {"reference_atoms": reference_atoms, "target_atoms": target_atoms}

    if "center_atoms" in config:
        center_atoms = list(config["center_atoms"])
    elif first_reference_atoms:
        center_atoms = first_reference_atoms
    else:
        center_atoms = sorted(reference_cube.coords)
    focus = _coords(reference_cube, center_atoms).mean(axis=0)

    if "view_matrix" in config:
        camera = _matrix4_from_config(config["view_matrix"])
    else:
        camera = _identity4()
    projection = config.get("projection", "Orthographic")
    if projection != "Orthographic":
        raise ComparisonConfigError("Only Orthographic projection is currently supported for validated shared-scale comparisons")

    projected_sets = []
    for label in panel_order:
        coords = np.vstack([cubes[label].coords[idx] for idx in sorted(cubes[label].coords)])
        transformed = _row_transform_coords(coords, globals_by_system.get(label, _identity4()))
        centered = transformed - focus
        projected = (camera[:3, :3] @ centered.T).T
        projected_sets.append(projected[:, :2])
    union = np.vstack(projected_sets)
    mins = union.min(axis=0)
    maxs = union.max(axis=0)
    spans = np.maximum(maxs - mins, 1e-9)
    screen_span = float(config.get("screen_span", 1.72))
    scale = screen_span / float(max(spans[0], spans[1]))
    center_xy = (mins + maxs) / 2.0
    translate = [-float(scale * center_xy[0]), -float(scale * center_xy[1]), 0.0]
    if "screen_translate_delta" in config:
        delta = list(config["screen_translate_delta"])
        if len(delta) != 3:
            raise ComparisonConfigError("screen_translate_delta must contain three numbers")
        translate = [translate[i] + float(delta[i]) for i in range(3)]
    if "screen_translate" in config:
        explicit_translate = list(config["screen_translate"])
        if len(explicit_translate) != 3:
            raise ComparisonConfigError("screen_translate must contain three numbers")
        translate = [float(value) for value in explicit_translate]

    per_panel_translate: dict[str, list[float]] = {label: list(translate) for label in panel_order}
    if "per_panel_screen_translate" in config:
        for label, values in config["per_panel_screen_translate"].items():
            if label not in per_panel_translate:
                raise ComparisonConfigError(f"per_panel_screen_translate references unknown system {label!r}")
            value_list = list(values)
            if len(value_list) != 3:
                raise ComparisonConfigError("per_panel_screen_translate entries must contain three numbers")
            per_panel_translate[label] = [float(value) for value in value_list]
    if "per_panel_screen_translate_delta" in config:
        for label, values in config["per_panel_screen_translate_delta"].items():
            if label not in per_panel_translate:
                raise ComparisonConfigError(f"per_panel_screen_translate_delta references unknown system {label!r}")
            value_list = list(values)
            if len(value_list) != 3:
                raise ComparisonConfigError("per_panel_screen_translate_delta entries must contain three numbers")
            per_panel_translate[label] = [per_panel_translate[label][i] + float(value_list[i]) for i in range(3)]

    transform_tcl = work_dir / "comparison_transform.tcl"
    transform_tcl.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"set comparison_projection {projection}",
        f"set comparison_camera {_matrix_tcl(camera)}",
        f"set comparison_center {_center_matrix_tcl(focus)}",
        f"set comparison_scale {scale:.10f}",
        "set comparison_translate {" + " ".join(f"{value:.10f}" for value in translate) + "}",
    ]
    panel_entries = []
    for label, values in per_panel_translate.items():
        panel_entries.append("{" + label + "} {" + " ".join(f"{value:.10f}" for value in values) + "}")
    lines.append("array set comparison_translate_panel {" + " ".join(panel_entries) + "}")
    for label, matrix in globals_by_system.items():
        lines.append(f"set comparison_global({label}) {_matrix_tcl(matrix)}")
    transform_tcl.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return ComparisonPlan(
        panel_order=panel_order,
        transform_tcl=transform_tcl,
        alignment_rmsd=rmsd,
        alignment_rmsd_A=rmsd_A,
        alignment_atom_mapping=atom_mapping,
        camera_matrix=camera.tolist(),
        scale=scale,
        projection=projection,
        screen_translate=per_panel_translate.get(reference, translate),
        per_panel_screen_translate=per_panel_translate,
        volume_registration="PASS: VMD molecule global_matrix is applied to the molecule containing atoms and both cube volumes",
    )
