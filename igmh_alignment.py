"""Generic active-site alignment helpers."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class KabschResult:
    rotation: list[list[float]]
    translation: list[float]
    rmsd: float


def kabsch_align(reference_xyz, target_xyz) -> KabschResult:
    """Return rigid transform that maps target coordinates onto reference coordinates."""
    ref = np.asarray(reference_xyz, dtype=float)
    target = np.asarray(target_xyz, dtype=float)
    if ref.shape != target.shape:
        raise ValueError(f"Reference and target shapes differ: {ref.shape} vs {target.shape}")
    if ref.ndim != 2 or ref.shape[1] != 3:
        raise ValueError("Coordinates must have shape (N, 3)")
    if ref.shape[0] < 3:
        raise ValueError("At least three corresponding atoms are required for Kabsch alignment")

    ref_centroid = ref.mean(axis=0)
    target_centroid = target.mean(axis=0)
    ref_centered = ref - ref_centroid
    target_centered = target - target_centroid
    covariance = target_centered.T @ ref_centered
    v, _s, wt = np.linalg.svd(covariance)
    determinant = np.linalg.det(v @ wt)
    correction = np.diag([1.0, 1.0, np.sign(determinant)])
    rotation = v @ correction @ wt
    aligned = target_centered @ rotation + ref_centroid
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - ref) ** 2, axis=1))))
    translation = ref_centroid - target_centroid @ rotation
    return KabschResult(rotation=rotation.tolist(), translation=translation.tolist(), rmsd=rmsd)


def apply_transform(xyz, rotation, translation):
    coords = np.asarray(xyz, dtype=float)
    rot = np.asarray(rotation, dtype=float)
    trans = np.asarray(translation, dtype=float)
    return coords @ rot + trans
