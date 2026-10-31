"""Geometric transformations and alignment metrics in Angstrom units."""
from __future__ import annotations

import numpy as np


def kabsch_rmsd(mobile: np.ndarray, reference: np.ndarray) -> float:
    """Return C-alpha RMSD (Angstrom) after optimal rigid-body superposition."""
    a = np.asarray(mobile, dtype=float)
    b = np.asarray(reference, dtype=float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 3 or len(a) < 3:
        raise ValueError("Expected two matching (N, 3) coordinate arrays, N >= 3")
    ac = a - a.mean(axis=0)
    bc = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(ac.T @ bc)
    correction = np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt))])
    aligned = ac @ (u @ correction @ vt)
    return float(np.sqrt(np.mean(np.sum((aligned - bc) ** 2, axis=1))))


def crankshaft_move(
    coords: np.ndarray,
    rng: np.random.Generator,
    max_angle_radians: float = np.pi / 3,
    max_span: int = 12,
) -> np.ndarray:
    """Rotate residues strictly between two anchors about their connecting axis.

    This preserves all *adjacent* C-alpha distances up to floating-point error,
    although it does not enforce bond angles, chirality or a physical force field.
    """
    points = np.asarray(coords, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 4:
        raise ValueError("crankshaft_move requires a (N,3) array with N >= 4")
    if not np.isfinite(max_angle_radians) or max_angle_radians <= 0 or max_span < 2:
        raise ValueError("max_angle_radians must be positive and max_span >= 2")
    i = int(rng.integers(0, len(points) - 2))
    j = int(rng.integers(i + 2, min(len(points) - 1, i + max_span) + 1))
    origin = points[i]
    axis = points[j] - origin
    axis_length = np.linalg.norm(axis)
    if axis_length < 1e-12:
        return points.copy()
    axis = axis / axis_length
    theta = float(rng.uniform(-max_angle_radians, max_angle_radians))
    v = points[i + 1 : j] - origin
    rotated = (
        v * np.cos(theta)
        + np.cross(axis, v) * np.sin(theta)
        + np.outer(v @ axis, axis) * (1 - np.cos(theta))
    )
    proposal = points.copy()
    proposal[i + 1 : j] = origin + rotated
    return proposal