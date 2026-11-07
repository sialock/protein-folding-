"""Descriptive C-alpha geometry and energy diagnostics for sampled structures."""
from __future__ import annotations

import numpy as np

from .energy import GoModel
from .geometry import kabsch_rmsd

EXPECTED_CA_BOND_ANGSTROM = 3.8


def structure_diagnostics(coords: np.ndarray, model: GoModel) -> dict[str, float | int]:
    """Summarize coarse C-alpha geometry and toy-model energy for one structure.

    The virtual angle is the angle between the two vectors from a central C-alpha
    bead to its adjacent beads. The 3.8-Angstrom bond reference is a descriptive
    sanity check, not a complete stereochemical validation.
    """
    xyz = np.asarray(coords, dtype=float)
    if xyz.shape != model.reference.shape or not np.isfinite(xyz).all():
        raise ValueError("Coordinates must be finite and match the model reference")

    bond_lengths = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    if np.any(bond_lengths <= 1e-12):
        raise ValueError("Cannot calculate geometry with zero-length C-alpha bonds")

    left = xyz[:-2] - xyz[1:-1]
    right = xyz[2:] - xyz[1:-1]
    denominator = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
    if np.any(denominator <= 1e-12):
        raise ValueError("Cannot calculate virtual angles for degenerate coordinates")
    cosine = np.sum(left * right, axis=1) / denominator
    virtual_angles = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))

    centered = xyz - xyz.mean(axis=0)
    radius_gyration = float(np.sqrt(np.mean(np.sum(centered**2, axis=1))))
    energy = model.energy_components(xyz)

    return {
        "n_ca_residues": int(len(xyz)),
        "radius_of_gyration_angstrom": radius_gyration,
        "ca_bond_length_mean_angstrom": float(np.mean(bond_lengths)),
        "ca_bond_length_sd_angstrom": float(np.std(bond_lengths)),
        "ca_bond_length_min_angstrom": float(np.min(bond_lengths)),
        "ca_bond_length_max_angstrom": float(np.max(bond_lengths)),
        "ca_bond_length_max_abs_deviation_from_3_8_angstrom": float(
            np.max(np.abs(bond_lengths - EXPECTED_CA_BOND_ANGSTROM))
        ),
        "virtual_ca_angle_mean_degrees": float(np.mean(virtual_angles)),
        "virtual_ca_angle_sd_degrees": float(np.std(virtual_angles)),
        "virtual_ca_angle_min_degrees": float(np.min(virtual_angles)),
        "virtual_ca_angle_max_degrees": float(np.max(virtual_angles)),
        "rmsd_to_reference_angstrom": float(kabsch_rmsd(xyz, model.reference)),
        "native_contact_fraction": float(model.native_contact_fraction(xyz)),
        **energy,
    }
