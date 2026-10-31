"""Dimensionless, *native-biased* C-alpha Go-like toy energy function."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class EnergyParameters:
    # All energies are reduced/dimensionless; temperature is k_B*T in these units.
    contact_strength: float = 1.0
    contact_width_angstrom: float = 1.25
    native_contact_cutoff_angstrom: float = 8.0
    excluded_radius_angstrom: float = 4.0
    repulsion_strength: float = 5.0
    min_sequence_separation: int = 4

    def __post_init__(self) -> None:
        values = (
            self.contact_strength,
            self.contact_width_angstrom,
            self.native_contact_cutoff_angstrom,
            self.excluded_radius_angstrom,
            self.repulsion_strength,
        )
        if any(not np.isfinite(v) or v <= 0 for v in values):
            raise ValueError("Energy parameters must be finite and strictly positive")
        if self.min_sequence_separation < 3:
            raise ValueError("min_sequence_separation must be >= 3")


class GoModel:
    """Reference-contact attraction + short-range nonlocal repulsion.

    IMPORTANT: contacts are *derived from the supplied reference structure*.
    This is an educational native-biased model, NOT de novo folding or MD.
    """

    def __init__(self, reference: np.ndarray, params: EnergyParameters | None = None):
        coords = np.asarray(reference, dtype=float)
        if coords.ndim != 2 or coords.shape[1] != 3 or len(coords) < 5:
            raise ValueError("Expected reference coordinates with shape (N,3), N >= 5")
        if not np.isfinite(coords).all():
            raise ValueError("Reference coordinates must be finite")
        self.reference = coords.copy()
        self.reference.setflags(write=False)
        self.params = params or EnergyParameters()
        self.pair_i, self.pair_j = np.triu_indices(
            len(coords), k=self.params.min_sequence_separation
        )
        native_distances = np.linalg.norm(
            coords[self.pair_i] - coords[self.pair_j], axis=1
        )
        mask = native_distances < self.params.native_contact_cutoff_angstrom
        self.contact_i = self.pair_i[mask]
        self.contact_j = self.pair_j[mask]
        self.contact_distances = native_distances[mask]
        if len(self.contact_i) == 0:
            raise ValueError("No native contacts detected; provide a compact reference structure")

    @property
    def n_contacts(self) -> int:
        return len(self.contact_i)

    def energy(self, coords: np.ndarray) -> float:
        xyz = np.asarray(coords, dtype=float)
        if xyz.shape != self.reference.shape or not np.isfinite(xyz).all():
            raise ValueError("Invalid coordinate shape or non-finite values")
        d = np.linalg.norm(xyz[self.pair_i] - xyz[self.pair_j], axis=1)
        overlap = np.maximum(self.params.excluded_radius_angstrom - d, 0.0)
        repulsive = 0.5 * self.params.repulsion_strength * np.sum(overlap**2)
        contact_d = np.linalg.norm(xyz[self.contact_i] - xyz[self.contact_j], axis=1)
        delta = (contact_d - self.contact_distances) / self.params.contact_width_angstrom
        attractive = -self.params.contact_strength * np.sum(np.exp(-0.5 * delta**2))
        return float(repulsive + attractive)

    def native_contact_fraction(self, coords: np.ndarray, tolerance_angstrom: float = 1.5) -> float:
        """Fraction of native contacts within an absolute distance tolerance."""
        xyz = np.asarray(coords, dtype=float)
        if xyz.shape != self.reference.shape or tolerance_angstrom <= 0:
            raise ValueError("Invalid coordinates or tolerance")
        distances = np.linalg.norm(xyz[self.contact_i] - xyz[self.contact_j], axis=1)
        return float(np.mean(np.abs(distances - self.contact_distances) <= tolerance_angstrom))