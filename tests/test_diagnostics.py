import numpy as np
import pytest

from protein_mc.diagnostics import structure_diagnostics
from protein_mc.energy import GoModel


def helix(n=20):
    t = np.deg2rad(np.arange(n) * 100)
    return np.column_stack([2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)])


def test_reference_geometry_and_energy_decomposition():
    reference = helix()
    model = GoModel(reference)
    metrics = structure_diagnostics(reference, model)
    components = model.energy_components(reference)

    assert metrics["n_ca_residues"] == len(reference)
    assert metrics["rmsd_to_reference_angstrom"] == pytest.approx(0.0, abs=1e-10)
    assert metrics["native_contact_fraction"] == pytest.approx(1.0)
    assert metrics["ca_bond_length_mean_angstrom"] == pytest.approx(3.8, abs=0.1)
    assert components["total_reduced"] == pytest.approx(model.energy(reference))
    assert components["repulsive_reduced"] >= 0.0
    assert components["attractive_reduced"] < 0.0
    assert metrics["radius_of_gyration_angstrom"] > 0.0


def test_diagnostics_reject_mismatched_or_nonfinite_coordinates():
    model = GoModel(helix())
    with pytest.raises(ValueError, match="match the model reference"):
        structure_diagnostics(np.zeros((5, 3)), model)
    invalid = helix()
    invalid[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        structure_diagnostics(invalid, model)


def test_energy_components_reject_invalid_coordinates():
    model = GoModel(helix())
    with pytest.raises(ValueError, match="Invalid coordinate"):
        model.energy_components(np.zeros((3, 3)))


def test_native_contact_fraction_rejects_nonfinite_inputs():
    model = GoModel(helix())
    invalid_coords = helix()
    invalid_coords[0, 0] = np.nan

    with pytest.raises(ValueError, match="Invalid coordinates or tolerance"):
        model.native_contact_fraction(invalid_coords)

    with pytest.raises(ValueError, match="Invalid coordinates or tolerance"):
        model.native_contact_fraction(helix(), tolerance_angstrom=np.nan)

    with pytest.raises(ValueError, match="Invalid coordinates or tolerance"):
        model.native_contact_fraction(helix(), tolerance_angstrom=np.inf)
