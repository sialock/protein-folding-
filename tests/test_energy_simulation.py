import numpy as np
import pytest

from protein_mc.energy import EnergyParameters, GoModel
from protein_mc.geometry import crankshaft_move
from protein_mc.simulation import RunConfig, simulate


def helix(n=20):
    t = np.deg2rad(np.arange(n) * 100)
    return np.column_stack([2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)])


def test_native_model_has_contacts_and_preferential_native_energy():
    ref = helix()
    model = GoModel(ref)
    perturbed = ref.copy()
    rng = np.random.default_rng(77)
    for _ in range(35):
        perturbed = crankshaft_move(perturbed, rng, max_angle_radians=2.0)
    assert model.n_contacts > 0
    assert model.native_contact_fraction(ref) == 1.0
    assert model.energy(ref) < model.energy(perturbed)


def test_reproducible_and_reference_not_mutated():
    ref = helix(18)
    original = ref.copy()
    model = GoModel(ref)
    config = RunConfig(steps=75, scramble_steps=30, seed=123, report_every=10)
    first, second = simulate(model, config), simulate(model, config)
    assert np.array_equal(first.final_coords, second.final_coords)
    assert first.trace == second.trace
    assert 0 <= first.acceptance_rate <= 1
    assert len(first.trace) == 1 + 8  # 10..70 plus 75
    assert np.array_equal(ref, original)
    assert np.array_equal(model.reference, original)
    assert first.best_energy <= first.trace[0]["energy_reduced"]


def test_parameter_validation():
    with pytest.raises(ValueError):
        EnergyParameters(contact_width_angstrom=0)
    with pytest.raises(ValueError):
        RunConfig(temperature=0)
    with pytest.raises(ValueError):
        GoModel(np.zeros((4, 3)))