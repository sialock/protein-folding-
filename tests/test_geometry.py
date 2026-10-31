import numpy as np
import pytest

from protein_mc.geometry import crankshaft_move, kabsch_rmsd


def test_kabsch_ignores_translation_and_rotation():
    rng = np.random.default_rng(7)
    coords = rng.normal(size=(15, 3))
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    moved = coords @ rotation + np.array([10, 5, -3])
    assert kabsch_rmsd(moved, coords) < 1e-10


def test_crankshaft_preserves_adjacent_ca_distances_and_input():
    coords = np.column_stack([np.cos(np.arange(20)), np.sin(np.arange(20)), np.arange(20)])
    original = coords.copy()
    proposal = crankshaft_move(coords, np.random.default_rng(3), max_angle_radians=1.2)
    assert np.allclose(np.linalg.norm(np.diff(coords, axis=0), axis=1),
                       np.linalg.norm(np.diff(proposal, axis=0), axis=1), atol=1e-10)
    assert np.array_equal(coords, original)
    assert not np.array_equal(coords, proposal)


def test_geometry_rejects_invalid_arrays():
    with pytest.raises(ValueError):
        kabsch_rmsd(np.ones((4, 3)), np.ones((5, 3)))
    with pytest.raises(ValueError):
        crankshaft_move(np.ones((3, 3)), np.random.default_rng(1))