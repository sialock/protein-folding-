import json
from pathlib import Path

import pytest

from protein_mc.benchmark import run_seed_sweep
from protein_mc.cli import main
from protein_mc.simulation import RunConfig


PDB = Path(__file__).parents[1] / "examples" / "synthetic_helix_ca.pdb"


def test_seed_sweep_reproducible(tmp_path):
    outdir = tmp_path / "a"
    args = [
        "benchmark", "--pdb", str(PDB), "--n-replicates", "3",
        "--steps", "30", "--scramble-steps", "12", "--first-seed", "123",
        "--outdir", str(outdir),
    ]
    assert main(args) == 0
    data = json.loads((outdir / "benchmark_summary.json").read_text())
    assert data["n_replicates"] == 3
    assert 0 <= data["metrics"]["acceptance_rate"]["mean"] <= 1
    assert len((outdir / "replicates.csv").read_text().splitlines()) == 4
    old = (outdir / "benchmark_summary.json").read_bytes()
    assert main(args) == 0
    assert (outdir / "benchmark_summary.json").read_bytes() == old


def test_invalid_n_replicates(tmp_path):
    with pytest.raises(ValueError):
        run_seed_sweep(PDB, tmp_path, n_replicates=1, config=RunConfig(steps=2))