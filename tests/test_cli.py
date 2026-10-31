import json
from pathlib import Path

from protein_mc.cli import main
from protein_mc.structure import load_ca


EXAMPLE = Path(__file__).parents[1] / "examples" / "synthetic_helix_ca.pdb"


def test_example_loads_and_cli_writes_reproducible_artifacts(tmp_path):
    protein = load_ca(EXAMPLE, chain_id="A")
    assert len(protein.coords) == 24
    output = tmp_path / "example"
    args = [
        "run", "--pdb", str(EXAMPLE), "--chain", "A", "--steps", "45",
        "--scramble-steps", "15", "--report-every", "7", "--seed", "5",
        "--outdir", str(output),
    ]
    assert main(args) == 0
    for filename in [
        "initial_ca.pdb", "best_ca.pdb", "final_ca.pdb", "metrics.csv",
        "summary.json", "energy_trace.svg", "rmsd_trace.svg",
    ]:
        assert (output / filename).is_file(), filename
    assert len(load_ca(output / "final_ca.pdb").coords) == 24
    info = json.loads((output / "summary.json").read_text())
    assert info["config"]["seed"] == 5
    assert info["n_native_contacts"] > 0
    assert "not de novo" in info["interpretation"]
    contents = (output / "summary.json").read_bytes()
    assert main(args) == 0
    assert (output / "summary.json").read_bytes() == contents


def test_missing_file_returns_failure(tmp_path):
    assert main(["run", "--pdb", str(tmp_path / "missing.pdb")]) == 2