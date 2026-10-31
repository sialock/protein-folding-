"""Seed-sweep diagnostic; NOT a claim of physical convergence or predictive power."""
from __future__ import annotations

import csv
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path

import numpy as np

from .energy import EnergyParameters, GoModel
from .geometry import kabsch_rmsd
from .simulation import RunConfig, simulate
from .structure import load_ca


def run_seed_sweep(
    pdb: Path,
    outdir: Path,
    n_replicates: int = 5,
    first_seed: int = 42,
    chain: str | None = None,
    config: RunConfig | None = None,
    params: EnergyParameters | None = None,
) -> dict:
    """Repeat independent Monte Carlo runs and report between-seed variability.

    The native reference itself is provided to the energy function, so this is
    strictly a model-behavior diagnostic rather than a folding benchmark.
    """
    if not 2 <= n_replicates <= 100:
        raise ValueError("n_replicates must be from 2 to 100")
    if first_seed < 0:
        raise ValueError("first_seed must be nonnegative")
    protein = load_ca(pdb, chain_id=chain)
    model = GoModel(protein.coords, params=params)
    config = config or RunConfig()
    outdir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for seed in range(first_seed, first_seed + n_replicates):
        result = simulate(model, replace(config, seed=seed))
        rows.append({
            "seed": seed,
            "acceptance_rate": result.acceptance_rate,
            "initial_rmsd_angstrom": kabsch_rmsd(result.initial_coords, model.reference),
            "best_rmsd_angstrom": kabsch_rmsd(result.best_coords, model.reference),
            "final_rmsd_angstrom": kabsch_rmsd(result.final_coords, model.reference),
            "best_energy_reduced": result.best_energy,
            "final_energy_reduced": model.energy(result.final_coords),
            "best_native_contact_fraction": model.native_contact_fraction(result.best_coords),
        })
    with (outdir / "replicates.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metric_names = [name for name in rows[0] if name != "seed"]
    aggregated = {
        name: {
            "mean": float(np.mean([r[name] for r in rows])),
            "sample_sd": float(np.std([r[name] for r in rows], ddof=1)),
            "min": float(np.min([r[name] for r in rows])),
            "max": float(np.max([r[name] for r in rows])),
        }
        for name in metric_names
    }
    summary = {
        "interpretation": "Independent seed-sweep diagnostic for a native-biased toy model; not a predictive folding benchmark or convergence proof",
        "pdb_sha256": hashlib.sha256(pdb.read_bytes()).hexdigest(),
        "reference_chain": protein.chain_id,
        "n_ca_residues": len(protein.coords),
        "n_contacts": model.n_contacts,
        "n_replicates": n_replicates,
        "first_seed": first_seed,
        "run_config": asdict(config),
        "energy_parameters": asdict(model.params),
        "metrics": aggregated,
    }
    (outdir / "benchmark_summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return summary