"""Entry point for a reproducible, local, native-biased MC experiment."""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

from . import __version__
from .benchmark import run_seed_sweep
from .energy import EnergyParameters, GoModel
from .geometry import kabsch_rmsd
from .report import write_trace_svg
from .simulation import RunConfig, simulate
from .structure import load_ca, write_ca_pdb


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="protein-mc",
        description="Native-biased C-alpha Monte Carlo demonstration (NOT de novo folding)",
    )
    sub = p.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Run a seeded native-biased C-alpha simulation")
    run.add_argument("--pdb", required=True, type=Path, help="Local PDB file; no network required")
    run.add_argument("--chain", default=None, help="PDB chain ID (default: first suitable chain)")
    run.add_argument("--outdir", type=Path, default=Path("results/demo"))
    run.add_argument("--steps", type=int, default=1500)
    run.add_argument("--scramble-steps", type=int, default=300)
    run.add_argument("--temperature", type=float, default=1.0, help="Reduced units, not Kelvin")
    run.add_argument("--max-angle", type=float, default=40.0, help="Degrees per MC proposal")
    run.add_argument("--report-every", type=int, default=25)
    run.add_argument("--seed", type=int, default=42)
    run.add_argument("--contact-cutoff", type=float, default=8.0, help="Angstrom")
    bench = sub.add_parser("benchmark", help="Assess seed-to-seed variability (not accuracy)")
    bench.add_argument("--pdb", type=Path, required=True)
    bench.add_argument("--chain", default=None)
    bench.add_argument("--outdir", type=Path, default=Path("results/seed-sweep"))
    bench.add_argument("--n-replicates", type=int, default=5)
    bench.add_argument("--first-seed", type=int, default=42)
    bench.add_argument("--steps", type=int, default=1000)
    bench.add_argument("--scramble-steps", type=int, default=200)
    bench.add_argument("--temperature", type=float, default=1.0)
    bench.add_argument("--max-angle", type=float, default=40.0)
    bench.add_argument("--contact-cutoff", type=float, default=8.0)
    return p


def run(args: argparse.Namespace) -> None:
    protein = load_ca(args.pdb, chain_id=args.chain)
    parameters = EnergyParameters(native_contact_cutoff_angstrom=args.contact_cutoff)
    model = GoModel(protein.coords, params=parameters)
    config = RunConfig(
        steps=args.steps,
        scramble_steps=args.scramble_steps,
        temperature=args.temperature,
        max_angle_degrees=args.max_angle,
        report_every=args.report_every,
        seed=args.seed,
    )
    result = simulate(model, config)
    args.outdir.mkdir(parents=True, exist_ok=True)
    for filename, coords, label in [
        ("initial_ca.pdb", result.initial_coords, "scrambled initial configuration"),
        ("best_ca.pdb", result.best_coords, "lowest toy-energy configuration"),
        ("final_ca.pdb", result.final_coords, "final accepted MC configuration"),
    ]:
        write_ca_pdb(args.outdir / filename, protein, coords, label)

    with (args.outdir / "metrics.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(result.trace[0]))
        writer.writeheader()
        writer.writerows(result.trace)

    summary = {
        "tool": "protein-mc-go",
        "version": __version__,
        "interpretation": "Native-biased educational Go-like model, not de novo folding or MD",
        "pdb_sha256": hashlib.sha256(args.pdb.read_bytes()).hexdigest(),
        "chain": protein.chain_id,
        "n_ca_residues": len(protein.coords),
        "n_native_contacts": model.n_contacts,
        "temperature_units": "dimensionless reduced energy (k_B*T)",
        "coordinate_units": "Angstrom",
        "config": asdict(config),
        "energy_parameters": asdict(parameters),
        "acceptance_rate": result.acceptance_rate,
        "initial_rmsd_angstrom": kabsch_rmsd(result.initial_coords, model.reference),
        "best_rmsd_angstrom": kabsch_rmsd(result.best_coords, model.reference),
        "final_rmsd_angstrom": kabsch_rmsd(result.final_coords, model.reference),
        "best_energy_reduced": result.best_energy,
        "final_energy_reduced": model.energy(result.final_coords),
    }
    (args.outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_trace_svg(
        args.outdir / "energy_trace.svg", result.trace, "energy_reduced", "Toy energy", "Reduced energy"
    )
    write_trace_svg(
        args.outdir / "rmsd_trace.svg", result.trace, "rmsd_angstrom", "CA RMSD", "Angstrom"
    )
    print(f"Wrote {args.outdir} | {len(protein.coords)} CA | {model.n_contacts} contacts")
    print(f"Acceptance {result.acceptance_rate:.2%} | best CA RMSD {summary['best_rmsd_angstrom']:.2f} A")
    print("Reminder: reference structure is used to construct the energy landscape.")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            run(args)
        elif args.command == "benchmark":
            summary = run_seed_sweep(
                pdb=args.pdb,
                outdir=args.outdir,
                n_replicates=args.n_replicates,
                first_seed=args.first_seed,
                chain=args.chain,
                config=RunConfig(
                    steps=args.steps,
                    scramble_steps=args.scramble_steps,
                    temperature=args.temperature,
                    max_angle_degrees=args.max_angle,
                    seed=args.first_seed,
                ),
                params=EnergyParameters(native_contact_cutoff_angstrom=args.contact_cutoff),
            )
            print(f"Wrote {args.outdir} | {summary['n_replicates']} independent seeds")
            print("Note: native-biased toy model; this is variability analysis, not a folding benchmark")
    except (ValueError, OSError, StopIteration) as exc:
        print(f"protein-mc: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())