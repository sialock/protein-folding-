# Development and release guide

This repository contains **Protein-MC-Go**, an installable Python package for reproducible, native-biased Cα Monte Carlo demonstrations.

## Local setup

```bash
git clone https://github.com/sialock/protein-folding-.git
cd protein-folding-
python -m pip install -e '.[dev]'
python -m pytest -q
python -m build
```

## Reproduce the bundled example

```bash
protein-mc run \
  --pdb examples/synthetic_helix_ca.pdb \
  --chain A --steps 600 --scramble-steps 120 \
  --seed 42 --report-every 10 --outdir results/demo
```

Run an independently seeded diagnostic sweep:

```bash
protein-mc benchmark \
  --pdb examples/synthetic_helix_ca.pdb \
  --chain A --steps 600 --n-replicates 5 \
  --outdir results/seed-sweep
```

The reference PDB is used to define native contacts. The results are demonstrations of a toy **native-biased** energy landscape, *not* de novo fold predictions or experimentally validated structures.

## Contribution workflow

1. Start from an up-to-date `main` and create a feature branch.
2. Make small, reviewed changes with deterministic fixtures; avoid including generated output in commits unless deliberately documenting an example.
3. Run the tests, package build, and CLI smoke test.
4. Open a pull request and verify GitHub Actions.
5. Keep methodology, limitations, and example provenance aligned with changes to the energy function or sampler.

## Release checklist

- Verify that CI passes on the proposed release commit.
- Update `CHANGELOG.md` and the package version.
- Build distributions with `python -m build`.
- Tag reviewed versions (for example `v0.1.0`) and publish release notes.
- Do not commit private data, credentials, or unpublished structures without permission.

## Repository organization

- `src/protein_mc/`: package implementation
- `tests/`: automated checks
- `examples/`: synthetic Cα input and diagnostics
- `assets/`: reproducible SVG plots
- `docs/`: scientific assumptions and historical migration notes

The earlier script is preserved in Git history. Donor-held-out single-cell benchmarking is developed in the **separate** [DonorShift-Bench repository](https://github.com/sialock/DonorShift-Bench).
