# Changelog

## 0.1.0 — 2026-10-31

Initial self-contained scientific rewrite of the older Monte Carlo example:

- Dimensionally consistent reduced-temperature Metropolis acceptance.
- Crankshaft geometry proposals; excluded volume and native-contact toy energy.
- Reproducible runs with independent NumPy random generators.
- Kabsch-aligned C-alpha RMSD and native-contact fraction.
- Single-run and multi-seed diagnostic CLI commands.
- PDB/CSV/JSON/SVG outputs; synthetic fixture, tests, and CI configuration.
- Explicit limitations: native-biased demonstration, **not de novo folding**.