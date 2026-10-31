# Migration / original code audit

Original repositories (retained; nothing here deletes them):

- https://github.com/sialock/protein-folding-
- https://github.com/sialock/protein-folding-using-Monte-Carlo-and-analyzing-the-structure

## Critical corrections

| Original behavior | Replacement |
| --- | --- |
| `k_B = 1.380649e-23` (J/K) applied to an arbitrary, uncalibrated energy | Explicit dimensionless toy-energy units with dimensionless temperature |
| Lennard-Jones pair potential used without a calibrated full force field | Explicit educational native-biased contact potential; limited claim |
| Pairwise division by zero on diagonal followed by filtering | Nonlocal unique pair indices, no diagonal pairs |
| Random independent atom displacements break covalent chain geometry | C-alpha crankshaft proposals preserve adjacent distances |
| Multiple concurrent simulations share mutable `Bio.PDB` structure | One immutable reference and one independently seeded RNG per run |
| Bare unversioned code with hard-coded downloaded PDB | Installable `src/` package, local synthetic fixture, CLI and tests |
| Files called `folded_structure*.pdb` | C-alpha-only `initial_ca`, `best_ca`, `final_ca` with limitations |

## Practical upload to existing GitHub repository

1. Extract the release archive into a fresh clone of `sialock/protein-folding-`.
2. Preserve its pre-existing `protein folding` file if you want the legacy script.
3. Run `python -m pip install -e ".[dev]"` and `python -m pytest -q`.
4. Commit the new files and publish. See the README's instructions.
5. Only archive the second duplicate repo after reviewing its uniqueness.

The original scripts are preserved online and should not be represented as a
validated folding predictor. The new model remains intentionally educational.