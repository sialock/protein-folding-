# Legacy migration and scientific corrections

The former `protein-folding-` repository contained a small exploratory Monte Carlo script. Protein-MC-Go v0.1 replaces it with a testable, installable, explicitly **native-biased Cα** demonstration. The original source remains available through Git history; no rewriting of historic commits is required.

A separate, related legacy experiment also exists at [protein-folding-using-Monte-Carlo-and-analyzing-the-structure](https://github.com/sialock/protein-folding-using-Monte-Carlo-and-analyzing-the-structure). It is not a dependency and none of its code is mixed into this repository.

## Main corrections

| Legacy approach | Current approach |
| --- | --- |
| `k_B = 1.380649e-23` in J/K applied to arbitrary numerical energies | Consistent dimensionless toy energy with reduced temperature |
| Lennard-Jones potential without a calibrated atomistic force field | Explicit native-contact attraction and repulsive toy potential |
| Pairwise divisions including zero-distance diagonal | Unique nonlocal pairs; no self-interaction |
| Independent atom displacements that can destroy chain continuity | Cα crankshaft proposals preserving neighbor distances |
| Shared mutable structure during threaded simulations | Separate seeded runs using immutable reference coordinates |
| Hard-coded remote PDB and a single untitled code file | Local synthetic example, CLI, packaging, documented tests |
| “Folded” outputs without scientific qualification | Explicit initial/best/final Cα models and toy-energy diagnostics |

## Interpretation limits

The native reference *defines* the energetic minimum. Therefore an apparent decrease in RMSD is not evidence of sequence-to-structure prediction. This code is not a force-field-validated molecular dynamics implementation or an independent enzyme-structure prediction method. See [scientific notes](SCIENTIFIC_NOTES.md).

## Reproducibility

Install and test the current package using the commands in [DEVELOPMENT.md](DEVELOPMENT.md). Earlier versions are inspectable via Git history; avoid reproducing their physical assumptions as established folding science.
