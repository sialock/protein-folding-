# Protein-MC-Go

Protein-MC-Go is a small Python project for trying out Monte Carlo moves on a protein's Cα chain. It reads a PDB file, changes the coordinates with crankshaft moves, and records how a simple energy function changes during the run.

**Please note:** this is not a protein-folding predictor or a molecular-dynamics package. The input structure is used to build the model's contact map, so the program cannot predict a structure from sequence. It is intended for learning and for exploring a basic, reference-based Monte Carlo model.

## Install

Clone the repository and install the package:

\`\`\`bash
git clone https://github.com/sialock/protein-folding-.git
cd protein-folding-

python -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate
python -m pip install -e .
\`\`\`

To run the tests, install the development dependencies as well:

\`\`\`bash
python -m pip install -e ".[dev]"
python -m pytest -q
\`\`\`

## Run an example

The repository includes a small synthetic helix. Use it to check that the installation works:

\`\`\`bash
protein-mc run \
  --pdb examples/synthetic_helix_ca.pdb \
  --chain A \
  --steps 1500 \
  --scramble-steps 300 \
  --seed 42 \
  --outdir results/helix
\`\`\`

You can also run the command through Python if the \`protein-mc\` entry point is not available:

\`\`\`bash
python -m protein_mc run --pdb examples/synthetic_helix_ca.pdb --chain A --outdir results/helix
\`\`\`

To use another structure, replace the example PDB path and set the chain you want to analyse.

## Input requirements

The input must be a local PDB file containing one standard protein chain with Cα atoms. The current reader supports chains with 5 to 300 residues and consecutive residue numbering. Structures with missing residues, numbering gaps, or insertion codes are not supported yet.

A Cα-only PDB can be used as input. The structures written by the program are also Cα-only; they are not complete all-atom protein models.

## Files written by a run

The output directory contains:

| File | What it contains |
| --- | --- |
| \`initial_ca.pdb\` | The structure after the initial scrambling moves |
| \`best_ca.pdb\` | The sampled structure with the lowest value of the model's energy function |
| \`final_ca.pdb\` | The structure at the end of the run |
| \`metrics.csv\` | Energy, RMSD, native-contact fraction and acceptance rate over the run |
| \`structure_diagnostics.csv\` | Geometry and energy measurements for the initial, best and final structures |
| \`summary.json\` | Run settings, input-file hash and summary measurements |
| \`energy_trace.svg\` | Energy during the run |
| \`rmsd_trace.svg\` | Cα RMSD to the input structure during the run |

The lowest-energy structure is not necessarily the one with the lowest RMSD. The program selects it using the toy energy function, not by choosing the structure closest to the reference.

## Run more than one seed

The \`benchmark\` command repeats the run with different random seeds. It writes \`replicates.csv\` and \`benchmark_summary.json\`, including the mean, standard deviation, minimum and maximum for each reported measurement.

\`\`\`bash
protein-mc benchmark \
  --pdb examples/synthetic_helix_ca.pdb \
  --chain A \
  --steps 600 \
  --n-replicates 5 \
  --outdir results/seed-sweep
\`\`\`

This is useful for seeing how much the results vary between runs. It does not show that the simulation has converged or that the program can predict an unknown protein structure.

## How the model works

The program keeps only the Cα coordinates. It builds a list of contacts from the supplied PDB, then uses a simple Go-like energy function with two parts: an attraction between those reference contacts and a penalty for nonlocal Cα pairs that get too close. Crankshaft moves change part of the chain while preserving the distances between neighbouring Cα atoms. The Metropolis rule decides whether each proposed move is accepted.

The energy and temperature are in reduced, dimensionless units. The temperature is not in kelvin, and Monte Carlo steps do not represent physical time.

## Limitations

This is a coarse-grained demonstration, not a validated physical model. It does not include side chains, solvent, electrostatics, hydrogen bonds or a complete description of backbone geometry. The geometry measurements in the output are useful checks, but they do not establish that a structure is physically realistic.

The contact map comes from the input structure itself. For that reason, results must not be presented as de novo folding predictions or as evidence that the method can recover a native structure from sequence.

More detail is available in [Scientific notes](docs/SCIENTIFIC_NOTES.md). The changes from the earlier script are described in [Migration notes](docs/MIGRATION.md).

## Development

Install the development dependencies, then run the tests and build the package:

\`\`\`bash
python -m pip install -e ".[dev]"
python -m pytest -q
python -m compileall -q src
python -m build
\`\`\`

GitHub Actions runs the test suite on Python 3.10 and 3.12.

## License

MIT. See [LICENSE](LICENSE).
