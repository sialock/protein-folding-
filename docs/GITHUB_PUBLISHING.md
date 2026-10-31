# Publishing this rewrite to the existing repository

This release is meant to **upgrade** `sialock/protein-folding-`. It does not
archive, delete, force-push or overwrite Git history. The original root file
named `protein folding` can be retained as a historical example.

## Local Git workflow (recommended)

1. Download and extract the `protein-mc-go-v0.1.0.zip` release archive.
2. Clone the existing project and move into it:

   ```bash
   git clone https://github.com/sialock/protein-folding-.git
   cd protein-folding-
   ```
3. Copy **the contents** of the extracted `protein-mc-go/` directory over
   the cloned project, including `.github/` and `.gitignore`. **Do not** copy
   the extracted directory itself as a nested subdirectory. Keep the cloned
   repository's `.git/` folder and existing `protein folding` file.
4. Run checks:

   ```bash
   python -m pip install -e ".[dev]"
   python -m pytest -q
   protein-mc run --pdb examples/synthetic_helix_ca.pdb --steps 100 \
      --scramble-steps 25 --outdir results/smoke
   ```
5. Publish after reviewing the exact diff:

   ```bash
   git status
   git diff --stat
   git add README.md pyproject.toml LICENSE .gitignore src tests examples assets docs .github
   git commit -m "feat: reproducible native-biased protein Monte Carlo toolkit"
   git push origin main
   ```

If the default branch has protection or permissions, create a feature branch
locally (`git switch -c feat/reproducible-protein-mc`), push that branch, and
open a pull request instead. **Do not force push**.

## Why this is a different scientific claim

The original script attempted a simplistic atomistic protein folding simulation.
This replacement explicitly uses a **known reference structure** to generate
its energy function and is only a native-biased, C-alpha-only teaching model.
Its outputs must not be described as novel protein predictions or experimental
folding results.