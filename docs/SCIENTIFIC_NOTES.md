# Scientific interpretation and limitations

## What is being modeled?

A **single continuous chain of C-alpha coordinates**, with a native-contact
potential defined using an already known reference PDB. The system is a
coarse-grained *teaching demonstration*, not an all-atom protein simulation.

1. Extract C-alpha points from one chain in the first model of a local PDB.
2. Build reference-native pairs with sequence separation >= 4 and reference
   distance < 8 Å (default).
3. Create a disturbed initial conformation via crankshaft rotations.
4. Run stochastic local crankshaft proposals with Metropolis acceptance.
5. Report sampled-energy changes, Kabsch-aligned C-alpha RMSD to the *known*
   reference, and reference-native-contact fraction.

## Energy and temperature

The energy function has an attractive Gaussian native-contact term and a
non-local repulsive quadratic overlap term. It omits side chains, solvation,
electrostatics, hydrogen bonding, disulfide bonds, torsional chemistry, bond
angle penalties, and explicit chirality constraints.

The energy parameters and temperature are **dimensionless reduced units**;
we NEVER multiply this energy by the physical Boltzmann constant in J/K.
Thus the Metropolis acceptance exponent is dimensionless by construction:
`exp(-delta_energy / reduced_temperature)`.

## What crankshaft moves preserve

Rotating internal C-alpha beads around the axis connecting two fixed bead
anchors preserves lengths between *consecutive* beads, to floating-point
precision. It does not guarantee physically realistic local backbone
geometry or freedom from nonlocal clashes; the repulsion term merely penalizes
the latter.

## Scientific limitations

- **Native-structure leakage by design.** Contact positions come from the
  reference PDB. This is not a sequence-to-structure predictor and must never
  be benchmarked as such.
- A lower toy energy does not prove a better experimental structure. The
  lowest-energy sample can have a larger RMSD than another sample.
- No solvent, protonation, atomistic dynamics or physical simulation time.
  Monte Carlo "steps" are *not* nanoseconds.
- Sampling can be slow to mix, strongly depends on move size, and may not
  approach equilibrium. One seed is not a convergence study.
- The energy function is illustrative, and its parameter values are **not
  calibrated** to experimental thermodynamics.
- The demo requires consecutive residue numbering and at least one native
  contact, and restricts input to <=300 CA sites for clarity and tractability.
- The output PDB is CA-only; it does not contain all atom coordinates or
  a properly rebuilt protein backbone.

## Required checks before a scientific study

Use multiple random seeds, independent chains and convergence diagnostics.
Inspect backbone geometry, structural ensembles and contact maps, document
parameter sensitivity, and validate against a credible physics-based or
experimentally benchmarked reference method.

## Useful primary references

- Metropolis N. et al. (1953). Equation of State Calculations by Fast
  Computing Machines. *J. Chem. Phys.* 21, 1087. DOI: 10.1063/1.1699114.
- Clementi C., Nymeyer H., Onuchic J.N. (2000). Topological and energetic
  factors: what determines the structural details of the transition state
  ensemble and "en-route" intermediates for protein folding? *J. Mol. Biol.*
  298, 937-953. DOI: 10.1006/jmbi.2000.3693.