---
title: On-shell assembly from one-loop vertex functions
status: derived
created: 2026-11-23
updated: 2026-11-23
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# On-shell assembly from one-loop vertex functions

The [tree amplitude](../../../decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-finite-mixing-three-meson-tree.md#result) consists of a four-point vertex and three cubic-exchange channels. The [loop estimate](alp-three-meson-one-loop-estimate.md#assembly-at-one-loop) gives the corresponding first variation. Use all incoming momenta and the [fixed-invariant Gram matrix](alp-three-meson-onshell-kinematics.md). Squared masses and invariants are in GeV$^2$.

For a complete orthonormal pole basis of the canonical LO mass matrix, the [executable assembly](../../../../../scripts/derive_alp_three_meson_onshell_assembly.py) performs the following operations.

1. Apply the inverse tree propagator at zero momentum to the one-point function and vary each vertex by the resulting stationary-point displacement.
2. Subtract the ALP bare-mass correction fixed by the ALP pole condition, and use the two-point function and its derivative to calculate the external mass and residue shifts. The sum over intermediate pole vectors includes the full spectrum.
3. Add the four-point correction and all three variations of each exchange channel, including the correction to its internal inverse propagator.
4. Vary each external pole vector and the daughter Gram matrix at fixed $s_{12}$ and $s_{23}$ in the full connected tree amplitude.

The program accepts the one-loop one-, two-, three- and four-point functions as callbacks in Minkowski conventions. The two-point callback is a symbolic function of the squared momentum and must retain its scalar-moment dependence under differentiation. The code assembles these inputs and does not itself derive their invariant source contractions.

Unit tests keep all physical parameters symbolic. They check a five-state basis with four selected external poles, the fixed ALP pole, residue and on-shell mass variations, and the cubic-exchange corrections in all three channels.

The external-pole variation omits an exactly degenerate intermediate pole when its projected off-diagonal two-point correction vanishes. A nonzero projection raises an error because the correction must first be diagonalized within the degenerate subspace. A symbolic test checks both cases. The code uses the supplied pole vectors as columns in the canonical field basis.

```sh
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_three_meson_onshell_assembly
```
