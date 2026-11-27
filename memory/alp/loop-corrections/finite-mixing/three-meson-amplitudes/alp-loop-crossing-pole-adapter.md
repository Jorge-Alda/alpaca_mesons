---
title: Dressed-source pole adapter for the three-meson crossing
status: mixed
created: 2026-11-23
updated: 2026-11-23
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Dressed-source pole adapter for the three-meson crossing

The [pole equations](../../../meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-exact-pole-equations.md) and [exact tree amplitude](../../../decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-finite-mixing-three-meson-tree.md) use the dressed-frame kinetic matrix $\mathbf1$. The [script](../../../../../scripts/alp_loop_crossing_poles.py) converts the dressed invariant matrices and topological coefficients into the full $9\times9$ meson squared-mass matrix and the mixing source $\langle\lambda_i\mathcal C\rangle$. It uses the nine generators with $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$. Inputs have GeV units and squared masses and pair invariants have GeV$^2$ units.

For a supplied $m_a^2$ and $\epsilon=f_\pi/f_a$, the fixed pole equation determines the ALP diagonal entry by the meson-block Schur complement. The script diagonalizes the complete real symmetric $10\times10$ mass matrix and assigns pole labels by maximum overlap with the unmixed basis. It computes the physical tree amplitude with the existing contact and exchange expressions. Daughter indices refer to the sorted unmixed meson eigenvectors. This assignment presupposes a continuous branch and isolated mass eigenstates; degeneracies and avoided crossings require additional care. The root-finding caller must verify the branch across its chosen interval.

Two symbolic tests verify the normalization of the invariant mixing source and the fixed-pole Schur equation. The numerical tree projection has not yet been independently validated for a complete physical channel.

```sh
PYTHONPATH=scripts .venv/bin/python -m unittest alp_loop_crossing_poles
```
