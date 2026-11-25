---
title: Symbolic check of the projected cubic and quintic vertices
status: derived
created: 2026-11-23
updated: 2026-11-23
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Projected cubic and quintic potential vertices

The [general dressed vertex](alp-dressed-vertex-power-counting-invariants.md#result-1-vertices-in-the-invariants-new-form-of-l1-verified) and [pole projection](alp-finite-mixing-three-meson-tree.md#result) provide the input. Units are GeV. The calculation uses diagonal symbolic entries for $\mathbb M_{m_a^2}^{-1}(\mathcal C)$ and $\widehat P_2$, a real symmetric $d$-$s$ source $\widehat S_1$, and symbolic $B_0$, $f_\pi$, $m_0^2$ and $\widehat\Theta_2$. These restrictions describe the symbolic check. The analytical formulas being checked remain valid for general noncommuting matrices.

The [script](../../../scripts/derive_alp_loop_cubic_quintic_projection.py) contracts the cubic and quintic derivatives with one or two ALP pole directions, leaving the common pole normalization outside. It compares each result with the corresponding symmetrized invariant trace, then forms both cubic–quintic products. The test input consists of symbolic source parameters and fixed generator directions. No physical masses or couplings are assigned values. This check does not establish the sum of all on-shell terms.

```sh
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_loop_cubic_quintic_projection
.venv/bin/python scripts/derive_alp_loop_cubic_quintic_projection.py
```
