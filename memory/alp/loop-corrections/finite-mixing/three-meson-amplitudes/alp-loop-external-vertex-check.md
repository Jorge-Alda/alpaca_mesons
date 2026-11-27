---
title: Symbolic checks of external-pole loop vertices
status: derived
created: 2026-11-23
updated: 2026-11-23
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Symbolic checks of external-pole loop vertices

The input Lagrangian and conventions are those of the [dressed-frame vertex counting](../../../decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-dressed-vertex-power-counting-invariants.md). Masses and momenta use GeV. The calculation uses the unnormalized ALP-like pole direction, whose meson and ALP components are respectively $2\epsilon\mathbb M_{m_a^2}^{-1}(\mathcal C)$ and $1$. Divide by $\sqrt{N(m_a^2)}$ for each such direction.

For this symbolic check, $\widehat S_0$ is diagonal, $\widehat S_1$ and $\widehat P_1$ have symbolic real $d$–$s$ entries, and $\mathbb M_{m_a^2}^{-1}(\mathcal C)$ is a real symmetric matrix with diagonal and $d$–$s$ entries. $\widehat S_2$ and $\widehat P_2$ are symbolic multiples of the identity. These restrictions apply to the check. The general expressions are given in Result 3 of the input record.

The script contracts the existing potential vertex with one or two ALP-like pole directions and independently constructs the invariant expressions in Result 3. It verifies the one-ALP cubic vertex, the closed-line cubic vertex, and the potential and kinetic parts of the closed-line quartic vertex. Two distinct symbolic meson matrices are used. The kinetic comparison uses independent formal external and loop momenta.

All four residuals vanish exactly. The one-ALP cubic vertices multiply to give the two-meson bubble numerator. The one- and two-ALP cubic vertices multiply to give the mixed bubble numerator. The closed-line cubic and quartic vertices supply the vacuum and external self-energy tadpole numerators, respectively. Scalar loop functions and their overall self-energy sign require the chosen propagator convention.

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_alp_loop_external_contractions.py
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_loop_external_contractions
```
