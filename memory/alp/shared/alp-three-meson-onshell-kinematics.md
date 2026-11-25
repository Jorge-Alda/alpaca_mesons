---
title: On-shell kinematic variation for three-meson amplitudes
status: derived
created: 2026-11-23
updated: 2026-11-23
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# On-shell kinematic variation

For the [three-meson amplitude](alp-finite-mixing-three-meson-tree.md#result), use all incoming momenta and choose the three daughter momenta as a basis. Their masses are $m_{n_1},m_{n_2},m_{n_3}$ and the ALP pole mass is $m_a$. Let $s_{12}=(p_{n_1}+p_{n_2})^2$ and $s_{23}=(p_{n_2}+p_{n_3})^2$. These squared masses and invariants are in GeV$^2$.

Momentum conservation and the on-shell conditions give
$$
(p_{n_1}+p_{n_3})^2
=m_a^2+m_{n_1}^2+m_{n_2}^2+m_{n_3}^2-s_{12}-s_{23}.
$$
The dot products follow from the three pair invariants. In particular,
$$
p_{n_1}\cdot p_{n_3}
=\frac12(m_a^2+m_{n_2}^2-s_{12}-s_{23}).
$$
For the one-loop mass shifts in the [on-shell assembly](alp-three-meson-one-loop-estimate.md#assembly-at-one-loop), hold $s_{12}$ and $s_{23}$ fixed. Vary each daughter squared momentum and each pair product by differentiating these expressions with respect to the four squared masses. The ALP momentum is minus the sum of the three daughter momenta.

The [SymPy calculation](../../../scripts/derive_alp_three_meson_onshell_kinematics.py) constructs the full Gram matrix and its linear variation. Its symbolic unit tests verify both fixed pair invariants and momentum conservation before and after the variation.

```sh
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_three_meson_onshell_kinematics
```
