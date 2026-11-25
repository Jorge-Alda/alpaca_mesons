---
title: Dependence of the radial loop-size weights on the field choice
status: derived
created: 2026-11-25
updated: 2026-11-25
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Dependence of the radial loop-size weights on the field choice

This record is task 6.1 of the OpenSpec change `phi-dependent-loop-tree-crossing-murayama`. It cross-checks the magnitude weights of the [size prescription](murayama-radial-loop-size-prescription.md) in polar and Cartesian fields, using exact projectors in the up-only model. It involves no loop: the check acts on the weights and on tree-level vertices.

## Result

- The normalization factors and the powers of $1/f_a$ in each weight are the same in the polar and the Cartesian fields.
- The vertex contractions that enter the weights are **not** the same:
  - the $\sigma PP$ vertex differs by a term proportional to $m_\sigma^2-q^2$, so the two fields agree only on shell;
  - the $\sigma\sigma PP$ vertex is $-2p_1\cdot p_2$ in polar fields and $-(m_\sigma^2-B_0m_u\epsilon^2)$ in Cartesian fields;
  - the four-pole contact vertex differs by a term that depends on $\lambda$.
- As a control, the on-shell four-pole amplitude is the same in both fields.

Therefore the magnitude weights depend on the field choice. The loop-to-tree estimate is a statement about magnitudes of a chosen parametrization, and its numerical value is parametrization dependent away from the mass shell. This is the reason for the polar choice of the [crossing record](murayama-radial-crossing-implicit.md), where the $\sigma$ line and its absence are separated without mixing with $a$.

## Limits

- Checked for the up-only model only.
- The contact term $\mathcal V^{(4)}_{a123}$ is a symbol in the script. It is zero for three distinct neutral daughters.
- Class 5 of the ALP-only record is not identified separately in this check.
- The check does not compute the loop. It shows that the weights of the same estimate change with the field choice, so the field choice must be stated wherever the estimate is quoted.

## Script

[Cartesian weights check](../../../scripts/derive_murayama_radial_cartesian_weights.py) and its [tests](../../../scripts/test_derive_murayama_radial_cartesian_weights.py).
