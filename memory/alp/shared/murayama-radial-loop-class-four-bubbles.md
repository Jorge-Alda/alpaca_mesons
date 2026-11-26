---
title: Radial and meson bubble diagram source-trace products
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Radial and meson bubble diagram products

## Inputs and assumptions

Use the conventions and compact vertex rules of the [generator-sum record](murayama-loop-generator-trace-products.md), following Alda:2026nrd and DiLuzio:2026poh, and the retained graphs and dimensional estimate of the [estimate record](murayama-radial-up-only-estimate-result.md). Each meson line is $\mathbf 1_{10}-v_av_a^T$ in the ten-field space of the nine meson fields and the bare ALP field: the completeness term, with $\sum_i\langle A\lambda_iB\lambda_i\rangle=2\langle A\rangle\langle B\rangle$ and $\sum_i\langle A\lambda_i\rangle\langle\lambda_iB\rangle=2\langle AB\rangle$ applied before any expansion, plus the bare angular term, minus the ALP pole term. Units are GeV. Loop-momentum products use the radial scale $m_\sigma^2$. No loop integral is evaluated.

## Method and result

The [script](../../../scripts/derive_radial_loop_class_four_bubbles.py) prints diagrams 4.1–4.6. The two vertices share the radial line and the meson line; the completeness term merges their traces into single traces of sums over orderings, and the bare angular and ALP pole terms are products of the two vertices. The charged-pion radial exchange is a shared factor in 4.5 and 4.6. The output gives each product without loop measure or symmetry factor. General sources may allow further graphs.

## Reproduction and verification

Run `PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_class_four_bubbles.py`. The test `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_class_four_bubbles` expands every printed product back into trace words and requires exact equality with the independent direct expansion of the Lagrangian described in the [generator-sum record](murayama-loop-generator-trace-products.md), in which each meson line is a pair of placeholder matrices removed by the completeness rewrites, plus the bare angular term, minus the ALP pole term. Words are compared up to cyclic rotation only. All coefficients are exact rationals times symbols. The check does not test diagrams outside the retained list, nor the sizing prescription as a loop integral.
