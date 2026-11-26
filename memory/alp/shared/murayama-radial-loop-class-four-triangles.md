---
title: Radial and two meson triangle diagram source-trace products
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Radial and two-meson triangle diagram products

## Inputs and assumptions

Use the pole-leg, source and loop conventions of the [class 4 bubble record](murayama-radial-loop-class-four-bubbles.md), with $m_\sigma^2$ in loop-momentum products. The diagrams are 4.7–4.14 of the [estimate record](murayama-radial-up-only-estimate-result.md), following Alda:2026nrd and DiLuzio:2026poh. Each meson line is $\mathbf 1_{10}-v_av_a^T$ in the ten-field space of the nine meson fields and the bare ALP field: the completeness term, with $\sum_i\langle A\lambda_iB\lambda_i\rangle=2\langle A\rangle\langle B\rangle$ and $\sum_i\langle A\lambda_i\rangle\langle\lambda_iB\rangle=2\langle AB\rangle$ applied before any expansion, plus the bare angular term, minus the ALP pole term. Units are GeV. The retained list does not establish completeness of the graph selection for general sources.

## Method and result

The [script](../../../scripts/derive_radial_loop_class_four_triangles.py) builds the three vertices; both meson lines meet the first vertex, and the radial line joins the other two. Each meson line takes its three terms, giving nine combinations. When both lines take the completeness term the two contractions are applied one after the other; the second contraction acts on a trace that already contains the sum over orderings produced by the first, so nested sums over orderings and nested commutators appear. The output lists the eight products in compact traces, including the radial exchange of 4.13 and 4.14 as a shared factor. Symmetry and loop factors are excluded. No loop integral is evaluated.

## Reproduction and verification

Run `PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_class_four_triangles.py`. The test `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_class_four_triangles` expands every printed product back into trace words and requires exact equality with the independent direct expansion of the Lagrangian described in the [generator-sum record](murayama-loop-generator-trace-products.md), in which each meson line is a pair of placeholder matrices removed by the completeness rewrites, plus the bare angular term, minus the ALP pole term. Words are compared up to cyclic rotation only. All coefficients are exact rationals times symbols. The check does not test diagrams outside the retained list, nor the sizing prescription as a loop integral.
