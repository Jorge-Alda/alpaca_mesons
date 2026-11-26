---
title: Radial and ALP triangle diagram source-trace products
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Radial and ALP triangle diagram products

## Inputs and assumptions

Use the conventions and dimensional prescription of the [bubble record](murayama-radial-loop-class-five-bubbles.md), following Alda:2026nrd and DiLuzio:2026poh, and the retained diagram set of the [estimate record](murayama-radial-up-only-estimate-result.md); the result does not establish that no other diagrams occur for general sources. Each meson line is $\mathbf 1_{10}-v_av_a^T$ in the ten-field space of the nine meson fields and the bare ALP field: the completeness term, with $\sum_i\langle A\lambda_iB\lambda_i\rangle=2\langle A\rangle\langle B\rangle$ and $\sum_i\langle A\lambda_i\rangle\langle\lambda_iB\rangle=2\langle AB\rangle$ applied before any expansion, plus the bare angular term, minus the ALP pole term. Units are GeV. Loop-momentum products are $m_a^2$. No loop integral is evaluated.

## Method and result

The [script](../../../scripts/derive_radial_loop_class_five_triangles.py) prints diagrams 5.5–5.12 with the vertex rules and the compact completeness contraction of the [vertex script](../../../scripts/derive_radial_loop_trace_contractions.py), described in the [generator-sum record](murayama-loop-generator-trace-products.md). The four-meson kinetic part enters in its commutator form, and the contraction keeps the commutators and the sums over orderings of the two vertices joined by the meson line. For each diagram the output gives the factor shared by the three meson-line terms (the remaining radial–ALP vertex and, for 5.11 and 5.12, the radial exchange), then the merged completeness trace, the bare angular product and the subtracted ALP pole product. No matrix component enters.

## Reproduction and verification

Run `PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_class_five_triangles.py`. The test `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_class_five_triangles` expands every printed product back into trace words and requires exact equality with the independent direct expansion of the Lagrangian described in the [generator-sum record](murayama-loop-generator-trace-products.md), in which each meson line is a pair of placeholder matrices removed by the completeness rewrites, plus the bare angular term, minus the ALP pole term. Words are compared up to cyclic rotation only. All coefficients are exact rationals times symbols. The check does not test diagrams outside the retained list, nor the sizing prescription as a loop integral.
