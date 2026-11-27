---
title: ALP and meson bubble diagram source-trace products
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# ALP and meson bubble diagram products

## Inputs and assumptions

Use the pole-leg, incoming-momentum and source conventions and the compact vertex rules of the [generator-sum record](murayama-loop-generator-trace-products.md), and the retained graph set of the [estimate record](../../../couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md), following Alda:2026nrd and DiLuzio:2026poh. Each meson line is $\mathbf 1_{10}-v_av_a^T$ in the ten-field space of the nine meson fields and the bare ALP field: the completeness term, with $\sum_i\langle A\lambda_iB\lambda_i\rangle=2\langle A\rangle\langle B\rangle$ and $\sum_i\langle A\lambda_i\rangle\langle\lambda_iB\rangle=2\langle AB\rangle$ applied before any expansion, plus the bare angular term, minus the ALP pole term. All dimensionful inputs are in GeV. Loop-momentum products use $m_a^2$. No loop integral is evaluated.

## Method and result

The [script](../../../../../../scripts/derive_radial_loop_class_two_bubbles.py) prints diagrams 2.1–2.5. Both vertices carry the internal ALP leg and one end of the meson line. In 2.1–2.3 both are four-point vertices, so their four-meson kinetic parts give commutators; after completeness these appear as nested commutators, for example $\langle[\lambda_a,[\lambda_a,\lambda_{\pi^+}]]\,[\lambda_{\pi^-},[\lambda_a,\lambda_{\pi^0}]]\rangle$, multiplied by products of $(k_x\cdot k_y-m_a^2)$. The radial exchanges of 2.4 and 2.5 are shared factors. Symmetry and loop factors are excluded. The retained diagrams need not exhaust the contributions for general sources.

## Reproduction and verification

Run `PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_class_two_bubbles.py`. The test `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_class_two_bubbles` expands every printed product back into trace words and requires exact equality with the independent direct expansion of the Lagrangian described in the [generator-sum record](murayama-loop-generator-trace-products.md), in which each meson line is a pair of placeholder matrices removed by the completeness rewrites, plus the bare angular term, minus the ALP pole term. Words are compared up to cyclic rotation only. All coefficients are exact rationals times symbols. The check does not test diagrams outside the retained list, nor the sizing prescription as a loop integral.
