---
title: Two-radial-line diagrams with internal mesons
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Two-radial-line diagrams with internal mesons

## Inputs and assumptions

Use the polar source and pole conventions and the compact vertex rules of the [generator-sum record](murayama-loop-generator-trace-products.md), following Alda:2026nrd and DiLuzio:2026poh. Each meson line is $\mathbf 1_{10}-v_av_a^T$ in the ten-field space of the nine meson fields and the bare ALP field: the completeness term, with $\sum_i\langle A\lambda_iB\lambda_i\rangle=2\langle A\rangle\langle B\rangle$ and $\sum_i\langle A\lambda_i\rangle\langle\lambda_iB\rangle=2\langle AB\rangle$ applied before any expansion, plus the bare angular term, minus the ALP pole term. Units are GeV. The radial cubic vertex is $V_{\sigma^3}=-6\lambda f_a+6B_0f_\pi^2\langle\widehat S_{30}\rangle/f_a^3$, the quartic radial potential giving the first term and the general source the second. The retained graphs are those of the [estimate record](murayama-radial-up-only-estimate-result.md); further graphs for general sources are not excluded.

## Method and result

The [script](../../../scripts/derive_radial_loop_class_six_meson.py) prints diagrams 6.2–6.6. In each, the vertices joined only by radial lines stay separate factors, and each meson line takes its three terms. The radial exchanges of 6.3 and 6.4, and $V_{\sigma^3}$, are shared factors. In the boxes 6.5 and 6.6 the two meson lines join different pairs of vertices, so each completeness term merges only the two vertices of its own line. Loop-momentum products use $m_\sigma^2$. Symmetry and loop factors are excluded; no loop integral is evaluated.

## Reproduction and verification

Run `PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_class_six_meson.py`. The test `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_class_six_meson` expands every printed product back into trace words and requires exact equality with the independent direct expansion of the Lagrangian described in the [generator-sum record](murayama-loop-generator-trace-products.md), in which each meson line is a pair of placeholder matrices removed by the completeness rewrites, plus the bare angular term, minus the ALP pole term. Words are compared up to cyclic rotation only. All coefficients are exact rationals times symbols. The checks keep $\langle\widehat S_{30}\rangle$ general. The check does not test diagrams outside the retained list, nor the sizing prescription as a loop integral.
