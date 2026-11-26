---
title: ALP tadpole diagram source-trace products
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# ALP tadpole diagram products

## Inputs and assumptions

Use the source, pole, momentum and loop-size conventions of the [estimate record](murayama-radial-up-only-estimate-result.md) and the compact vertex rules of the [generator-sum record](murayama-loop-generator-trace-products.md), following Alda:2026nrd and DiLuzio:2026poh. The retained graphs need not exhaust the contributions for general sources. Units are GeV. Loop-momentum products are $m_a^2$, including the product of the two internal momenta of the closed line; no loop integral is evaluated.

## Method and result

The [script](../../../scripts/derive_radial_loop_class_one_tadpoles.py) prints diagrams 1.1–1.3, each a single vertex with a closed ALP line, times the radial exchange for 1.2 and 1.3. The six-pole vertex of 1.1 contains the six-meson kinetic part from $\mathcal L_6=\langle[\Phi,[\Phi,\partial_\mu\Phi]][\Phi,[\Phi,\partial^\mu\Phi]]\rangle/(1440f_\pi^4)$: for each pair of legs $c,f$ carrying the derivatives, the coefficient $-k_c\cdot k_f/(720f_\pi^4)$ multiplies the sum over the orderings of the other four meson matrices in the open positions of $\langle[\cdot,[\cdot,\lambda_c]][\cdot,[\cdot,\lambda_f]]\rangle$. Pairs giving the same trace are combined. Symmetry and loop factors are omitted.

## Reproduction and verification

Run `PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_class_one_tadpoles.py`. The test `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_class_one_tadpoles` expands every printed product back into trace words and requires exact equality with the independent direct expansion of the Lagrangian described in the [generator-sum record](murayama-loop-generator-trace-products.md), in which each meson line is a pair of placeholder matrices removed by the completeness rewrites, plus the bare angular term, minus the ALP pole term. Words are compared up to cyclic rotation only. All coefficients are exact rationals times symbols. The direct expansion derives the six-meson term from $\frac{f_\pi^2}{4}\langle\partial_\mu U\partial^\mu U^\dagger\rangle$ with $U=u^2$, so the check also fixes the coefficient $1/1440$. The check does not test diagrams outside the retained list, nor the sizing prescription as a loop integral.
