---
title: Quartic vertex with two ALP-like pole legs
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Quartic vertex with two ALP-like pole legs

## Inputs and assumptions

Use polar connection-free fields, incoming momenta, source coefficients and pole legs of the [estimate record](../../../couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md), and the compact vertex rules of the [generator-sum record](murayama-loop-generator-trace-products.md). The four-meson kinetic interaction is the leading-order one of Alda:2026nrd; the source conventions follow DiLuzio:2026poh. Masses and decay constants are in GeV. The two other legs $x=(\lambda_x,\xi_x)$ and $y=(\lambda_y,\xi_y)$ are general pole legs with incoming momenta $k_x$, $k_y$, not necessarily the physical pions. Products with the loop momentum are $m_a^2$. No loop integral is calculated.

## Method and result

The [script](../../../../../../scripts/derive_radial_loop_three_one.py) evaluates $V(a,a;x,y)$ with both ALP pole legs on the internal line. The potential part is the sum over the assignments of each leg to its meson matrix or angular component, with the sum over orderings of the meson matrices. In the kinetic part only the split of the four legs into the pairs $(a,x)$ and $(a,y)$ survives, since $[\lambda_a,\lambda_a]=0$, and the two ways of forming these pairs add:

$$
V_{\rm kin}(a,a;x,y)=-\frac{k_x\cdot k_y-m_a^2}{12f_\pi^2}\,\big\langle[\lambda_a,\lambda_x][\lambda_a,\lambda_y]\big\rangle .
$$

The script also prints the class-3 products 3.1–3.3 with the rules of the [generator-sum record](murayama-loop-generator-trace-products.md): two ALP lines, with $(x,y)=(a,\pi^0)$ and $(\pi^+,\pi^-)$, and the radial exchanges of 3.2 and 3.3 as shared factors. Symmetry factors ($1/2$) and the loop factor are excluded.

## Reproduction and verification

Run `PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_three_one.py`. The test `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_three_one` checks exactly that the printed $V(a,a;x,y)$ equals the direct expansion of the Lagrangian of the [generator-sum record](murayama-loop-generator-trace-products.md) for the general legs $x$, $y$; that its coefficient of $k_x\cdot k_y$ is the single commutator trace displayed above; and that the products 3.1–3.3, expanded into trace words, equal the direct expansion with the ALP lines inserted. Words are compared up to cyclic rotation only, with exact coefficients. These checks do not establish completeness of the graph selection for general sources.
