---
title: Equivalence of the Cartesian and polar tree amplitudes of sigma to pi+ pi-
status: derived
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Equivalence of the Cartesian and polar trees of $\sigma\to\pi^+\pi^-$

Task 2.2 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama-cartesian`. This is the only comparison with the polar results used in the derivation of the Cartesian estimate. The Cartesian estimate is an independent order-of-magnitude estimate. No loop integral, scalar loop function, loop coefficient or loop amplitude is evaluated.

## Conventions and inputs

- Cartesian tree $M_0=V(s;\pi^+,\pi^-)=-\frac{B_0\epsilon}{2f_\pi}\langle\widehat S^C_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle$, derived from the Cartesian Lagrangian: [Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md), Step 6.
- Polar tree, taken as a stated input (no polar vertex code is used): $M_0^{\rm polar}=-\frac{B_0\epsilon}{2f_\pi}\langle\widehat S_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle$ ([polar tree and legs record](../one-loop-estimate/murayama-sigma-pipi-loop-tree-and-legs.md)), with $\widehat M(\sigma,a)=\sum(\widehat S_{r,s}+i\widehat P_{r,s})(\sigma/f_a)^r(a/f_a)^s$.
- Exact field map $f_a+s+ip=(f_a+\sigma)e^{ia/f_a}$ ([general-source radial record](../../../couplings-and-decays/finite-mixing/radial-decay-polar-cartesian-general-sources.md)). It is used here, and only here, to relate the coefficients.

Units GeV; conventions of the frame record.

## Method

1. **Relation between the coefficients.** $\widehat M$ is one function on field space. With $x=\sigma/f_a$, $y=a/f_a$, $X=s/f_a$, $Y=p/f_a$, the map $1+X+iY=(1+x)e^{iy}$ gives $x=\sqrt{(1+X)^2+Y^2}-1$ and $y=\arctan\big(Y/(1+X)\big)$. The script inserts the Taylor series of $x$ and $y$ into $\sum(\widehat S_{r,s}+i\widehat P_{r,s})x^ry^s$ and collects $X^RY^Q$. The Taylor coefficients are real, so $\widehat S^C$ ($\widehat P^C$) is the same real combination of the $\widehat S$ ($\widehat P$). The tree needs only the first order:
$$
\widehat S^C_{0,0}=\widehat S_{0,0},\quad \widehat S^C_{1,0}=\widehat S_{1,0},\quad \widehat S^C_{0,1}=\widehat S_{0,1},\qquad
\widehat P^C_{0,0}=\widehat P_{0,0},\quad \widehat P^C_{1,0}=\widehat P_{1,0},\quad \widehat P^C_{0,1}=\widehat P_{0,1}.
$$
The script checks that $(1+x)e^{iy}-(1+X+iY)$ vanishes to the order of the series.
2. **On-shell tree.** The external $\sigma$ is the pole of $s$, which has no meson or $p$ component and a canonical kinetic term. The charged legs are $(\lambda_{\pi^\pm},0)$. So the on-shell Cartesian amplitude is the vertex $V(s;\pi^+,\pi^-)$ with no momentum dependence. Inserting $\widehat S^C_{1,0}=\widehat S_{1,0}$ gives
$$
V(s;\pi^+,\pi^-)-M_0^{\rm polar}=0 \qquad\text{(exact).}
$$

**Remark.** The charged legs have no ALP component (charge conservation, a fact). Therefore neither the stationarity condition of $V_\varphi$ nor the ALP-like pole equation enters this comparison. The general on-shell identity of the [general-source radial record](../../../couplings-and-decays/finite-mixing/radial-decay-polar-cartesian-general-sources.md) needs them only for final states with an ALP component. This is a remark; the requirement of the OpenSpec change is not modified.

No off-shell, quartic or quintic vertex is compared with the polar frame.

## Verification

Exact symbolic checks: the map series to orders 1 and 2; the first-order relation; the difference of the trees is zero after the relation is inserted; without the relation the Cartesian and polar letters are different, so the comparison is not empty; the charged legs have no $p$ component.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_tree_equivalence.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_tree_equivalence
```

[Script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_tree_equivalence.py), [tests](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_tree_equivalence.py).
