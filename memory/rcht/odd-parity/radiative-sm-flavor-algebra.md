---
title: Flavor traces for Standard Model vector radiative decays
status: derived
created: 2026-10-04
updated: 2026-10-13
bibtex_keys: [Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg, Alda:2026nrd]
provenance: analytical-derivation
---

# Flavor traces for $V\leftrightarrow P\gamma$

The action-derived implementation and complete channel amplitude are recorded in [radiative-sm-amplitude.md](radiative-sm-amplitude.md). This file defines the flavor traces and conventions used there; the numerical fit inputs remain explicit.

## Field and spurion conventions

For canonical RChT vectors, $V=\lambda^aV^a/\sqrt2$ with $\operatorname{tr}(\lambda^a\lambda^b)=2\delta^{ab}$. Thus $T_{\rho^0}=\operatorname{diag}(1,-1,0)/\sqrt2$, $T_\omega=\operatorname{diag}(1,1,0)/\sqrt2$, and $T_\phi=\operatorname{diag}(0,0,1)$. The electromagnetic matrix is $Q=\operatorname{diag}(2/3,-1/3,-1/3)$. The pseudoscalar coefficient matrix $T_P$ is supplied in the RChT field convention; if it is expressed in Alda:2026nrd's field matrix, first use $\Phi_R=F\Phi_A/(\sqrt2f_\pi)$ as derived in the [SM amplitude record](radiative-sm-amplitude.md).

Let $\chi=2B_0m_q$ be the real diagonal mass-squared spurion. For external vector coefficient $T_V$ and pseudoscalar coefficient $T_P$, the direct VJP traces are

$$
F_Q=\operatorname{tr}(\{T_V,Q\}T_P),\qquad
F_H=\operatorname{tr}(\{T_V,Q\}\{T_P,\chi\}),\qquad
F_K=\operatorname{tr}(T_V[[T_P,Q],\chi]).
$$

$F_Q$ is dimensionless and $F_H,F_K$ have dimension GeV$^2$. The last trace comes from the one-meson term in $f_-$ and supplies the $c_4$ contribution for charged channels. In $K^{*+}\to K^+\gamma$, $F_K=-(\chi_s-\chi_u)$ with the charged-leg conjugation convention; for neutral diagonal pseudoscalars it vanishes.

## Complete internal-photon flavor sum

The photon couples to each neutral vector $T_\alpha$ with overlap $w_\alpha=\operatorname{tr}(T_\alpha Q)$. Define

$$
F_{Q,\alpha}=\operatorname{tr}(\{T_V,T_\alpha\}T_P),\qquad
F_{H,\alpha}=\operatorname{tr}(\{T_V,T_\alpha\}\{T_P,\chi\}).
$$

The exchange amplitude sums $w_\alpha F_{Q,\alpha}$ and $w_\alpha F_{H,\alpha}$ over $\rho^0,\omega,\phi$, each divided by its pole mass squared. With a common pole, $\sum_\alpha w_\alpha T_\alpha=Q$, and the exchange traces reduce to the direct charge and mass traces. A rho-only exchange is therefore not the complete neutral-photon response. The full Lorentz coefficients and normalization are given in [radiative-sm-amplitude.md](radiative-sm-amplitude.md).

## Inputs and reproduction

No ALP spurion enters these Standard Model channels. Physical pseudoscalar mixing, decay constants, residues, and mass-spurion traces are not inferred from meson masses. The exploratory channel matrices and rank assumptions are listed in [radiative-rank-saturation.md](radiative-rank-saturation.md).

```text
.venv/bin/python scripts/derive_radiative_sm_amplitude.py
```
