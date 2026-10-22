---
title: Neutral-kaon selection rules in the odd ALP photon-pseudoscalar-pair amplitude
status: derived
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and inputs

This record gives the neutral-kaon selection rules for the retained odd-resonance contribution to $a\to\gamma P_1P_2$. The ALP axial and pseudoscalar sources are already in the no-mixing basis. We use the amplitude and width conventions in the [odd $a\to\gamma PP$ amplitude record](alp-photon-pp-amplitude.md) and [its Dalitz-width record](alp-photon-pp-width.md), together with the [neutral-kaon projection convention](../../kaon/neutral-kaon-amplitude-projection.md) [@Alda:2026nrd; @ParticleDataGroup:2024cfk; @Ecker:1988te; @Ecker:1989yg; @Ruiz-Femenia:2003jdx]. The statement is for the specified tree-level odd amplitude and coherent physical-kaon production amplitudes; it does not add local anomalous or even radiative terms.

# Flavor and exchange symmetry

For fixed momenta, let $\mathcal A_{0\bar 0}(s,t,u)$ denote the scalar coefficient multiplying
$$
\epsilon_{\mu\nu\rho\sigma}\varepsilon_\gamma^\mu p_1^\nu p_2^\rho k^\sigma
$$
for the ordered state $(K^0(p_1),\bar K^0(p_2))$. The flavor traces of the retained odd amplitude give
$$
\mathcal A_{\bar 0 0}(s,t,u)=-\mathcal A_{0\bar 0}(s,t,u),
$$
when the labels $p_1,p_2$ are held fixed and the flavor assignments alone are exchanged. The full amplitude is Bose symmetric: exchanging the two equal-mass kaon momenta also changes the sign of the Levi-Civita tensor, so this antisymmetry of the scalar coefficient is required and does not imply a vanishing physical $K_LK_S$ amplitude.

Use
$$
|K_S\rangle=p|K^0\rangle+q|\bar K^0\rangle,\qquad
|K_L\rangle=p|K^0\rangle-q|\bar K^0\rangle,
$$
where $p=(1+\epsilon_K)/\sqrt{2(1+|\epsilon_K|^2)}$ and $q=(1-\epsilon_K)/\sqrt{2(1+|\epsilon_K|^2)}$. As final-state bras, the coefficients are conjugated. Therefore
$$
\begin{aligned}
\mathcal A(K_L(p_1)K_S(p_2))
&=p^*q^*\mathcal A_{0\bar0}+(-q^*)p^*\mathcal A_{\bar00}\\
&=p^*q^*(\mathcal A_{0\bar0}-\mathcal A_{\bar00})
=2p^*q^*\mathcal A_{0\bar0}.
\end{aligned}
$$
Thus $K_LK_S$ is allowed and remains nonzero for arbitrary finite $\epsilon_K$ (unless the underlying source amplitude itself vanishes). Its rate carries the coherent factor $|2p^*q^*|^2$ relative to the definite-flavor ordered scalar amplitude, before the common phase-space integration.

For equal physical kaon masses, the equal-state projections instead give
$$
\begin{aligned}
\mathcal A(K_S(p_1)K_S(p_2))
&=p^*q^*(\mathcal A_{0\bar0}+\mathcal A_{\bar00})=0,\\
\mathcal A(K_L(p_1)K_L(p_2))
&=-p^*q^*(\mathcal A_{0\bar0}+\mathcal A_{\bar00})=0.
\end{aligned}
$$
The definite-flavor $K^0K^0$ and $\bar K^0\bar K^0$ amplitudes vanish for the diagonal neutral ALP sources considered here by strangeness conservation. Consequently, among the four ordered neutral physical pairs, only $K_LK_S$ and its reversed ordering $K_SK_L$ are allowed. Reversing the ordered labels changes the scalar coefficient sign and exchanges the Dalitz variables; it describes the same physical channel when treated with a consistent phase-space convention.

# Numerical implementation and verification

The selected-channel implementation is [`scripts/implementation/alp_photon_pp_width.py`](../../../scripts/implementation/alp_photon_pp_width.py). It expands each physical kaon bra into $K^0,\bar K^0$ components, applies the conjugated mixing coefficients before assembling flavor traces, and groups equal invariant routings and vector poles before quadrature. Flavor coefficients that cancel at the scale-relative floating-point criterion are set to exact zero. The implementation thus returns exact empty amplitudes for $K_SK_S$, $K_LK_L$, $K^0K^0$, and $\bar K^0\bar K^0$ before integration; it retains the five grouped terms for $K_LK_S$.

Checks are executable with
```text
.venv/bin/python -m unittest scripts.implementation.test_alp_photon_pp_width
.venv/bin/python scripts/derive_alp_photon_two_pseudoscalar.py --test
```
The tests verify the complex-$\epsilon_K$ bra projection, exact pre-integration zeros, and the nonzero $K_LK_S$ amplitude. Representative allowed-channel integrations converge without warnings. These checks support the stated odd-resonance tree truncation; they do not establish the same cancellation for additional amplitudes with different flavor or kinematic structures.
