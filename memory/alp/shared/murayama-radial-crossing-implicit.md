---
title: Loop-to-tree crossing criterion with the radial mode, theta-hat zero
status: mixed
created: 2026-11-23
updated: 2026-11-25
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh, Manohar:1983md]
provenance: assumption
---

# Loop-to-tree crossing criterion with the radial mode

This record states the crossing criterion for the ALP-sourced three-meson amplitude in the Murayama model and its generalizations, with the radial mode $\sigma$ included and $\widehat\Theta_1=0$. The criterion is $R(\kappa_u)=R_{\max}$, with $R=|M_1^{\rm est}|/|M_0|$ and $R_{\max}$ symbolic. The numerical evaluation is in [the estimate record](murayama-radial-up-only-estimate-result.md).

## Scope

- Model: $\widehat\Theta_1=0$, so $\mathcal C=B_0\widehat P_{01}$ ([radial-decay record](murayama-sigma-two-pseudoscalar-decays.md)).
- Final state: three pseudoscalars, with the fixed independent Dalitz invariants of the [ALP-only crossing record](alp-loop-crossing-exact-tree-condition.md).
- $\kappa_u$ and $\lambda$ are analytical parameters. $R_{\max}$ is symbolic.

## The criterion

The exact tree is
$$
M_0=V(a,\pi^+,\pi^-,\pi^0)+\frac{V(\sigma,a,\pi^0)\,V(\sigma,\pi^+,\pi^-)}{m_\sigma^2-s_0},
$$
where $V$ denotes the tree vertex of the [radial-exchange record](../three-pseudoscalar/alp-three-pseudoscalar-radial-exchange.md#step-3-off-shell-sigma-vertex) and $s_0$ the fixed independent invariant. The explicit vertices are given in the [estimate record](murayama-radial-up-only-estimate-result.md). The $\sigma$ exchange in the channels $(a\pi^\pm)$ and $(\pi^0\pi^\pm)$ vanishes by charge.

The one-loop estimate is
$$
M_1^{\rm est}=\sum_d S_d\,\frac{(m_d^2)^{2-r_d}}{16\pi^2}\,\operatorname{tr}\big[V(B_1)W_{L_1}\cdots V(B_{r_d})W_{L_{r_d}}\big],
$$
with the cycles, line weights and scales as defined in the estimate record. The normalization $(m_d^2)^{2-r_d}/(16\pi^2)$ is dimensionally consistent, and the result is independent of the unit system.

The crossing criterion is
$$
R(\kappa_u)\equiv\frac{|M_1^{\rm est}(\kappa_u)|}{|M_0(\kappa_u)|}=R_{\max}.
$$

## Caveats

These are stated, not computed.

- **Non-degenerate poles are assumed.** The pole legs use simple roots $m_a^2$ and $m_{\pi^0}^2$ of the pole function, with $N=F'$ finite.
- **Near-degeneracy with the $\eta'$.** The ALP-like mass passes $m_{\eta'}$ near $\kappa_u\approx0.4$. This region is not computed.
- **Chiral perturbation theory at large $m_a$.** Its validity for $m_a$ from about 1 to 2.3 GeV ($\kappa_u\gtrsim0.43$) is not assessed.
- **Truncation.** Cycles with three or more A or S lines, and the sigma self-energies with pseudoscalar lines, are excluded; they are listed in the estimate record.
- **Fixed kinematic point.** The Dalitz point is the symmetric point. Other points are not computed.
- **Fixed coupling.** $\lambda=10^{-3}$ is the default value used in the numerical evaluation.

## Frame invariance

The sources enter only through the dressed invariants $\widehat S_{rs}$, $\widehat P_{rs}$ of a frame, and the pole data built from them. These are unchanged under a constant axial redefinition $U=\xi U'\xi$ (see the [frame-invariance record](murayama-radial-frame-invariance.md)). The estimate was checked with the polar frame and a redefined frame with $\kappa_{\rm ax}={\rm diag}(0.2,-0.2,0)$, which changes the raw sources by more than $10^{-4}$. The relative differences were:

| quantity | $\kappa_u=0.3$ | $\kappa_u=1.0$ |
|---|---|---|
| $M_1^{\rm est}/M_0$ | $7.6\times10^{-16}$ | $1.3\times10^{-15}$ |
| largest single contribution divided by $M_0$ | $3.6\times10^{-14}$ | $5.7\times10^{-12}$ |
| $M_0$ | $\le1.7\times10^{-15}$ | $\le1.7\times10^{-15}$ |
