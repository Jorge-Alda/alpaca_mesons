---
title: DiLuzio:2026poh effective chiral theory after integrating out the PQ scalars
status: literature
created: 2026-11-13
updated: 2026-11-13
bibtex_keys: [DiLuzio:2026poh]
provenance: literature
---

# Heavy-scalar effective chiral theory

Regime: $n_s=0$, $n_{u,d}\le1$, with $m_{\sigma,a}\gtrsim m_{\pi,K,\eta}$. The PQ scalars are integrated out at the matching scale $m_\varphi$ and the light degrees of freedom are $\Pi$ of [@DiLuzio:2026poh, Eq. (8)], including $\eta'$.

Effective Lagrangian [@DiLuzio:2026poh, Eq. (39)]:
$$
\mathscr L_{\rm eff}=\mathscr L_{\rm QCD}(\tilde M_q)+\sum_iC_i\mathcal O_i^{\rm PQ},
$$
with spurions $\tilde M_q$, $I_{\rm PQ}$ of Eqs. (36)-(38) (see [angular-mode-couplings.md](angular-mode-couplings.md)). Operators for $n_q=0,1$:
$$
\mathcal O_1^{\rm PQ}=\mathrm{tr}[I_{\rm PQ}^\dagger U]\,\mathrm{tr}[I_{\rm PQ}U^\dagger]\ (p^4)\ [\text{Eq. (40)}],
$$
$$
\mathcal O_2^{\rm PQ}=\mathrm{tr}[I_{\rm PQ}^\dagger\partial_\mu U]\,\mathrm{tr}[I_{\rm PQ}\partial^\mu U^\dagger],\quad
\mathcal O_3^{\rm PQ}=\mathrm{tr}[I_{\rm PQ}^\dagger(\partial_\mu U)(\partial^\mu U)^\dagger I_{\rm PQ}]\ (p^6)\ [\text{Eq. (41)}].
$$
Tree-level matching to the model [@DiLuzio:2026poh, Eq. (42)]:
$$
C_1=\frac14\frac{f_\pi^4B_0^2}{m_\varphi^2},\qquad C_2=\frac14\frac{f_\pi^4B_0^2}{m_\varphi^4},\qquad C_3=0.
$$
Neutral and charged pion masses [@DiLuzio:2026poh, Eqs. (43)-(44)] (for $m_s\gg m_{u,d}$, $\tilde m=m_u(1-n_u)+m_d(1-n_d)$):
$$
m_{\pi^0}^2\approx B_0\tilde m\Big(1+\frac{x_{\rm KM}m_s}2\Big)+\frac{8C_1\kappa_u\kappa_dn_un_d}{f_\pi^2},\qquad
m_{\pi^\pm}^2\approx B_0\tilde m\Big(1+\frac{x_{\rm KM}m_s}2\Big)+\frac{2C_1(\kappa_un_u+\kappa_dn_d)^2}{f_\pi^2},
$$
and $m_{\pi^0}^2-m_{\pi^\pm}^2=-\frac{2C_1}{f_\pi^2}(\kappa_un_u-\kappa_dn_d)^2$ [@DiLuzio:2026poh, Eq. (45)].

Explicit pion-only form of $\mathcal O_1^{\rm PQ}$ [@DiLuzio:2026poh, Eq. (48)]:
$$
\mathcal O_1^{\rm PQ}=(n_u\kappa_u+n_d\kappa_d)^2\cos^2\pi+(n_u\kappa_u-n_d\kappa_d)^2\frac{(\pi^0)^2}{\pi^2}\sin^2\pi .
$$

## Assumptions and caveats stated by the paper

- Tree-level matching only; loop corrections are expected to be large because $\kappa\sim1$ and $C_1$ breaks the usual $p^2$ chiral counting.
- $\mathcal O_{2,3}^{\rm PQ}$ act as wave-function renormalizations and give corrections suppressed by $1/m_\varphi^2$ to the masses.
- The restriction to $n_q=0,1$ excludes operators with higher powers of quarks.
- Eq. (45) is quoted by the paper as consistent with its tree-level result Eq. (24) for $n_{u,d}\le1$.
