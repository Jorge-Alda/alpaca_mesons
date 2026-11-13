---
title: DiLuzio:2026poh couplings of the angular and radial modes to quarks and mesons
status: mixed
created: 2026-11-13
updated: 2026-11-13
bibtex_keys: [DiLuzio:2026poh]
provenance: literature
---

# Couplings of $a$ and $\sigma$ to quarks and mesons

The paper does not expand the couplings in $a$ or $\sigma$: it keeps the exact dependence through $M_q(\varphi)$. Statements marked "extractor" are elementary rewritings of the paper's equations, not equations printed in the paper.

## Exact quark-mass dependence

Substituting [@DiLuzio:2026poh, Eq. (9)] into $(\sqrt2\varphi)^{n_q}$ in Eq. (3) and Eq. (7) gives, to all orders in $a$ and $\sigma$ (extractor):
$$
(\sqrt2\varphi)^{n_q}=(1+\sigma)^{n_q}e^{in_qa},\qquad
-m_q(\sqrt2\varphi)^{n_q}\bar q_Lq_R+\text{h.c.}=-m_q(1+\sigma)^{n_q}\big[\cos(n_qa)\,\bar qq+i\sin(n_qa)\,\bar q\gamma_5q\big].
$$
The second form uses $\bar q_Lq_R=\bar q\frac{1+\gamma_5}2q$. At $a=\sigma=0$ the quark mass is $m_q$. The paper states only the linear scalar coupling [@DiLuzio:2026poh, Sec. 5]: with $\mathscr L\supset-y_q\sigma_{\rm can}\bar qq$ it quotes
$$
y_q=\frac{n_qm_q}{f_\varphi},
$$
using Eq. (3). This is the term linear in $\sigma_{\rm can}=f_\varphi\sigma$ of the expansion above (extractor), valid for any $n_q$ (coefficient $n_q$).

## Exact potential for the CP-odd fields

Putting $\varphi$ in Eq. (6) and keeping $\sigma=0$, the CP-odd Lagrangian is [@DiLuzio:2026poh, Eq. (14)]:
$$
\begin{aligned}
\mathscr L={}&\frac{f_\varphi^2}2(\partial a)^2+\frac{f_\pi^2}2\big[(\partial\pi^0)^2+(\partial\eta)^2+(\partial\eta')^2\big]-\frac{m_1^2}2f_\pi^2(\eta'-\theta)^2\\
&-B_0f_\pi^2\Big[m_s\cos\!\Big(n_sa+\frac{2\eta-\sqrt2\eta'}{\sqrt3}\Big)
+m_u\cos\!\Big(n_ua-\frac\eta{\sqrt3}-\sqrt{\tfrac23}\eta'-\pi^0\Big)
+m_d\cos\!\Big(n_da-\frac\eta{\sqrt3}-\sqrt{\tfrac23}\eta'+\pi^0\Big)\Big].
\end{aligned}
$$
Each cosine contains the complete phase $n_qa$ plus the meson phase from $U^\dagger$, so the $a$-dependence is exact. The potential is invariant under the PQ shift of Eq. (15) (see [angular-mode-parametrization.md](angular-mode-parametrization.md)). The $\theta$-term normalization issue is noted in [model-lagrangian-and-field-normalization.md](model-lagrangian-and-field-normalization.md).

The $\sigma$-dependence at the same level follows from $M_q\propto(1+\sigma)^{n_q}$: each $m_q$ in the bracket above multiplies $(1+\sigma)^{n_q}$ (extractor). The paper uses this only through the vacuum condition and $m_\sigma^2$, Eqs. (10), (13).

## Quadratic mixing of $a$ with mesons

The Hessian at the minimum in the basis $(a,\pi^0,\eta,\eta')$ is the symmetric matrix of [@DiLuzio:2026poh, Eq. (16)]:
$$
\mathcal M^2=B_0\begin{pmatrix}
\frac{f_\pi^2}{f_\varphi^2}\sum_qm_qn_q^2&\frac{f_\pi}{f_\varphi}(m_dn_d-m_un_u)&\frac{-f_\pi}{\sqrt3f_\varphi}(m_un_u+m_dn_d-2m_sn_s)&\frac{-\sqrt2f_\pi}{\sqrt3f_\varphi}\sum_qm_qn_q\\
\cdot&m_u+m_d&\frac{m_u-m_d}{\sqrt3}&\sqrt{\tfrac23}(m_u-m_d)\\
\cdot&\cdot&\frac13(m_d+m_u+4m_s)&\frac{\sqrt2}3(m_d+m_u-2m_s)\\
\cdot&\cdot&\cdot&\frac{m_1^2}{B_0}+\frac23(m_u+m_d+m_s)
\end{pmatrix}.
$$
Its determinant is $\det\mathcal M^2=\frac{4f_\pi^2}{3f_\varphi^2}B_0^3m_um_dm_s(n_u+n_d+n_s)^2m_1^2$ [@DiLuzio:2026poh, Eq. (18)]. Eq. (16) is the mass matrix of the canonically normalized fields $(a_{\rm can},\pi^0_{\rm can},\eta_{\rm can},\eta'_{\rm can})$: dividing the Hessian of Eq. (14) in the dimensionless fields by the kinetic matrix $\mathrm{diag}(f_\varphi^2,f_\pi^2,f_\pi^2,f_\pi^2)$ gives the $(1,1)$ entry $B_0\frac{f_\pi^2}{f_\varphi^2}\sum_qm_qn_q^2$ and the first-row entries proportional to $f_\pi/f_\varphi$ (extractor's check). The paper does not state this step explicitly.

## Spurion form in the heavy-scalar effective theory

For $n_s=0$, $n_{u,d}\le1$ the paper defines [@DiLuzio:2026poh, Eqs. (36)-(38)]
$$
\tilde M_q=\mathrm{diag}\big(m_u(1-n_u),m_d(1-n_d),m_s\big),\quad
I_{\rm PQ}=\mathrm{diag}(\kappa_un_u,\kappa_dn_d,0),\quad
I_{\rm PQ}\to I_{\rm PQ}e^{-i\alpha Q_{\rm PQ}}\ (Q_{\rm PQ}=\mathrm{diag}(n_u,n_d,0)),
$$
with $I_{\rm PQ}\to LI_{\rm PQ}R^\dagger$ under the chiral group; see [heavy-scalar-eft-operators.md](heavy-scalar-eft-operators.md).

## Couplings to the QCD scalar and pseudoscalar operators

In the appendix on the effective pion potential, $\kappa_q\varphi$ is upgraded to a source $\Phi\to L\Phi R^\dagger$ for the operators $\mathcal O_{LR}^{ij}=\bar q_L^iq_R^j$; see [scalar-pseudoscalar-source-operators.md](scalar-pseudoscalar-source-operators.md). The paper's explicit effective Lagrangian in the background of the neutral pion, including the dependence on $\pi$ to all orders, is [@DiLuzio:2026poh, Eq. (75)]
$$
\begin{aligned}
\mathscr L_{\rm eff}={}&\varphi(p)\varphi^*(-p)\Big\{(p^2-m_\varphi^2)+\Pi_A(p)(\kappa_u^2n_u^2+\kappa_d^2n_d^2)\\
&+\Pi_B(p)\Big[(\kappa_un_u+\kappa_dn_d)^2\cos^2\pi+(\kappa_un_u-\kappa_dn_d)^2\frac{(\pi^0)^2}{\pi^2}\sin^2\pi\Big]\Big\}\\
&+\Pi_0\,\varphi(p)\Big[(\kappa_un_u+\kappa_dn_d)\cos\pi+i(\kappa_un_u-\kappa_dn_d)\frac{\pi^0}\pi\sin\pi\Big]+\text{h.c.},
\end{aligned}
$$
with $\pi^2\equiv\pi^a\pi^a=(\pi^0)^2+2\pi^-\pi^+$ and $\Pi_0=f_\pi^2B_0/2$. Here $\varphi(p)$ is the Fourier transform of the complex field and is not split into $\sigma$ and $a$.

## Ambiguities

- The paper never writes the $\sigma$ and $a$ couplings to $\bar qq$ and $\bar q\gamma_5q$ explicitly; the rewriting above is the extractor's.
- Eq. (16) is given for canonically normalized fields although Eq. (14) is written in dimensionless fields; the rescaling is not stated in the paper.
- No expansion parameter is used for the exact phase; $\sigma$ and $a$ couplings are polynomial and trigonometric in the fields. The only expansion in $f_\pi/f_\varphi$ appears in the light-axion limit, Eq. (52).
