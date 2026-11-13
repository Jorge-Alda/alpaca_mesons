---
title: DiLuzio:2026poh parametrization of varphi, angular mode and its normalization
status: literature
created: 2026-11-13
updated: 2026-11-13
bibtex_keys: [DiLuzio:2026poh]
provenance: literature
---

# Parametrization of $\varphi$ and the angular mode

## Polar decomposition

The complex field is expanded in a CP-even radial mode $\sigma$ and a CP-odd phase $a$ [@DiLuzio:2026poh, Eq. (9)]:
$$
\varphi=\frac{1+\sigma}{\sqrt2}\,e^{ia}.
$$
Both $\sigma$ and $a$ are dimensionless here, because $\varphi$ is dimensionless. The radial mode is a scalar and $a$ is a pseudoscalar. For $\theta=0$ (always attainable by a PQ transformation) the CP-odd vacuum expectation values vanish by CP invariance.

## Normalization

- Inserting Eq. (9) into $f_\varphi^2|\partial\varphi|^2$ gives (with $\sigma=0$) the kinetic term $\frac{f_\varphi^2}2(\partial a)^2$, shown in [@DiLuzio:2026poh, Eq. (14)]. Hence the canonically normalized fields are $a_{\rm can}=f_\varphi a$ and $\sigma_{\rm can}=f_\varphi\sigma$. The paper states $\sigma_{\rm can}$ explicitly only in Sec. 5, in $\mathscr L\supset-y_q\sigma_{\rm can}\bar qq$, and gives the meson rule $\pi^0_{\rm can}=f_\pi\pi^0$ under Eq. (8). The rule $a_{\rm can}=f_\varphi a$ is the one implied by Eq. (14) (the extractor's reading, not a stated sentence).
- Pion kinetic term: $\frac{f_\pi^2}2[(\partial\pi^0)^2+(\partial\eta)^2+(\partial\eta')^2]$ [@DiLuzio:2026poh, Eq. (14)].
- Axion decay constant in the light-axion limit $f_\varphi\gg f_\pi$ [@DiLuzio:2026poh, Eq. (53)]:
$$
f_a\equiv\frac{f_\varphi}{\sum_qn_q}.
$$
So the canonical axion field is $f_\varphi a$ and the coupling to the anomaly is $a\sum_qn_q$, i.e. $a_{\rm can}/f_a$.

## Vacuum conditions and radial-mode mass

With the vacuum at $\sigma=0$, $a=0$ [@DiLuzio:2026poh, Eqs. (10), (13)]:
$$
m_\varphi^2=\frac{B_0f_\pi^2}{f_\varphi^2}\sum_qm_qn_q-\lambda f_\varphi^2,\qquad
m_\sigma^2=m_\varphi^2+3\lambda f_\varphi^2=\frac{B_0f_\pi^2}{f_\varphi^2}\sum_qm_qn_q+2\lambda f_\varphi^2 .
$$
For $n_q\le1$ the $(1,1)$ entry of the CP-odd mass matrix is $\mathcal M^2_{11}=m_\varphi^2+\lambda f_\varphi^2=m_\sigma^2-2\lambda f_\varphi^2$ [@DiLuzio:2026poh, Eq. (17)]. Eq. (13) is quoted for general charges, Eq. (17) only for $n_q\le1$.

## PQ transformation of the angular mode

Under $\bar q_Lq_R\to\bar q_Lq_Re^{i\alpha Q_{\rm PQ}}$ with $Q_{\rm PQ}=\mathrm{diag}(n_u,n_d,n_s)$ [@DiLuzio:2026poh, Eq. (15)]:
$$
a\to a-\alpha,\quad U\to e^{-i\alpha Q_{\rm PQ}}U,\quad
\pi^0\to\pi^0-\tfrac\alpha2(n_u-n_d),\quad
\eta\to\eta-\tfrac{\alpha}{2\sqrt3}(n_u+n_d-2n_s),\quad
\eta'\to\eta'-\tfrac{\alpha}{\sqrt6}(n_u+n_d+n_s).
$$

## Other parametrization in the same paper

In the two-scalar appendix model the phase is written with an explicit decay constant, $\varphi_{u,d}=\frac{f_{u,d}+\sigma_{u,d}}{\sqrt2}e^{ia_{u,d}/f_{u,d}}$, with $\langle\varphi_{u,d}\rangle=f_{u,d}/\sqrt2$ and $f_\varphi^2\equiv f_u^2n_d^2+f_d^2n_u^2$ (equations of the appendix "Renormalizable models with more PQ scalars", between Eqs. (58) and (62)). There the fields $a_{u,d}$ carry mass dimension and $\varphi_{u,d}$ are canonically normalized (kinetic term $|\partial\varphi_u|^2$, Eq. (58)).

## Ambiguities

- Main text: $a$, $\sigma$ dimensionless, $\varphi$ with $f_\varphi^2|\partial\varphi|^2$. Appendix: $a_{u,d}$ dimensionful, kinetic terms canonical. The two conventions must not be mixed.
- $a_{\rm can}$ is never defined by an equation in the paper.
