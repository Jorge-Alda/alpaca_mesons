---
title: Polar and Cartesian radial-mode vertices for general sources
status: mixed
created: 2026-11-18
updated: 2026-11-18
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Polar and Cartesian $\sigma\to P_iP_j$ vertices for general sources

The equality is derived analytically with SymPy and checked numerically with random sources.

## Conventions and units

Natural units; all fields, masses, $B_0$, $f_\pi$, $f_a$, $\widehat S_{rs}$, $\widehat P_{rs}$ and $g_{ij}$ are in GeV. The conventions are those of the draft section on radial-mode decays and of Alda:2026nrd [@Alda:2026nrd]. The polar fields and model follow DiLuzio:2026poh [@DiLuzio:2026poh].

- Partonic source $\widehat M(a,\sigma)=\widehat S+i\widehat P$, with $\widehat S=\sum_{r,s}\widehat S_{rs}(\sigma/f_a)^r(a/f_a)^s$ and likewise for $\widehat P$. The coefficients are general Hermitian $3\times3$ matrices of mass dimension one, as in [the only-up record](murayama-sigma-two-pseudoscalar-decays.md). No flavour structure is assumed.
- $\Phi=\sum_kP_k\lambda_k$ with $\langle\lambda_k\lambda_l\rangle=2\delta_{kl}$ and $\lambda_0=\sqrt{2/3}\,\mathbf 1$. The fields are $u=e^{i\Phi/(2f_\pi)}$ and $U=u^2$, and $\chi=2B_0\widehat M$. The Lagrangian is $\tfrac{f_\pi^2}{4}\langle u_\mu u^\mu+\chi_+\rangle-\tfrac{m_1^2}{12}\langle\Phi\rangle^2$, which gives $\langle\chi_+\rangle=4B_0\langle\widehat S+\widehat P\Phi/f_\pi-\widehat S\Phi^2/(2f_\pi^2)\rangle+\mathcal O(\Phi^3)$.
- Polar kinetic term $\tfrac12(\partial\sigma)^2+\tfrac12(1+\sigma/f_a)^2(\partial a)^2$ and radial potential $V_\varphi=\tfrac{m_\varphi^2}{2}(f_a+\sigma)^2+\tfrac{\lambda}{4}(f_a+\sigma)^4$.
- Cartesian fields from the exact map $f_a+s+ip=(f_a+\sigma)e^{ia/f_a}$, with kinetic term $\tfrac12(\partial s)^2+\tfrac12(\partial p)^2$ and potential $V(s,p)=V_\varphi(s,p)-B_0f_\pi^2\langle\widehat S(\sigma(s,p),a(s,p))\rangle$ at $U=\mathbf1$.
- $\mathbb M_s(X)=B_0\{\widehat S_{00},X\}+\tfrac{m_1^2}{3}\langle X\rangle\mathbf 1-sX$ on Hermitian $3\times3$ matrices.
- $\mathcal L\supset-\tfrac12\sigma\sum_{ij}g_{ij}P_iP_j$, with outgoing on-shell momenta and positive angular components of the pole directions.

## Assumptions

1. No tadpoles: $\widehat P_{00}=0$ and $\langle\widehat S_{01}\rangle=0$. The radial tadpole is removed by the stationarity condition below.
2. No radial mixing at quadratic order: $\widehat P_{10}=0$ and $\langle\widehat S_{11}\rangle=0$. In Cartesian coordinates the $sp$ mixing coefficient is $B_0f_\pi^2(\langle\widehat S_{11}\rangle-\langle\widehat S_{01}\rangle)/f_a^2$, so assumption 1 is also used here.
3. The anomalous singlet mass term does not depend on $\sigma$ or $a$; any gluonic ALP coupling is included in the dressed source.
4. Work at tree level and leading $U(3)$ chiral order, with no expansion in $f_\pi/f_a$. Use the original representative without a radial connection.
5. Use regular simple poles with a nonzero angular component, so that $\mathbb M_{m_r^2}$ is invertible. For general $\widehat S_{00}$ the eigenvalues $B_0(m_a+m_b)$ of $B_0\{\widehat S_{00},\cdot\}$ are doubly degenerate for $a\neq b$. One combination in each pair has no ALP component and is an exact meson pole, as for the $K^0$ in [the $K^0$ record](../../../meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-k0-doublet-exact-pole.md). The trace forms below do not apply to such poles. The direct comparison of the two Lagrangians covers them.

The random-source check uses the full $U(3)$ basis $\lambda_0,\dots,\lambda_8$. Non-commuting sources violate charge conservation, so the neutral and charged fields mix.

## Inputs

- Pole equation, $\mathbb M_s$ and normalization: [only-up record](murayama-sigma-two-pseudoscalar-decays.md) and [operator inverse](../../../meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-finite-mixing-operator-inverse.md).
- Exact coordinate map and general on-shell equivalence: [polar and Cartesian record](../../../decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-polar-cartesian-amplitude-equivalence.md).
- Radial potential: DiLuzio:2026poh, [model extract](../../../../extracts/DiLuzio:2026poh/model-lagrangian-and-field-normalization.md).

## Step 1: pole equation, direction and angular row

The quadratic Lagrangian gives $\mathbb M_s(\Phi)=(2B_0f_\pi/f_a)\widehat P_{01}a$. For $a=1$,
$$
\Omega_r\equiv\sum_kw_k^{(r)}\lambda_k=\frac{2B_0f_\pi}{f_a}\mathbb M_{m_r^2}^{-1}(\widehat P_{01}),\qquad
N(m_r^2)=1+\tfrac12\langle\Omega_r^2\rangle=1+\frac{2B_0^2f_\pi^2}{f_a^2}\left\langle\mathbb M_{m_r^2}^{-1}(\widehat P_{01})^2\right\rangle .
$$
The angular row is
$$
m_r^2+\frac{2B_0f_\pi^2}{f_a^2}\langle\widehat S_{02}\rangle=-\frac{B_0f_\pi}{f_a}\langle\widehat P_{01}\Omega_r\rangle .
$$
Inserting $\Omega_r$ turns it into the pole equation $s+\tfrac{2B_0f_\pi^2}{f_a^2}[\langle\widehat S_{02}\rangle+B_0\langle\widehat P_{01}\mathbb M_s^{-1}(\widehat P_{01})\rangle]=0$ at $s=m_r^2$. This holds identically as a function of $s$, for general sources, and is not a special case.

## Step 2: radial vacuum

The coordinate map gives $\sigma=s+p^2/(2f_a)-sp^2/(2f_a^2)+\dots$ and $a=p-sp/f_a+\dots$. Stationarity in $s$ and the second derivative give
$$
m_\varphi^2=\frac{B_0f_\pi^2\langle\widehat S_{10}\rangle}{f_a^2}-\lambda f_a^2,\qquad
m_\sigma^2=\left.\frac{\partial^2V}{\partial s^2}\right|_0=2\lambda f_a^2+\frac{B_0f_\pi^2}{f_a^2}\big(\langle\widehat S_{10}\rangle-2\langle\widehat S_{20}\rangle\big).
$$
The same $m_\sigma^2$ follows in polar coordinates. Only $\langle\widehat S_{10}\rangle$ enters the stationarity condition. $\langle\widehat S_{02}\rangle$ enters the masses of the pseudoscalars, not $m_\varphi^2$.

## Step 3: cubic Lagrangians

Polar:
$$
\mathcal L_{\sigma PP}=\frac{\sigma}{f_a}\left[-\frac{B_0}{2}\langle\widehat S_{10}\Phi^2\rangle+\frac{B_0f_\pi}{f_a}a\langle\widehat P_{11}\Phi\rangle+\frac{B_0f_\pi^2}{f_a^2}a^2\langle\widehat S_{12}\rangle+(\partial a)^2\right].
$$
The Cartesian Lagrangian is obtained by substituting the exact map. The polar terms $\langle\widehat P_{01}\Phi\rangle a$, $\langle\widehat S_{02}\rangle a^2$, $\langle\widehat S_{10}\rangle\sigma$, $\langle\widehat S_{20}\rangle\sigma^2$ and $V_\varphi$ then produce cubic terms:
$$
\mathcal L_{sPP}=-\frac{B_0s}{2f_a}\langle\widehat S_{10}\Phi^2\rangle+\frac{B_0f_\pi}{f_a^2}sp\langle(\widehat P_{11}-\widehat P_{01})\Phi\rangle
+\left[\frac{B_0f_\pi^2}{f_a^3}\Big(-\tfrac12\langle\widehat S_{10}\rangle+\langle\widehat S_{20}\rangle-2\langle\widehat S_{02}\rangle+\langle\widehat S_{12}\rangle\Big)-\lambda f_a\right]sp^2 .
$$
The Cartesian kinetic term has no cubic vertex. In the only-up model $\widehat P_{11}=\widehat P_{01}$, and the trace combination in round brackets vanishes. Only $-\lambda f_asp^2$ is then left, which is the result of [the only-up record](murayama-sigma-two-pseudoscalar-decays.md).

## Step 4: on-shell vertices

Insert $\Phi=\sum_r\Omega_rP_r/\sqrt{N(m_r^2)}$, $a=p=\sum_rP_r/\sqrt{N(m_r^2)}$ and $2p_i\cdot p_j=m_\sigma^2-m_i^2-m_j^2$:
$$
\begin{aligned}
g^{\rm polar}_{ij}={}&\frac1{f_a}\Big\{\frac{2B_0^3f_\pi^2}{f_a^2}\big\langle\widehat S_{10}\{\mathbb M_{m_i^2}^{-1}(\widehat P_{01}),\mathbb M_{m_j^2}^{-1}(\widehat P_{01})\}\big\rangle
-\frac{2B_0^2f_\pi^2}{f_a^2}\big\langle\widehat P_{11}\big(\mathbb M_{m_i^2}^{-1}(\widehat P_{01})+\mathbb M_{m_j^2}^{-1}(\widehat P_{01})\big)\big\rangle\\
&-\frac{2B_0f_\pi^2}{f_a^2}\langle\widehat S_{12}\rangle+m_\sigma^2-m_i^2-m_j^2\Big\}\big[N(m_i^2)N(m_j^2)\big]^{-1/2},\\[1mm]
g^{\rm Cart}_{ij}={}&\frac1{f_a}\Big\{2\lambda f_a^2+\frac{2B_0^3f_\pi^2}{f_a^2}\big\langle\widehat S_{10}\{\mathbb M_{m_i^2}^{-1}(\widehat P_{01}),\mathbb M_{m_j^2}^{-1}(\widehat P_{01})\}\big\rangle
-\frac{2B_0^2f_\pi^2}{f_a^2}\big\langle(\widehat P_{11}-\widehat P_{01})\big(\mathbb M_{m_i^2}^{-1}(\widehat P_{01})+\mathbb M_{m_j^2}^{-1}(\widehat P_{01})\big)\big\rangle\\
&+\frac{B_0f_\pi^2}{f_a^2}\big(\langle\widehat S_{10}\rangle-2\langle\widehat S_{20}\rangle+4\langle\widehat S_{02}\rangle-2\langle\widehat S_{12}\rangle\big)\Big\}\big[N(m_i^2)N(m_j^2)\big]^{-1/2}.
\end{aligned}
$$

## Result: general on-shell identity

Subtracting the two expressions,
$$
f_a\sqrt{N(m_i^2)N(m_j^2)}\,\big(g^{\rm polar}_{ij}-g^{\rm Cart}_{ij}\big)
=\Big[m_\sigma^2-2\lambda f_a^2-\frac{B_0f_\pi^2}{f_a^2}\big(\langle\widehat S_{10}\rangle-2\langle\widehat S_{20}\rangle\big)\Big]
-\sum_{r=i,j}\Big[m_r^2+\frac{2B_0f_\pi^2}{f_a^2}\Big(\langle\widehat S_{02}\rangle+B_0\big\langle\widehat P_{01}\,\mathbb M_{m_r^2}^{-1}(\widehat P_{01})\big\rangle\Big)\Big].
$$
The first bracket is zero because the radial potential is stationary (Step 2). Each term of the sum is zero because it is the pole equation at $s=m_r^2$ (Step 1). Hence $g^{\rm polar}_{ij}=g^{\rm Cart}_{ij}$ on shell for general Hermitian sources, under assumptions 1–5, including the ALP-like state among the daughters. The $\widehat S_{10}$ anticommutator terms are identical in the two forms, and so are the $\widehat P_{11}$ and $\widehat S_{12}$ terms. Only the $\widehat P_{01}$, $\widehat S_{02}$, $\widehat S_{10}$, $\widehat S_{20}$ and $\lambda$ terms need the two relations.

This agrees with the general argument of [the polar and Cartesian record](../../../decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-polar-cartesian-amplitude-equivalence.md). The map is the identity at linear order and the tadpoles vanish, so the change of the cubic vertex is proportional to an inverse propagator acting on an external pole vector. If $m_\varphi^2$ violates stationarity, the radial tadpole multiplies the cubic part of $\sigma(s,p)$, and the two vertices differ. The script checks this.

## Verification

[The script](../../../../../scripts/derive_radial_decay_polar_cartesian_general.py) performs the following steps.

- SymPy expands the map. It checks that the pulled-back Cartesian metric equals the polar metric. It derives the stationarity condition, $m_\sigma^2$ and both cubic Lagrangians in trace symbols, and reduces the vertex difference to zero by using the two relations.
- Numerically, it builds both Lagrangians from the exact map with random non-commuting Hermitian sources. It uses all coefficients with $r+s\le3$ other than those fixed by the assumptions, and fixed seeds 2026, 7 and 11. The expansion is exact to cubic order: it is a polynomial in three formal variables with vanishing squares. The script then computes the poles from the Hessian and contracts the cubic terms, including the polar derivative term, with the pole vectors.

Residuals for seed 2026, with 7 of the 10 poles carrying an ALP component:

| check | residual |
|---|---|
| tadpole; polar minus Cartesian Hessian; radial mixing | $3\times10^{-18}$; $2\times10^{-19}$; $0$ |
| $m_\sigma^2$ trace form | $0$ |
| pole equation and angular row at the numerical poles (GeV$^2$) | $6\times10^{-12}$ |
| $\Omega_r$ against the eigenvectors, and $N$ (relative) | $3\times10^{-12}$, $6\times10^{-12}$ |
| direct $g^{\rm polar}$ against direct $g^{\rm Cart}$ (all poles, relative) | $3\times10^{-16}$ |
| direct against trace form, polar and Cartesian (relative) | $3\times10^{-13}$ |

The only-up model gives a direct polar–Cartesian residual of $10^{-17}$ and a residual against the trace forms of $\le2\times10^{-12}$. The tests also show that the only-up Cartesian form fails for general sources, with a relative difference above $10^{-3}$, and that the two direct vertices differ once $m_\varphi^2$ is shifted away from stationarity.

```sh
.venv/bin/python scripts/derive_radial_decay_polar_cartesian_general.py
.venv/bin/python -m unittest scripts.test_derive_radial_decay_polar_cartesian_general
```
