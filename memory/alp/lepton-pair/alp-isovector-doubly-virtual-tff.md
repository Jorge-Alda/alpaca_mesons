---
title: Isovector part of the doubly virtual no-mixing ALP TFF
status: derived
created: 2026-10-23
updated: 2026-10-23
bibtex_keys: [Alda:2026nrd, Bauer:2020jbp, Ecker:1988te, Ecker:1989yg, Messerli:2025rnv, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and conventions

This record isolates every retained tree contribution with at least one
neutral isovector $\rho^0$ response in the doubly virtual ALP TFF. It is the
low-energy rho-pole part needed as an input to the two-photon loop for
$a\to\ell^+\ell^-$. The full TFF convention is
$$
i\mathcal M=i e^2 R F_{a\gamma^*\gamma^*}(s_1,s_2),\qquad
R=\epsilon_{\mu\nu\rho\sigma}\varepsilon_1^\mu\varepsilon_2^\nu
q_1^\rho q_2^\sigma,\qquad s_i=q_i^2,
$$
with metric $(+---)$ and $\epsilon^{0123}=+1$. Use stable tree poles,
leading electromagnetic order, the original seven VJP and four VVP operators,
and diagonal ALP sources already in the no-mixing basis. Hadronic cuts,
widths, additional odd operators and short-distance matching are not included.

The sources and flavor conventions are those of the
[source-resolved TFF](alp-doubly-virtual-tff.md) and
[LO no-mixing source assembly](../../rcht/odd-parity/odd-alp-source-lsz-assembly.md).
In particular, the no-mixing coefficients obey $c_q^A=-S$ and
$c_q^P=-(c_q^{AP}+S)$, where $S=\mathbb M_a^{-1}(\mathcal C)$. Thus
$$
U_a=\frac{2c_q^A}{f_a}=-\frac{2S}{f_a},\qquad
Y_a=-\frac{2\{c_q^P,\chi\}}{f_a}
=\frac{2\{c_q^{AP}+S,\chi\}}{f_a},\qquad
\chi=\operatorname{diag}(\chi_u,\chi_d,\chi_s),\quad \chi_q=2B_0m_q.
$$
No second field rotation is made.

Assume an isospin-diagonal even vector sector without $\rho$–isoscalar
mixing. Set
$$
T_\rho=\frac{\operatorname{diag}(1,-1,0)}{\sqrt2},\quad
T_{\omega_I}=\frac{\operatorname{diag}(1,1,0)}{\sqrt2},\quad
T_{\phi_I}=\operatorname{diag}(0,0,1).
$$
For physical isoscalars use
$T_\omega=\cos\theta_VT_{\omega_I}-\sin\theta_VT_{\phi_I}$ and
$T_\phi=\sin\theta_VT_{\omega_I}+\cos\theta_VT_{\phi_I}$. Define
the signed photon responses $W_i(s)=R_i/(M_i^2-s)$. The stage-4
single-trace source gives
$$
R_\rho=\operatorname{tr}(T_\rho S_\gamma)
=\frac{F_V+2x_lz_V}{\sqrt2},\qquad x_l=\chi_+=2\chi_u=2\chi_d,
$$
in the isospin limit. The double-trace coefficient $y_V$ has zero rho trace.
More generally retain the chosen signed $R_\rho$ as an input if another
declared photon-source fit is used.

# Flavor traces and ordered numerators

For diagonal no-mixing couplings $c_q^A=\operatorname{diag}
(c_u^A,c_d^A,c_s^A)$ and $c_q^P=\operatorname{diag}
(c_u^P,c_d^P,c_s^P)$, the required traces are
$$
\begin{aligned}
A_\rho&\equiv\operatorname{tr}(T_\rho Q U_a)
=\frac{\sqrt2(2c_u^A+c_d^A)}{3f_a},&
B_\rho&\equiv\operatorname{tr}(T_\rho QY_a)
=-\frac{2\sqrt2(2\chi_uc_u^P+\chi_dc_d^P)}{3f_a},\\
C_{\rho\rho}&\equiv\operatorname{tr}(T_\rho^2U_a)
=\frac{c_u^A+c_d^A}{f_a},&
D_{\rho\rho}&\equiv\operatorname{tr}(T_\rho^2Y_a)
=-\frac{2(\chi_uc_u^P+\chi_dc_d^P)}{f_a}.
\end{aligned}
$$
For any isoscalar generator $T_I=\operatorname{diag}(t_l,t_l,t_s)$,
$$
E_{\rho I}\equiv\operatorname{tr}(T_\rho T_IU_a)
=\frac{\sqrt2t_l(c_u^A-c_d^A)}{f_a},\qquad
H_{\rho I}\equiv\operatorname{tr}(T_\rho T_IY_a)
=-\frac{2\sqrt2t_l(\chi_uc_u^P-\chi_dc_d^P)}{f_a}.
$$
The strange entries cancel because $T_\rho$ has no strange component.

Let $m_a$ be the ALP mass and let $M_V$ be the odd-operator normalization
mass. For vector leg 1 with the other leg a photon, define
$$
L(s,t)=(c_1+c_2-c_5)m_a^2
+(-c_1+c_2+c_5-2c_6)s+(c_1-c_2+c_5)t.
$$
The ordered VJP numerator and the ordered VVP numerator for vector legs $i,j$
are
$$
V_{\rho\gamma}(s,t)=\frac{2}{M_V}
\left[A_\rho L(s,t)-4c_3B_\rho\right],
$$
$$
P_{\rho j}(s,t)=4\left[
\operatorname{tr}(T_\rho T_jU_a)
\big((d_1-d_3)m_a^2+d_3(s+t)\big)
-4d_2\operatorname{tr}(T_\rho T_jY_a)\right].
$$
These are the contractions in the canonical $c_i,d_i$ convention. The
neutral $c_4,c_7,d_4$ structures do not contribute here.

# Isovector TFF

The contribution with at least one rho response is
$$
\boxed{\begin{aligned}
F^{(I=1)}_{a\gamma^*\gamma^*}(s_1,s_2)={}&
2\sqrt2\left[W_\rho(s_1)V_{\rho\gamma}(s_1,s_2)
+W_\rho(s_2)V_{\rho\gamma}(s_2,s_1)\right]\\
&-2W_\rho(s_1)W_\rho(s_2)P_{\rho\rho}(s_1,s_2)\\
&-2\sum_{I=\omega,\phi}\left[
W_\rho(s_1)W_I(s_2)P_{\rho I}(s_1,s_2)
+W_I(s_1)W_\rho(s_2)P_{I\rho}(s_1,s_2)
\right].
\end{aligned}}
$$
The rho-rho term appears once; the mixed terms retain both ordered vector
assignments. In exact isospin, $c_u^A=c_d^A$ and
$\chi_uc_u^P=\chi_dc_d^P$, so $E_{\rho I}=H_{\rho I}=0$ and the mixed
rho-isoscalar terms vanish. They must be retained for general diagonal
no-mixing sources, which need not be isospin symmetric. The VJP rho terms
remain because the electromagnetic charge matrix is not isoscalar.

The derivation is executable with
`.venv/bin/python scripts/derive_alp_doubly_virtual_tff.py`. It contracts the
rho flavor matrices, derives $R_\rho$, and assembles VJP, rho-rho and mixed
rho-isoscalar terms with the correct photon placements and ordered vector
assignments. In the loop, use the unnormalized TFF; do not divide by
$F_{a\gamma\gamma}$ when the real-photon source contributions may cancel.

# Limits

This is a tree-level rho-pole component rather than a complete isovector
spectral function. In particular, it does not implement the two-pion cut or
the short-distance completion described in the
[lepton-pair plan](alp-lepton-pair-calculation-plan.md). The stable-pole
approximation is not a precision representation near the rho resonance.
