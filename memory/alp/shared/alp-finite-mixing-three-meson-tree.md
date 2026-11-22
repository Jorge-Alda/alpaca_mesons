---
title: Exact tree amplitude of A -> P P P at finite f_pi/f_a in the dressed frame, in M^{-1}(C) and the invariants
status: mixed
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Exact tree amplitude of $A\to P_{n_1}P_{n_2}P_{n_3}$ at finite $f_\pi/f_a$

This record is task 2.1 of the OpenSpec change `estimate-loop-tree-crossing-finite-mixing`. It gives the tree term of the observable class chosen in the [selection record](alp-loop-observable-selection-a-to-3pi0.md) (section "Observable class used for the estimate"): the on-shell amplitude of the ALP-like pole $A$ into three pseudoscalar poles, exact in $\epsilon=f_\pi/f_a$, in matrix notation. The formula is verified exactly against the [pole solver](alp-lo-pole-solver.md). No loop is evaluated.

## Expansion order kept

The amplitude is exact in $\epsilon$: no expansion is used. The pole masses are the exact roots of $D_a(s)=0$, or the unmixed masses for poles with zero ALP component. An expansion is used only for comparison (section "Difference from the leading expansion").

## Conventions and units

- Units GeV, squared masses and $(p_{n_i}+p_{n_j})^2$ in GeV$^2$. Notation of [@Alda:2026nrd]: $U=\exp(i\Phi/f_\pi)$, $\Phi=P_i\lambda_i$, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, $\langle\cdot\rangle$ the flavour trace, $N_f=3$, $x=a/f_a$, $\epsilon=f_\pi/f_a$, $f\equiv f_\pi$ in the formulas below.
- Dressed frame ([derivative-coupling removal](alp-dressed-frame-derivative-coupling-removal.md)): $K=\mathbf 1$, $\chi=2B_0\widehat M(x)$, $\widehat M(x)=\sum_n\frac{x^n}{n!}(\widehat S_n+i\widehat P_n)$, $\widehat\Theta(x)=\sum_n\frac{x^n}{n!}\widehat\Theta_n$, and
$$
\mathcal L=\frac{f^2}{4}\langle u_\mu u^\mu\rangle+\frac{f^2}{4}\langle\chi_+\rangle-\frac{f^2m_0^2}{4N_f}\Big(\frac{\langle\Phi\rangle}{f}+2\widehat\Theta(x)\Big)^2+\frac12(\partial a)^2-\frac12m_{a,0}^2a^2 .
$$
- Vertices are derivatives of $\mathcal L$ without the factor $i$, all momenta incoming. The amplitude is the contact term minus the three exchanges, with no symmetry factor for identical daughters ([pole-projected record](alp-lo-pole-projected-tree-amplitudes.md)); the $1/3!$ for $3\pi^0$ belongs in the rate.
- $\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1$; $\mathbb M_s(X)=\frac12\{\chi_0,X\}+\frac{m_0^2}{N_f}\langle X\rangle\mathbf 1-sX$ with $\chi_0=2B_0\widehat S_0$; $\mathbb M_s^{-1}$ its inverse ([operator-inverse record](alp-lo-finite-mixing-operator-inverse.md), with $s$ the squared-mass argument); $Y=\mathbb M^{-1}_{m_A^2}(\mathcal C)$ as in the [candidate record](alp-loop-observable-candidates-comparison.md).
- Pole vectors $v^{(n)}$ with $v^{(m)\mathsf T}v^{(n)}=\delta_{mn}$ ([pole solver](alp-lo-pole-solver.md)); $v^{(n)}_a$ is the ALP component.
- **Symbol pending approval:** $\Phi_\ell\equiv\sum_kv^{(\ell)}_k\lambda_k$, the meson part of the pole vector of leg $\ell\in\{A,n_1,n_2,n_3\}$ in matrix form. It is only a display shorthand; every formula below holds with the sum written out. $\Phi_{n,0}\equiv\sum_kv^{(n)}_k|_{\epsilon=0}\lambda_k$ is the unmixed matrix.

## Assumptions

1. LO $U(3)$ Lagrangian above, tree level, stable real poles, non-degenerate poles. $V_{\rm ext}=\frac12m_{a,0}^2a^2$ only (no other pure-ALP potential).
2. Basis conditions: $\widehat S_0=m_q$ diagonal (three different entries allowed, $m_u\neq m_d$ not excluded), $\widehat P_0=0$, $\widehat\Theta_0=0$.
3. $\langle\widehat S_1\rangle=0$ (no ALP tadpole; vacuum alignment). **Flagged assumption.** $\widehat S_1$ itself is kept general.
4. All invariants are Hermitian and commute with the quark charge matrix (a $d$-$s$ entry is allowed, a $u$-$d$ entry is not). No value, representative or relation among them is assumed; in particular $\langle\lambda_3\widehat P_1\rangle$ is not set to zero.
5. Sign of each leg (a convention, not physics): $A$ has positive ALP component (solver rule); each meson-like daughter is signed by continuity with its unmixed matrix, i.e. $\langle\Phi_n\Phi_{n,0}\rangle>0$, with $\Phi_{n,0}$ signed by the solver rule (largest component in $(\lambda_3,\lambda_8,\lambda_0)$ positive). Outgoing $\pi^\pm$ use $\sqrt2E_{21}$ and $\sqrt2E_{12}$.
6. Exceptions: a pole with zero ALP component (charged mesons, the combination orthogonal to a $d$-$s$ entry of $\mathcal C$, or an accidental zero overlap of $\mathcal C$ with a neutral unmixed direction) uses its unmixed matrix, $\Phi_n=\Phi_{n,0}$, $v^{(n)}_a=0$. Exchanged momenta $(p_{n_j}+p_{n_k})^2$ must stay away from the poles.

## Inputs

- Operator $\mathbb M_s$, its explicit inverse, $D_a(s)$, $N(s)=D_a'(s)$: [operator-inverse record](alp-lo-finite-mixing-operator-inverse.md), steps 2–7, specialised to $K=\mathbf 1$ ($c_q^A\to0$).
- Dressed frame and on-shell equality of tree amplitudes between frames: [derivative-coupling removal](alp-dressed-frame-derivative-coupling-removal.md).
- Contact plus exchange construction: [pole-projected record](alp-lo-pole-projected-tree-amplitudes.md).
- Leading term of $A\to3\pi^0$ through $\epsilon^5$ for $m_u=m_d$: [candidate record](alp-loop-observable-candidates-comparison.md), candidate 4.

## Result

**Quadratic part (verified from the script's quadratic blocks).** $M^2_{ka}=-\epsilon\langle\mathcal C\lambda_k\rangle$ and
$$
M^2_{aa}=m_{a,0}^2-\epsilon^2B_0\langle\widehat S_2\rangle+2\epsilon^2\frac{m_0^2}{N_f}\widehat\Theta_1^2,\qquad
D_a(s)=s-M^2_{aa}+2\epsilon^2\big\langle\mathcal C\,\mathbb M^{-1}_s(\mathcal C)\big\rangle,\qquad
N(s)=1+2\epsilon^2\big\langle\mathbb M^{-1}_s(\mathcal C)^2\big\rangle .
$$

**Pole vectors.** Every pole $n$ with nonzero ALP component is a root of $D_a$, and
$$
\Phi_A=\frac{2\epsilon Y}{\sqrt{N(m_A^2)}},\quad v^{(A)}_a=\frac{1}{\sqrt{N(m_A^2)}};\qquad
\Phi_n=\operatorname{sgn}\big\langle\epsilon\,\mathbb M^{-1}_{m_n^2}(\mathcal C)\Phi_{n,0}\big\rangle\frac{2\epsilon\,\mathbb M^{-1}_{m_n^2}(\mathcal C)}{\sqrt{N(m_n^2)}},\quad
v^{(n)}_a=\operatorname{sgn}\big\langle\epsilon\,\mathbb M^{-1}_{m_n^2}(\mathcal C)\Phi_{n,0}\big\rangle\frac{1}{\sqrt{N(m_n^2)}} .
$$

**Amplitude.** With $(i,j,k)$ running over the cyclic orderings of the daughters $(1,2,3)$, $q^2=(p_{n_j}+p_{n_k})^2$, and sums over legs $\ell,m\in\{A,n_1,n_2,n_3\}$:
$$
\mathcal M(A\to P_{n_1}P_{n_2}P_{n_3})=\mathcal V_{\rm kin}+\mathcal V_{PPPP}+\mathcal V_{PPPa}+\mathcal V_{PPaa}+\mathcal V_{Paaa}+\mathcal V_{aaaa}-\sum_{i=1}^3E_i ,
$$
$$
\mathcal V_{\rm kin}=-\frac{1}{24f^2}\sum_{i}\Big[(p_{n_j}+p_{n_k})^2-\tfrac12\big(m_A^2+m_{n_1}^2+m_{n_2}^2+m_{n_3}^2\big)\Big]\Big(\big\langle[\Phi_{n_j},\Phi_A][\Phi_{n_k},\Phi_{n_i}]\big\rangle+\big\langle[\Phi_{n_k},\Phi_A][\Phi_{n_j},\Phi_{n_i}]\big\rangle\Big),
$$
$$
\mathcal V_{PPPP}=\frac{B_0}{24f^2}\sum_{\sigma\in S_4}\big\langle\widehat S_0\Phi_{\sigma_1}\Phi_{\sigma_2}\Phi_{\sigma_3}\Phi_{\sigma_4}\big\rangle,\qquad
\mathcal V_{PPPa}=-\frac{\epsilon B_0}{6f^2}\sum_\ell v^{(\ell)}_a\sum_{S_3\ (\text{other three legs})}\big\langle\widehat P_1\Phi\Phi\Phi\big\rangle,
$$
$$
\mathcal V_{PPaa}=-\frac{\epsilon^2B_0}{2f^2}\sum_{\ell<m}v^{(\ell)}_av^{(m)}_a\big\langle\widehat S_2\{\Phi_{k'},\Phi_{l'}\}\big\rangle_{\{k',l'\}=\text{other two legs}},\qquad
\mathcal V_{Paaa}=\frac{\epsilon^3}{f^2}\sum_\ell\Big(\prod_{m\neq\ell}v^{(m)}_a\Big)\Big\langle\Big(B_0\widehat P_3-\frac{m_0^2}{N_f}\widehat\Theta_3\Big)\Phi_\ell\Big\rangle,
$$
$$
\mathcal V_{aaaa}=\frac{\epsilon^4}{f^2}\Big[B_0\langle\widehat S_4\rangle-\frac{m_0^2}{N_f}\big(8\widehat\Theta_1\widehat\Theta_3+6\widehat\Theta_2^2\big)\Big]\prod_\ell v^{(\ell)}_a .
$$
The exchanges use the cubic vertex with two external legs $\ell,m$ and one internal line, in matrix form for the meson line and as a number for the ALP line:
$$
\sum_k\mathcal V_{\ell mk}\lambda_k=-\frac{\epsilon B_0}{f}\big(v^{(\ell)}_a\{\widehat S_1,\Phi_m\}+v^{(m)}_a\{\widehat S_1,\Phi_\ell\}\big)+\frac{2\epsilon^2}{f}v^{(\ell)}_av^{(m)}_a\Big(B_0\widehat P_2-\frac{m_0^2}{N_f}\widehat\Theta_2\mathbf 1\Big),
$$
$$
\mathcal V_{\ell ma}=-\frac{\epsilon B_0}{2f}\big\langle\widehat S_1\{\Phi_\ell,\Phi_m\}\big\rangle+\frac{\epsilon^2}{f}\Big(v^{(\ell)}_a\Big\langle\Big(B_0\widehat P_2-\frac{m_0^2}{N_f}\widehat\Theta_2\Big)\Phi_m\Big\rangle+(\ell\leftrightarrow m)\Big)+\frac{\epsilon^3}{f}\Big(B_0\langle\widehat S_3\rangle-6\frac{m_0^2}{N_f}\widehat\Theta_1\widehat\Theta_2\Big)v^{(\ell)}_av^{(m)}_a ,
$$
and, with the exact propagator in Schur form (meson block $-\mathbb M^{-1}_{q^2}$, plus the direction $(2\epsilon\,\mathbb M^{-1}_{q^2}(\mathcal C),1)$ divided by $D_a(q^2)$),
$$
E_i=-\frac12\Big\langle\Big(\sum_k\mathcal V_{An_ik}\lambda_k\Big)\mathbb M^{-1}_{q^2}\Big(\sum_k\mathcal V_{n_jn_kk}\lambda_k\Big)\Big\rangle
+\frac{\Big(\epsilon\big\langle\sum_k\mathcal V_{An_ik}\lambda_k\,\mathbb M^{-1}_{q^2}(\mathcal C)\big\rangle+\mathcal V_{An_ia}\Big)\Big(\epsilon\big\langle\sum_k\mathcal V_{n_jn_kk}\lambda_k\,\mathbb M^{-1}_{q^2}(\mathcal C)\big\rangle+\mathcal V_{n_jn_ka}\Big)}{D_a(q^2)} .
$$
There is no three-meson vertex at LO. All cubic vertices contain CP-odd invariants ($\widehat S_1$, $\widehat P_2$, $\widehat\Theta_2$, $\widehat S_3$), so for CP-conserving sources the exchanges vanish. For flavour-diagonal invariants the neutral matrices commute and $\mathcal V_{\rm kin}=0$; it is nonzero for charged daughters or a $d$-$s$ entry.

## Difference from the leading expansion

Expansion at fixed $m_A^2$ (then $m_{a,0}^2$ adjusts with $\epsilon$), fixed invariants and fixed $(p_{n_1}+p_{n_2})^2$, $(p_{n_2}+p_{n_3})^2$; the third invariant follows from the exact daughter masses.

- **Oddness (verified).** $a$ enters only through $x=\epsilon a/f$ and $V_{\rm ext}$ is even in $a$, so with the sign convention of assumption 5 $\mathcal M(-\epsilon)=-\mathcal M(\epsilon)$. The coefficients of $\epsilon^0,\epsilon^2,\epsilon^4$ vanish exactly at the verification inputs.
- **Leading term, coefficient of $\epsilon$ (verified).** $\mathcal V_{\rm kin}+\mathcal V_{PPPP}$ evaluated at $\Phi_A\to2Y$, $\Phi_n\to\Phi_{n,0}$, $m_n\to m_{n,0}$, minus $\frac{B_0}{6f^2}\sum_{S_3}\langle\widehat P_1\Phi_{n_1,0}\Phi_{n_2,0}\Phi_{n_3,0}\rangle$. For flavour-diagonal invariants and three equal daughters $\pi$ it is $\frac{B_0}{f^2}\big\langle(2Y\widehat S_0-\widehat P_1)\Phi_{\pi,0}^3\big\rangle$, the general form of the candidate-record result (which assumed $m_u=m_d$).
- **Exact minus leading is $O(\epsilon^3)$.** It collects: the normalisations $N(m_n^2)$, the ALP components of the daughters, the shift of the daughter poles from $m_{n,0}^2$, the $\epsilon$ dependence of $\mathbb M^{-1}_{m_n^2}(\mathcal C)$ for the daughters, the classes $\mathcal V_{PPaa}$ ($\epsilon^3$), $\mathcal V_{Paaa}$ ($\epsilon^5$), $\mathcal V_{aaaa}$ ($\epsilon^7$), and the exchanges (first nonzero at $\epsilon^3$). Its $\epsilon^3$ coefficient is given only at the verification inputs; no general closed form of it is derived here.
- **Reference point of $Y$ (verified).** Evaluating $Y$ at $m_{a,0}^2$ (as in the candidate record) instead of $m_A^2$ changes the leading term at $O(\epsilon^3)$ (coefficient $-1.50469$ at the inputs below).

## Verification

Exact rational inputs (checks only, not physical values): $B_0=2$, $\widehat S_0=m_q=\operatorname{diag}(1,6,27)/200$, $m_0^2=12/25$, $f=1/10$; $\widehat S_1=\operatorname{diag}(\frac13,-\frac15,-\frac2{15})$, $\widehat S_2=\operatorname{diag}(\frac12,-\frac14,\frac17)$, $\widehat S_3=\operatorname{diag}(-\frac13,\frac25,\frac16)$, $\widehat S_4=\operatorname{diag}(\frac14,\frac19,-\frac12)$; $\widehat P_1=\operatorname{diag}(\frac15,-\frac13,\frac12)$, $\widehat P_2=\operatorname{diag}(\frac16,\frac14,-\frac15)$, $\widehat P_3=\operatorname{diag}(-\frac12,\frac13,\frac18)$; $\widehat\Theta_{1,2,3}=\frac13,-\frac14,\frac15$. Unmixed neutral masses $3/50$, $7/25$, $41/50$ GeV$^2$; $m_A^2=2$ GeV$^2$; $(p_{n_1}+p_{n_2})^2=1/2$, $(p_{n_2}+p_{n_3})^2=3/5$ GeV$^2$. All contact classes and all exchanges are nonzero at these inputs (except $\mathcal V_{\rm kin}$ for neutral daughters).

Method: at fixed rational sources and $m_A^2$, $D_a(m_A^2)=0$ fixes $m_{a,0}^2$ linearly, and $D_a(m_A^2)-D_a(m_\pi^2)=0$ for a chosen rational pion-like pole $m_\pi^2$ fixes $\epsilon^2$ linearly ($M^2_{aa}$ cancels). Then $\det D$ has the rational roots $m_A^2$, $m_\pi^2$ and [`solve_spectrum_exact`](../../../scripts/solve_alp_meson_spectrum.py) works unchanged; $\epsilon=\sqrt{\epsilon^2}$. The solver-side amplitude uses the solver pole vectors, the CCWZ vertices of the existing scripts and the direct $D(q^2)^{-1}$; it shares no formula with the closed expression. Equality is exact ($\mathcal M\sqrt{N(m_A^2)N(m_\pi^2)^3}$ is rational at a tuned point).

$A\to\pi\pi\pi$ (three equal pion-like daughters), closed formula $=$ solver exactly at all three points:

| $m_\pi^2$ (GeV$^2$) | $\epsilon^2$ | $\epsilon$ | $\mathcal M$ | coefficient-of-$\epsilon$ term | (exact $-$ leading)$/\epsilon^3$ |
|---|---|---|---|---|---|
| 59/1000 | 372487446981/149661968053600 | 0.0499 | $-5.002868$ | $-4.946784$ | $-451.6876$ |
| 17/300 | 11327230581/1401437592208 | 0.0899 | $-9.244556$ | $-8.914517$ | $-454.1920$ |
| 1/20 | 3922412571/173274434128 | 0.1505 | $-16.484429$ | $-14.918743$ | $-459.7022$ |

- Coefficient of $\epsilon$: $-156401432320\sqrt{418}/32248240827=-99.15688$, equal to the series coefficient of the closed formula exactly.
- Coefficient of $\epsilon^3$ from the series of the closed formula: $-450.53685$ (exact value printed by the script). The column above approaches it; the remainder after subtracting it, divided by $\epsilon^2$, is about $462$, $452$, $405$, consistent with $O(\epsilon^5)$.
- Relative difference exact$-$leading: $1.1\%$, $3.6\%$, $9.5\%$ of the exact amplitude at the three points.
- $A\to\pi^+\pi^-\pi^0$ at $\epsilon^2=11327230581/1401437592208$: $\mathcal M=-1.618615$, exact agreement; here $\mathcal V_{\rm kin}$ and all three exchanges are nonzero.
- $d$-$s$ entry in the invariants, at $\epsilon^2=77557547788107/10006444089585776$: $\mathcal M=-9.299677$, exact agreement. The two tuned poles are rational; their pole vectors are the $K$-normalised null vectors of $D$ (the solver algorithm), since the other poles are not found in closed form by `solve_spectrum_exact`.
- Also checked: the contact and cubic vertex lists and the absence of a three-meson vertex against the CCWZ vertices; $M^2_{aa}$ and $M^2_{ka}$; $N=D_a'$; oddness in $\epsilon$; vanishing even coefficients.

Not verified: two different neutral daughters in one amplitude (the tuning has no free input left; the formula is unchanged); poles near degeneracy.

[`scripts/derive_alp_finite_mixing_three_meson_tree.py`](../../../scripts/derive_alp_finite_mixing_three_meson_tree.py); tests in [`scripts/test_derive_alp_finite_mixing_three_meson_tree.py`](../../../scripts/test_derive_alp_finite_mixing_three_meson_tree.py) (15 tests). The existing suites `test_solve_alp_meson_spectrum.py` (10) and `test_derive_alp_dressed_frame_unitarity.py` (27) pass.

```sh
.venv/bin/python scripts/derive_alp_finite_mixing_three_meson_tree.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_finite_mixing_three_meson_tree.py
```
