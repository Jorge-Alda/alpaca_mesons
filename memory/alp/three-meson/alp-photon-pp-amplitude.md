---
title: Odd vector-exchange amplitude for an ALP to a photon and two pseudoscalars
status: derived
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Inputs and scope

This record derives the odd-parity vector-exchange contribution to
$a(p_a)\to\gamma(k,\varepsilon_\gamma)P_1(p_1)P_2(p_2)$ in natural units,
with metric $(+---)$ and $\epsilon^{0123}=+1$. The ALP sources are already in
Alda:2026nrd's no-mixing basis [@Alda:2026nrd]. We retain the original seven
$VJP$ and four $VVP$ operators [@Ruiz-Femenia:2003jdx], the leading even
$G_V$ $VPP$ vertex [@Ecker:1988te], and the tree vector-photon source
diagonalization. We assume electromagnetically compatible sources
$[Q,c_q^A]=[Q,c_q^P]=0$. The calculation is leading order in $e$ and sums a
complete orthonormal Hermitian vector nonet. Pseudoscalars use Alda:2026nrd's
normalization
$$
\Phi=\sum_P\lambda_P P,\qquad
\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij},\qquad
u=\exp\!\left(\frac{i\Phi}{2f_\pi}\right),\qquad
f_\pi=0.093\,{\rm GeV}.
$$
The RChT vertices are converted with $T_P=\lambda_P/\sqrt2$ and
$F=f_\pi$ before writing the amplitude in these generators. Its executable
derivation is
[`scripts/derive_alp_photon_two_pseudoscalar.py`](../../../scripts/derive_alp_photon_two_pseudoscalar.py):

```text
.venv/bin/python scripts/derive_alp_photon_two_pseudoscalar.py
.venv/bin/python scripts/derive_alp_photon_two_pseudoscalar.py --test
```

This is the odd resonance-exchange piece, not a complete physical amplitude:
independent local WZW/anomalous terms and the even scalar/vector-source
radiative contributions are not included. The full antisymmetric-tensor
propagator is used, including the local part required by representation
equivalence [@Ecker:1989yg].

# Kinematics and radiative source kernels

Set
$$
s=(p_1+p_2)^2,\qquad t=(k+p_2)^2,\qquad u=(k+p_1)^2,
\qquad p_a=k+p_1+p_2,\qquad k^2=0.
$$
For any pseudoscalar source $X$, define
$$
\mathcal U_a=\frac{2c_q^A}{f_a},\quad
\mathcal Y_a=-\frac{2}{f_a}\{c_q^P,\chi\},\qquad
\mathcal U_i=-\frac{1}{f_\pi}\lambda_i^\dagger,\quad
\mathcal Y_i=+\frac{1}{f_\pi}\{\lambda_i^\dagger,\chi\},\quad
\chi=2B_0m_q.
$$
The real-photon vector response in the fitted diagonal source convention is
$$
\mathcal W_\gamma(0)=\sum_i\frac{R_i}{M_i^2}T_i,\qquad
V_{\mu\nu}=\sqrt2e\,\mathcal W_\gamma(0)F_{\gamma\mu\nu},
$$
where $R_i$ are the signed photon-source residues and $M_i$ the physical
vector pole masses, as in the [photon-vector diagonalization](../../rcht/foundations/rcht-photon-vector-diagonalization.md).

For an exchanged vector generator $T_R$ of mass $M_R$, the radiative odd
numerator is the sum of VVP conversion, direct VJP, and (for a pseudoscalar
source) the $c_4$ VJP term:
$$
\mathcal H_R(X;r^2)=\sqrt2e\,\mathcal N_R(X;0,r^2)
-2e\,\mathcal V_{R Q}(X;r^2,0)
-\delta_{X=P_i}\frac{8e c_4}{f_\pi M_V}
\operatorname{tr}\!\left(T_R[[\lambda_i^\dagger,Q],\chi]\right).
$$
Here $r^2$ is the exchanged-vector virtuality, $M_V$ is the odd-operator
normalization scale, and
$$
\begin{aligned}
D_R(X)&=\frac12\operatorname{tr}\!\left(\{T_R,\mathcal W_\gamma(0)\}\mathcal U_X\right),&
P_R(X)&=\frac12\operatorname{tr}\!\left(\{T_R,\mathcal W_\gamma(0)\}\mathcal Y_X\right),\\
D_{RQ}(X)&=\frac12\operatorname{tr}\!\left(\{T_R,Q\}\mathcal U_X\right),&
P_{RQ}(X)&=\frac12\operatorname{tr}\!\left(\{T_R,Q\}\mathcal Y_X\right),
\end{aligned}
$$
$$
\begin{aligned}
\mathcal N_R(X;0,r^2)&=4\left[d_1m_X^2D_R(X)
-d_3(m_X^2-r^2)D_R(X)-4d_2P_R(X)\right],\\
\mathcal V_{RQ}(X;r^2,0)&=\frac4{M_V}\left[
D_{RQ}(X)\left\{\frac{m_X^2-r^2}{2}c_1
+\left(\frac{m_X^2+r^2}{2}\right)c_2
-\frac{m_X^2-r^2}{2}c_5-r^2c_6\right\}
-2c_3P_{RQ}(X)\right].
\end{aligned}
$$
The direct photon leg has $q_\gamma^2=0$. The $c_4$ term vanishes for the
ALP insertion at the vacuum but contributes to the crossed pseudoscalar
vertices when $[Q,\lambda_i^\dagger]\ne0$ and the quark masses break the
corresponding flavor symmetry. Omitting it would miss a charged-channel
term. The $c_7$ and $d_4$ contractions vanish for these routed transverse
currents.

# Amplitude

With $\mathcal M=\epsilon_{\mu\nu\rho\sigma}
\varepsilon_\gamma^\mu p_1^\nu p_2^\rho k^\sigma\,\mathcal A_{\gamma PP}$,
the result is
$$
\begin{aligned}
\mathcal A_{\gamma PP}={}&\sum_{T_R}\Bigg[
-\frac{G_V\,\mathcal H_R(a;s)
\operatorname{tr}(T_R[\lambda_1^\dagger,\lambda_2^\dagger])}
{\sqrt2 f_\pi^2(s-M_R^2+i0)}\\
&+\frac{\sqrt2G_V\,\mathcal H_R(P_2;t)
\operatorname{tr}(T_R[\lambda_1^\dagger,c_q^A])}
{f_\pi f_a(t-M_R^2+i0)}
-\frac{\sqrt2G_V\,\mathcal H_R(P_1;u)
\operatorname{tr}(T_R[\lambda_2^\dagger,c_q^A])}
{f_\pi f_a(u-M_R^2+i0)}\Bigg].
\end{aligned}
$$
The first graph is $a\to T_R^*\gamma$ followed by $T_R^*\to P_1P_2$;
its $\mathcal H_R(a;s)$ contains both direct VJP photon insertion and VVP
conversion of one vector leg. The other two graphs are the crossed
$a\to P_iT_R^*$ and $T_R^*\to P_j\gamma$ exchanges. The $\mathcal H_R(P_i)$
include the same VJP/VVP terms and the charged-source $c_4$ term. The sum is
over every vector nonet state selected by its traces, with tree denominators
shown; near a pole, replace these by one explicitly declared dressed
propagator prescription throughout the amplitude.

# Checks and limitations

The derivation obtains the stationary response
$r_i=\sqrt2eR_iF_\gamma/M_i^2$, verifies the photon Ward identity for the
operator contractions, checks the four original-fit coefficient null
directions, and tests the charged pion/kaon $c_4$ mass-difference structure.
The four fit-null projections vanish at $q_\gamma^2=0$ in this retained
amplitude. These are symbolic checks of the stated tree truncation, not a
complete higher-order uncertainty analysis. Widths, loops, rescattering,
additional odd operators, local anomalous contacts, and the even
$c_q^{VS}$-dependent radiative amplitude remain outside this record.
