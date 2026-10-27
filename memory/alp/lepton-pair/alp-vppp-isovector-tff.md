---
title: Local VPPP contribution to the ALP isovector TFF production wave
status: derived
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Alda:2026nrd, GomezDumm:2010ra, Ecker:1988te, Ecker:1989yg, Messerli:2025rnv]
provenance: analytical-derivation
---

# Inputs and assumptions

Use the five $g_i\mathcal O^i_{VPPP}/M_V$ operators in the
[source extract](../../extracts/GomezDumm:2010ra/vppp-basis-and-omega-three-pion.md)
[@GomezDumm:2010ra], and the canonical epsilon and odd-coupling convention in
[the convention record](../../rcht/odd-parity/levi-civita-convention.md).
ALP couplings are the supplied no-mixing $c_q^A,c_q^P$ of
[Alda:2026nrd's expansion](../../extracts/Alda:2026nrd/appendix-a-expansion.md).
They are independent inputs. This calculation does not replace them by pion
couplings. Natural units, tree level, one ALP, two on-shell charged pions,
and diagonal axial and pseudoscalar sources are assumed. No equality of the
light-quark spurions is imposed.
It uses no numerical $g_i$ prescription or resonance widths.

Specifically, every $c_q^A,c_q^P$ below denotes the **already transformed**
$\widehat c_q^A,\widehat c_q^P$; hats are suppressed exactly as in the
existing pion-production record. In terms of the original UV sources and
the profile $S=\mathbb M_a^{-1}(C)$ in
[the LSZ assembly record](../../rcht/odd-parity/odd-alp-source-lsz-assembly.md),
substitute $c_q^A=-S$ and $c_q^P=-c_q^{AP}-S$. Then the source matrices
below are precisely $\mathcal U_a=-2S/f_a$ and
$\mathcal Y_a=2\{c_q^{AP}+S,\chi\}/f_a$. The script also prints this
substituted form. No additional meson LSZ insertion is to be added. Using
untransformed UV $c_q^A,c_q^P$ in the displayed formula would instead give
only the direct contribution and would omit meson mixing.

The photon response is the full neutral-vector matrix
$$
W_\gamma(y)=\sum_i\frac{R_iT_i}{M_i^2-y-i0},
$$
with the signed residues and normalization of
[the photon-source calculation](../../rcht/odd-parity/isoscalar-vppp-photon-response.md)
[@Ecker:1988te; @Ecker:1989yg]. In the explicit charged-channel evaluation
The photon response may contain rho and isoscalar vectors. The isovector
label refers to the charged two-pion cut of the other current.

# Source insertion and derivation

The executable [SymPy calculation](../../../scripts/derive_alp_vppp_isovector_tff.py)
uses $u_\mu=i p_\mu\mathcal U$ and $\chi_-=-i\mathcal Y$ on each
outgoing field. Its source matrices are
$$
\mathcal U_a=\frac{2c_q^A}{f_a},\quad
\mathcal Y_a=-\frac{2\{c_q^P,\chi\}}{f_a},\quad
\mathcal U_{\pi^\pm}=-\frac{\lambda_\pm^\dagger}{f_\pi},\quad
\mathcal Y_{\pi^\pm}=\frac{\{\lambda_\pm^\dagger,\chi\}}{f_\pi}.
$$
Here $\chi=\operatorname{diag}(\chi_u,\chi_d,\chi_s)$; no isospin equality
is imposed.
The six field assignments in each operator include an ALP in either ordinary
vielbein or in $h_{\mu\nu}=\nabla_\mu u_\nu+\nabla_\nu u_\mu$.
For $g_4,g_5$ they also include the ALP in $\chi_-$, with the two
vielbeins producing the charged pions.

The $c_q^V$ contribution begins as
$-i(\partial_\mu a)[\Phi,c_q^V]/(f_\pi f_a)$ in $u_\mu$.
It already contains an ALP and a meson. Each of the other two building blocks
of a $VPPP$ operator contains at least one pseudoscalar field, so this term
cannot contribute at $a\pi\pi$ order. This field count also excludes
nonlinear vielbeins, chiral connections, and nonlinear $\chi_-$ terms at
this order. It does not require setting $c_q^V$ to zero.

Take outgoing momenta $(p_+,p_-,-P)$ on the three pseudoscalar slots and
$k=P-p_+-p_-$ on the photon, so
$$
s=(p_++p_-)^2,\quad t=(k+p_-)^2,\quad u=(k+p_+)^2,\quad k^2=y,
\qquad s+t+u=m_a^2+y+2m_{\pi^\pm}^2.
$$
The script contracts the explicit charged generators and all derivative
slots before imposing this last relation. It then inserts
$V_i^{\mu\nu}=\sqrt2 eR_i F_\gamma^{\mu\nu}/(M_i^2-y-i0)$.
For plane waves $\partial_\mu\to+i p_\mu$,
$F_\gamma^{\mu\nu}=i(k^\mu\varepsilon^\nu-k^\nu\varepsilon^\mu)$;
relative to the intermediate $q^\mu\varepsilon^\nu-q^\nu\varepsilon^\mu$
with $q=-k$, this supplies $-i$. It removes the $i$ multiplying the
external-vector epsilon scalar. Crossing gives
$\epsilon(p_+,p_-,-P,\varepsilon)=\epsilon(\varepsilon,p_+,p_-,k)$.
The local external-vector vertex retains its $1/M_i$ tensor normalization;
the photon-source insertion cancels it in the response.

# Result: local vector vertex and photon conversion

Use the charged-pion generators $\lambda_\pm^\dagger$ and the supplied
no-mixing ALP source matrices $\mathcal U_a,\mathcal Y_a$. Let $T_i$ be an
arbitrary diagonal neutral-vector generator. No ideal-mixing, tracelessness,
or $\chi_u=\chi_d$ condition is imposed. Define
$$
\mathcal U_a=\frac{2c_q^A}{f_a},\qquad
\mathcal Y_a=-\frac{2\{c_q^P,\chi\}}{f_a},\qquad
\chi=\operatorname{diag}(\chi_u,\chi_d,\chi_s).
$$

For the crossed production process $a(p_a)+V_i(q_i)\to
\pi^+(p_+)\pi^-(p_-)$, define $s=(p_++p_-)^2=(p_a+q_i)^2$ and
$y=q_i^2$. In the real-epsilon convention of the linked pion-production
record, write
$$
\mathcal M(aV_i\to\pi^+\pi^-)
=i\epsilon_{\mu\nu\rho\sigma}
\varepsilon_i^\mu p_+^\nu p_-^\rho q_i^\sigma\,
\mathcal F^{VPPP}_{aV_i\pi\pi}(s,y).
$$
The local scalar is
$$
\begin{aligned}
\mathcal F^{VPPP}_{aV_i\pi\pi}(s,y)
=\frac{2}{M_V M_i f_\pi^2}\Big\{&
\Big[(g_1-g_3)(m_a^2+4m_{\pi^\pm}^2-s-y)
+2g_2(y-m_a^2-s)\Big]
\operatorname{tr}\!\left(T_i\mathcal U_a
[\lambda_+^\dagger,\lambda_-^\dagger]\right)\\
&+\Big[(g_1-g_3)(y-m_a^2-s)
+2g_2(2m_{\pi^\pm}^2-s)\Big]
\operatorname{tr}\!\left(T_i
(\lambda_+^\dagger\mathcal U_a\lambda_-^\dagger
-\lambda_-^\dagger\mathcal U_a\lambda_+^\dagger)\right)\\
&+g_4\operatorname{tr}\!\Big(
\{T_i,[\lambda_+^\dagger,\lambda_-^\dagger]\}\mathcal Y_a
+\{T_i,[\lambda_+^\dagger,\mathcal U_a]\}
 \{\lambda_-^\dagger,\chi\}
-\{T_i,[\lambda_-^\dagger,\mathcal U_a]\}
 \{\lambda_+^\dagger,\chi\}\Big)\\
&+g_5\operatorname{tr}\!\Big(
(\lambda_+^\dagger T_i\lambda_-^\dagger
-\lambda_-^\dagger T_i\lambda_+^\dagger)\mathcal Y_a\\
&\hspace{14mm}
+(\lambda_+^\dagger T_i\mathcal U_a
-\mathcal U_aT_i\lambda_+^\dagger)
 \{\lambda_-^\dagger,\chi\}\\
&\hspace{14mm}
-(\lambda_-^\dagger T_i\mathcal U_a
-\mathcal U_aT_i\lambda_-^\dagger)
 \{\lambda_+^\dagger,\chi\}\Big)\Big\}.
\end{aligned}
$$
The factor $1/M_i$ is the external antisymmetric-tensor vector
normalization; set $y=M_i^2$ for an on-shell vector. The first two terms
come from the ALP in the vielbeins, including $h_{\mu\nu}$. The last two
terms include the ALP insertion in $\chi_-$ and the pion mass-source
insertions.

Insert vector--photon mixing only after forming the local vector vertex:
replace $T_i/M_i$ by $\sqrt2W_\gamma(y)$, with
$$
W_\gamma(y)=\sum_i\frac{R_iT_i}{M_i^2-y-i0}.
$$
Thus the charge-stripped current partial wave is
$$
\begin{aligned}
\mathcal F^{VPPP}_{a\pi\pi\gamma^*}(s,y)
=\frac{2\sqrt2}{M_V f_\pi^2}\Big\{&
\Big[(g_1-g_3)(m_a^2+4m_{\pi^\pm}^2-s-y)
+2g_2(y-m_a^2-s)\Big]
\operatorname{tr}\!\left(W_\gamma\mathcal U_a
[\lambda_+^\dagger,\lambda_-^\dagger]\right)\\
&+\Big[(g_1-g_3)(y-m_a^2-s)
+2g_2(2m_{\pi^\pm}^2-s)\Big]
\operatorname{tr}\!\left(W_\gamma
(\lambda_+^\dagger\mathcal U_a\lambda_-^\dagger
-\lambda_-^\dagger\mathcal U_a\lambda_+^\dagger)\right)\\
&+g_4\operatorname{tr}\!\Big(
\{W_\gamma,[\lambda_+^\dagger,\lambda_-^\dagger]\}\mathcal Y_a
+\{W_\gamma,[\lambda_+^\dagger,\mathcal U_a]\}
 \{\lambda_-^\dagger,\chi\}
-\{W_\gamma,[\lambda_-^\dagger,\mathcal U_a]\}
 \{\lambda_+^\dagger,\chi\}\Big)\\
&+g_5\operatorname{tr}\!\Big(
(\lambda_+^\dagger W_\gamma\lambda_-^\dagger
-\lambda_-^\dagger W_\gamma\lambda_+^\dagger)\mathcal Y_a\\
&\hspace{14mm}
+(\lambda_+^\dagger W_\gamma\mathcal U_a
-\mathcal U_aW_\gamma\lambda_+^\dagger)
 \{\lambda_-^\dagger,\chi\}\\
&\hspace{14mm}
-(\lambda_-^\dagger W_\gamma\mathcal U_a
-\mathcal U_aW_\gamma\lambda_-^\dagger)
 \{\lambda_+^\dagger,\chi\}\Big)\Big\}.
\end{aligned}
$$
The expression is independent of the pion scattering angle, so its reduced
$P$ wave equals this scalar. Unequal $\chi_u,\chi_d$ remain explicit in the
mass-source anticommutators; the pion-mass relation above is only the
on-shell kinematic input.

# Contribution to the TFF

This is a local production vertex followed by photon--vector mixing.
It enters the isovector TFF through the charged-pion unitarity cut:
$$
\operatorname{Im}_s\Delta F_{a\gamma^*\gamma^*}(s,y)
=\frac{s\sigma_{\pi^\pm}^3(s)}{96\pi}
[F_\pi^V(s)]^*\mathcal F^{VPPP}_{a\pi\pi\gamma^*}(s,y),
\qquad \sigma_{\pi^\pm}(s)=\sqrt{1-4m_{\pi^\pm}^2/s}.
$$
The production wave must be completed consistently with pion rescattering
before this is interpreted as a physical unitarity relation. It is not a
separate tree $a\gamma^*\gamma^*$ vertex.
In the per-ordering double-spectral convention of
[the spectral record](alp-isovector-tff-spectral-construction.md)
[@Messerli:2025rnv], its contribution is
$$
\Delta\rho_a^{I=1}(s,y)=\frac{s\sigma_{\pi^\pm}^3(s)}{192\pi}
\operatorname{Im}\!\left\{[F_\pi^V(s)]^*
\mathcal F^{VPPP}_{a\pi\pi\gamma^*}(s,y)\right\}.
$$
Insert this into that record's double dispersion and add the exchanged photon
ordering once. Local polynomial terms obtained when a numerator proportional
to $y$ is divided by $M_i^2-y$ can also need single-spectral/subtraction
terms. Thus the double-spectral integral alone need not reconstruct the
complete contribution. Stable poles, rescattered spectral functions, and
subtractions must be matched before being combined.

No numerical TFF follows without $g_i$, phase-compatible production-wave
rescattering, vector residues and a subtraction/high-energy prescription.
The earlier eleven-operator $c_i,d_i$ fit does not determine the independent
$g_i$ couplings.

# Verification

Run:

    .venv/bin/python scripts/derive_alp_vppp_isovector_tff.py
    .venv/bin/python scripts/derive_alp_vppp_isovector_tff.py --test

Five focused tests pass. They compare the compact vector and photon-mixed
trace formulas with the full six-assignment expansion for arbitrary diagonal
$T_i,\mathcal U_a,\mathcal Y_a,\chi$, including $\chi_u\ne\chi_d$, and verify
the linear vector-response insertion and angular independence. They verify
the local tree calculation and its conventions; they do not supply physical
rescattering or a numerical spectral density.
