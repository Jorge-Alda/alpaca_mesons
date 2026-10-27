---
title: Isoscalar ALP production through local vector three-pion contacts
status: derived
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Alda:2026nrd, Ecker:1989yg, GomezDumm:2010ra, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope, conventions, and inputs

This derives the charged isoscalar $V_i\to\pi^+\pi^-\pi^0$ endpoint from
the five local operators in the [source extract](../../extracts/GomezDumm:2010ra/vppp-basis-and-omega-three-pion.md),
and contracts it with the source-resolved ALP--photon--vector vertex. The
operator normalization is $\mathcal L_{VPPP}=\sum_{r=1}^5g_r\mathcal
O^r_{VPPP}/M_V$, with independent real $g_r$. Use the repository conventions
in [Alda:2026nrd](../../extracts/Alda:2026nrd/INDEX.md), metric $(+---)$,
$\epsilon^{0123}=+1$, and the antisymmetric-tensor propagator derived in the
[pole-factorization record](alp-three-pion-vector-pole-factorization.md).
The derivation scripts use a fixed all-incoming routing and calibrate the
production-source orientation by reproducing the established wedge numerator
$N_i(s,y)$; the physical crossing/LSZ phase is not inferred anew here.

Inputs are the complete five-operator source basis in
[@GomezDumm:2010ra], the ordered source-resolved $N_i(s,y)$ and its
VJP/VVP coupling conventions in the
[pole-factorization record](alp-three-pion-vector-pole-factorization.md),
and the full tensor response in that record. The retained sequential
$\rho$-exchange term is a separate contribution, given in the
[off-shell exchange record](alp-isoscalar-vector-three-pion-offshell.md).

# Local $VPPP$ endpoint

Let $s=(p_++p_-+p_0)^2$, $s_{ab}=(p_a+p_b)^2$, and
$s_{+-}+s_{+0}+s_{-0}=s+3m_\pi^2$. For an isoscalar generator
$T_i=\operatorname{diag}(t_n,t_n,t_s)$, the five coefficients multiplying
the geometric charged-pion epsilon tensor in the Lagrangian, before the
common Feynman-vertex factor $i$, are
$$
\left(\frac{16it_n(s-3m_\pi^2)}{M_i f_\pi^3},
-\frac{16it_n(s-3m_\pi^2)}{M_i f_\pi^3},
-\frac{16it_n(s-3m_\pi^2)}{M_i f_\pi^3},
\frac{96it_nm_\pi^2}{M_i f_\pi^3},
\frac{48it_nm_\pi^2}{M_i f_\pi^3}\right).
$$
Thus the polarized local scalar is
$$
\mathcal F^{\rm local}_{i,\rm pol}(s)=
\frac{16it_n}{M_iM_Vf_\pi^3}
\left[(g_1-g_2-g_3)(s-3m_\pi^2)
+3m_\pi^2(2g_4+g_5)\right].
$$
It has no Dalitz-variable dependence and is therefore entirely in total
spin one. At the ideal-$\omega$ pole, $t_n=1/\sqrt2$ and $s=M_\omega^2$;
the result agrees with the local term in Appendix B.1 of
[@GomezDumm:2010ra], after the declared tensor and coupling convention map.
The channel probes $g_1-g_2-g_3$ and $2g_4+g_5$. The paper's quoted
$2g_4+g_5=-0.60\pm0.02$ follows only after its resonance-model constraints
and branching-fraction input; it is not imposed here as a universal ALP
input.

The unpolarized endpoint is not determined by this scalar alone. Define
$$
D_{ab}^{\mu\nu}=\epsilon^{\mu\nu\alpha\beta}p_{a,\alpha}p_{b,\beta},
\qquad (ab)=(+-),(+0),(-0),
$$
and let $B_{\rm contact}^{\mu\nu}$ be the sum of these three two-forms
weighted by the five rows of the executable coefficient matrix. In the
three-pion rest frame, for rows $r=1,\ldots,5$, those weights are
$$
\frac{4it_n}{f_\pi^3}\begin{pmatrix}
-(4m_\pi^2-s_{+0}-s_{-0}) & 4m_\pi^2-s_{+-}-s_{-0} & -(4m_\pi^2-s_{+-}-s_{+0})\\
4m_\pi^2-s_{+0}-s_{-0} & -(4m_\pi^2-s_{+-}-s_{-0}) & 4m_\pi^2-s_{+-}-s_{+0}\\
4m_\pi^2-s_{+0}-s_{-0} & -(4m_\pi^2-s_{+-}-s_{-0}) & 4m_\pi^2-s_{+-}-s_{+0}\\
4m_\pi^2 & -4m_\pi^2 & 4m_\pi^2\\
2m_\pi^2 & -2m_\pi^2 & 2m_\pi^2
\end{pmatrix}.
$$
Let $K_{r,ab}$ denote the displayed matrix, and define the lower-index
two-form $D_{ab,\mu\nu}=\epsilon_{\mu\nu\alpha\beta}
p_a^\alpha p_b^\beta$.
The endpoint components used in the propagator contraction are
$$
B_{\rm contact}^{01}=-\sum_{r=1}^5\frac{g_r}{M_V}
\sum_{ab}K_{r,ab}D_{ab,01},\qquad
B_{\rm contact}^{13}=+\sum_{r=1}^5\frac{g_r}{M_V}
\sum_{ab}K_{r,ab}D_{ab,13}.
$$
The relative signs follow from raising the $01$ versus $13$ components of
the operator source. The common Fourier factors are fixed by the script.
The matrix preserves all three independent Dalitz structures; the polarized
scalar is its spin-one projection.

# Full tensor propagation into the ALP production source

The ALP production source is derived before restricting the internal vector
to a spin-one polarization. In the rest frame of $k=p+q$, take
$q=(q^0,0,0,Q)$, photon polarization $\varepsilon_\gamma=(0,0,1,0)$,
$q^0=(m_a^2-s-y)/(2\sqrt{s})$, and $y=q^2$. With
$\mathcal L=V_{\mu\nu}A_i^{\mu\nu}$, the only nonzero upper-source
components are
$$
A_i^{01}=-\frac{Q N_i(s,y)}{2},\qquad
A_i^{13}=\frac{\mathcal C_{i,13}(s,y)}{2}.
$$
The first equality recovers the strict spin-one source and the previously
derived numerator. The second is the complementary tensor component, not
fixed by $N_i$ alone.

The contact-relevant production coefficient is
$$
\begin{aligned}
\mathcal C_{i,13}(s,y)={}&-\frac{2e}{M_V\sqrt{s}}\Big\{
\operatorname{tr}(T_iQ U_a)
\big[c_1(m_a^4-(s-y)^2)
+c_2((m_a^2-s-y)^2-4sy)
+c_5(s^2-(m_a^2-y)^2)
-2c_7s(m_a^2-s+y)\big]
 +4c_3\operatorname{tr}(T_iQY_a)(-m_a^2+s+y)
\Big\}\\
&+\frac{\sqrt2e}{\sqrt{s}}\sum_jW_j(y)\Big\{
\operatorname{tr}(T_iT_jU_a)
\big[2d_1m_a^2(m_a^2-s-y)
-2d_3(m_a^4-m_a^2s-2m_a^2y-sy+y^2)
-2d_4s(m_a^2-s+y)\big]
 +8d_2\operatorname{tr}(T_iT_jY_a)(-m_a^2+s+y)
\Big\}.
\end{aligned}
$$
Here $U_a=2c_q^A/f_a$, $Y_a=-2\{c_q^P,\chi\}/f_a$, and $W_j(y)$ is
the photon response defined in the linked source-resolved TFF record. The
matrix traces retain arbitrary diagonal ALP sources; the isoscalar choice is
made by $T_i=\omega,\phi$ on the three-pion leg. The complementary source
contains $c_7$ and $d_4$, which do not appear in the wedge numerator $N_i$.

For any decay endpoint $B_i^{\mu\nu}$, the full free-tensor contraction in
this polarization channel is
$$
\mathcal A_{i,\rm tensor}=\frac{2QN_i(s,y)B_{i,01}}{s-M_i^2+i0}
-\frac{2\mathcal C_{i,13}(s,y)B_{i,13}}{M_i^2}.
$$
For the direct local endpoint, use $B_i=B_{\rm contact}$ from the previous
section. The first term is the spin-one pole; the second is the local tensor
projector contribution. It generally survives because the full local
$VPPP$ endpoint has spatial tensor components in the $k$ rest frame. Thus
the strict pole residue is not the complete contribution of the local
contact operators. This formula is a fixed free-propagator tree
continuation; widths, rescattering, field redefinitions, and matching to a
dispersive amplitude are not included.

# Inputs not fixed by the derivation

No values for the independent $g_i$ are selected. The on-shell ideal-$\omega$
constraint reported by [@GomezDumm:2010ra] is conditional on that paper's
short-distance model, coupling assumptions, and branching-ratio input. An
ALP phenomenology fit needs a compatible matching prescription or independent
data for these couplings, plus complex rescattered decay functions and a
common phase convention. The strict pole and regular tensor term must be
matched together before they are combined with a dispersive $V_i\to3\pi$
amplitude.

# Reproducibility

Run the charged endpoint calculation with
`.venv/bin/python scripts/derive_vector_three_pion_contacts.py` and the
unprojected ALP-source calculation with
`.venv/bin/python scripts/derive_alp_contact_production_source.py`. The first
derives the five polarized coefficients, the complete $5\times3$ two-form
coefficient matrix, and the free-tensor response. The second differentiates
the retained VJP/VVP operators with respect to all six antisymmetric vector
components, reproduces $N_i$ under the wedge restriction, and prints the
complementary $c_i,d_i$ source, the five-coupling contact endpoint, and both
propagator contractions. Both use natural units and deterministic SymPy
expressions.
