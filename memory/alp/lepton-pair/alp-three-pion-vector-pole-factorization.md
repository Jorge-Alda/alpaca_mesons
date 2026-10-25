---
title: Spin-one pole factorization for the ALP three-pion current
status: derived
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Alda:2026nrd, Ecker:1989yg, Ruiz-Femenia:2003jdx, Niecknig:2012sj, Hoferichter:2025lcz]
provenance: analytical-derivation
---

# Scope and conventions

This calculation derives the strict spin-one pole contribution to
$a(p)+\gamma^*(q)\to\pi^+(p_+)\pi^-(p_-)\pi^0(p_0)$ for a neutral vector
$V_i=\omega,\phi$ on the three-pion leg. Let $k=p+q=p_++p_-+p_0$,
$s=k^2$, and $y=q^2$. Use natural units, metric $(+---)$, and
$\epsilon^{0123}=+1$. The ALP source matrices $c_q^A,c_q^P$ are independent
inputs in the no-mixing convention of Alda:2026nrd [@Alda:2026nrd]. The
ordered VJP/VVP amplitudes use the canonical couplings and tensor convention
of the [source-resolved doubly virtual TFF](alp-doubly-virtual-tff.md).

The three-pion decay endpoint is restricted to the transverse spin-one lift
of an on-shell vector decay amplitude. This defines a pole residue; it does
not determine regular off-shell terms or direct $V3\pi$ contact operators.
The available $V\to3\pi$ RChT extract supplies only a sequential on-shell
contribution and notes that direct terms are missing
[@Ruiz-Femenia:2003jdx]. The dispersive Khuri--Treiman amplitude supplies a
spin-one decay subamplitude, not the full ALP production amplitude
[@Niecknig:2012sj; @Hoferichter:2025lcz].

# Inputs and ordered production numerator

Use
$$
U_a=\frac{2c_q^A}{f_a},\qquad
Y_a=-\frac{2\{c_q^P,\chi\}}{f_a},\qquad
W_j(y)=\frac{R_j}{M_j^2-y-i0},
$$
where $R_j$ and $M_j$ are the signed photon residue and pole mass of each
neutral intermediate vector. For $j=\rho,\omega,\phi$, the ordered numerators
are
$$
\begin{aligned}
V_{i\gamma}(s,y)={}&\frac{2}{M_V}\Big[
\operatorname{tr}(T_i Q U_a)\big((c_1+c_2-c_5)m_a^2
+(-c_1+c_2+c_5-2c_6)s+(c_1-c_2+c_5)y\big)
-4c_3\operatorname{tr}(T_i QY_a)\Big],\\
P_{ij}(s,y)={}&4\Big[
\operatorname{tr}(T_iT_jU_a)\big((d_1-d_3)m_a^2+d_3(s+y)\big)
-4d_2\operatorname{tr}(T_iT_jY_a)\Big].
\end{aligned}
$$
The direct VJP and photon-converted VVP terms combine into the wedge-projected
production numerator
$$
N_i(s,y)=-2eV_{i\gamma}(s,y)
+\sqrt2e\sum_{j=\rho,\omega,\phi}W_j(y)P_{ij}(s,y).
$$
The executable obtains $V_{i\gamma}$ and $P_{ij}$ from the ordered Lorentz
contractions and substitutes the ALP mass $m_a$ in both kernels. The full
charge matrix remains on the photon leg; selecting $i=\omega,\phi$ does not
project that leg onto the light-quark isoscalar current.

# Antisymmetric-tensor propagator

The six independent components of the canonical free tensor field are
formed from
$$
\mathcal L_0=-\frac12(\partial^\mu V_{\mu\nu})
(\partial_\rho V^{\rho\nu})
+\frac{M_i^2}{4}V_{\mu\nu}V^{\mu\nu}
$$
with the source term normalized as $V_{\mu\nu}B^{\mu\nu}$. SymPy forms the
six-by-six Hessian and inverts it before imposing any source restriction. If
$B^{\mu\nu}$ is a general antisymmetric source, define its spin-one
projection by
$$
B_{(1)}^{\mu\nu}=\frac{1}{s}\left[
k^\mu B^{\rho\nu}k_\rho-k^\nu B^{\rho\mu}k_\rho
\right].
$$
The executed inverse-Hessian decomposition is
$$
V^{\mu\nu}=\frac{2B_{(1)}^{\mu\nu}}{s-M_i^2}
-\frac{2(B^{\mu\nu}-B_{(1)}^{\mu\nu})}{M_i^2}.
$$
The second term is the local tensor projector contribution. It is retained
until the decay source is specified.

Write the geometric three-pion current as
$$
J^\beta_{3\pi}=\epsilon^{\beta\lambda\kappa\delta}
p_{+\lambda}p_{-\kappa}p_{0\delta},\qquad
\mathcal J^\beta_{V_i\to3\pi}=F_i^{\rm on}(t,u)J^\beta_{3\pi},
$$
where $F_i^{\rm on}$ denotes the on-shell decay function in the chosen
Khuri--Treiman convention. Momentum conservation gives
$k_\beta J^\beta_{3\pi}=0$. For the declared pole-residue lift
$$
B_{3\pi}^{\mu\nu}=k^\mu J^\nu_{3\pi}-k^\nu J^\mu_{3\pi},
$$
SymPy finds $B_{(1)}=B_{3\pi}$ and the local projector term vanishes. The
tensor response on this source is therefore purely the spin-one pole. This
does not determine other allowed off-shell decay sources.

# Pole-factorized ALP production amplitude

Let
$$
\mathcal E_i^{\alpha\mu}=\epsilon^{\alpha\mu\rho\sigma}
k_\rho q_\sigma,\qquad
\Pi^{(i)}_{\alpha\beta}(k)=-g_{\alpha\beta}
+\frac{k_\alpha k_\beta}{M_i^2}.
$$
The production coupling at the vector pole is the on-shell residue of the
ordered numerator divided by the external-vector tensor normalization $M_i$.
Keep $\zeta_i$ explicit for the common crossing/LSZ phase between the
production and decay amplitudes:
$$
\mathcal M^\mu_{i,\rm pole}=
\frac{i\,\zeta_i F_i^{\rm on}(t,u)}{M_i(M_i^2-s-i0)}
N_i(M_i^2,y)
\left[\mathcal E_i^{T}g\Pi^{(i)}gJ_{3\pi}\right]^\mu.
$$
The phase convention must be fixed by matching the real-photon on-shell
production amplitude to the existing $a\to V_i\gamma$ record. The on-shell
decay amplitude alone does not set that relative phase. The script outputs
the contracted four-vector and verifies both Ward contractions:
$$
q_\mu\mathcal E_i^{\alpha\mu}=0,\qquad
k_\alpha\mathcal E_i^{\alpha\mu}=0,\qquad
k_\beta J_{3\pi}^\beta=0.
$$
Consequently the $k_\alpha k_\beta/M_i^2$ term in $\Pi^{(i)}$ drops. The
projector contraction reduces to the transverse contraction of the two
epsilon tensors, producing the ordinary-parity vector structure expected on
the ALP side of the three-pion cut.

The scalar numerator can be separated algebraically as
$$
\frac{N_i(s,y)}{M_i^2-s}
=\frac{N_i(M_i^2,y)}{M_i^2-s}
+\frac{N_i(s,y)-N_i(M_i^2,y)}{M_i^2-s}.
$$
Only the first term is fixed by the on-shell spin-one residue. This record
does not promote the regular continuation to a physical prediction. It also
does not supply the complex $F_i^{\rm on}$ grid, its covariance, additional
three-pion partial waves, local ALP-current terms, or a double-spectral TFF.

# Reproducibility and checks

Run:
```text
.venv/bin/python scripts/derive_alp_three_pion_vector_pole.py
```
The executed SymPy calculation inverts the full six-component Hessian and
prints its exact agreement with the spin-one/local-projector decomposition.
For the pure transverse decay lift it returns a zero local term and the
transversality identity $k\cdot J_{3\pi}=0$. It also derives the ordered
VJP/VVP numerators with $m_P$ replaced by the ALP mass, contracts the
production epsilon tensor with the spin-one projector and decay current, and
returns zero for the projector-reduction difference and both production Ward
contractions. No tests, fit, finite-width model, or numerical decay function
were run or selected.
