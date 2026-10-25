---
title: Retained odd-resonance virtual-photon kaon production amplitude
status: derived
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx, Messerli:2025rnv]
provenance: analytical-derivation
---

# Scope and inputs

This extends the retained odd-resonance
[$a\to\gamma PP$ amplitude](../photon-two-pseudoscalar/alp-photon-pp-amplitude.md)
to $a\gamma^*\to K^+K^-$ and $a\gamma^*\to K^0\bar K^0$.
The original seven VJP and four VVP operators and the leading $G_V$ vertex
are retained, at leading order in $e$ [@Ruiz-Femenia:2003jdx;
@Ecker:1988te; @Ecker:1989yg]. The independent diagonal ALP axial and
pseudoscalar sources already use the no-mixing convention of
Alda:2026nrd [@Alda:2026nrd]. No ALP source is reconstructed from a
physical-meson amplitude.

The input Lorentz contractions are those of the
[virtual-pion derivation](alp-pion-vector-and-production-form-factors.md).
The [current decomposition](alp-isoscalar-tff-source-decomposition.md)
distinguishes projection of the electromagnetic current from selection of
an omega/phi pole. The electromagnetic charge is removed from current
matrix elements, as in the pion record. Natural units, masses in GeV,
metric $(+---)$, and $\epsilon^{0123}=+1$ are used.

This is an explicitly restricted contribution. Independent local anomalous
terms, additional operators, even scalar/vector-source radiative terms,
loops and rescattering are absent. It is not a complex physical kaon
production amplitude or a complete kaon spectral contribution.

# Kinematics and executable derivation

Cross the decay photon with $k=-q$, $k^2=q^2=y$, and retain
$$
s=(p_1+p_2)^2,\qquad t=(k+p_2)^2,\qquad
u=(k+p_1)^2,\qquad p_a=k+p_1+p_2.
$$
For either equal-mass kaon channel,
$s+t+u=m_a^2+2M_K^2+y$. Define the scalar coefficient by the existing
decay tensor,
$$
\mathcal M^{\rm odd}=\epsilon_{\mu\nu\rho\sigma}
\varepsilon^{*\mu}p_1^\nu p_2^\rho k^\sigma
\mathcal A^{\rm odd}_{aK\bar K\gamma^*}(s,t,u;y).
$$
The [script](../../../scripts/derive_alp_virtual_photon_kaon_amplitude.py)
derives the following result by re-evaluating the original contractions at
nonzero $y$, rather than attaching a form factor to the real-photon result.

1. Keep the source matrices $\mathcal U_X,\mathcal Y_X$ and all four
   trace contractions $D_R(X),P_R(X),D_{RQ}(X),P_{RQ}(X)$ from the
   real-photon record, replacing the response there by
   $$
   \mathcal W_\gamma(y)=\sum_{i=\rho,\omega,\phi}
   \frac{R_iT_i}{M_i^2-y-i0}.
   $$
   The full neutral response remains on the photon leg. Both isovector
   and isoscalar terms are required for arbitrary diagonal ALP sources.

2. Substitute $s_1=y,s_2=r^2$ in the VVP contractions, and
   $s_1=r^2,s_2=y$ in the ordered VJP contractions. SymPy obtains
   $$
   \begin{aligned}
   \mathcal N_R(X;y,r^2)&=4\big(
   [(d_1-d_3)m_X^2+d_3(r^2+y)]D_R(X)-4d_2P_R(X)\big),\\
   \mathcal V_{RQ}(X;r^2,y)&=\frac{2}{M_V}\big(
   [(c_1+c_2-c_5)m_X^2+(-c_1+c_2+c_5-2c_6)r^2
   +(c_1-c_2+c_5)y]D_{RQ}(X)-4c_3P_{RQ}(X)\big).
   \end{aligned}
   $$
   The term proportional to $y$ in the direct vertex and the $d_3y$
   conversion term are therefore retained independently of
   $\mathcal W_\gamma(y)$.

3. Retain the original $c_4$ contraction at charged pseudoscalar
   insertions. The complete retained radiative numerator is
   $$
   \mathcal H_R(X;r^2,y)=\sqrt2e\,\mathcal N_R(X;y,r^2)
   -2e\,\mathcal V_{RQ}(X;r^2,y)
   -\delta_{X=P_i}\frac{8ec_4}{f_\pi M_V}
   \operatorname{tr}\big(T_R[[\lambda_i^\dagger,Q],\chi]\big).
   $$
   This term does not acquire an additional $y$ factor. The charge
   commutator vanishes for $K^0,\bar K^0$ because $Q_d=Q_s$.
   For $K^\pm$ the subsequent mass commutator is proportional to
   $\chi_s-\chi_u$, so setting it to zero would omit a retained
   charged-kaon contribution. The $c_7,d_4$ contractions still vanish
   for these routed transverse tensor currents.

4. Insert these numerators into the three original exchange graphs:
   $$
   \begin{aligned}
   \mathcal A^{\rm odd}_{aK\bar K\gamma^*}
   =\sum_R\Bigg[&
   -\frac{G_V\mathcal H_R(a;s,y)
   \operatorname{tr}(T_R[\lambda_1^\dagger,\lambda_2^\dagger])}
   {\sqrt2f_\pi^2(s-M_R^2+i0)}\\
   &+\frac{\sqrt2G_V\mathcal H_R(P_2;t,y)
   \operatorname{tr}(T_R[\lambda_1^\dagger,c_q^A])}
   {f_\pi f_a(t-M_R^2+i0)}\\
   &-\frac{\sqrt2G_V\mathcal H_R(P_1;u,y)
   \operatorname{tr}(T_R[\lambda_2^\dagger,c_q^A])}
   {f_\pi f_a(u-M_R^2+i0)}\Bigg].
   \end{aligned}
   $$
   The sum runs over an orthonormal Hermitian vector nonet, with the same
   mass for the two Hermitian components of a given charged vector.
   Evaluating the supplied kaon generators gives rho, omega and phi
   exchange in $s$ for both channels; $K^{*\pm}$ exchange in $t,u$ for
   $K^+K^-$; and $K^{*0},\bar K^{*0}$ exchange for $K^0\bar K^0$.
   No charged-rho exchange survives these flavor traces. Symbolic
   charged/neutral kaon and vector masses are kept distinct; this does
   not supply a full calculation of isospin-breaking interactions.

# Current form factors and kaon cut

Define the two external complex form factors by
$$
\begin{aligned}
\langle K^+(p_1)K^-(p_2)|J^\mu_{\rm em}|0\rangle
&=(p_1-p_2)^\mu F_{K^+}^V(s),\\
\langle K^0(p_1)\bar K^0(p_2)|J^\mu_{\rm em}|0\rangle
&=(p_1-p_2)^\mu F_{K^0}^V(s).
\end{aligned}
$$
Charge normalization requires $F_{K^+}^V(0)=1$ and $F_{K^0}^V(0)=0$;
it does not determine their complex timelike values or correlations.
In exact strong isospin, the isoscalar part of either current matrix
element is their half-sum; their half-difference gives the isovector
part, with opposite signs for the charged and neutral channel.
The arbitrary ALP production source should be projected on the hadronic
state separately. Using physical mass splittings while claiming exact
isospin in all other pieces is an approximation that must be stated.

The existing distinct-particle pion-cut derivation applies to each
equal-mass kaon pair with its own threshold and current production wave.
For the charge-stripped coefficient
$f_{1,\mathrm{current}}=(3/4)\int_{-1}^{1}dz(1-z^2)
\mathcal A^{\rm odd}/e$ in the same real-epsilon crossing convention,
the contribution to the single discontinuity is
$$
\frac{\operatorname{Disc}_s F_{a\gamma^*\gamma^*}(s,y)}{2i}
\bigg|_{K\bar K}
=\sum_{K=K^+,K^0}
\theta(s-4M_K^2)\frac{s(1-4M_K^2/s)^{3/2}}{96\pi}
[F_K^V(s)]^*f_{1,\mathrm{current}}(s,y).
$$
Here the cut product is conditional on the complete compatible current
amplitudes, rather than a prediction obtained by inserting the retained
tree term alone. The script reproduces the $1/(96\pi)$ coefficient
from the angular measure and the $2\operatorname{Im}T$ normalization in
the [pion record](alp-pion-vector-and-production-form-factors.md).
Charged and neutral kaons are distinct intermediate particles, so no
$1/2!$ is inserted. The two channels are both summed once. For arbitrary
complex $y$, this denotes the specified cut discontinuity, not an
unqualified imaginary part. The three-flavor isoscalar-current projection
requires the corresponding half-sum on the current side, with consistent
isospin and threshold conventions.

# Phi region, validity, and remaining inputs

The user specifies $2m_e<m_a<2m_\tau$ and the unbounded virtualities of
the complete two-photon lepton loop. The distinct thresholds are
$4M_{K^+}^2$ and $4M_{K^0}^2$. A narrow stable phi pole is not a uniform
description of that scope. Once physical complex kaon form factors and
rescattered production amplitudes are used, their phi strength must replace
the same explicit phi contribution in a common coupled-channel or matched
description. Independently adding a phi pole, a phi-to-three-pion term,
and a phi-to-kaon term without shared propagator/self-energy and subtraction
bookkeeping counts overlapping resonance strength. This derivation selects
no finite-width or continuum prescription.

To evaluate the kaon contribution one still needs compatible complex
charged/neutral kaon current form factors, the completed source-dependent
production amplitudes, their phases and correlated uncertainties, and a
three-pion/kaon treatment of the phi region. Higher channels and a controlled
high-energy continuation remain necessary for the full loop
[@Messerli:2025rnv].

# Reproducibility and checks

Run:

```text
.venv/bin/python scripts/derive_alp_virtual_photon_kaon_amplitude.py
.venv/bin/python scripts/derive_alp_virtual_photon_kaon_amplitude.py --test
```

The script prints the original contractions, both virtual radiative
vertices, the $c_4$ term, the three exchange contributions, the selected
kaon flavor sums, and the distinct-particle cut coefficient. Four tests
compare each $y=0$ vertex piece to the existing real-photon derivation,
verify the photon Ward identity at the original-operator and assembled
tensor levels, check allowed vector exchanges and charged/neutral $c_4$,
and verify source linearity and the angular normalization. These checks
validate the stated retained tree contribution, not the missing physical
spectral inputs. The derivation and all four tests were executed successfully.
