---
title: Conditional finite-width isoscalar vector response and pole replacement
status: derived
created: 2026-10-26
updated: 2026-10-26
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx, Niecknig:2012sj]
provenance: analytical-derivation
---

# Inputs, conventions, and scope

This derives the general absorptive and dispersive bookkeeping needed to
replace a stable isoscalar vector response by a finite-width response. It
uses the same transverse vector exchange that appears in the
[source-resolved tree TFF](alp-doubly-virtual-tff.md) and keeps the
three-pion and kaon channels discussed in the [isoscalar
plan](alp-isoscalar-tff-plan.md) and [virtual-kaon
record](alp-virtual-photon-kaon-amplitude.md). The decay amplitudes and
unitarity conventions are those of the [three-pion
records](alp-three-pion-vector-partial-wave.md) and
[discontinuity derivation](alp-three-pion-discontinuity.md)
[@Alda:2026nrd; @Messerli:2025rnv; @Ecker:1988te; @Ecker:1989yg;
@Ruiz-Femenia:2003jdx; @Niecknig:2012sj].

Use natural units, invariant masses in GeV, relativistic state normalization,
metric $(+---)$, $S=1+iT$, and the upper boundary $s+i0$. For this formal
derivation take a single diagonal transverse response, with inverse
propagator
$$
D_i(s)=M_{0i}^2-s-\Delta_i^{\mathrm{self}}(s).
$$
The notation $\Delta_i^{\mathrm{self}}$ denotes the transverse self-energy
correction of vector $i$. The width is defined conditionally by the chosen
off-shell transverse self-energy and the complete set of open states:
$$
\Gamma_i(s)=\frac{1}{2\sqrt{s}}\sum_n\int d\Phi_n\,
\overline{|\mathcal M_{i\to n}|^2},
$$
where the bar includes the average over the three transverse vector
polarizations. This is not a selected channel model: no numerical decay
amplitudes, thresholds beyond those in their channel records, subtraction
constants, or high-energy continuation are supplied here.

# Optical theorem and spectral response

With the stated $S$-matrix convention, the cut of the transverse
self-energy gives
$$
2\operatorname{Im}\Delta_i^{\mathrm{self}}(s+i0)
=\sum_n\int d\Phi_n\,\overline{|\mathcal M_{i\to n}|^2}.
$$
Combining this with the width definition fixes the sign for the chosen
inverse propagator:
$$
\operatorname{Im}\Delta_i^{\mathrm{self}}(s+i0)=+\sqrt{s}\,\Gamma_i(s).
$$
Thus the denominator on the upper edge is
$$
M_{0i}^2-s-\operatorname{Re}\Delta_i^{\mathrm{self}}(s)
-i\sqrt{s}\Gamma_i(s),
$$
which yields a positive spectral density for the dressed response
$G_i(s)=1/D_i(s)$:
$$
\rho_i^{\mathrm{vec}}(s)
=\frac{1}{\pi}\operatorname{Im}G_i(s+i0)
=\frac{\sqrt{s}\Gamma_i(s)}
{\pi\left[
\left(M_{0i}^2-s-\operatorname{Re}\Delta_i^{\mathrm{self}}(s)\right)^2
+s\Gamma_i^2(s)\right]}.
$$
This response-level density does not supply the ALP production vertices,
regular terms, or direct multiparticle contacts. In particular, the
strict-pole factorization in the three-pion record fixes a residue, not a
complete off-shell production amplitude.

# Subtracted self-energy and stable-pole limit

If $N$ subtractions at $s=0$ are sufficient, and the subtraction point lies
below all included thresholds, Cauchy's formula gives
$$
\Delta_i^{\mathrm{self}}(s)
=\sum_{k=0}^{N-1}c_{ik}s^k
+\frac{s^N}{\pi}\int_{s_{\mathrm{th}}}^{\infty}
\frac{\operatorname{Im}\Delta_i^{\mathrm{self}}(x+i0)}
{x^N(x-s-i0)}\,dx.
$$
The real coefficients $c_{ik}$ are independent subtraction inputs. The
number of subtractions, thresholds included, and convergence of the integral
must be established for a chosen channel and high-energy description; they
are not fixed by this formal expression.

# On-shell three-pion partial width

The existing on-shell decay amplitude is
$$
\mathcal M(V_i\to\pi^+\pi^-\pi^0)
=iF_i^{\mathrm{on}}(t,u)\,n_\mu J_{3\pi}^{\mu},
\qquad
J_{3\pi}^{\mu}=\epsilon^{\mu\nu\alpha\beta}
p_{+\nu}p_{-\alpha}p_{0\beta}.
$$
For a vector of invariant mass $\sqrt{s}$, the polarization-averaged
contraction is $-J_{3\pi}^2/3$. In the charged-pair variables
$u=(p_++p_-)^2$ and angle $z$, the established phase-space and Gram
identities give
$$
\frac{d\Phi_3}{du\,dz}
=\frac{\sqrt{\lambda(s,u,m_0^2)\lambda(u,m_+^2,m_-^2)}}
{256\pi^3 s u},
\qquad
-J_{3\pi}^2
=\frac{(1-z^2)\lambda(s,u,m_0^2)\lambda(u,m_+^2,m_-^2)}
{16u}.
$$
Substituting these into the width integral derives the on-shell charged
channel width
$$
\Gamma_{i\to\pi^+\pi^-\pi^0}(M_i^2)
=\frac{1}{24576\pi^3M_i^3}
\int_{(m_++m_-)^2}^{(M_i-m_0)^2}\!du
\int_{-1}^{1}\!dz\,
\frac{(1-z^2)
[\lambda(M_i^2,u,m_0^2)\lambda(u,m_+^2,m_-^2)]^{3/2}}{u^2}
\left|F_i^{\mathrm{on}}(u,t(u,z;M_i^2))\right|^2.
$$
The charged particles are distinct, so no $1/3!$ appears. This formula
expresses the on-shell width in terms of the complex decay solution; the
overall phase drops out of its modulus. The available
$F_i^{\mathrm{on}}$ is defined at $s=M_i^2$. Promoting this expression to an
energy-dependent $\Gamma_i(s)$ in a self-energy requires an off-shell
$V_i\to3\pi$ amplitude, which is not fixed by the on-shell partial-wave
solution. The physical phase relative to the signed ALP production residue
also remains independent of the width.

To identify the overlap with a stable pole, regulate its absorptive part by
a positive infinitesimal $\varepsilon$. SymPy derives the normalized
distribution
$$
\rho_i^{\mathrm{stable}}(s)
=\frac{\varepsilon}{\pi[(M_i^2-s)^2+\varepsilon^2]}
\xrightarrow{\varepsilon\to0^+}\delta(s-M_i^2).
$$
The once-subtracted pole reconstruction is
$$
\frac{R_i}{M_i^2-s}
=\frac{R_i}{M_i^2}
+\frac{R_i s}{M_i^2(M_i^2-s)}.
$$
Consequently the stable pole and its spectral delta are two representations
of the same strength. Once a finite-width response is represented through
its spectral density, the same explicit pole contribution must be removed
or matched; retaining both unchanged counts it twice. In a coupled
three-pion/kaon treatment of the $\phi$ region, the shared resonance strength
must likewise be represented once.

# Limits and inputs still required

- The derivation is diagonal in vector space. Common decay channels can
  induce off-diagonal $\omega$--$\phi$ self-energies; those require a
  coupled propagator matrix and are not calculated here.
- A numerical $\Delta_i^{\mathrm{self}}$ needs compatible complex decay
  amplitudes and their phase-space integrals for every retained open
  channel. The current records do not yet provide all such amplitudes and
  correlated uncertainties.
- Subtraction conditions must set the physical mass and residue convention.
  No conditions or channel truncation are selected.
- The finite-width response must be used consistently in current and ALP
  production amplitudes. It does not by itself complete the source-dependent
  $a\gamma^*\to3\pi$ or $a\gamma^*\to K\bar K$ amplitudes.

# Reproducibility

The executable
[derive_alp_isoscalar_vector_self_energy.py](../../../scripts/derive_alp_isoscalar_vector_self_energy.py)
uses SymPy to derive the optical-theorem sign, spectral density, general
$N$-subtracted dispersion relation, on-shell three-pion width, and
stable-pole distribution and replacement identities. Run:

    .venv/bin/python scripts/derive_alp_isoscalar_vector_self_energy.py

The script was executed successfully. Its four focused tests also pass;
they check the optical-theorem sign and positivity, the subtracted Cauchy
decomposition, the on-shell three-pion width prefactor, and the unit-weight
stable-pole limit and reconstruction. These checks validate the formal
identities, not a numerical vector self-energy or a physical continuum/pole
fit.
