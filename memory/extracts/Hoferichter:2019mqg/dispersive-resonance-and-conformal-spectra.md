---
title: Resonance and conformal spectral functions for the three-pion normalization
status: literature
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Hoferichter:2019mqg, Hoferichter:2025lcz, Niecknig:2012sj]
provenance: literature
---

# Dispersive normalization

[@Hoferichter:2019mqg, Eqs. (13)--(17)] parameterizes the subtraction
normalization by
$$
a(q^2)=\alpha_A+
\frac{q^2}{\pi}\int_{s_{\mathrm{thr}}}^{\infty}
\frac{\operatorname{Im}\mathcal A(s')}{s'(s'-q^2)}\,ds'
+C_p(q^2),
$$
with $s_{\mathrm{thr}}=M_{\pi^0}^2$ because of $\omega\to\pi^0\gamma$.
The resonance function is
$$
\mathcal A(q^2)=\sum_V
\frac{c_V}{M_V^2-q^2-i\sqrt{q^2}\,\Gamma_V(q^2)},
\qquad
V=\omega,\phi,\omega',\omega''.
$$
For real $s>0$, real $c_V$, and the stated width convention, a per-unit-residue
spectral template is proportional to
$$
\frac{\sqrt{s}\,\Gamma_V(s)}
{(M_V^2-s)^2+s\Gamma_V(s)^2}.
$$
The paper specifies that the $\omega$ and $\phi$ widths include the principal
decay modes; the $3\pi$ partial widths include final-state rescattering as in
[@Niecknig:2012sj]. The remaining roughly 2% missing channels are handled by
rescaling partial widths. The excited $\omega'$ and $\omega''$ are treated as
100% $3\pi$ decays in the fit model [@Hoferichter:2019mqg].

# Conformal spectral basis and sum rule

The inelastic contribution is
$$
C_p(q^2)=\sum_{n=1}^{p}c_n\big[z(q^2)^n-z(0)^n\big],
\qquad
z(q^2)=\frac{\sqrt{s_{\mathrm{inel}}-s_1}-\sqrt{s_{\mathrm{inel}}-q^2}}
{\sqrt{s_{\mathrm{inel}}-s_1}+\sqrt{s_{\mathrm{inel}}-q^2}},
$$
with $s_{\mathrm{inel}}=1\,\mathrm{GeV}^2$ and
$s_1=-1\,\mathrm{GeV}^2$. The spectral template for $c_n$ is the upper-lip
imaginary part of $z(s+i0)^n-z(0)^n$ for $s\ge s_{\mathrm{inel}}$.
Additional constraints remove the induced $S$-wave cusp and impose
$$
\alpha_A=\frac{1}{\pi}\int_{s_{\mathrm{thr}}}^{\infty}
\frac{\operatorname{Im}\mathcal A(s')}{s'}\,ds'
+\frac{1}{\pi}\int_{s_{\mathrm{inel}}}^{\infty}
\frac{\operatorname{Im}C_p(s')}{s'}\,ds'
$$
exactly [@Hoferichter:2019mqg].

# Numerical input limitation

The source gives the spectral parameterization but not reusable numerical
grids for $\Gamma_\omega(s)$ or $\Gamma_\phi(s)$. In particular, obtaining
their $3\pi$ partial widths requires the normalized, rescattered
Khuri--Treiman amplitude; other decay-channel widths and the chosen high-
energy continuation are also needed. Therefore the resonance spectral
templates must take explicit width functions or arrays as input. Substituting
constant widths or the zero-width limit without labeling that approximation
would not reproduce the published fit basis. [@Hoferichter:2025lcz] updates
the global-fit data and line-shape treatment but likewise does not publish
these width arrays.
