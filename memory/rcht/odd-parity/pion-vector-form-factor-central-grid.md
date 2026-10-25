---
title: Central dispersive pion vector form-factor grid
status: numerical
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [GomezDumm:2013sib, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# Scope and inputs

This calculates the complex isospin-symmetric pion vector form factor needed
as one input to the low-energy pion TFF spectral density. The dispersive
definition and phase model follow GomezDumm:2013sib, and the selected fit is
its exact-isospin fit I [@GomezDumm:2013sib, Eqs. (10), (15), Table 1]. The
central fit parameters are $M_\rho=0.8430\,\mathrm{GeV}$,
$F_\pi=0.0901\,\mathrm{GeV}$, $\alpha_1=1.87\,\mathrm{GeV}^{-2}$, and
$\alpha_2=4.29\,\mathrm{GeV}^{-4}$. The common isospin-limit charged masses
used are $m_\pi=0.13957039\,\mathrm{GeV}$ and
$m_K=0.493677\,\mathrm{GeV}$ [@ParticleDataGroup:2024cfk].

# Phase and form factor

For $P=\pi,K$, the one-loop function on the upper rim is
$$
A_P(s)=\log\frac{m_P^2}{M_\rho^2}+8\frac{m_P^2}{s}-\frac53
+\sigma_P^3(s)\log\frac{\sigma_P(s)+1}{\sigma_P(s)-1},\qquad
\sigma_P(s)=\sqrt{1-\frac{4m_P^2}{s}},
$$
with $\operatorname{Im}A_P(s+i0)=-\pi\sigma_P^3(s)$ above threshold.
The resummed form used to define the elastic phase is
$$
F_V^{\pi(0)}(s)=
\frac{M_\rho^2}
{M_\rho^2\left[1+\frac{s}{96\pi^2F_\pi^2}
\left(A_\pi(s)+\frac12 A_K(s)\right)\right]-s},
\qquad
\delta_1^1(s)=\arg F_V^{\pi(0)}(s).
$$
We use the source's central phase prescription: the elastic phase is used
through $s_1=0.98\,\mathrm{GeV}^2$, it is continued linearly to $\pi$ at
$s_2=M_\tau^2=(1.77686\,\mathrm{GeV})^2$, and it remains $\pi$ above
$s_2$. The three-subtracted Omnès form factor is
$$
F_\pi^V(s+i0)=\exp\left[
\alpha_1s+\frac{\alpha_2}{2}s^2+
\frac{s^3}{\pi}\operatorname{PV}\int_{4m_\pi^2}^{\infty}
\frac{\delta_1^1(s')}{s'^3(s'-s)}\,ds'
+i\delta_1^1(s)\right].
$$
The upper-rim imaginary part is the elastic phase; the normalization is
$F_\pi^V(0)=1$. The central prescription $s_1\simeq0.98\,\mathrm{GeV}^2$
and $s_2\simeq M_\tau^2$ follows the stated source choices. The constant
phase tail and integration to infinity are explicit model assumptions,
consistent with the source's listed variation ranges
[@GomezDumm:2013sib].

# Numerical result

The calculation stores 401 complex values from $x=4m_\pi^2$ to
$1.69\,\mathrm{GeV}^2$ in
[`pion_vector_form_factor_central.json`](../../../data/odd_parity/pion_vector_form_factor_central.json).
Selected values are
$$
\begin{array}{c|rr}
x\,[\mathrm{GeV}^2]&\operatorname{Re}F_\pi^V&\operatorname{Im}F_\pi^V\\\hline
0.10&1.24583&0.006504\\
0.30&2.20676&0.38007\\
0.60&-0.93774&5.08536\\
0.98&-1.24462&0.68683\\
1.50&-0.63005&0.25444
\end{array}
$$
The grid's density scale and complex phase follow directly from the
three-subtracted Omnès construction. The Python calculation reports a maximum
absolute change of $1.6\times10^{-8}$ at the selected points when tightening
the Cauchy-principal-value quadrature tolerances from $2\times10^{-10}$ to
$2\times10^{-12}$. This is a numerical integration diagnostic, not an
uncertainty on the physical form factor.

# Assumptions and limits

- Exact isospin is used, and the charged-pion/charged-kaon masses are used as
  common masses in the loop functions.
- The source fit is only described as reliable through the low-energy fit
  region, about $s\lesssim1.5\,\mathrm{GeV}^2$. Grid entries above that
  value are the declared low-energy phase continuation, not a separately
  validated inelastic form factor.
- No electromagnetic corrections, $\rho$--$\omega$ mixing, or excited
  vector contributions are included. Fit covariance and correlated replicas
  are not available here.
- This grid supplies $F_\pi^V$ only. The complex $f_1(x,y)$ input is still
  missing, so no double-spectral density or correlated spectral moments are
  determined yet.

# Reproducibility

The paired executable is
[`derive_pion_vector_form_factor.py`](../../../scripts/derive_pion_vector_form_factor.py).
Run `.venv/bin/python scripts/derive_pion_vector_form_factor.py` to regenerate
the grid. This command was executed successfully and wrote the saved JSON
using SciPy Cauchy-principal-value quadrature. The selected values and
tolerance diagnostic above are taken from that Python output.

# Inputs and provenance

- The Omnès representation, loop-derived phase, fit values and continuation
  ranges come from [GomezDumm:2013sib](../../extracts/GomezDumm:2013sib/INDEX.md).
- Particle-mass and tau-mass conventions use
  [ParticleDataGroup:2024cfk](../../extracts/ParticleDataGroup:2024cfk/INDEX.md).
- The downstream array shape and spectral-moment definitions are in the
  [pion TFF input contract](pion-tff-spectral-input-contract.md).
