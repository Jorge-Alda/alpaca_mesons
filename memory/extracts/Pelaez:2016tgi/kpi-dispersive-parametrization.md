---
title: Dispersively constrained pi K amplitude through 1.6 GeV
status: literature
created: 2026-09-23
updated: 2026-09-23
bibtex_keys:
  - Pelaez:2016tgi
provenance: literature
---

The paper provides parametrizations of $\pi K$ partial waves from threshold through $1.6\,\mathrm{GeV}$ constrained by forward dispersion relations and threshold sum rules. The complex partial waves encode the phase shift and, where the parametrization is inelastic, the departure from elastic unitarity. [Pelaez:2016tgi]

The authors restrict their claimed simultaneous data and dispersive consistency to energies at or below $1.6\,\mathrm{GeV}$. [Pelaez:2016tgi]

For the elastic $I=1/2$ scalar wave, the conformal parametrization is

$$
\cot\delta_0^{1/2}(s)
=\frac{\sqrt{s}}{2q(s-s_A)}\left(B_0+B_1\omega(s)\right),
$$

with

$$
\omega(s)=
\frac{\sqrt{y(s)}-\alpha\sqrt{y_0-y(s)}}
{\sqrt{y(s)}+\alpha\sqrt{y_0-y(s)}},
\qquad
y(s)=\left(\frac{s-\Delta_{K\pi}}{s+\Delta_{K\pi}}\right)^2.
$$

The paper uses $m_\pi=139.57\,\mathrm{MeV}$, $m_K=496\,\mathrm{MeV}$, $m_\eta=547\,\mathrm{MeV}$, $\alpha=1.15$, and conformal endpoint $s_0=(1.1\,\mathrm{GeV})^2$. Its printed Adler-zero equation is

$$
s_A=\frac{\Sigma_{K\pi}+2\sqrt{\Delta_{K\pi}^2+m_K^2m_\pi^2}}{5}.
$$

With the paper's definition $\Delta_{K\pi}=m_K^2-m_\pi^2$, this expression does not reproduce its quoted $s_A\simeq0.236\,\mathrm{GeV}^2$. The later review `Pelaez:2020gnd` prints $4\Delta_{K\pi}^2$ under the square root, which does reproduce the quoted value. The numerical implementation uses the later expression and treats this as resolution of an apparent typographical inconsistency. [@Pelaez:2020gnd]

The central constrained-fit parameters are $B_0=0.411$ and $B_1=0.162$; the unconstrained values are $B_0=0.411$ and $B_1=0.181$. This elastic expression is used only through the $K\eta$ threshold, above which the paper switches to an inelastic parametrization. [Pelaez:2016tgi, Eqs. for $q_{12}$, the conformal variable, and $\cot\delta_0^{1/2}$; elastic $S^{1/2}$ parameter table]

For the elastic $I=3/2$ scalar wave,

$$
\cot\delta_0^{3/2}(s)
=\frac{\sqrt{s}}{2q(s_A-s)}
\left(B_0+B_1\omega(s)+B_2\omega(s)^2\right),
$$

with $s_A=\Sigma_{K\pi}$, $\alpha=1.4$, and conformal endpoint $s_0=(1.84\,\mathrm{GeV})^2$. The constrained central values are $(B_0,B_1,B_2)=(2.27,3.94,3.36)$ and the unconstrained values are $(2.25,4.21,2.45)$. The analysis reports no evidence for inelasticity in this wave up to approximately $1.8\,\mathrm{GeV}$. [Pelaez:2016tgi, Eq. for $\cot\delta_0^{3/2}$ and $S^{3/2}$ parameter table]
