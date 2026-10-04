---
title: BaBar ISR Born cross-section definition and its radiative conventions
status: literature
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [BaBar:2004ytv]
provenance: literature
---

# Definition of the Tabulated `sigma` Column

The arXiv source archive is cached at `memory/sources/BaBar:2004ytv/` (`paper.tex`
plus `arxiv-source.tar`; the extracted `paper.tex` SHA256 matches the archive
member). Line numbers below refer to that `paper.tex`.

## Kinematic quantity

The tabulated variable is $\sqrt{s'}$, the invariant mass of the three-pion
system **after** initial-state radiation, not the beam energy $\sqrt s$. With
$x\equiv 2E_\gamma/\sqrt s$, one has $s'=s(1-x)$ (`paper.tex` Eq. (1), l. 786).

## The cross section is Born

The analysis fits the resolution-corrected spectrum with

$$
\frac{{\rm d}N}{{\rm d}m}=\sigma_{3\pi}(m)\frac{{\rm d}L}{{\rm d}m}R\,\varepsilon
\qquad (m\equiv\sqrt{s'}),
$$

where $\varepsilon$ is the mass-dependent detection efficiency and $R$ is a
radiative correction factor for multi-photon emission (l. 1643--1660). The
published cross section is obtained as

$$
\sigma_{3\pi}(m)=\frac{({\rm d}N/{\rm d}m)_{\rm corr}}{\varepsilon\,R\,{\rm d}L/{\rm d}m}.
$$

$R$ is **not** fitted per bin: it is a single generator-level structure-function
result, $R=0.9994$ for the spectrum fit (l. 1706) and $R=1.005\pm0.002\pm0.010$
for the $J/\psi$ sample (l. 2075). The ISR differential luminosity uses
$\theta_0=20^\circ$ (l. 1690).

So `sigma_nb` is the **Born** cross section with ISR already unfolded, with no
initial-state radiative correction remaining beyond $\approx 0.06\%$.

## Vacuum polarization is absorbed, not applied

This is the convention that matters for predicting the cross section. Verbatim
(l. 1711--1714):

> The radiative correction factor does not include the corrections due to
> leptonic and hadronic vacuum polarization. Here we follow the generally
> accepted practice of including the vacuum polarization correction in the
> resonance electronic width.

Consequently a prediction of this column must dress the vector coupling with
$\Pi(s=m_V^2)^{-1}$, equivalently multiplying $\Gamma(V\to e^+e^-)$ by
$|\Delta_V|^2=1/\Pi_V$ at the resonance pole. A *bare* ChPT coupling, or one
dressed with leptonic vacuum polarization only, does not reproduce the published
convention. This is a modeling requirement, not a cosmetic correction: at
$m_V=m_\phi$ the hadronic $\Pi$ is $O(10^3)$.

## Binning and coverage

Table `sumtab` (l. 1879--1932), caption verbatim:

> $e^+e^-$ c.m. energy ($\sqrt{s'}$), number of selected events after
> $3\pi$-mass-resolution correction ($N_{corr}$), detection efficiency
> ($\varepsilon$), differential ISR luminosity ($L$), and measured cross section
> ($\sigma$) for $e^+e^-\to \pi^+\pi^-\pi^0$. All values are calculated for 25
> MeV bin size. The quoted errors are statistical and systematic.

The 78 tabulated rows are bin centres on a uniform 25 MeV grid spanning
$1.0625$--$2.9875$ GeV, i.e. the full 1.05--3.00 GeV range.

**Neither the $\omega$ (782.7 MeV) nor the $\phi$ (1019.5 MeV) peak lies inside
this range.** The measured cross section is maximal at $\sqrt{s'}=1.3125$ GeV,
in the $\omega'(1350)$ region, not at a one-multiplet resonance.

## Quoted uncertainties

$\sigma$ errors combine event statistics, detection efficiency, and the
radiative correction (l. 1991--1998). The paper states that the systematic
errors for different mass bins are **fully correlated**, and that the resolution
deconvolution itself raises the errors by $\approx20\%$ and introduces
correlation. The $1.2\%$ luminosity term is treated separately. This matches the
[existing extract](three-pion-isr-spectrum.md) and the sidecar likelihood block;
no $\pi^0\pi^0/\pi^+\pi^-$ isospin factor is applied anywhere in this analysis, so
the column refers to the charged mode $\pi^+\pi^-\pi^0$ as written.