---
title: Source audit for physical pion TFF spectral amplitudes
status: mixed
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Hoferichter:2018kwz, Hoferichter:2023bjm, Hoferichter:2025lcz, GomezDumm:2013sib, Niecknig:2012sj, Hanhart:2012wi]
provenance: literature
---

# Scope

This is a source-availability audit for the physical pion double-spectral
input. The numerical interface and correlated moments are specified in the
[spectral input contract](pion-tff-spectral-input-contract.md). The required
amplitudes are the complex pion vector form factor $F_\pi^V(x)$ and the
complex $P$-wave $f_1(x,y)$ for $\gamma^*(y)\to3\pi$, in one consistent
normalization and phase convention.

# Available local material

The arXiv source for Hoferichter:2018kwz gives the Khuri--Treiman equation
for the $\gamma^*\to3\pi$ partial wave, the dispersive construction of
$F_\pi^V$, and the fit strategy for the subtraction function $a(y)$. It does
not include a numerical complex $f_1(x,y)$ grid, complete paired fit
parameter sets with covariance, or a numerical solver. The publisher's
electronic supplement contains `TFF_full_band_ind_half.dat` and its
readme. That file tabulates the final spacelike doubly-virtual pion TFF and
separate uncertainty components; it is a validation benchmark, not either
of the two spectral amplitudes.

The local Niecknig:2012sj archive supplies the dispersive $V\to3\pi$
formalism and discusses Dalitz data for vector-meson decays. It does not
provide the off-shell $\gamma^*(y)\to3\pi$ partial-wave grid required here.
The local Hanhart:2012wi extract specifies the elastic Omnès structure for
$F_\pi^V$, but not the phase-shift tables and fitted variants needed to
generate numerical complex replicas. A central $F_\pi^V$ grid is now
calculated using the explicit phase model and exact-isospin fit I in
GomezDumm:2013sib; see
[the numerical derivation](pion-vector-form-factor-central-grid.md).
This does not reproduce the separate Bern/Madrid phase variants used in
Hoferichter:2018kwz, and it has no correlated uncertainty replicas. The
repository's existing
`PelaezRodas2016ElasticPhase` implementation is for $K\pi$ S waves and is
not a valid substitute for the $\pi\pi$ P wave.

Hoferichter:2023bjm supplies the dispersive $\gamma^*\to3\pi$ amplitude,
Khuri--Treiman equation, and form of the fitted normalization $a(y)$. Its
ancillary $\bar\eta_{3\pi}$ data are radiative-correction inputs, not
partial-wave values. Hoferichter:2025lcz updates the global SND, CMD-2, and
BaBar fit and reports variations with conformal-polynomial order and the
phase of $\rho$--$\omega$ mixing. Together these papers give a source-backed
central prescription to reproduce, but not a ready complex $f_1(x,y)$ grid
or complete correlated parameter replicas. The separate Belle II fits in
Hoferichter:2025lcz have poorer fit quality and resonance/anomaly tensions,
so they should remain a separate data scenario. See the
[source extract](../../extracts/Hoferichter:2025lcz/gamma-star-to-three-pion-fit.md)
[@Hoferichter:2025lcz].

# Consequence and next input

The central $F_\pi^V$ is available, but no full spectral density or moments
can be computed until the complex $f_1$ input is reconstructed. Do not infer
$f_1$ from a total $e^+e^-\to3\pi$ cross section alone, or combine
independent midpoints of published parameter ranges as if they were a
correlated fit replica.

The next executable milestone is a fixed-prescription central reconstruction
of $f_1(x,y)$ using the global SND/CMD-2/BaBar fit in Hoferichter:2025lcz,
with one declared $\pi\pi$ P-wave phase and continuation. It needs the
complex-contour prescription for crossed-channel integrals, numerical
Khuri--Treiman solution, full normalization parameter set, and energy-
dependent widths. Reuse the calculated central $F_\pi^V(x)$ grid. First
reproduce selected $e^+e^-\to3\pi$ cross sections and on-shell
$\omega,\phi\to3\pi$ normalizations and demonstrate grid convergence. Keep
the fixed central calculation, complete discrete model variants, and
statistical uncertainty propagation distinct. Do not call independent
marginal-error draws correlated replicas. Then populate the existing
[input contract](pion-tff-spectral-input-contract.md), calculate the four
moments, and compare the assembled spacelike TFF with the
Hoferichter:2018kwz supplementary table before matching the effective-pole
and asymptotic pieces. Use a common cutoff prescription to avoid counting
the same resonance strength twice.

# Provenance

- Hoferichter:2018kwz, source equations for the dispersive pion TFF and
  publisher supplementary data, indexed under
  [the source archive](../../sources/Hoferichter:2018kwz/INDEX.md).
- Hoferichter:2023bjm, dispersive $\gamma^*\to3\pi$ amplitude and fitted
  normalization, indexed under
  [the source archive](../../sources/Hoferichter:2023bjm/INDEX.md).
- Niecknig:2012sj, dispersive three-pion amplitude framework and its stated
  data scope, indexed under
  [the source extracts](../../extracts/Niecknig:2012sj/INDEX.md).
- Hanhart:2012wi, elastic Omnès representation for the pion vector form
  factor, indexed under
  [the source extracts](../../extracts/Hanhart:2012wi/INDEX.md).
