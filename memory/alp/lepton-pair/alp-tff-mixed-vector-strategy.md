---
title: Mixed dispersive rho and narrow-pole omega/phi strategy for the ALP TFF
status: mixed
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv]
provenance: assumption
---

# Decision and scope

The first TFF deliverable uses a dispersive two-pion treatment for the broad
isovector $\rho$ and narrow-pole terms for the isoscalar $\omega$ and $\phi$.
This is a declared approximation strategy, not an exact decomposition of
the amplitude. In general, vector resonances occur within hadronic
continua: poles in an idealized stable-particle limit and cuts from on-shell
states are related aspects of the analytic structure. The $\rho$ is
especially tied to the $\pi\pi$ cut, whose phase and rescattering produce
its broad line shape [@Messerli:2025rnv].

The stable $\omega$ and $\phi$ poles approximate their spectral strength.
Their narrow widths motivate this treatment when the observable does not
resolve their timelike line shapes and does not require precise phases.
Narrow widths alone do not guarantee accuracy near the corresponding
timelike peaks, nor do pole terms reproduce the full spectral tails. No
quantitative error estimate is part of this first deliverable.

# Matching rule

The selected low-energy isovector piece is the source-specific two-pion
double-spectral contribution. The modeled isoscalar resonance piece
contains the pure narrow-pole $\omega,\phi$ terms and the source-dependent
mixed $\rho I$ VVP contribution, $I=\omega,\phi$. Retain the mixed term;
dropping it omits a generally nonzero contribution. In that term use a
propagator-level dressed $\rho$ response on the $\rho$ leg and the selected
narrow-pole $\omega$ or $\phi$ response on the other leg. The dressed
response has a $\pi\pi$ self-energy whose imaginary part is
$\sqrt{s}\,\Gamma_{\rho\to\pi\pi}(s)$ and whose real part follows from a
subtracted dispersion relation. It is neither a constant-width denominator
nor the complete source-specific two-pion double-spectral TFF. Match the
mixed term to the two-pion continuum so the same $\rho$ strength is counted
once; that numerical matching remains open. See the
[mixed VVP record](alp-tff-mixed-rho-isoscalar-vvp.md).

The [two-pion spectral construction](alp-isovector-tff-spectral-construction.md)
defines the isovector low-energy component. The formal
[finite-width self-energy record](alp-isoscalar-vector-self-energy.md)
shows how open channels generate an absorptive self-energy and spectral
response, but does not yet supply a numerical common prescription. In the
selected approximation, add the dressed mixed $\rho I$ term to the pure
$\omega,\phi$ pole contribution; do not subtract it from the pure pole
term, which does not contain it.
The didactic explanation of cuts, thresholds, widths, and the pole picture
is recorded in [the TFF dispersive-treatment note](alp-tff-dispersive-treatment-explained.md).

# Deferred extensions

Quantifying the error of the narrow $\omega,\phi$ pole terms and replacing
them with dispersive finite-width treatments are future extensions. A full
dispersive treatment of all vectors also requires channel amplitudes,
subtractions, phases, and a pole/continuum prescription that avoids double
counting. The present decision establishes the first deliverable's modeling
scope; it does not claim that the mixed strategy is numerically validated.
