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

Include the $\rho$ dynamics once. If the dispersive two-pion production
amplitude and pion vector form factor contain the $\rho$, do not also add a
stable $\rho$ pole for the same contribution. For mixed $\rho$--isoscalar
terms, use the dispersive response on the $\rho$ leg and specify whether the
other leg is the $\omega$ or $\phi$ narrow-pole response. Keep any residual
terms separate until their overlap has been established.

The stable-pole TFF expression in the
[neutral-rho pole record](alp-isovector-doubly-virtual-tff.md) is a tree-level
component with widths and hadronic cuts omitted; the
[two-pion spectral construction](alp-isovector-tff-spectral-construction.md)
is a separate dispersive component. Combining them requires the matching
rule above rather than direct addition. The formal
[finite-width self-energy record](alp-isoscalar-vector-self-energy.md)
shows how open channels generate an absorptive self-energy and spectral
response, but does not yet supply a numerical common prescription.
The didactic explanation of cuts, thresholds, widths, and the pole picture
is recorded in [the TFF dispersive-treatment note](alp-tff-dispersive-treatment-explained.md).

# Deferred extensions

Quantifying the error of the narrow $\omega,\phi$ pole terms and replacing
them with dispersive finite-width treatments are future extensions. A full
dispersive treatment of all vectors also requires channel amplitudes,
subtractions, phases, and a pole/continuum prescription that avoids double
counting. The present decision establishes the first deliverable's modeling
scope; it does not claim that the mixed strategy is numerically validated.
