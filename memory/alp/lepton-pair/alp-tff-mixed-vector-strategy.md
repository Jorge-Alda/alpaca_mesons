---
title: Mixed dispersive rho and narrow-pole omega/phi strategy for the ALP TFF
status: mixed
created: 2026-10-27
updated: 2026-11-06
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

The first approximation retains the mixed $\rho I$ VVP term in its
cut-subtracted form. Let $F^{(\rho I),\mathrm{dress}}$ be the ordered mixed
term with a propagator-level dressed $\rho$ response on the $\rho$ leg and
the narrow-pole $\omega$ or $\phi$ response on the other leg. Define
$F^{(\rho I),\pi\pi,[N]}(\Lambda)$ by replacing each dressed $\rho$ response
with its $N$-subtracted $\pi\pi$ dispersive projection up to the same
$\Lambda$ used in the two-pion construction. Define
$F^{(\rho I),\mathrm{tail},[N]}(\Lambda)$ by the same replacement using the
projection above $\Lambda$. The retained mixed remainder is
$$
\widetilde F^{(\rho I),[N]}(\Lambda)=
F^{(\rho I),\mathrm{dress}}
-F^{(\rho I),\pi\pi,[N]}(\Lambda)
-F^{(\rho I),\mathrm{tail},[N]}(\Lambda).
$$
This removes the mixed term's full $\rho\to\pi\pi$ cut from its explicit
addition. The physical two-pion single and double dispersive contributions
remain and carry the $\pi\pi$ cut once. The remainder contains the
narrow-isoscalar-pole and subtraction part; assign it once. In the
low-energy first approximation omit the above-$\Lambda$ tail, or include
it only as a separately matched higher-energy component. Do not add the
unmodified $F^{(\rho I),\mathrm{dress}}$ on top of a two-pion construction
that already carries its $\rho$ cut.

The dressed response uses
$W_\rho^{\mathrm{dress}}(s)=R_\rho/[M_{0\rho}^2-s-\Delta_\rho^{\mathrm{self}}(s)]$,
with $\operatorname{Im}\Delta_\rho^{\mathrm{self}}(s+i0)
=+\sqrt{s}\,\Gamma_{\rho\to\pi\pi}(s)$ and a subtracted real part. It is
neither a constant-width denominator nor the complete source-specific
two-pion TFF. The pure narrow-pole $\omega,\phi$ terms remain separate and
are not subtracted by this $\rho I$ prescription. See the
[symbolic cut-replacement derivation](alp-tff-rho-isoscalar-overlap.md)
and [mixed VVP record](alp-tff-mixed-rho-isoscalar-vvp.md).

The [two-pion spectral construction](alp-isovector-tff-spectral-construction.md)
defines the physical low-energy cut. Its present script is a symbolic
scaffold, not a numerical first-pass assembly. The source-resolved pion
current and production inputs, subtraction functions, and cutoff are not
available for a numerical TFF.
The didactic explanation of cuts, thresholds, widths, and the pole picture
is recorded in [the TFF dispersive-treatment note](alp-tff-dispersive-treatment-explained.md).

# Deferred extensions

Quantifying the error of the narrow $\omega,\phi$ pole terms and replacing
them with dispersive finite-width treatments are future extensions. A full
dispersive treatment of all vectors also requires channel amplitudes,
subtractions, phases, and a pole/continuum prescription that avoids double
counting. The present decision establishes the first deliverable's modeling
scope; it does not claim that the mixed strategy is numerically validated.
