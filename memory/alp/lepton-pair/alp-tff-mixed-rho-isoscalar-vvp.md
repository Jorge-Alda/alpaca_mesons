---
title: Mixed rho-isoscalar VVP term in the ALP TFF mixed vector strategy
status: mixed
created: 2026-10-28
updated: 2026-10-28
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: assumption
---

# Scope and conclusion

In the first-deliverable mixed strategy, the low-energy isovector TFF is
represented by the source-specific two-pion double-spectral contribution.
The low-energy isoscalar contribution includes the pure narrow-pole
$\omega,\phi$ VJP/VVP terms and, for general diagonal ALP sources, the
mixed $\rho I$ VVP term, where $I=\omega,\phi$. The mixed term is generally
nonzero and must not simply be dropped: doing so omits retained source
contributions. It vanishes in the exact strong-isospin limit for the ALP
source traces specified in the [neutral-rho pole record](alp-isovector-doubly-virtual-tff.md).

In this selected approximation, there is no separate stable-pole
$F^{(\rho)}$ contribution. The isovector low-energy piece is the two-pion
double-spectral term. The modeled isoscalar resonance piece is the pure
$\omega,\phi$ pole selection plus the mixed $\rho I$ VVP contribution,
with the $\rho$ leg dressed as below. Do not subtract the mixed term from
the pure $\omega,\phi$ selection, because it is not contained in that
selection. The mixed term is part of the model and must be retained.

The dressed mixed term and the two-pion double-spectral contribution can
share $\rho$-mediated strength. Their quantitative matching must ensure
that this strength is represented once. The needed source-specific matching
is not yet available, so the displayed organization defines the intended
approximation but not a completed numerical sum.

# Rho response in the mixed term

Do not use the stable $W_\rho(s)=R_\rho/(M_\rho^2-s)$ in the mixed term if
the broad rho is being treated with its $\pi\pi$ dynamics. Also, a constant
on-shell width in the denominator is only a near-resonance line-shape
approximation; it is not the dispersive treatment selected for the rho.

At the propagator-response level, use
$$
W_\rho^{\mathrm{dress}}(s)=
\frac{R_\rho}{M_{0\rho}^2-s-\Delta_\rho^{\mathrm{self}}(s)},
\qquad
\operatorname{Im}\Delta_\rho^{\mathrm{self}}(s+i0)
=\sqrt{s}\,\Gamma_{\rho\to\pi\pi}(s),
$$
with the real part of $\Delta_\rho^{\mathrm{self}}$ fixed by a subtracted
dispersion relation and subtraction conditions fixing the chosen physical
mass and residue convention. Use this same response in each rho leg of the
ordered mixed term, while the other vector leg uses the selected narrow-pole
$\omega$ or $\phi$ response. The response-level self-energy dressing is a
model for the rho leg; it is not identical to the full source-specific
two-pion double-spectral TFF, whose density also depends on the pion vector
form factor and the ALP production amplitude [@Messerli:2025rnv].

# Selected-model bookkeeping

Using the notation of the draft, write the modeled isoscalar resonance
piece as
$$
F_{a\gamma^*\gamma^*}^{I=0,\,\mathrm{res}}(s_1,s_2)
=F_{a\gamma^*\gamma^*}^{(\omega,\phi)}(s_1,s_2)
+F_{a\gamma^*\gamma^*}^{(\rho I),\,\mathrm{dress}}(s_1,s_2),
$$
where
$$
F_{a\gamma^*\gamma^*}^{(\rho I),\,\mathrm{dress}}
=-2\sum_{I=\omega,\phi}\left[
W_\rho^{\mathrm{dress}}(s_1)W_I^{\mathrm{pole}}(s_2)P_{\rho I}(s_1,s_2)
+W_I^{\mathrm{pole}}(s_1)W_\rho^{\mathrm{dress}}(s_2)P_{I\rho}(s_1,s_2)
\right].
$$
Here $F^{(\omega,\phi)}$ is the pure $\omega,\phi$ term defined in the
draft; it does not include the mixed term. This model contribution is only
the resonance part of the state-defined isoscalar TFF; the $3\pi$,
$K\bar K$, higher-state and subtraction contributions remain separate.
