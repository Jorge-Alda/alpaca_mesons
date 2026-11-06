---
title: Mixed rho-isoscalar VVP term in the ALP TFF mixed vector strategy
status: mixed
created: 2026-10-28
updated: 2026-11-06
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: assumption
---

# Scope and conclusion

The first approximation uses the physical source-resolved two-pion
single- and double-cut reconstruction for the broad isovector $\rho$,
the pure narrow-pole $\omega,\phi$ terms, and a retained mixed $\rho I$
VVP remainder for general diagonal ALP sources. This mixed source
contribution is generally nonzero; it vanishes in the exact strong-isospin
limit for the source traces specified in the
[neutral-rho pole record](alp-isovector-doubly-virtual-tff.md).

The full dressed mixed term is not added on top of the physical two-pion
cut. Retain it as
$$
\widetilde F^{(\rho I),[N]}(\Lambda)
=F^{(\rho I),\mathrm{dress}}
-F^{(\rho I),\pi\pi,[N]}(\Lambda)
-F^{(\rho I),\mathrm{tail},[N]}(\Lambda),
$$
with both ordered photon placements. The low projection replaces each
dressed rho response by its $N$-subtracted $\pi\pi$ spectral integral up to
the same $\Lambda$ used for the two-pion reconstruction. The tail
projection uses the integral above $\Lambda$. Omit that tail from the
low-energy first approximation or retain it only as a separately matched
higher-energy term. The physical two-pion cut remains once; the mixed
remainder is an isoscalar-pole/subtraction contribution assigned once.
There is no separate stable-rho pole. The pure $\omega,\phi$ selection is
not part of the $\rho I$ subtraction.

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

Write the selected first-approximation TFF as
$$
F_{a\gamma^*\gamma^*}^{\mathrm{first}}(s_1,s_2)
=F_{a\gamma^*\gamma^*}^{2\pi,\mathrm{cut-only},[N]}(s_1,s_2;\Lambda)
+F_{a\gamma^*\gamma^*}^{(\omega,\phi),\mathrm{pole}}(s_1,s_2)
+\widetilde F_{a\gamma^*\gamma^*}^{(\rho I),[N]}(s_1,s_2;\Lambda),
$$
where
$$
F_{a\gamma^*\gamma^*}^{(\rho I),\,\mathrm{dress}}
=-2\sum_{I=\omega,\phi}\left[
W_\rho^{\mathrm{dress}}(s_1)W_I^{\mathrm{pole}}(s_2)P_{\rho I}(s_1,s_2)
+W_I^{\mathrm{pole}}(s_1)W_\rho^{\mathrm{dress}}(s_2)P_{I\rho}(s_1,s_2)
\right].
$$
The cut-only two-pion term contains the physical single- and double-cut
reconstruction with the stated subtractions. The pure narrow-pole
$\omega,\phi$ term is separate. The dressed mixed term in the second
display is an intermediate expression; the contribution in the first
display is its cut-subtracted remainder defined above. Standalone $3\pi$,
$K\bar K$, higher-state and unmatched short-distance contributions remain
outside this first approximation.
