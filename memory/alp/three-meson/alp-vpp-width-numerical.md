---
title: Numerical implementation of selected ALP vector-two-pseudoscalar widths
status: numerical
created: 2026-10-21
updated: 2026-10-22
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: numerical-computation
---

# Scope and inputs

The selected-channel calculator is
[`scripts/implementation/alp_vpp_width.py`](../../../scripts/implementation/alp_vpp_width.py).
It implements the amplitude and width in the [strong amplitude record](alp-vpp-amplitude.md)
and [invariant width derivation](alp-vpp-width.md), with flavor support from
the [selection-rule derivation](alp-vpp-selection-rules.md). Each call evaluates
one ordered state $(V,P_1,P_2)$, with $V=\omega$ or $\phi$, and returns a
partial width in GeV. Couplings $c_q^A$, $c_q^P$, and the optional $c_q^{VS}$
are supplied directly in the no-mixing basis. $m_a$, $f_a$, and $G_V$ are in
GeV; $G_V$ is required explicitly because the provisional odd-sector fit does
not determine it.

The original-eleven provisional fit supplies $d_1,d_2,d_3$, the fitted neutral
vector generators and masses, charged-vector mass parameters, and shared
total widths. The external vector uses its physical mass in phase space. The
pseudoscalar fields and generators now follow Alda:2026nrd directly:
$f_\pi=93\,\mathrm{MeV}$,
$\Phi=\sum_i\lambda_i\phi_i$,
$\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, and
$u=\exp(i\Phi/(2f_\pi))$. Consequently the charged and neutral external and
internal pseudoscalar generators passed to the amplitude all have the same
Alda normalization; the implementation does not divide the charged
generators by $\sqrt2$ or use the RChT-normalized $T_P=\lambda_P/\sqrt2$ as
though it were an Alda generator. Any conversion to the RChT vertex basis is
made explicitly in the vertex factors. Charged and neutral kaons share the
stated LO mass.

# Numerical method

The script builds coherent coefficients $F_{\rm odd},F_1,F_2$ for the
polarization-summed invariant amplitude recorded in [the width derivation](alp-vpp-width.md).
It evaluates the two-dimensional Dalitz integral with nested adaptive
quadrature, including internal vector-pole locations as integration
subdivision points. Identical pseudoscalars receive the $1/2!$ factor after
Bose-symmetric exchange terms have been assembled.

The internal vector denominators use the constant-width prescription
$r^2-M_R^2+iM_R\Gamma_R$. This is an implementation choice; the amplitude
record leaves the line shape to the calculation. No rescattering, loops,
electromagnetic conversion, or additional vector-mass-spurion interaction is
included. The external vector is treated as a narrow on-shell state.

Charge and threshold checks return zero before integration. Flavor terms that
share the same coefficient, invariant, pole mass, and width are combined
before quadrature. Their sum is set to zero only when it satisfies the
scale-relative $64\epsilon_{\rm mach}\sum_i|c_i|$ cancellation criterion
already used for the $\pi^0K_LK_S$ vector-width slowdown. This operates on the
coherent physical $K_S/K_L$ generators, so no channel-name-specific kaon
shortcut is used. $K_S$ and $K_L$ follow the repository's $\epsilon_K$
convention and equal-mass approximation.

# Alda normalization update (2026-10-22)

The width phase-space expression and polarization sum are unchanged. Their
amplitude input is now built with the corrected Alda-normalized pseudoscalar
generators and the outgoing source
$i\chi_-|_{P_i}=+\{\lambda_i^\dagger,\chi\}/f_\pi$ from the
[Appendix A derivation](chi-minus-meson-source-convention.md). This correction
changes the crossed odd-exchange mass-source terms as well as the generator
prefactors; it must be propagated through both the standalone implementation
and the library function. See the [strong amplitude record](alp-vpp-amplitude.md)
for the complete coefficients and the explicit RChT conversion.

# Verification and limits

The focused convention tests in
`scripts/implementation/test_alp_vpp_width.py` pass. They check that every
external and internal pseudoscalar generator has trace norm two in both
implementations and compare a pure crossed $d_2$ contribution against an
independent trace evaluation in both the standalone calculator and library.
For that cross-check, the library's charged-vector width is set to the
standalone calculator's fitted width; the two interfaces otherwise retain
their existing distinct width inputs (PDG widths in the library and shared fit
widths in the standalone calculator). The invariant width derivation remains
unchanged. Earlier checks covered forbidden states, the CP-conserving
$K_SK_L$ projection, exchange of $P_1,P_2$, and a direct rest-frame
polarization contraction, but were not rerun as part of this convention update.
A production width has not been numerically benchmarked after the update. The
result inherits the tree truncation and line-shape limitations stated above.
