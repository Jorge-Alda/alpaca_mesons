---
title: Numerical implementation of selected ALP vector-two-pseudoscalar widths
status: numerical
created: 2026-10-21
updated: 2026-10-21
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
pseudoscalar exchanges use tree masses and normalized generators from the
Alda:2026nrd inputs; charged and neutral kaons share the stated LO mass. The
generator conversion follows the [Alda:2026nrd normalization record](../../extracts/Alda:2026nrd/neutral-meson-generators.md):
the charged generators inherited from `three_meson_width` are divided by
$\sqrt2$, while the neutral $T_P$ matrices already include this conversion.

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

# Verification and limits

Four focused checks were run during implementation: forbidden states avoid
quadrature, the CP-conserving $K_S K_L$ projection agrees with its flavor
amplitude, exchanging $P_1,P_2$ preserves the squared amplitude with the
corresponding invariant relabeling, and the invariant polarization sum agrees
with a direct rest-frame contraction. The temporary test file was removed;
the final input-validation edits were not rerun through those checks. A
production width has therefore not been numerically benchmarked after the
final edits. The result inherits the tree truncation and line-shape limitations
stated above.
