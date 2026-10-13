---
title: Stage-3 vector mass and electromagnetic-residue calibration
status: numerical
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [ParticleDataGroup:2024cfk, Alda:2026nrd, Ecker:1988te]
provenance: numerical-computation
---

# Stage-3 calibration and adequacy tests

## Inputs, scope and conventions

Use the [PDG 2024 snapshot and provenance](../../extracts/ParticleDataGroup:2024cfk/vector-sector-calibration-inputs.md),
the [arXiv-v1 Alda inputs](../../extracts/Alda:2026nrd/neutral-meson-generators.md),
the [minimal mass action](vector-nonideal-mass-action.md), the
[additional mass invariant](vector-mass-extension.md), and the
[leptonic-width/source derivation](vector-leptonic-width.md).
All generators, rotations and source phases retain the positive-strange phi
convention; angles below are departures from the ideal basis in that convention.

**Assumptions:** tree level, exact strong isospin, canonical tensor kinetic terms,
fixed on-shell alpha and electron mass, and diagonal Gaussian summary errors.
PDG listing fitted masses are used as tree masses without electromagnetic or
complex-pole conversion. Rho and K-star charged/neutral scenarios remain separate;
the alternative phi electronic width is also separate. No radiative or three-pion
data are used. Direct widths are listing summaries, not three independent raw
measurements; their full cross-covariance is unavailable. Chi-squared values below
are conditional adequacy diagnostics, not precision global-fit significances.

Set $x_l=4B_0\widehat m=0.036616$ GeV squared and
$x_s=4B_0m_s=0.943684$ GeV squared from Alda's central parameters.
Spurion errors are fixed in this diagnostic; quoted propagated errors below
include PDG input errors only. They do not include chiral truncation, omitted
self-energies, source corrections or correlations. No model error is silently
added to an input variance.

## Reproduction and likelihood

Run [the calibration script](../../../scripts/calibrate_vector_stage3.py):

```text
.venv/bin/python scripts/derive_vector_mass_extension.py
.venv/bin/python scripts/derive_vector_leptonic_width.py
.venv/bin/python scripts/calibrate_vector_stage3.py
```

The script verifies the four source PDF hashes, uses SymPy-generated mass-matrix
and determinant kernels, minimizes the mass residuals with SciPy, evaluates both
exact extended-action branches, and fits the electronic widths. The deterministic
[results JSON](../../../data/vector_sector/stage3_results.json) pins its input and
script hashes. Natural units are GeV internally; output masses and widths use
MeV and keV.

For each declared block, minimize

$$
\chi^2=\sum_i\left[\frac{f_i(\vartheta)-y_i}{\sigma_i}\right]^2.
$$

The minimal mass parameters are $(m_\rho,m_{K^*},\Delta_V)$ and the four
predictions are those two masses plus the square roots of the neutral eigenvalues.
Multiple starting singlet shifts give the retained minimum. The extended mass
action solves the trace/determinant equations exactly. Local mass-input errors
are propagated to each branch with a central-difference delta method.

At fixed angle, the universal-source prediction is linear in $F_V^2$; its positive
weighted optimum is calculated directly. Then fit positive $(F_l,F_s)$ for the
single-trace spurion extension, and separately test positive light/negative
strange source strengths. The common overall sign fixes the light strength
positive; no measured width fixes the relative strange sign. Finally interpolate three signed residues with
$(F_V,z_V,y_V)$ at fixed angle. The last step fixes the signs to those of the
universal source at the chosen angle: widths alone cannot establish those signs.

## Minimal mass action fails the quoted-error diagnostic

| Rho category | K-star category | Conditional mass chi-squared (one nominal degree of freedom) | Fitted angle |
|---|---|---:|---:|
| neutral | charged | 1640.21 | $1.925^\circ$ |
| neutral | neutral | 1097.03 | $1.905^\circ$ |
| charged | charged | 1232.34 | $2.418^\circ$ |
| charged | neutral | 744.45 | $2.263^\circ$ |

For neutral rho/charged K-star, the least-squares masses differ from the inputs
by $(-4.739,+8.989,+0.564,-0.019)$ MeV for rho, K-star, omega, phi.
For neutral rho/neutral K-star they differ by
$(-4.623,+5.191,+0.549,-0.019)$ MeV. These fitted angles are not adopted as
physical estimates: the action does not describe the selected inputs at their
quoted errors. This does not establish failure of a general vector EFT or
determine the omitted model uncertainty.

## One additional mass invariant gives two branches

The assumed $\kappa_V\langle V\rangle\langle V\chi_+\rangle/4$ operator
can interpolate all four central masses, but gives two branches. With neutral rho:

| K-star input | Near-ideal branch | Other branch |
|---|---:|---:|
| charged | $5.174\pm0.086^\circ$ | $-75.703\pm0.086^\circ$ |
| neutral | $3.912\pm0.068^\circ$ | $-74.441\pm0.068^\circ$ |

Errors are conditional PDG-input propagation only. The near-ideal branch has
$(\kappa_V,\Delta_V)=(0.043224,0.017751\ {\rm GeV}^2)$ or
$(0.030389,0.016935\ {\rm GeV}^2)$ for charged or neutral K-star respectively.
The large-angle branches require roughly $\kappa_V=-0.61$ and
$\Delta_V=0.68$ GeV squared, so treating them as small OZI corrections would be
an additional and questionable assumption. Both branches are retained.
The user subsequently selected proximity to ideal mixing as the tie-breaker
when the other criteria are equivalent; both branches remain recorded and
available, and this preference does not establish a measured angle.
The 1.26-degree charged/neutral change
remains separate from the much smaller conditional input errors.

## Electromagnetic residues require more than a common source

With the primary direct widths, the near-ideal extended mass branches give:

| K-star mass scenario (neutral rho) | Universal $F_V$ chi-squared (two nominal degrees of freedom) | Single-trace source-breaking chi-squared (one nominal degree of freedom) |
|---|---:|---:|
| charged | 356.67 | 341.23 |
| neutral | 269.07 | 253.12 |

At the charged scenario's mass angle, independent effective $F_V$ values inferred
from rho, omega, phi widths are $(0.15642,0.12252,0.17338)$ GeV; for the neutral
scenario they are $(0.15642,0.12580,0.17022)$ GeV. Hence the common-source hypothesis
does not fit those fixed angles, and the single-trace spurion correction alone
does not repair the positive light/strange source-strength fit. The other mass
branches likewise fail that restricted source diagnostic; all results remain
in the JSON.

The opposite strange-source sign must not be silently excluded. On the
near-ideal mass branches, allowing $F_l>0,F_s<0$ yields chi-squared 0.28 and
2.18 for charged and neutral K-star respectively (one nominal degree of freedom).
The charged scenario gives $(F_l,F_s)=(0.15652,-0.15207)$ GeV,
$F_V=0.16898$ GeV and $z_V=-0.17010$ inverse GeV; the neutral scenario gives
$(0.15614,-0.15484)$ GeV, $F_V=0.16870$ GeV and
$z_V=-0.17142$ inverse GeV. Such a reversal requires an order-one strange
spurion correction rather than a small perturbation around the universal
positive source. The widths cannot reject this branch on sign grounds.
Consequently the double-trace source below is a small-breaking candidate,
not an operator proven necessary by widths without a sign/truncation assumption.

If the angle instead floats independently of the mass action within
$[-20^\circ,20^\circ]$, the primary widths alone give
$F_V=0.15616\pm0.00063$ GeV and $\theta_V=-4.45\pm0.53^\circ$, with conditional
chi-squared 1.41 for one nominal degree of freedom. The alternative phi average
gives $F_V=0.15561$ GeV, $\theta_V=-3.86^\circ$, and chi-squared 5.85.
This is a restricted width-only diagnostic, not a consistent mass/source fit or
a new physical mixing-angle prior. The opposite-sign angle requirement exposes
the incompatibility of these restricted action choices under the adopted inputs.

Adding the separately declared spurion double-trace source operator supplies
the missing residue direction. At the charged/neutral K-star near-ideal branches,
the primary-width interpolations give respectively

$$
(F_V,z_V,y_V)=(0.156692,-0.003672,0.022020),\qquad
(0.156661,-0.003245,0.019121),
$$

with $F_V$ in GeV and $z_V,y_V$ in inverse GeV. This is exact interpolation of
three widths by three coefficients at fixed mass angle, with zero remaining
degrees of freedom. It supplies no independent validation of the angle.
The unbroken singlet source trace alone cannot do this because
$\operatorname{Tr}Q=0$. Relative residue signs remain unmeasured.

## Stage-3 outcome and stage-4 requirements

The minimal universal model is inadequate for this snapshot. Explicit additional
mass and source invariants can accommodate the central inputs within a declared
larger truncation, but produce a conditional prescription rather than a unique
precision omega–phi mixing determination. These are assumed extensions, not
literature-established fitted coefficients or a complete operator basis.

Stage 4 should retain the mass branches, both tested relative source-sign
branches, and source coefficients explicitly,
derive their contributions to every photon/vector slot, and verify the baseline
and Ward identities before recomputing odd-parity sensitivity. Do not substitute
independent effective $F_V$ values into a universal action without its operators.
Branch/sign discrimination requires additional amplitude data and the declared
interaction-breaking model; a general kinetic/self-energy treatment remains
outside this stage. No pooled charged/neutral estimate, decay-fitted prior,
physical-angle confidence interval or odd-parity fit is established here.
