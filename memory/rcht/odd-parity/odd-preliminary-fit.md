---
title: Focused preliminary odd-parity fit and input sensitivities
status: mixed
created: 2026-10-19
updated: 2026-10-22
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, Ecker:1988te, Ecker:1989yg]
provenance: numerical-computation
---

# Scope and input choices

The production fit now exports only the central negative single-trace result.
The working calibration uses the original seven VJP and four VVP operators
[@Ruiz-Femenia:2003jdx], the twelve-summary primary dataset, charged-rho and
charged-K-star mass inputs, the evaluated phi electronic width, and the
near-ideal mass root. These are explicit pragmatic input choices for a
conditional preliminary result, not an experimentally determined
isospin-limit mass prescription. The charged and neutral categories are
documented in the [PDG input extract](../../extracts/ParticleDataGroup:2024cfk/vector-sector-calibration-inputs.md).
Other source hypotheses and sensitivity cases remain in the research archives;
they are not exported in the production fit JSON. The relative strange-source
sign is not inferred from electronic widths alone.

The original eleven-operator model is the smallest existing unsaturated
model. Its measured image agrees with the soft-leading model at the
archived fitted nuisance points, as established in the
[calibration record](odd-low-energy-calibration.md). The saturated extensions
cannot select their operator content from these twelve rows. The present
choice therefore gives identifiable low-energy combinations in the original
operator notation, without treating the larger bounded operator enumeration
as a request to fit every coefficient.

# Inputs, assumptions and reproducibility

Use the [calibration archive](odd-low-energy-calibration.md),
[finite profiles](odd-low-energy-profiles.md),
[normalized kernels](odd-bounded-sm-numerical-kernels.md),
[measurement Jacobians](odd-low-energy-measurement-jacobians.md), and
[canonical tensor conversion](levi-civita-convention.md).
The source-backed [PDG measurement extract](../../extracts/ParticleDataGroup:2024cfk/odd-parity-low-energy-measurements.md)
specifies every fitted observable. Alda:2026nrd's central neutral inputs are retained
[@Alda:2026nrd]. Bibliography records were inspected before using these
existing research results.

The inherited assumptions are tree level, fixed Alda:2026nrd neutral inputs,
diagonal quoted-summary covariance where the experimental covariance is
unavailable, shared total widths included once, LO source-compatible
vertex masses with physical phase-space masses, and the documented partial
isospin and broad-rho approximations. There is no OPE, naturalness or OZI
coefficient prior. Units are GeV, with dimensionless original $c_i,d_i$
and the fixed operator normalization $M_V=0.770$ GeV. Canonical displayed
odd coefficients use $\epsilon^{0123}=+1$ and the recorded sign adapter for
archived opposite-convention coefficients.

The focused artifact reuses the best stored valid result for each requested
case. Reusing an archived optimum does not constitute a fresh minimization,
an independent global-minimum proof, or a new experimental determination.
Focused input changes retain their separate scenario labels. Removing
broad-rho or constrained-fit K-star rows is a sensitivity diagnostic; it
does not reconstruct their missing experimental covariance.

# Interpretation and next step

Regenerate the conditional $G_V$ value and export the central fit:

```text
.venv/bin/python scripts/fit_vector_gv.py data/vector_sector/gv_measurements.json data/vector_sector/gv_conditional_fit.json
.venv/bin/python scripts/summarize_odd_preliminary_fit.py
```

The [production JSON](../../../data/vector_sector/odd_preliminary_fit.json)
contains only `metadata.bibtex_references` and `parameters` for
`central_single_trace_negative`. The parameter map includes the even-sector
parameters, $c_i,d_i$, the $\omega$--$\phi$ mixing angle, fitted total widths,
and $G_V$. The standalone [`G_V` output](../../../data/vector_sector/gv_conditional_fit.json)
has the same minimal shape. The $G_V$ value remains conditional on the
provisional point. No parameter correlation matrix was calculated, so no
correlation field is stored.

The generator maps the selected archived amplitude vector into the original
eleven-operator image and checks its reconstruction numerically. It performs
no additional provisional minimization. The largest relative
amplitude-reconstruction defect is $2.38\times10^{-15}$ and the largest
coefficient-null residual norm is $1.55\times10^{-14}$. The minimal fit
artifacts were regenerated on 2026-10-22; the stored provisional parameter
values were reused and $G_V$ was recomputed. These checks do not independently
repeat the nonlinear provisional fit.

At the chosen fixed setup, the achieved total $\chi^2$ values are
$362.14547$, $346.29371$, $5.5675233$ and $54.216008$ for the universal,
positive-strange single-trace, negative-strange single-trace and double-trace
source hypotheses, respectively. The archived double-trace number imposed an
additional requirement that each signed vector-current residue have the same
sign as its universal-source value. That requirement does not follow from the
double-trace Lagrangian and incorrectly removes the $y_V=0$ embedding of the
negative-strange single-trace fit. The corrected double-trace fit imposes no
such residue-sign restriction; starting at that exact embedding gives
$\chi^2=5.5675233$, and the local refit reaches $5.5411233$ with
$y_V=-2.19\times10^{-4}\,\mathrm{GeV}^{-1}$. The corrected fit and its
reproduction script are in
[the seeded double-trace result](../../../data/vector_sector/odd_double_trace_seeded_refit.json)
and [its script](../../../scripts/refit_double_trace_from_single_trace.py).
The older $54.216008$ value is not comparable to this fit because its
parameter domain was artificially restricted. The negative-strange case is an adequate
conditional working example with six nominal residual degrees of freedom;
the nominal count assumes the inherited diagonal-summary likelihood.
Electronic widths alone do not determine that relative source sign.

For this working example, replacing only the rho mass input by its neutral
summary gives $\chi^2=5.6466225$; replacing only the K-star mass input gives
$12.158901$; replacing only the phi electronic evaluation by its average
gives $5.6276328$. The K-star replacement is a material sensitivity to the
mass prescription in the exact-isospin even sector. It is not an extra
independent dataset to be averaged with the charged input. Removing broad-rho
rows gives $0.37741791$; removing constrained-fit K-star rows gives
$3.5374270$. These totals have different observation counts and are not
likelihood-ratio model comparisons.

Seven measured odd combinations do not determine eleven individual
coefficients. Any displayed $c_i,d_i$ representative must state its
four-dimensional unmeasured coefficient space. A chosen zero-null
representative has no measured uncertainty in the omitted directions;
individual coefficient errors must not be manufactured from that choice.
The source hypothesis, amplitude-sign solution and all fixed-input
assumptions remain attached to every preliminary number.

The immediate fit task is to inspect the compact focused artifact and use
the adequate original-operator case as a conditional working result.
Additional optimization is warranted only by a specific unresolved
observable or material input sensitivity. Completing normalized ALP
predictors and propagating neutral-input uncertainty remain separate
scientific tasks; they do not require repeating the 1792-case fit grid.
