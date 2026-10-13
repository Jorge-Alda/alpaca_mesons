---
title: Stage-5 joint decay likelihood specification and readiness
status: mixed
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, CLEO:1997fho, BaBar:2009rrj, Niecknig:2012sj]
provenance: assumption
---

# Joint calibration specification

## Subsequent selected calibration scope

On 2026-10-13 the user selected an
[odd flavor-breaking single/double-trace extension plus low-energy calibration](odd-breaking-low-energy-calibration-plan.md).
That plan supersedes the all-channel dataset as the next fit objective.
The retained hard-TFF matching and dispersive three-pion gates below apply
to a later expanded fit, not to the selected low-energy dataset. Its first
fit has no inherited OPE prior and requires neither a demonstrated ultraviolet
completion nor an unimplemented three-pion predictor. Adequacy, normalized
amplitudes and identifiability remain mandatory for that narrower scope.

## Decision and scope

User instruction on 2026-10-13: retain tree-level ALP accuracy and defer
momentum-dependent self-energies. No loop subtraction scheme is selected;
no dynamical mixing parameter is fitted. This closes the precision-target
decision without claiming loop effects are numerically negligible.

The user now requests stage 5; the earlier ALP-audit task-order hold is
superseded. Missing normalized ALP predictors remain separate work.

## Inputs and truncation

Use [stage-4 amplitudes](vector-stage4-amplitudes.md),
[Alda arXiv inputs](alda-neutral-amplitude-inputs.md),
[mass/source calibration](../foundations/vector-stage3-calibration.md),
[validated likelihood data](odd-parity-data-implementation.md), and
[fit infrastructure](correlated-odd-parity-fit-infrastructure.md).
All cited sources already have records in `bibliography.bib`; no new source
is introduced here. Natural units and the existing positive-strange phi phase
apply; masses are in GeV. Preserve $F=0.093$ GeV and the independent
operator scale $M_V=0.770$ GeV in the central-input diagnostics.

Retain the declared extended mass action and the universal, single-trace
positive, single-trace negative, and single-plus-double-trace source hypotheses.
Retain both mass branches, four charged/neutral rho/K-star selections, and
both phi electronic-width summaries separately. The near-ideal preference
is only a tie-breaker after adequacy/truncation assessment.

Odd vertices remain the existing single-trace one-pseudoscalar action.
Additional odd flavor-breaking and multi-pseudoscalar contact operators,
kinetic breaking, rho–isoscalar mixing and loops are omitted. The even block
is isospin projected while odd spurions and neutral pseudoscalar generators
retain Alda isospin breaking. This partial-isospin truncation is an assumption,
not a complete symmetry-breaking action. Do not reinterpret the three-pion
exchange contribution as a complete amplitude.

## Likelihood and source treatment

The intended joint residual includes mass, leptonic, radiative, TFF and
validated dispersive three-pion blocks. Each block uses a fixed covariance
and Cholesky whitening. Model alternatives and omitted-operator uncertainty
remain separate from experimental variances.

In the eventual joint fit, float the mass/source coefficients of the selected
action alongside observable odd combinations. Evaluate masses, generators
and signed residues from those same coefficients on every call. Electronic
widths constrain squared residues; retain sign hypotheses explicitly.
Stage-3 estimates are starting values, not independent Gaussian priors:
reusing them as priors while retaining their input widths would double-count
data. Fixed central sectors are allowed only for labelled conditional diagnostics.
Do not infer a physical mixing confidence interval from such a diagnostic.

Covariance assumptions for the held data are:

| Block | Specification and limitation |
|---|---|
| Mass/electronic summaries | Stage-3 diagonal Gaussian summary errors; alternative summaries are separate scenarios, never simultaneous measurements. Unavailable PDG cross-covariances preclude a precision significance claim. |
| Omega radiative branching fraction | Held PDG value $0.0835\pm0.0027$; predict the branching fraction directly. A supplied total width is conditional until its uncertainty/correlation is represented. Do not manufacture an independent partial-width datum from it. |
| Other radiative decays | Ten amplitude kernels exist, but a provenance-checked measurement/covariance snapshot is not yet connected. They are diagnostic rows only. |
| CLEO TFFs | Existing per-bin statistical plus systematic diagonal approximation, with zero inter-bin/inter-meson covariance explicitly assumed. No invented source-wide nuisance. This policy does not establish high-energy model adequacy. |
| BaBar pion TFF | Existing mean asymmetric error conversion and separate standardized 2.3% multiplicative normalization nuisance. No extra cross-section normalization nuisance for this observable. |
| Three pions | No accepted likelihood yet: raw spectrum/Dalitz predictions and their dispersive matching are absent. The held KLOE amplitude decomposition is model-dependent, not independent raw Dalitz bins. |

For a future assembled likelihood, assume zero cross-covariance between
different experiments and mass/leptonic/radiative summary blocks only as an
explicit conditional approximation. Identified shared normalizations must
be represented once, not both in a diagonal variance and a nuisance prior.
PDG summary overlap requires an audit before claiming precision coverage.
Do not enable inherited affine OPE targets: their enlarged-action matching
is absent. No fit or numerical coupling prior is authorized by rank alone.

## Completion gates

1. The [enlarged-action matching calculation](vector-physical-tff-matching.md)
   now establishes a source-breaking obstruction: none of the 48 fixed-central
   source-breaking scenarios cancels all physical TFF hard-limit constants.
   The 16 universal-source scenarios satisfy necessary cancellations but do
   not establish QCD matching or leptonic adequacy. The
   [chiral-octet calculation](vector-chiral-octet-matching.md) also exposes the
   incompatibility between the old correlator OPE and single-virtual falloff.
   Specify and derive a sufficient interaction/resonance extension before
   fitting the held high-energy TFF tables. No arbitrary low-$Q^2$ cutoff or
   inherited affine precision prior bypasses this requirement.
2. Connect the remaining radiative measurements with source-level provenance,
   total-width treatment and covariance bookkeeping.
3. Complete the [scoped dispersive three-pion treatment](dispersive-three-pion-treatment-plan.md):
   establish the multi-pseudoscalar contact action, subtraction-to-RChT matching,
   crossed-channel solution, normalized Dalitz integration and validation.
   The existing Omnès kernel alone does not satisfy these gates. A chi-squared
   difference between published shape fits is not a covariance or a quantitative
   width systematic for an uncomputed predictor.
4. Propagate Alda input uncertainty from documented correlations (or label
   explicitly conditional input variations when covariance is unavailable),
   then minimize separate scenarios and report data/theory components,
   identifiable combinations, source-marginalized rank and model adequacy.

Stage 5's specification and precision decision are established. Full joint
calibration remains incomplete at the concrete gates above; the first matching
attempt now has an explicit failure result for the declared truncation. No physical
odd-coefficient fit or dynamical mixing determination has been produced.
