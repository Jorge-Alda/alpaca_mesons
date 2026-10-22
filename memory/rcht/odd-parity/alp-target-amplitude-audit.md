---
title: Inventory of ALP odd-parity amplitudes and relevance of SM blind directions
status: mixed
created: 2026-10-13
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg, Bai:2025fvl]
provenance: analytical-derivation
---

# ALP target audit before the next fit stage

## Scope and meaning of the current fit

The current repository has [stage-4 sensitivity/calibration infrastructure](vector-stage4-rank.md),
not a completed physical fit of the eleven odd coefficients. Its structural
null directions, in the order $(c_1,\ldots,c_7,d_1,\ldots,d_4)$, are
$c_7$, $d_4$ and $\delta c_2=\delta c_5=\delta c_6$.
Uncertainties in the other directions, even-sector nuisance parameters and
discrete residue signs have not been established by that rank calculation.

Use Alda:2026nrd's archived arXiv v1 for ALP LSZ/source conventions, as indexed in
[the neutral input record](../../extracts/Alda:2026nrd/neutral-meson-generators.md).
Its invariant definition and no-mixing formulas were checked in the downloaded
source at `eq:invariants`, `eq:mixing_coeff`, `eq:ALPcouplings_nomixing` and
`eq:mixing_validity`. The old source/operator derivations agree with these
definitions. No local manuscript is substituted for the arXiv source.

## What existed before this audit

| Target | Existing calculation | Missing completion |
|---|---|---|
| $a\to3P$ | LO chiral and even-parity $G_V$ exchange plus the representation-induced local contact | No old explicit odd-sector contribution; fresh tree audit below establishes its zero in the retained vacuum action |
| $a\to VV$ | Physical-ALP VVP operator insertions and a two-background-source anomaly construction | Complete normalized physical amplitude/width in the calibrated vector action |
| $a\to V\gamma$ | Direct VJP source insertion and the direct $c_7$ on-shell zero | Full direct-plus-VVP/photon-response amplitude, anomaly/induced-term bookkeeping and width |
| $a\to\gamma\gamma$ | Derived no-vector WZW/contact invariant in Bauer's photon-operator normalization and the retained leading-nonet WZW approximation, plus unreduced resonance source insertions | Normalized resonance-corrected ALP amplitude/width with both photon assignments |

Evidence is the [even three-P calculation](../../alp/three-pseudoscalar/a3p-vector-exchange-amplitude.md),
[two-source ALP derivation](two-source-alp-rcht-odd-parity-extension.md),
[direct VJP derivation](odd-parity-coupling-prescriptions-and-direct-alp-vgamma.md),
[operator-level numerical kernel](odd-parity-rcht-numerical-kernel.md), and
[no-vector photon invariant](../../alp/photon/alp-photon-invariant-from-mixing.md).
SM amplitudes in stage 4 are not automatically physical ALP amplitudes:
the ALP derivative and pseudoscalar mass spurions are independent.
The generic anomaly bilinear congruence is not a completed physical predictor.

## Fresh blindness checks

The [two-boson tensor audit](alp-two-boson-blind-direction-audit.md) and
[real-photon ALP audit](../../alp/photon/alp-photon-rcht-blind-direction-audit.md)
verify that the projected VJP/VVP sensitivities annihilate all three structural
blind directions, without imposing OPE identities or identifying the ALP
derivative and mass spurions. Signed spectral/source weights do not resurrect
an identically zero contraction. The [three-P audit](../../alp/three-pseudoscalar/odd-parity-three-pseudoscalar-audit.md)
establishes zero single-odd tree sensitivity for the four-scalar vacuum amplitude.

For the two-boson result let $N$ contain the three blind columns and $r$ be a
projected coefficient row. The symbolic audits find $rN=0$. Because the current
SM sensitivity matrix has rank 8 with this three-dimensional kernel, its row
space contains every such row. This is a linear-algebra coverage statement,
not an uncertainty or normalization calculation.

The independent [span check](../../../scripts/audit_alp_blind_direction_span.py)
uses arbitrary independent derivative/mass weights and the stored stage-4 rows.
Across all 64 scenarios, the SM observable rows alone have rank 8; the sampled
target-row residual from that span is below $2.3\times10^{-15}$, and below
$6.5\times10^{-15}$ when inherited OPE rows are also included.
The [results JSON](../../../data/vector_sector/alp_blind_direction_audit.json)
pins scripts/input hashes. This blindness conclusion for the retained projected
vertices does not depend on accepting the old OPE targets in the enlarged action.

## Decision and limits

The three structural blind directions do not enter the checked leading
on-shell two-boson resonance projections; there is no new unconstrained
combination of the retained eleven coefficients in those projections.
The retained tree odd sector gives no $a\to3P$ contribution. These are newly
verified limited results, not claims that all physical ALP predictors already
exist or that every higher-order contribution is blind.

Full assembly still needs correct normalization, contact/LSZ bookkeeping,
Ward identities and reparametrization invariance. Generic off-shell tensors,
two-odd-insertion loops, additional singlet/OZI operators and independently
matched ALP local terms are outside the checked scope. UV ALP coefficients
are external inputs, not quantities constrained by SM meson calibration.
Unstable-vector daughter amplitudes must be audited before extending a narrow
on-shell conclusion to inclusive spectra.

Pause advancement of the previous vector/joint-fit plan while the ALP-target
completion is assessed. The missing calculations have indexed plans:
[aVV](alp-vv-completion-plan.md), [aVgamma](alp-vgamma-completion-plan.md),
and [a gamma gamma](../../alp/photon/alp-photon-rcht-completion-plan.md).
For three P the tree zero is sufficient at this truncation; a conditional
loop-extension plan is contained in its audit record. No new fit is run here.
