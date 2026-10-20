---
title: Coupling sensitivity after stage-4 vector amplitude integration
status: numerical
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# Stage-4 conditional coupling ranks

## Inputs and approximations

Use the [integrated action kernels](vector-stage4-amplitudes.md),
[arXiv-derived neutral inputs](alda-neutral-amplitude-inputs.md),
[stage-3 calibration](../foundations/vector-stage3-calibration.md), and
its [PDG snapshot](../../extracts/ParticleDataGroup:2024cfk/vector-sector-calibration-inputs.md).
The [results JSON](../../../data/vector_sector/stage4_rank_results.json)
pins the input/kernel/script hashes and every retained row and singular value.

Assume the declared tree action and fixed central source/mass inputs. Retain
all four separate rho/K-star choices, both phi electronic-width summaries,
both mass branches and four source hypotheses: universal, positive single-trace,
opposite-sign single-trace, and single-plus-double-trace interpolation.
This gives 64 separate scenarios, including inadequate calibration hypotheses
as sensitivity controls. They are not 64 independent datasets or equally valid
physical models. The user-selected near-ideal preference applies only when the
other criteria are equivalent; it is not included as a likelihood or rank row.

Calibrated external vector masses follow the same exact-isospin mass action:
charged rho uses the selected rho mass and both K-star charges use the selected
K-star mass. No independent charged/neutral vector mass insertion is made into
these action scenarios. The old baseline controls retain their previous physical
kinematic inputs and are explicitly labelled as controls.

The even-vector block remains isospin projected. The odd-vertex quark spurion
and physical neutral generators retain Alda:2026nrd's fitted isospin breaking. Rho–isoscalar
mixing is omitted, so this is a declared partial-isospin approximation. Physical
pseudoscalar kinematics also remain distinct from the LO eigenmasses, especially
for eta. See the input record for the convention conversion and limitations.

**Short-distance limitation:** the five rows used here are inherited homogeneous
sensitivities of the old OPE relations, not a rematched OPE prescription for the
enlarged mass/source action. The ranks below describe this specified matrix;
they do not justify imposing the old affine targets as a precision prior.

## Method and reproduction

Run:

```text
.venv/bin/python scripts/diagnose_vector_stage4.py
```

The [script](../../../scripts/diagnose_vector_stage4.py) evaluates each retained
coupling basis vector and subtracts the zero-coupling amplitude. Since the
amplitudes are affine in the eleven odd coefficients, these differences give
the exact linear sensitivity rows, including removal of the WZW offset.
The ten radiative channels follow the existing list; the neutral TFF rows cover
pi0, eta and eta-prime at $(s_1,s_2)=(0,0),(-1,0),(-5,0)$ GeV squared.
$F=0.093$ GeV and the independent operator normalization is $M_V=0.770$ GeV.

Normalize every nonzero row to unit Euclidean norm, then use singular-value
decomposition with absolute tolerance $10^{-11}$. No measurement errors,
likelihood weights or source-parameter marginalization enter this calculation.

## Results

For every retained calibrated scenario:

| Matrix | Rank of 11 |
|---|---:|
| inherited OPE rows | 5 |
| inherited OPE plus ten radiative rows | 8 |
| above plus three real-photon TFF rows | 8 |
| above plus six single-virtual TFF rows | 8 |

The exact persistent null directions are $c_7$, $d_4$ and the common shift
$\delta c_2=\delta c_5=\delta c_6$, with all other entries unchanged.
The executable verifies their annihilation for every individual amplitude row.
Neutral TFFs continue to be insensitive to $c_4$; the charged-kaon radiative
channel provides that direction in the combined matrix.

The eighth normalized singular value spans about 0.00288–0.04879 for radiative
rows, 0.04228–0.06118 after real-photon rows, and 0.05677–0.07135 after virtual
rows. These row-normalized diagnostics carry no experimental information scale.
The extra TFF rows can alter redundancy/conditioning without adding a new
direction; unchanged rank does not mean unchanged fit uncertainty.

Two controls pass. The old illustrative common/split-pole inputs reproduce
radiative ranks 7/8 and combined TFF rank 8. With the new full Alda:2026nrd generators
and independent quark spurion but a universal common-pole sector, radiative
rank remains 7 and rises to 8 with TFFs. Thus merely updating the eta inputs
does not explain the calibrated radiative rank 8; the spectral/source response
changes its relative flavor weights.

## Interpretation

Stage 4 establishes amplitude implementation and a reproducible conditional
sensitivity result for the declared retained action. It does not establish
short-distance matching, statistical usefulness, marginal identifiability
when the even-vector coefficients float, a physical choice of relative signs,
or an odd-parity fit. Further odd symmetry-breaking operators could change
the parameter space. Both mass branches remain recorded even where a later
adequacy/truncation assessment favors the near-ideal one.
