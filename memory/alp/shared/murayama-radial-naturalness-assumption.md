---
title: Naturalness assumption for the loop estimate with the radial mode
status: conjectured
created: 2026-11-23
updated: 2026-11-23
bibtex_keys: [Alda:2026nrd, Manohar:1983md, DiLuzio:2026poh]
provenance: assumption
---

# Naturalness assumption for the loop estimate with the radial mode


This record is task 3.3 of the OpenSpec change `phi-dependent-loop-tree-crossing-murayama`. It extends the [ALP-only naturalness assumption](alp-loop-naturalness-assumption.md) to the classes that contain the radial mode $\sigma$. It is an assumption. It is not a result, and it is not checked against any loop calculation.

## Assumption

The neglected next-order correction to an on-shell ALP-sourced amplitude is taken to have the magnitude of the loop-size estimate of the [radial counting record](murayama-radial-loop-class-counting.md), with the weights and reference scales of the [size prescription record](murayama-radial-loop-size-prescription.md). The assumption covers every part of that correction, including:

- the part that depends on the renormalization scale $\mu$;
- the counterterm part, that is, the local contributions whose coefficients are not fixed by the notes.

The estimate contains no $\mu$. The assumption does not analyse what the counterterms or the $L_i$ absorb. That question depends on the operator basis, which the notes do not fix, and it is not needed for the estimate. This follows the choice made in design D3 of the ALP-only change ([naturalness record](alp-loop-naturalness-assumption.md)).

## Conditions under which the assumption is taken

- Large logarithms, threshold enhancements and infrared singularities do not control the net correction.
- The radial scale $m_\sigma$ and the ALP-like scale $m_a$ enter the loop estimate with the scales of the size prescription.

## Limits: cancellations for which the estimate overestimates the correction

The estimate overestimates the correction when the following cancellations occur. Each is a limit of the assumption, and none is excluded by it.

1. **Equal-power sum.** The ALP-only classes 2, 3, 5 and 6 and the radial diagrams with one radial line (classes 4 and 6 of the radial estimate, see [the estimate record](murayama-radial-up-only-estimate-result.md)) are all of relative order $\epsilon^2$. If their signed, mass-weighted sum at that power vanishes at fixed kinematics, the estimate based on the uncancelled size does not apply. This is the ALP-only condition extended to the $\sigma$ classes ([cancellation record](alp-loop-estimate-cancellations.md)).
2. **Zeros of the $\sigma$ contractions.** The ALP-only zeros of the invariant contractions do not transfer to contractions built from $\sigma$ couplings. Any such zero is a cancellation that the estimate does not capture.
3. **Relative sign with the $\sigma$-exchange tree.** The tree $\sigma$ exchange is of the same relative order $\epsilon^2$ and is kept exact in the crossing denominator. If the loop estimate and this exchange partly cancel, the estimate of the ratio overstates the loop-to-tree size. The prescription does not fix the relative sign.
4. **Non-generic cases.** Where the vertex scale is set by a kinematic invariant that vanishes at the fixed Dalitz point, the scale of the class is not the one assumed, and the estimate does not apply there.

## Status

- The assumption is stated for all classes of the estimate, including the radial ones, and for the $\mu$-dependent and counterterm parts.
- Its validity is not checked. The cancellation limits are listed, not evaluated.
- No value of $\mu$, no counterterm coefficient and no loop value appears in this record.
