---
title: Size prescription and cancellation limits for the loop estimate with the radial mode
status: mixed
created: 2026-11-23
updated: 2026-11-23
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: assumption
---

# Size prescription and cancellation limits with the radial mode


This record is task 3.2 of the OpenSpec change `phi-dependent-loop-tree-crossing-murayama`. It fixes the weight given to each loop contribution of the estimate, and the limits within which the estimate may be used, when the radial mode $\sigma$ is included. The classes are those of [the radial counting record](../diagram-traces/murayama-radial-loop-class-counting.md).

**No loop moment is evaluated.** The prescription assigns sizes to scalar moments and to their coefficients. It does not compute any scalar loop function, coefficient, or loop value. The weights below are estimates, not results.

## Inputs

- Scalar-moment reduction, size per moment and simplex weights: [momentum-reduction record](../../../../loop-corrections/finite-mixing/three-meson-amplitudes/alp-loop-momentum-reduction.md).
- Class-dependent mass scales and the naturalness assumption: [naturalness record](../../../../loop-corrections/finite-mixing/power-counting-and-validity/loop-counting/alp-loop-naturalness-assumption.md).
- Class numbering, ALP-only cancellation conditions and their symbolic checks: [cancellation record](../../../../loop-corrections/finite-mixing/power-counting-and-validity/crossing-and-validity/alp-loop-estimate-cancellations.md).
- Classes and powers with the radial mode: [radial counting record](../diagram-traces/murayama-radial-loop-class-counting.md).

## Weight of each contribution

1. **Reduce each loop to scalar moments.** Shift the loop momentum and reduce the numerator with the rule of the momentum-reduction record. Each scalar moment with $r$ remaining denominators carries a size $(m^2)^{2-r}/(16\pi^2)$ at the reference scale $m^2$ of that loop, with its polynomial coefficient kept.
2. **Apply the derivative rule before sizing.** The external-momentum derivative needed for a residue, and any double propagator from the mixing expansion, are differentiated first, as the momentum-reduction record requires.
3. **Triangles.** Parameter integrals over the triangle use the simplex weight $2\,i!\,j!\,k!/(i+j+k+2)!$ for each monomial $z_1^iz_2^j(1-z_1-z_2)^k$.
4. **Combine before taking magnitudes.** Coefficients of equal scalar moments, at the same masses and invariants, are combined first, including their weights from the one-particle-irreducible, exchange, vacuum and pole terms. The magnitudes of the distinct moments are then added.
5. **Reference scale.** The scale $m$ is that of the class in the naturalness record for ALP-like classes. For the radial classes:

| class | reference scale $m$ |
|---|---|
| closed radial line, one-point tadpole (class 8) | $m_\sigma$ |
| radial bubble or triangle (class 4 and class 6) | $s$ or $m_\sigma$, whichever is the scale of the vertex, with the choice stated for each term |
| radial line with an ALP line (class 5), or two radial lines (class 6) | $s$ |

The naturalness assumption is applied with these scales. They are labelled choices. In the up-only model $m_\sigma$ and $m_a$ are close (only-up radial record), so the two scales coincide in practice.

## Limits of the prescription

- **Relative signs and phases are not given.** The prescription yields magnitudes only. Two distinct moments can cancel, and the estimate does not detect this.
- **Valid when the loop scales are comparable.** The prescription assumes no enhancement from a threshold, an infrared singularity, or a large logarithm. Where a threshold lies in the physical region, the estimate is not applicable and this must be stated for the channel.
- **Radial exchange in the tree is not a loop.** The $\sigma$-exchange tree term is of the same relative order $\epsilon^2$ as the ALP-only classes 2, 3, 5, 6 and the radial diagrams with one radial line (counting record). It is kept exact in the denominator of the crossing. Its relative sign with respect to the loop estimate is not fixed by the prescription.

## Cancellation limits with the radial classes

The ALP-only cancellation conditions still apply to the ALP-only classes. With the radial mode they are extended as follows.

1. **Equal-power sum.** The ALP-only classes 2, 3, 5 and 6, and the radial diagrams with one radial line (classes 4 and 6 of the radial estimate), share the leading power $\epsilon^2$. The leading cancellation condition of the cancellation record, that the signed mass-weighted sum at that power vanishes, applies to their sum. The $\sigma$ classes enter this sum with their own scale and contraction.
2. **Zeros of the contractions.** The ALP-only zeros of the cancellation record are conditions on the ALP invariants. They do not apply to $\sigma$ contractions, which are built from $\widehat S_{10}$, $\widehat S_{20}$, $\widehat P_{11}$ and the $\sigma PP$ couplings. Each $\sigma$ zero must be stated separately, and none is assumed here.
3. **Limit of use.** Where the leading sum vanishes at fixed kinematics, the crossing estimate based on its uncancelled size does not apply, as in the ALP-only case. The estimate cannot show that such a cancellation happens: that needs the loop coefficients, which this change does not compute.

## Status

- Weights, scales and cancellation limits are stated. No moment, coefficient or loop value is evaluated.
- The relative signs of the classes are not determined. Any use of the estimate must state this limit.
