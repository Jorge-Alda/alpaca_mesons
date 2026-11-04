# Design

## Context

See [proposal.md](./proposal.md) for the measured tolerance failures. The current
`phase_points` routine constructs the complete energy-bin, pair, and angular
grid for each channel in memory. The saved certification stops at pair and
angular order 32; any higher-order run must retain all five saved points and
all comparisons.

## Goals / Non-Goals

**Goals:**

- Evaluate the independent and combined qualification comparisons from pair
  and angle order 64 to 128 while keeping peak grid memory bounded.
- Determine pass/failure from the existing per-bin prediction, standardized
  residual, residual-Jacobian, and total-$\chi^2$ comparisons.
- Keep the saved point vectors, physical Dalitz limits, pole maps, derivative
  settings, and input hashes reproducible.

**Non-Goals:**

- Changing the likelihood, experimental errors, or the four declared
  tolerances.
- Refitting parameters or replacing the saved reference and displacement
  points.
- Treating a lower $\chi^2$ as convergence evidence without the successive
  quadrature comparisons.

## Decisions

1. **Evaluate the grid in energy-bin blocks.** Split each channel's energy
   bins into deterministic blocks, build the same pair and angular maps for
   each block, and immediately accumulate the per-bin predictions and
   derivative rows. This bounds peak array size while preserving the physical
   domain, pole locations, weights, and likelihood. A single fully materialized
   high-order grid was considered, but it scales as the product of all three
   orders and all bins and risks exhausting memory before order 128.

2. **Use a distinct high-order qualification baseline.** The previous result
   already records the failed $16\to32$ comparisons. Keep that file unchanged
   as historical diagnostic evidence; do not mix its failures into the new
   qualification status. Compare energy orders $2\to4$ with pair and angle
   held at 64; compare pair orders $64\to128$ with energy 2 and angle 64; and
   compare angular orders $64\to128$ with energy 2 and pair 64. The combined
   comparison is $(2,64,64)\to(4,128,128)$. The new result reports these
   comparisons and applies the unchanged rule that every comparison in this
   high-order run must pass.

3. **Keep the certification contract unchanged.** Reuse the saved five point
   vectors and the existing prediction, residual, Jacobian, and $\chi^2$
   tolerances. Every reported comparison in the new high-order run must pass
   for its status to be `certified`. The earlier $16\to32$ result remains
   separately linked historical evidence. If order 128 still fails, preserve
   the observed differences and report `uncertified`; do not change thresholds
   or discard points.

4. **Verify chunking against the current implementation at supported orders.**
   For both channels, compare predictions and derivatives from block evaluation
   with the existing full-grid calculation at low and moderate orders. Include
   edge bins spanning the $\phi$ pole-map boundaries and both K* angular maps.
   Then run the full five-point certification and replay its saved invocation
   against the recorded input and source hashes.

5. **Keep Beads acceptance tied to the saved evidence.** Leave `alp-uzh.2`
   open if any point or comparison fails. Update it only after a replayable
   result shows that every selected point satisfies the unchanged tolerances.

## Risks / Trade-offs

- [Risk] Block reductions could change floating-point summation order → Use a
  fixed block and reduction order; compare against full-grid results at the
  supported low and moderate orders.
- [Risk] Order 128 may still fail a tolerance → Preserve the failure data and
  keep `alp-uzh.2` open rather than relaxing the contract.
- [Risk] The five full Jacobian evaluations at high orders may be slow → Reuse
  only parameter-independent geometry within each order and point; keep the
  recorded derivative steps fixed.

## Migration Plan

No migration is required. Keep the existing certification file as prior
evidence until the new result is complete, then store the higher-order result
with its own invocation and provenance. If chunked output disagrees with the
current implementation at supported orders, stop the high-order run and retain
the old evidence while correcting the discrepancy.
