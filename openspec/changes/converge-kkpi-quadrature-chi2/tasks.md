# Tasks

## 1. Enable bounded-memory higher-order evaluation

- [x] 1.1 Add a high-order baseline option and deterministic energy-bin block evaluation through pair/angular order 128; verify predictions, residuals, and Jacobians match the full-grid path at supported low and moderate orders for both channels.
- [x] 1.2 Add focused tests for block boundaries, the physical Dalitz limits, the ϕ and both K* pole maps, and deterministic reductions; verify the quadrature certification and block-equivalence suites pass.

## 2. Retry and record the certification

- [x] 2.1 Run all five saved points with energy orders 2,4, pair orders 64,128, and angular orders 64,128; store every independent and combined high-order comparison separately from the historical 16→32 result, and verify no point or failed comparison is omitted.
- [x] 2.2 Replay the recorded result and verify all point vectors, comparisons, outcomes, and input/source hashes reproduce; add a memory record and index entry, and update `alp-uzh.2` only if every point passes the unchanged tolerances.
