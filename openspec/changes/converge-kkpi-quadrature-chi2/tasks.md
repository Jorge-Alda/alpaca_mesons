# Tasks

## 1. Enable bounded-memory higher-order evaluation

- [ ] 1.1 Implement deterministic energy-bin block evaluation for pair/angular orders through 128, and verify predictions, residuals, and Jacobians match the full-grid path at supported low and moderate orders for both channels.
- [ ] 1.2 Add focused tests for block boundaries, the physical Dalitz limits, the ϕ and both K* pole maps, and deterministic reductions; verify the quadrature certification and block-equivalence suites pass.

## 2. Retry and record the certification

- [ ] 2.1 Run all five saved points with energy orders 2,4 and pair/angular orders 16,32,64,128; store every independent and combined comparison in a machine-readable result and verify no point or failed comparison is omitted.
- [ ] 2.2 Replay the recorded result and verify all point vectors, comparisons, outcomes, and input/source hashes reproduce; add a memory record and index entry, and update `alp-uzh.2` only if every point passes the unchanged tolerances.
