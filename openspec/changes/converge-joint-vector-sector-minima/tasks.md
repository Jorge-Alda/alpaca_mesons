# Tasks

## 1. Define deterministic starts and scaling

- [x] 1.1 Extract and record the existing coordinate bounds, branch checks, and scale factors; define named deterministic starts inside the permitted domain and verify every start passes the likelihood's physical-domain checks.
- [x] 1.2 Add the multi-start attempt ledger and CLI/configuration for optimizer settings; verify serialization retains each start, scaling, bounds, quadrature, inputs, and complete termination fields, including unsuccessful attempts.
- [x] 1.3 Add focused tests for start validation, reproducibility, bound/branch preservation, and distinct endpoint retention; verify the tests pass without changing existing single-start outputs.

## 2. Run and qualify candidate minima

- [x] 2.1 Run all named starts against the unchanged raw likelihood using the existing propagator reuse, channel-specific Jacobian path, and bounded-bin evaluation; verify fixed-candidate predictions, residuals, and Jacobians against the established path within recorded tolerances.
- [x] 2.2 Preserve per-attempt residual and block diagnostics, objective, active bounds, optimizer status/message, evaluation counts, optimality/gradient, wall time, software versions, and input/source hashes; verify each ledger entry is independently interpretable.
- [x] 2.3 Re-evaluate successful promising candidates at $(2,64,64)$ and $(4,128,128)$ with four-bin blocks; verify prediction, standardized residual, scaled-Jacobian, and $\chi^2$ changes against the certified tolerances and reject candidates that fail. (No attempts succeeded; each was explicitly recorded as not eligible for refinement.)
- [x] 2.4 Write the separate machine-readable results and indexed memory record documenting inputs, assumptions, exact invocation, outcomes, and limitations; verify index links and replay the recorded invocation.

## 3. Continue the best endpoint and assess convergence

- [x] 3.1 Compare the four recorded endpoints using their scaled coordinates and objectives, select the physically valid lowest-objective endpoint as the continuation seed, and retain its unsuccessful status; verify the chosen vector matches the source ledger exactly.
- [x] 3.2 Add a `--resume-from` path that loads a named endpoint from a prior ledger, revalidates its bounds and physical branch, and records progress checkpoints containing coordinates, objective, evaluation count, and available optimality; verify restart-vector preservation, checkpoint serialization, and source-ledger immutability with focused tests.
- [x] 3.3 Continue from `stored_negative_c4` at the default quadrature with `max_nfev=5000`, the existing bounds, and the original fixed scale vector; write a separate result ledger and verify it records the seed provenance, per-run and cumulative evaluation counts, full optimizer diagnostics, and endpoint validity.
- [x] 3.4 If the continuation succeeds, refit sequentially from its endpoint at $(2,64,64)$ and $(4,128,128)$ using four-bin blocks, the same bounds/scales, and a 5000-evaluation cap; verify successful termination, physical validity, and stationarity information at each quadrature level, retaining each fit separately.
- [x] 3.5 If the certified-order continuation succeeds, define and record three deterministic starts that perturb separate parameter blocks from that endpoint using the saved scale vector; validate each start before fitting at default quadrature, retain all outcomes, and refine any successful competitive endpoint with the certified quadrature checks.
- [x] 3.6 If continuation or a required higher-order refit terminates unsuccessfully, analyze saved progress, optimality, coordinate movement, active bounds, and residual blocks; record that no minimum was accepted without claiming that no minimum exists, and verify the diagnosis is supported by the attempt ledger. (The continuation and both required higher-order refits succeeded, so this failure branch was not triggered.)
- [x] 3.7 Compare all successful and competing endpoints with declared scaled-coordinate and objective tolerances, accepting a minimum only if optimizer, domain, residual, and quadrature gates pass; confirm historical fit artifacts are unchanged and run the focused fit, Jacobian, bounded-block, and quadrature-certification checks, recording commands and outcomes.
