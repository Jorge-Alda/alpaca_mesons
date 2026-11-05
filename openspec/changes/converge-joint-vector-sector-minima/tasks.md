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

## 3. Assess integrated result

- [ ] 3.1 Compare successful and competing endpoints with declared scaled-coordinate and objective tolerances, retain raw attempts, and mark accepted minima only when every optimizer, domain, residual, and refinement gate passes; verify summary labels match the recorded evidence.
- [ ] 3.2 Confirm the original baseline and optimized artifacts remain unchanged and run the focused fit, Jacobian, bounded-block, and quadrature-certification checks; record commands and outcomes in the result and memory record.
