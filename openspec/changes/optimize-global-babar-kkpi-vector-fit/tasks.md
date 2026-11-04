# Tasks

## 1. Establish the baseline and profile

- [x] 1.1 Run the current default fit with a preserved output snapshot and record exact command, input hashes, Python/SciPy/NumPy versions, machine context, elapsed wall time, `nfev`/`njev`, success/message, chi-squared, parameter values, and predictions; verify the output matches the stored baseline at the same candidate and settings.
- [x] 1.2 Profile representative residual, Jacobian, and full-fit execution without including profiler overhead in timing claims; document cumulative and per-call costs and verify the selected optimization target is among the measured dominant costs.

## 2. Optimize the measured bottleneck

- [x] 2.1 Implement the smallest bottleneck-focused optimization while preserving residual definitions, quadrature orders, parameterization, optimizer settings, output schema, and diagnostics; verify affected fixed-candidate residuals and predictions against the baseline within documented tolerances.
- [x] 2.2 Add or update focused regression coverage for the changed evaluation path, including cache invalidation or derivative behavior when applicable; verify the relevant tests pass and compare Jacobians or other affected quantities against the baseline.

## 3. Verify speed and numerical agreement

- [x] 3.1 Run optimized and baseline implementations with identical inputs, initialization, and fit settings under comparable conditions; report repeated unprofiled elapsed times and median/range, and verify the measured runtime improves.
- [x] 3.2 Compare fixed-point per-bin predictions, residual vectors, Jacobians where affected, objective values, parameter vectors, and full fit diagnostics with explicit tolerances; record any optimizer-trajectory differences without claiming distinct candidates are identical.
- [x] 3.3 Update the focused memory/reproducibility record with profiler findings, implementation summary, before/after measurements, numerical comparisons, exact invocation, and environment metadata; verify links and reproduction instructions are valid.
