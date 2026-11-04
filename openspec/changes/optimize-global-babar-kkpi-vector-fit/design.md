# Design

## Context

See `proposal.md` for motivation. The current script evaluates 75 residual entries for 32 coordinates, constructs two-channel Dalitz quadrature grids, caches propagator/current terms, computes selected analytic derivatives and finite-differences the remaining columns, and runs SciPy `least_squares`. The stored baseline reaches `max_nfev=1500` without convergence; the existing memory record documents its input snapshot and fit state.

## Goals / Non-Goals

**Goals:**

- Measure the default fit's runtime by stage and identify the actual dominant cost before choosing an optimization.
- Reduce wall time per comparable run without altering the mathematical likelihood or numerical precision needed to reproduce fit diagnostics.
- Leave a reproducible comparison of environment, command, inputs, wall time, objective, parameter values, predictions, and convergence diagnostics.

**Non-Goals:**

- Change the data, residual definitions, quadrature prescription, parameterization, optimizer method, stopping criteria, or physical assumptions as a route to apparent speedup.
- Claim global convergence or revise the interpretation of the baseline fit.
- Add a runtime dependency unless profiling establishes a compelling need and it is separately justified.

## Decisions

- **Profile before selecting a target.** Use a deterministic, representative run with the existing default arguments and capture both whole-run timing and function-level cost. This avoids optimizing visually complex code that may not dominate total runtime. Keep profiler overhead separate from the comparable wall-time measurements.
- **Preserve a baseline reference.** Capture the current revision's output and measured observables before modifying the numerical path. Record the command, input hashes, Python/library versions, hardware context, elapsed time, `nfev`/`njev`, success/message, chi-squared, parameter vector, and predictions. Repeat timed runs when practical and report the spread.
- **Make the smallest bottleneck-focused change.** Reuse existing cache and bounded-block pathways where appropriate. Any derivative, array layout, caching, or loop changes must retain the same objective and output fields; avoid reducing quadrature orders or convergence requirements to improve timing.
- **Compare numerical behavior at fixed points and end-to-end.** Compare residuals, per-bin predictions, Jacobian columns where affected, chi-squared, and optimizer diagnostics against the baseline using explicit tolerances justified by existing quadrature certification or floating-point equivalence. Then run both versions with matching settings and initial conditions. If optimization changes the path taken by the optimizer, document that and separately compare predictions/objective at the resulting candidates; do not describe different minima as identical.
- **Keep profiling artifacts auditable.** Store concise measurements and comparison results with the existing fit artifacts or a focused memory record, including exact invocation and environment metadata; do not commit bulky profiler traces unless needed to reproduce a finding.

## Risks / Trade-offs

- **Timing noise can hide or exaggerate gains** → use the same machine, environment, inputs, command, warm-up policy, and repeated measurements; report median and range.
- **A faster approximation may bias the likelihood** → retain quadrature settings and require prediction/residual/Jacobian/objective comparisons before accepting the change.
- **Caching may return stale values as parameters change** → include all dependencies in cache keys and compare cached and uncached evaluations on fixed candidates.
- **Optimizer trajectories are sensitive to small numerical changes** → separate fixed-point equivalence from end-to-end convergence comparison and preserve complete diagnostics.
- **Profiling instrumentation adds overhead** → use it only to locate costs; do not use profiled runtime as the before/after performance claim.

