# Proposal

## Why

The global BaBar $KK\pi$ vector-sector fit currently reaches its 1500-evaluation limit and takes close to an hour, which makes iteration and reproducibility checks costly. Profiling and optimizing the dominant numerical work should reduce wall time while preserving the likelihood, fit diagnostics, and numerical results.

## What Changes

- Profile the default fit under documented, reproducible conditions and identify the dominant runtime costs.
- Optimize the measured bottleneck while preserving all residual blocks, parameter definitions, convergence settings, and reported diagnostics unless a change is shown to be numerically equivalent.
- Record comparable before/after runtime and verify optimized predictions, residuals, objective values, and fit outputs against the current implementation within stated tolerances.

## Capabilities

### New Capabilities

None.

### Modified Capabilities

None. This is an internal performance and reproducibility change with no spec-level behavior change; `.openspec.yaml` declares `skip_specs: true`.

## Impact

- `scripts/fit_babar_kkpi_global_vector_sector.py` and, if profiling identifies a shared numerical bottleneck, its directly used implementation modules.
- Fit output and profiling/reproducibility records under the existing vector-sector data and memory conventions.
- No intended changes to the likelihood, input data, fitted parameterization, or interpretation of the current nonconverged baseline.
