# Proposal

## Why

The certified joint vector-sector likelihood is reproducible and now evaluates substantially faster, but all four documented starts still terminate at the 1500-evaluation limit. The lowest-objective physically valid endpoint has $\chi^2=226.195$ but unsuccessful optimizer status, so continuation and a clear failure diagnosis are needed before any minimum can be reported as accepted.

## What Changes

- Add a reproducible multi-start fit procedure for the same raw 73-summary likelihood, preserving the established shared nuisances, measurement accounting, and physical branch constraints.
- Record every attempted start, parameter scaling and bounds, quadrature settings, residual diagnostics, competing endpoint, and complete optimizer termination information.
- Continue from the lowest-objective valid endpoint in a separate run with `max_nfev=5000`, preserving its coordinates, likelihood, bounds, scaling, and default quadrature.
- If continuation converges, refit at the certified quadrature orders and then run three additional deterministic starts to probe competing basins. If it still fails, diagnose the optimizer trajectory and report that no minimum was accepted without claiming that no minimum exists.
- Use the established runtime improvements—shared propagator evaluation, cached channel terms, and bounded energy-bin blocks—so repeated starts and high-order refinement remain practical without changing the likelihood.
- Accept a minimum only after successful optimizer termination, physical-domain and residual checks, and stability under the certified quadrature refinement; preserve failures and competing solutions.

## Capabilities

### New Capabilities

- `joint-vector-fit-minima`: Reproducible, physically constrained multi-start minimization and evidence-based reporting of converged joint vector-sector minima.

### Modified Capabilities

None.

## Impact

- Fit driver and focused regression coverage in `scripts/fit_babar_kkpi_global_vector_sector.py` and related scripts.
- Machine-readable per-start results under `data/vector_sector/` and an indexed numerical reproducibility record under `memory/rcht/odd-parity/`.
- Existing likelihood inputs, definitions, parameterization, shared even-sector audit, and certified quadrature tolerances remain authoritative. The existing propagator reuse, channel-specific Jacobian evaluation, and four-bin bounded-memory path should be retained and exercised by the multi-start workflow.
