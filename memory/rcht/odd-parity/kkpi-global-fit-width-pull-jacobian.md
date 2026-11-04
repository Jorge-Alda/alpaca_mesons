---
title: Width-pull derivatives in the global KKπ fit
status: numerical
created: 2026-11-04
updated: 2026-11-04
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Width-pull derivatives in the global KKπ fit

## Inputs, conventions, and assumptions

The [global fit baseline](babar-kkpi-global-vector-sector-fit-baseline.md) gives the observations, fit coordinates, conventions, and assumptions. The [running-propagator record](kkpi-global-fit-running-propagator-reuse.md) describes the preceding numerical improvement. Invariants are in GeV squared, widths and masses are in GeV, predictions are in nb, and the residual vector and pull coordinates are dimensionless. The stored reference candidate is not a converged fit.

The $\rho^+$ and $\eta'$ total-width pulls affect the low-energy and shared-width residuals, but the BaBar predictions and normalization constraints do not depend on them in the baseline likelihood. The calculation uses that existing dependency; it introduces no new physical assumption or approximation.

## Method

The executable is [`scripts/fit_babar_kkpi_global_vector_sector.py`](../../../scripts/fit_babar_kkpi_global_vector_sector.py). The Jacobian retains the baseline forward-difference steps and proceeds as follows for these two pull coordinates:

1. Evaluate the full baseline residual and the original odd-coordinate derivatives.
2. Shift the pull by its unchanged forward-difference step, unpack the shifted physical quantities, and evaluate the first 28 residuals.
3. Divide the first 28 shifted-minus-baseline differences by the same step. For the remaining entries, subtract each unchanged baseline entry from itself before dividing, retaining the floating-point behavior for either sign of step.
4. If the baseline is non-finite or is the penalty vector, or if the shifted input raises an error handled by the full residual, use the original full-residual evaluation for that column.

At this stage all other columns used the preceding calculation, and the number of full residual evaluations for a regular Jacobian fell from 21 to 19. The later [channel-specific calculation](kkpi-global-fit-channel-jacobian.md) reduces that count further. The dependency must be reconsidered if future versions add either width to the BaBar predictions.

## Verification and timing

The regression script [`scripts/test_global_kkpi_running_propagator.py`](../../../scripts/test_global_kkpi_running_propagator.py) contains the preceding full-residual Jacobian calculation for comparison. At the stored candidate, the complete Jacobian agrees bit-for-bit with that calculation using both the default full grid and a smaller blocked grid. The same suite compares residuals and predictions after displacing all 32 coordinates. Additional checks compare an invalid shifted width, a negative difference step including the signs of zero entries, and a penalty baseline against the full-residual calculation. The penalty check also confirms that all 21 full residual evaluations are retained in that case.

The following focused and block tests passed:

```text
.venv/bin/python -m unittest scripts.test_global_kkpi_running_propagator scripts.test_kkpi_quadrature_blocks
```

Seven alternating unprofiled evaluations of the default-grid Jacobian at the stored candidate gave median times of $0.819$ s before and $0.760$ s after this change, a ratio of $0.928$. This local comparison followed the running-propagator improvement and does not measure the full fit. The machine and software are identified in the [baseline profile](../../../data/vector_sector/babar_kkpi_global_vector_sector_fit_baseline_20261104_profile.json). After a warm-up, the comparison used:

```bash
.venv/bin/python -c 'import json,sys,time,statistics,numpy as np; sys.path.insert(0,"scripts"); import fit_babar_kkpi_global_vector_sector as fit; from scripts.test_global_kkpi_running_propagator import original_jacobian; values=np.array(json.load(open("data/vector_sector/babar_kkpi_global_vector_sector_fit.json"))["parameter_values"]); model=fit.GlobalLikelihood(); model.jacobian(values); original=[]; improved=[]; exec("for index in range(7):\n start=time.perf_counter(); original_jacobian(model,values); original.append(time.perf_counter()-start)\n start=time.perf_counter(); model.jacobian(values); improved.append(time.perf_counter()-start)"); print("previous Jacobian seconds",original); print("improved Jacobian seconds",improved); print("median ratio improved/previous",statistics.median(improved)/statistics.median(original))'
```
