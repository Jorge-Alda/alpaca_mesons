---
title: Reuse of the running propagator in the global KKπ fit
status: numerical
created: 2026-11-04
updated: 2026-11-04
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Reuse of the running propagator in the global KKπ fit

## Inputs, conventions, and assumptions

The likelihood, its inputs and assumptions, and the stored nonconverged candidate are described in the [global fit baseline](babar-kkpi-global-vector-sector-fit-baseline.md). The quadrature prescription and its high-order qualification are described in the [high-order quadrature record](kkpi-high-order-quadrature-certification.md). The conventions are those of Alda:2026nrd as specified in the baseline record. Invariant masses are in GeV squared, masses and widths are in GeV, propagators are in inverse GeV squared, and predictions are in nb.

The calculation retains the original running-width prescription, fit parameters, data, numerical quadrature, and arithmetic order. No new physical assumption or approximation is introduced. The saved candidate is a nonconverged snapshot; this change does not establish a global minimum.

## Method

The executable is [`scripts/fit_babar_kkpi_global_vector_sector.py`](../../../scripts/fit_babar_kkpi_global_vector_sector.py). Its running propagator is evaluated as follows:

1. Calculate the reference propagator with both excited-state coefficients set to zero at the requested invariant mass.
2. Invert that value to obtain the imaginary denominator and rescale its running width to the fitted shared on-shell width, as in the original implementation.
3. Evaluate the reference propagator with the requested excited-state coefficients and subtract the value from step 1 to obtain the excited-state contribution.
4. Add that contribution to the rescaled ground-state propagator. The value from step 1 is reused only within this single calculation; a changed invariant, fit parameter, or quadrature grid is evaluated again.

Previously the value in step 1 was calculated twice. This change removes one of four calls to the reference propagator per running-propagator calculation while retaining the original subtraction and addition order. The output remains the same 75-component residual, 45 BaBar predictions, and 75-by-32 Jacobian used by the global fit.

## Verification and timing

The regression script [`scripts/test_global_kkpi_running_propagator.py`](../../../scripts/test_global_kkpi_running_propagator.py) reproduces the preceding running-propagator calculation independently. It compares the full residual vector, both prediction arrays, and every Jacobian entry with the preceding calculation at the stored candidate. It also compares residuals and predictions after displacing each of the 32 coordinates, returning to the candidate, and using both the default full grid and a smaller blocked grid. These comparisons pass bit-for-bit. The default-grid residuals and predictions also agree with the stored baseline at relative tolerance $10^{-12}$ and absolute tolerance $10^{-11}$. A direct check confirms that the revised calculation calls the reference propagator three times and evaluates changed inputs afresh.

The focused regression and existing block tests both passed, six tests in total:

```text
.venv/bin/python -m unittest scripts.test_global_kkpi_running_propagator scripts.test_kkpi_quadrature_blocks
```

On the machine and software environment identified in the [baseline profile](../../../data/vector_sector/babar_kkpi_global_vector_sector_fit_baseline_20261104_profile.json), five alternating unprofiled Jacobian evaluations at the stored candidate gave median times of $1.286$ s for the preceding calculation and $0.937$ s for the revised calculation, a ratio of $0.728$. This is a local timing measurement, not a measurement of full-fit runtime. It was obtained after a warm-up with:

```bash
.venv/bin/python -c 'import json,sys,time,statistics,numpy as np; from unittest.mock import patch; sys.path.insert(0,"scripts"); import fit_babar_kkpi_global_vector_sector as fit; from scripts.test_global_kkpi_running_propagator import original_running_propagator; x=np.array(json.load(open("data/vector_sector/babar_kkpi_global_vector_sector_fit.json"))["parameter_values"]); m=fit.GlobalLikelihood(); m.jacobian(x); original=[]; improved=[]; exec("for index in range(5):\n with patch.object(fit,\"running_propagator\",side_effect=original_running_propagator):\n  start=time.perf_counter(); m.jacobian(x); original.append(time.perf_counter()-start)\n start=time.perf_counter(); m.jacobian(x); improved.append(time.perf_counter()-start)"); print("original seconds",original); print("improved seconds",improved); print("median ratio improved/original",statistics.median(improved)/statistics.median(original))'
```
