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

1. Evaluate the reference propagator at the requested invariant mass with both excited-state coefficients set to zero and with their requested values. An additional array dimension selects between the two coefficient choices, so the invariant-dependent widths are calculated once.
2. Invert that value to obtain the imaginary denominator and rescale its running width to the fitted shared on-shell width, as in the original implementation.
3. Subtract the zero-coefficient result from the requested-coefficient result to obtain the excited-state contribution.
4. Add that contribution to the rescaled ground-state propagator. Both values from step 1 are used only within this single calculation; a changed invariant, fit parameter, or quadrature grid is evaluated again.

The preceding calculation called the reference propagator four times for each running propagator; the revised calculation calls it twice, including once for the on-shell normalization. The same combined evaluation replaces the two separate reference calls for each neutral excited-state contribution. The original subtraction and addition order is retained. The output remains the same 75-component residual, 45 BaBar predictions, and 75-by-32 Jacobian used by the global fit.

## Verification and timing

The regression script [`scripts/test_global_kkpi_running_propagator.py`](../../../scripts/test_global_kkpi_running_propagator.py) reproduces the preceding running-propagator and separate reference-propagator calculations independently. It compares the full residual vector, both prediction arrays, and every Jacobian entry with the preceding calculation at the stored candidate. It also compares residuals and predictions after displacing each of the 32 coordinates, returning to the candidate, and using both the default full grid and a smaller blocked grid. These comparisons pass bit-for-bit. The default-grid residuals and predictions also agree with the stored baseline at relative tolerance $10^{-12}$ and absolute tolerance $10^{-11}$. Direct checks compare the combined and separate evaluations for all four vector species, scalar and array invariants, and zero and nonzero excited-state coefficients. They confirm that the revised running propagator calls the reference propagator twice and evaluates changed inputs afresh.

The focused regression and existing block tests both passed after the further [width-pull Jacobian improvement](kkpi-global-fit-width-pull-jacobian.md) and [channel-specific calculation](kkpi-global-fit-channel-jacobian.md):

```text
.venv/bin/python -m unittest scripts.test_global_kkpi_running_propagator scripts.test_kkpi_quadrature_blocks
```

On the machine and software environment identified in the [baseline profile](../../../data/vector_sector/babar_kkpi_global_vector_sector_fit_baseline_20261104_profile.json), five alternating unprofiled Jacobian evaluations at the stored candidate gave median times of $1.276$ s for the preceding calculation and $0.667$ s for the revised calculation, a ratio of $0.523$. This local timing predates the width-pull improvement; it is not a measurement of full-fit runtime. It was obtained after a warm-up with:

```bash
.venv/bin/python -c 'import json,sys,time,statistics,numpy as np; from unittest.mock import patch; sys.path.insert(0,"scripts"); import fit_babar_kkpi_global_vector_sector as fit; from scripts.test_global_kkpi_running_propagator import original_running_propagator,original_reference_propagator_pair; x=np.array(json.load(open("data/vector_sector/babar_kkpi_global_vector_sector_fit.json"))["parameter_values"]); m=fit.GlobalLikelihood(); m.jacobian(x); original=[]; improved=[]; exec("for index in range(5):\n with patch.object(fit,\"running_propagator\",side_effect=original_running_propagator),patch.object(fit,\"reference_propagator_pair\",side_effect=original_reference_propagator_pair):\n  start=time.perf_counter(); m.jacobian(x); original.append(time.perf_counter()-start)\n start=time.perf_counter(); m.jacobian(x); improved.append(time.perf_counter()-start)"); print("original seconds",original); print("improved seconds",improved); print("median ratio improved/original",statistics.median(improved)/statistics.median(original))'
```
