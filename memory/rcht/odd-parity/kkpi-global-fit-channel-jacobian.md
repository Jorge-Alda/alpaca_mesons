---
title: Channel-specific finite differences in the global KKπ fit
status: numerical
created: 2026-11-04
updated: 2026-11-04
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Channel-specific finite differences in the global KKπ fit

## Inputs, conventions, and assumptions

The [global fit baseline](babar-kkpi-global-vector-sector-fit-baseline.md) defines the data, coordinates, likelihood, conventions and physical assumptions. The [running-propagator record](kkpi-global-fit-running-propagator-reuse.md) and [width-pull Jacobian record](kkpi-global-fit-width-pull-jacobian.md) describe earlier improvements. Invariants are in GeV squared, masses and widths in GeV, and cross-section predictions in nb. Pulls and residuals are dimensionless. The stored reference candidate remains a nonconverged fit snapshot.

In the existing model, the two tensor coefficients enter only the $K_S K\pi$ channel. Each pair of excited-vector coefficients and each normalization pull is specific to one of the two BaBar channels. None of these eight coordinates changes the first 28 residuals. This calculation uses those existing dependencies and introduces no additional physical approximation.

## Method

The executable is [`scripts/fit_babar_kkpi_global_vector_sector.py`](../../../scripts/fit_babar_kkpi_global_vector_sector.py). For the finite differences of the tensor and excited-vector coordinates, it evaluates only the affected channel at the shifted parameters. It copies the unchanged residual entries through subtraction of the baseline from itself, divides by the original step, and computes the affected entries with the same normalization and residual arithmetic as the full likelihood. Shifted physical quantities are still unpacked so that the existing parameter-domain checks apply.

For each normalization pull, the channel predictions from the baseline residual are reused. The shifted channel residual is formed from those predictions with the same multiplication, subtraction, and division as the full likelihood; the two normalization constraint entries are also differenced directly. Thus the forward-difference steps and their floating-point arithmetic are retained. A non-finite input or baseline, a penalty baseline, or an exception in the selected evaluation falls back to the full residual. Both full-grid and bounded-block integration support selecting a channel.

The regular Jacobian now calls the full residual 11 times rather than the 19 calls after the width-pull improvement. These channel dependencies must be revisited if the model is extended to couple the two channels or these coordinates to other residual blocks.

## Verification and timing

The regression script [`scripts/test_global_kkpi_running_propagator.py`](../../../scripts/test_global_kkpi_running_propagator.py) compares IEEE-754 bit patterns against the preceding full-residual Jacobian. The saved candidate passes at the default quadrature and on a smaller blocked grid. Jacobians at all 32 individual coordinate displacements pass on the smaller grid; the eight tensor, excited-vector, and normalization displacements also pass at default quadrature. Selected-channel predictions and residuals agree bit-for-bit with full integration in both grid modes. Additional checks cover a beta shift outside its permitted domain, a negative finite-difference step, an invalid width shift, and penalty and non-finite baselines.

All 32 focused, bounded-block, and quadrature-certification tests passed:

```text
.venv/bin/python -m unittest scripts.test_global_kkpi_running_propagator scripts.test_kkpi_quadrature_blocks scripts.test_kkpi_quadrature_certification scripts.test_certify_kkpi_quadrature
```

At the stored candidate on the machine and software identified in the [baseline profile](../../../data/vector_sector/babar_kkpi_global_vector_sector_fit_baseline_20261104_profile.json), seven alternating unprofiled default-grid Jacobian evaluations gave median times of $0.618$ s before and $0.482$ s after this channel-specific calculation, a ratio of $0.781$. This local timing does not measure the full fit. Following a warm-up, it used:

```bash
.venv/bin/python -c 'import json,sys,time,statistics,numpy as np; sys.path.insert(0,"scripts"); import fit_babar_kkpi_global_vector_sector as fit; from scripts.test_global_kkpi_running_propagator import width_optimized_jacobian; values=np.array(json.load(open("data/vector_sector/babar_kkpi_global_vector_sector_fit.json"))["parameter_values"]); model=fit.GlobalLikelihood(); model.jacobian(values); original=[]; improved=[]; exec("for index in range(7):\n start=time.perf_counter(); width_optimized_jacobian(model,values); original.append(time.perf_counter()-start)\n start=time.perf_counter(); model.jacobian(values); improved.append(time.perf_counter()-start)"); print("previous Jacobian seconds",original); print("improved Jacobian seconds",improved); print("median ratio improved/previous",statistics.median(improved)/statistics.median(original))'
```

The final implementation was run twice end to end and compared with the preserved baseline rerun. The baseline and both optimized runs used identical input hashes, initialization, quadrature and optimizer settings. One baseline run took $1955.96\,\mathrm{s}$; the optimized runs took $783.24\,\mathrm{s}$ and $779.23\,\mathrm{s}$, with median $781.235\,\mathrm{s}$ and range $[779.23,783.24]$ s. Relative to the single baseline timing, the optimized median is a factor $2.50$ faster and reduces wall time by $60.06\%$. The baseline repetition was omitted because each baseline run takes over 32 minutes. The saved baseline and both optimized outputs agree bit-for-bit in all 32 parameter values, 75 residuals, both BaBar prediction arrays, $\chi^2=307.984932373217$, the `nfev`/`njev` counts $(1500,1453)$, and the nonconverged status and message. Exact commands, machine details, input hashes and comparison fields are in the [baseline profile](../../../data/vector_sector/babar_kkpi_global_vector_sector_fit_baseline_20261104_profile.json); the optimized outputs are [first run](../../../data/vector_sector/babar_kkpi_global_vector_sector_fit_optimized_final_20261104.json) and [repeat](../../../data/vector_sector/babar_kkpi_global_vector_sector_fit_optimized_repeat_20261104.json).
