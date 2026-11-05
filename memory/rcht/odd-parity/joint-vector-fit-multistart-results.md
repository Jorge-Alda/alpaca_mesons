---
title: Joint vector fit multi-start results
status: numerical
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Joint vector fit multi-start results

## Inputs, conventions, and method

The [joint likelihood baseline](babar-kkpi-global-vector-sector-fit-baseline.md) defines the 32 coordinates, 75 residuals, datasets, covariance approximations, and physical conventions. The [shared even-sector audit](shared-even-sector-gv-constraint-audit.md) records the $G_V$ constraints and branch assumptions. The [quadrature certification](kkpi-high-order-quadrature-certification.md) supplies the accepted comparison tolerances. This calculation does not change any likelihood input or assumption.

The executable is [`scripts/fit_joint_vector_multistart.py`](../../../scripts/fit_joint_vector_multistart.py). It used four deterministic starts—`stored`, `preliminary`, `stored_negative_c4`, and `stored_negative_odd`—with the existing rectangular bounds and `GlobalLikelihood.unpack` domain checks. The fixed optimizer-coordinate scale is the reciprocal residual-Jacobian column norm at the stored start; the six even physical coordinates retain the existing `EVEN_SCALES`. The likelihood was evaluated at energy order 2, pair order 16, angular order 16, in four-bin blocks, with trust-region reflective least squares and `max_nfev=1500`.

The exact command, bounds, starts, all coordinate scales, optimizer settings, Python/NumPy/SciPy/platform versions, likelihood input hashes, and loaded-source hashes are in the [attempt ledger](../../../data/vector_sector/joint_vector_multistart.json). Per-attempt records retain the final 32 coordinates, all 75 residuals, block contributions, both BaBar prediction arrays, active bounds, gradient and optimality, termination message, evaluation counts, elapsed time, and an independent endpoint-domain check.

## Results

Every attempt reached the 1,500-evaluation limit. All four endpoints passed the independent physical-domain and finite-residual check, but SciPy returned `success=false` and status 0 for each. They are retained as unsuccessful endpoints and none is an accepted minimum.

| Start | $\chi^2$ | $n_{\rm fev}$ | $n_{\rm jev}$ | Optimality | Termination |
| --- | ---: | ---: | ---: | ---: | --- |
| `stored` | 260.311519 | 1500 | 1490 | 3709.98 | Evaluation limit |
| `preliminary` | 337.143665 | 1500 | 1475 | 2791.93 | Evaluation limit |
| `stored_negative_c4` | 226.194810 | 1500 | 1430 | 9.84851 | Evaluation limit |
| `stored_negative_odd` | 311.497137 | 1500 | 1470 | 642.048 | Evaluation limit |

The lowest recorded objective is from `stored_negative_c4`, but its unsuccessful optimizer status excludes it from acceptance. The blockwise residual sums, predictions, final vectors, and full diagnostics for every endpoint remain in the JSON ledger; the four attempts have not been merged or discarded.

## Refinement and verification

The separate [refinement record](../../../data/vector_sector/joint_vector_multistart_refinement.json) verifies the saved input and source hashes. It reports `not_evaluated` for all four endpoints because the optimizer-success prerequisite was not met. Consequently this run supplies no high-order stability result and no accepted minimum. The certified $(2,64,64)\to(4,128,128)$ comparison remains valid only for the historical points documented in the quadrature-certification record.

The focused evaluation-path and ledger regressions passed, 22 tests total:

```text
.venv/bin/python -m unittest scripts.test_fit_joint_vector_multistart scripts.test_global_kkpi_running_propagator scripts.test_kkpi_quadrature_blocks
```

For the four named starts at small test quadrature, blocked predictions, residuals, and all Jacobian columns agree bit-for-bit with full-grid evaluation using the preceding separate-propagator and full-residual finite-difference paths. This verifies the runtime-improved calculation path; it does not turn the failed optimizer endpoints into fitted minima.

The four-start campaign was replayed with the same start vectors, scales, bounds, quadrature, and optimizer settings, writing only to `/tmp/joint_vector_multistart_replay.json`. All four final coordinate vectors and 75-component residual vectors matched the primary ledger bit-for-bit; each $\chi^2$, `nfev`, `njev`, and success status also matched exactly. The replay confirms reproducibility of these unsuccessful optimizer outcomes, not convergence or acceptance.
