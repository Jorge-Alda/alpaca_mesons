---
title: Prepared joint vector fit with the isospin-breaking VPP option (plan and dry run)
status: numerical
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh]
provenance: numerical-computation
---

# Joint vector fit with the isospin-breaking $VPP$ option: plan and dry run

**The fit has not been run.** This record fixes the starting points, bounds, quadrature gates, acceptance and adequacy criteria and degrees of freedom of the joint vector fit with the coordinates $r_1,r_2,r_3$ switched on, before any optimizer call (as required by `openspec/changes/extend-joint-vector-model-diagnostics/design.md`, "Added coordinates get explicit bounds and starting values recorded before running"). The only numbers below come from a dry run: one residual and one Jacobian evaluation at each starting point.

Driver: [`scripts/fit_joint_vector_isospin_breaking_vpp.py`](../../../scripts/fit_joint_vector_isospin_breaking_vpp.py); tests: [`scripts/test_fit_joint_vector_isospin_breaking_vpp.py`](../../../scripts/test_fit_joint_vector_isospin_breaking_vpp.py).

Command to run the fit later (writes `data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json`, not yet created):

```bash
.venv/bin/python scripts/fit_joint_vector_isospin_breaking_vpp.py --output data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json
```

Dry run and tests (the only commands executed):

```bash
.venv/bin/python scripts/fit_joint_vector_isospin_breaking_vpp.py --dry-run
.venv/bin/python -m unittest scripts.test_fit_joint_vector_isospin_breaking_vpp
```

## Conventions and units

Conventions of Alda:2026nrd and of the baseline likelihood [@Alda:2026nrd]. GeV, nb; $r_i$ in GeV$^{-2}$; residuals, $\chi^2$ and scaled distances dimensionless. The notation $r_i$ is the one proposed (awaiting approval) in [joint-vector-isospin-breaking-vpp-operators.md](joint-vector-isospin-breaking-vpp-operators.md).

## Inputs

- Likelihood option `IsospinBreakingVppLikelihood` with `isospin_breaking_vpp=True` and the default switches `vpp_kstar_excitation_scaling='ground_only'`, `vpp_s_channel_factors=True`: [joint-vector-isospin-breaking-vpp-likelihood-option.md](joint-vector-isospin-breaking-vpp-likelihood-option.md), $s$-channel factors in [joint-vector-isospin-breaking-vpp-s-channel-factors.md](joint-vector-isospin-breaking-vpp-s-channel-factors.md).
- Scan point $r=(-0.350,\,1.855,\,0.729)$ GeV$^{-2}$: `scripts/scan_isospin_breaking_vpp_simultaneous.py`, quoted in the option record.
- Accepted endpoints (best: `basin_probes/even_shared_widths`) and the 32 fixed optimizer scales: [joint-vector-fit-endpoint-comparison.md](joint-vector-fit-endpoint-comparison.md).
- Baseline $(4,128,128)$ $\chi^2=225.727$ and local errors: [joint-vector-fit-adequacy-assessment-method.md](joint-vector-fit-adequacy-assessment-method.md), [joint-vector-fit-goodness-of-fit.md](joint-vector-fit-goodness-of-fit.md).
- Refit procedure, gates and ledger format of alp-uzh.6: [joint-vector-block-exclusion-and-variation-refits.md](joint-vector-block-exclusion-and-variation-refits.md), [joint-vector-fit-sensitivity-variation-results.md](joint-vector-fit-sensitivity-variation-results.md).
- Data: PDG rows [@ParticleDataGroup:2024cfk], BaBar $KK\pi$ [@BaBar:2007ceh], through the baseline likelihood [babar-kkpi-global-vector-sector-fit-baseline.md](babar-kkpi-global-vector-sector-fit-baseline.md).

## Coordinates, bounds and degrees of freedom

- 35 coordinates: the 32 baseline coordinates with their baseline bounds (`fit_babar_kkpi_global_vector_sector.parameter_bounds`), followed by $r_1,r_2,r_3$ with bounds $[-5,5]$ GeV$^{-2}$ each (the scan range).
- 75 residual rows (73 measurements and 2 normalization constraints). Nominal degrees of freedom $75-35=40$ (baseline 43). Per-block $p$-values use the row count as degrees of freedom (assumption, as in alp-uzh.6).

## Starting points

The 32 baseline coordinates of both starts are those of the accepted endpoint `basin_probes/even_shared_widths`.

| start | $(r_1,r_2,r_3)$ [GeV$^{-2}$] | $r$ | free coordinates | dof |
|---|---|---|---|---|
| `scan_point` | $(-0.350,\,1.855,\,0.729)$ | free | 35 | 40 |
| `zero_control` | $(0,0,0)$ | held fixed | 32 | 43 |

**Choice (labelled):** in `zero_control` the $r_i$ are held fixed at zero. With $r=0$ the likelihood equals the baseline bitwise, so this refit must reproduce the baseline endpoint: scaled distance to the best endpoint $\le0.01$ and $|\chi^2_{(4,128,128)}-225.727|\le0.01$ (the tolerances of the endpoint comparison). A start at $r=0$ with $r$ free is not part of the plan.

## Method (gates, identical to the alp-uzh.6 variation refits)

1. Source-hash gate before anything else: the six source ledgers of the endpoint comparison must have the hashes recorded in it and pass `compare_joint_vector_endpoints.validate_ledgers`; all 42 input and source hashes recorded in `data/vector_sector/joint_vector_sensitivity_variations.json` must match the current files.
2. The ledger with model options, parameter names, bounds, starts, optimizer scales, gates and input/source hashes is written (state `prepared`) before any optimizer call.
3. Refit with scipy `least_squares` (trf, `max_nfev=5000`, `ftol=1e-8`, `xtol=1e-9`, `gtol=1e-7`) at quadrature $(2,16,16)$ with four-bin blocks, fixed optimizer scales. Optimizer status, message, `nfev`, `njev`, optimality, active bounds and near-singular directions (relative singular value $<10^{-6}$ of the scaled Jacobian) are recorded; failed refits are kept.
4. Optimizer scales: the 32 scales of the endpoint comparison; for $r_1,r_2,r_3$ the reciprocal Jacobian column norms at the `zero_control` start at $(2,16,16)$ (choice; same definition as the baseline scales). Dry-run values: $9.76\times10^{-4}$, $5.83\times10^{-2}$, $1.82\times10^{-3}$ GeV$^{-2}$.
5. Re-evaluation at $(4,128,128)$: $\chi^2$, dof, nominal $p$, block $\chi^2$.
6. Local errors $\sqrt{\mathrm{diag}((J^TJ)^{-1})}$ of the free coordinates from the $(4,128,128)$ Jacobian (linearized, not profile intervals).
7. Quadrature ladder of the endpoint comparison, from $(2,64,64)$ with energy orders $\{2,4\}$, pair and angle orders $\{64,128\}$. A run is accepted when the optimizer succeeds and the ladder is certified.
8. Multistart: each start is refit from the best accepted endpoint; if the scaled distance of its endpoint to that endpoint exceeds 1, it is refit from the other three accepted endpoints with the same $r$ start. The selected run has the lowest $(4,128,128)$ $\chi^2$ among accepted runs. The scaled distance and the parameter shifts (in baseline local errors) use the 32 baseline coordinates only; $r$ is reported in GeV$^{-2}$ with its own local errors.

## Acceptance and adequacy criteria

- Run accepted: optimizer success and certified quadrature ladder (step 7).
- Adequacy (criterion of `openspec/changes/assess-joint-vector-fit-adequacy/design.md`): total nominal $p>0.01$ and no block with nominal $p<10^{-3}$, at $(4,128,128)$, with 40 dof for the free-$r$ fit.

## Dry-run result (numerical; starting points only, not a fit)

Orders $(2,16,16)$, four-bin blocks; 75 rows, 35 Jacobian columns at both starts; residuals and Jacobians finite; no penalty value. Option on with default switches, 35 coordinates and the bounds above checked; source-hash gate passed.

| block (rows) | $\chi^2$ at `scan_point` | $\chi^2$ at `zero_control` |
|---|---|---|
| odd (12) | 374.65 | 35.18 |
| masses and leptonic widths (7) | 12.93 | 12.93 |
| shared widths (7) | 3581.20 | 106.84 |
| $\phi\to K\bar K$ (2) | 649.28 | 17.88 |
| BaBar (45) | 67.77 | 52.77 |
| normalization constraints (2) | 0.06 | 0.06 |
| total | 4685.90 | 225.67 |

At `zero_control` the residual and the first 32 Jacobian columns are bitwise equal to those of `GlobalLikelihood(2,16,16)` at the accepted endpoint. At `scan_point`, with all baseline coordinates held at the accepted endpoint, both $K^*$ widths drop by about half (option record), which gives the large shared-width $\chi^2$; this says nothing about the minimum of the fit, where $G_V$ and the other coordinates readjust.

## Verification

`scripts.test_fit_joint_vector_isospin_breaking_vpp` (19 tests, run 2026-11-09, pass; cheap orders $(1,3,4)$, optimizer always mocked): argument parsing; refit settings equal to those of `refit_joint_vector_sensitivity.py`; model check and rejection of non-default switches; bounds; start construction and rejection of $|r_i|>5$; `zero_control` residual bitwise equal to the baseline; prepared-ledger schema and JSON round trip; refit record keys with a mocked optimizer and kept failure records; `--dry-run` never calls `least_squares`, the run functions or the process pool, and `main --dry-run` writes no ledger; adequacy and control-reproduction rules; the source-hash gate detects a changed hash.

Not checked: any optimizer run, the $(4,128,128)$ evaluation and the quadrature ladder with 35 coordinates (exercised only through the shared baseline code with renamed column labels).
