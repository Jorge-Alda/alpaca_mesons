---
title: Joint vector fit with the double-trace singlet-source coordinate (bounds, starts and dry run)
status: numerical
created: 2026-11-17
updated: 2026-11-17
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh]
provenance: numerical-computation
---

# Joint vector fit with the double-trace singlet-source coordinate $y_V$: bounds, starts and dry run

This record fixes the bounds, starting values, optimizer scales, gates and degrees of freedom of the refit of task 4.2 of `extend-joint-vector-model-diagnostics`, **before any optimizer call**, as required by the change design ("Added coordinates get explicit bounds and starting values recorded before running"). The only numbers below come from a dry run: one residual and one Jacobian evaluation at each start. The result of the refit is in the [outcome record](joint-vector-double-trace-singlet-source-fit-outcome.md).

Driver: [`scripts/fit_joint_vector_double_trace_singlet_source.py`](../../../scripts/fit_joint_vector_double_trace_singlet_source.py); likelihood option: [joint-vector-double-trace-singlet-source-likelihood-option.md](joint-vector-double-trace-singlet-source-likelihood-option.md).

```bash
.venv/bin/python scripts/fit_joint_vector_double_trace_singlet_source.py --dry-run
.venv/bin/python scripts/fit_joint_vector_double_trace_singlet_source.py --output data/vector_sector/joint_vector_double_trace_singlet_source_fit.json --workers 2
```

## Conventions and units

Conventions of Alda:2026nrd and of the baseline likelihood [@Alda:2026nrd]. GeV, nb; $x_q$ in GeV$^2$; $z_V$ and $y_V$ in GeV$^{-1}$; residuals, $\chi^2$ and scaled distances dimensionless. $y_V$ is the coefficient of the double-trace singlet-source structure $y_V\langle\chi_+Q\rangle\mathbf 1$ of the photon source matrix ([4.1 record](joint-vector-double-trace-source-mass-and-kkpi.md)); the notation was approved by the user.

## Inputs

- Likelihood option `DoubleTraceSingletSourceLikelihood(double_trace_singlet_source=True)`: [option record](joint-vector-double-trace-singlet-source-likelihood-option.md).
- Accepted endpoints of alp-uzh.6 and the 32 fixed optimizer scales: [joint-vector-fit-endpoint-comparison.md](joint-vector-fit-endpoint-comparison.md), [joint-vector-fit-multistart-results.md](joint-vector-fit-multistart-results.md).
- Baseline $(4,128,128)$ $\chi^2=225.727$ for 43 dof and local errors: [joint-vector-fit-adequacy-assessment-method.md](joint-vector-fit-adequacy-assessment-method.md), [joint-vector-fit-goodness-of-fit.md](joint-vector-fit-goodness-of-fit.md).
- Refit gates and ledger format: [joint-vector-block-exclusion-and-variation-refits.md](joint-vector-block-exclusion-and-variation-refits.md); template of task 3.2: [joint-vector-isospin-breaking-vpp-fit-plan.md](joint-vector-isospin-breaking-vpp-fit-plan.md).
- Stage-3 double-trace values of $y_V$: `data/vector_sector/stage3_results.json`, scenario $\rho^\pm$, $K^{*\pm}$, $\phi$ evaluation, case `extended_mass_branch_1` (near-ideal): $y_V=0.0220$ GeV$^{-1}$ with $F_V=0.157$ GeV, $z_V=-0.0037$ GeV$^{-1}$. That solution has the residue signs of the universal source, so it lies on a different source sign branch than the baseline ($F_s<0$); it is used here only as a scale.
- Source parameters at the best accepted endpoint `basin_probes/even_shared_widths`: $F_V=0.1690$ GeV, $z_V=-0.1710$ GeV$^{-1}$, $F_l=0.1564$ GeV, $F_s=-0.1537$ GeV, with $x_l=0.03662$ GeV$^2$, $x_s=0.94368$ GeV$^2$ (`stage3_results.json` assumptions).

## Coordinates, bounds and degrees of freedom

- 33 coordinates: the 32 baseline coordinates with their baseline bounds, followed by $y_V$ (`double_trace_singlet_source_y_V_GeV_inverse`).
- **Bounds (choice): $-0.5\le y_V\le0.5$ GeV$^{-1}$.** Reasons:
  1. it is about three times $|z_V|=0.171$ GeV$^{-1}$, the single-trace source coefficient of the same dimension at the accepted endpoint;
  2. it is about 23 times the stage-3 scale $|y_V|\approx0.022$ GeV$^{-1}$;
  3. at the bound the shifts of the isoscalar insertion factors are $2|y_V|(x_s-x_l)=0.91$ GeV $\approx5.8F_l$ for $\omega_I$ and $|y_V|(x_s-x_l)=0.45$ GeV $\approx3|F_s|$ for $\phi_I$, so both factors can change sign and the range covers the baseline's source scale several times over;
  4. it is about $10^3$ optimizer scales ($4.6\times10^{-4}$ GeV$^{-1}$, below) from $y_V=0$.
  Outside the bounds the residual is the baseline penalty vector. The baseline source-branch check ($F_l>0$, $F_s<0$) is kept and acts on $F_V,z_V$ only (assumption, stated in the option record).
- 75 residual rows. Nominal degrees of freedom $75-33=42$ (baseline 43); the control with $y_V$ fixed has 43. Per-block $p$-values use the row count as degrees of freedom (assumption, as in alp-uzh.6).

## Starting points

| start | endpoint of the other 32 coordinates | $y_V$ start [GeV$^{-1}$] | $y_V$ |
|---|---|---|---|
| `y_zero_from_basin_probes__even_shared_widths` | `basin_probes/even_shared_widths` (best) | $0$ | free |
| `y_zero_from_basin_probes__babar` | `basin_probes/babar` | $0$ | free |
| `y_zero_from_continuation_4128128__stored_negative_c4` | `continuation_4128128/stored_negative_c4` | $0$ | free |
| `y_zero_from_basin_probes__odd` | `basin_probes/odd` | $0$ | free |
| `y_plus_stage3_scale_from_best` | best | $+0.022$ | free |
| `y_minus_stage3_scale_from_best` | best | $-0.022$ | free |
| `zero_control` | best | $0$ | fixed |

**Choices (labelled):** $y_V=0$ is the baseline, so the four accepted minima with $y_V=0$ are the natural multistart; all four are used, so the basin-switch rerun of the template is not needed. The two offset starts use the stage-3 magnitude with both signs, because the stage-3 solution lies on another source sign branch and its sign carries no information here. `zero_control` must reproduce the baseline endpoint: scaled distance $\le0.01$ and $|\chi^2_{(4,128,128)}-225.727|\le0.01$.

## Method (gates, identical to the alp-uzh.6 refits)

1. Source-hash gate (the accepted-endpoint ledgers and every file hashed in `data/vector_sector/joint_vector_sensitivity_variations.json` unchanged).
2. Ledger with options, names, bounds, starts, scales, gates and hashes written (state `prepared`) before any optimizer call.
3. scipy `least_squares` (trf, `max_nfev=5000`, `ftol=1e-8`, `xtol=1e-9`, `gtol=1e-7`) at quadrature $(2,16,16)$ with four-bin blocks; status, message, `nfev`, `njev`, optimality, active bounds and near-singular directions recorded; failed runs kept.
4. Optimizer scales: the 32 scales of the endpoint comparison, and for $y_V$ the reciprocal Jacobian column norm at the `zero_control` start at $(2,16,16)$ (choice; same definition as the baseline scales): $4.640\times10^{-4}$ GeV$^{-1}$.
5. Re-evaluation at $(4,128,128)$: $\chi^2$, dof, nominal $p$, block $\chi^2$; local errors $\sqrt{\mathrm{diag}((J^TJ)^{-1})}$ from the $(4,128,128)$ Jacobian (linearized).
6. Quadrature ladder of the endpoint comparison from $(2,64,64)$. A run is accepted when the optimizer succeeds and the ladder is certified.
7. Selected run: lowest $(4,128,128)$ $\chi^2$ among accepted runs with $y_V$ free.
8. Adequacy: total nominal $p>0.01$ and no block with nominal $p<10^{-3}$ at $(4,128,128)$, with 42 dof.

## Dry-run result (starting points only, not a fit)

Orders $(2,16,16)$, four-bin blocks; 75 rows, 33 Jacobian columns; residuals and Jacobians finite; source-hash gate passed. At every $y_V=0$ start the residual and the first 32 Jacobian columns are bitwise equal to those of `GlobalLikelihood(2,16,16)` at the same endpoint.

| block (rows) | $y_V=0$ (best endpoint) | $y_V=+0.022$ | $y_V=-0.022$ |
|---|---|---|---|
| odd (12) | 35.18 | 3147.16 | 994.81 |
| masses and leptonic widths (7) | 12.93 | 287.30 | 403.72 |
| shared widths (7) | 106.84 | 106.84 | 106.84 |
| $\phi\to K\bar K$ (2) | 17.88 | 17.88 | 17.88 |
| BaBar (45) | 52.77 | 80.05 | 98.56 |
| normalization constraints (2) | 0.06 | 0.06 | 0.06 |
| total | 225.67 | 3639.29 | 1621.87 |

The other three $y_V=0$ starts give totals $225.668$, $225.669$, $225.668$. With the other coordinates held at the accepted endpoint, $y_V=\pm0.022$ GeV$^{-1}$ moves the $\omega,\phi$ photon residues by 15–30% (option record), which the odd radiative rows and the $\omega,\phi$ leptonic widths do not tolerate. This says nothing about the minimum, where $F_V$, $z_V$ and the odd coefficients readjust. The shared-width and $\phi\to K\bar K$ blocks do not depend on $y_V$.

## Verification

Option tests: [option record](joint-vector-double-trace-singlet-source-likelihood-option.md). Driver tests: [`scripts/test_fit_joint_vector_double_trace_singlet_source.py`](../../../scripts/test_fit_joint_vector_double_trace_singlet_source.py), 10 tests, run 2026-11-17, pass (orders $(1,3,4)$; optimizer mocked or forbidden): start definitions (four $y_V=0$ starts on the four accepted endpoints, $\pm0.022$, one fixed control); model check and 42 dof; every $y_V=0$ start bitwise equal to the baseline residual; point validation; optimizer scales; refit record keys with the fixed coordinate kept and failures kept; selection rule; prepared-ledger schema and JSON round trip; the dry run never calls the optimizer.

```bash
.venv/bin/python -m unittest scripts.test_fit_joint_vector_double_trace_singlet_source
``` Not checked by the dry run: any optimizer run, the $(4,128,128)$ evaluation and the quadrature ladder with 33 coordinates.
