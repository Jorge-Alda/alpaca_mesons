---
title: Outcome of the joint vector fit with the double-trace singlet-source coordinate (task 4.2)
status: numerical
created: 2026-11-17
updated: 2026-11-17
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh]
provenance: numerical-computation
---

# Outcome of the joint vector fit with the double-trace singlet-source coordinate $y_V$

This record gives the result of task 4.2 of `extend-joint-vector-model-diagnostics`: the joint likelihood with the coefficient $y_V$ of the double-trace singlet-source structure free ([option record](joint-vector-double-trace-singlet-source-likelihood-option.md)), refitted with the bounds, starts and gates fixed beforehand in the [fit plan](joint-vector-double-trace-singlet-source-fit-plan.md). **The adequacy criterion is not met.** $y_V$ lowers $\chi^2$ by $1.08$ and is compatible with zero.

## Ledger

`data/vector_sector/joint_vector_double_trace_singlet_source_fit.json` (per-run parts in `...fit.json.parts/`, driver log `...fit.log`). It records the command, optimizer settings, bounds, starts, scales, source-hash gate (passed), input and source SHA-256 values and all seven runs, including the two that did not converge. Command:

```bash
.venv/bin/python scripts/fit_joint_vector_double_trace_singlet_source.py --output data/vector_sector/joint_vector_double_trace_singlet_source_fit.json --workers 2
```

Settings: scipy `least_squares` trf, `max_nfev=5000`, `ftol=1e-8`, `xtol=1e-9`, `gtol=1e-7`, quadrature $(2,16,16)$ with four-bin blocks, re-evaluation at $(4,128,128)$, quadrature ladder from $(2,64,64)$. Wall time 2.6 h with two workers.

## Runs

Conventions of Alda:2026nrd [@Alda:2026nrd]; $y_V$ in GeV$^{-1}$; $\chi^2$ at $(4,128,128)$; "distance" is the scaled distance of the 32 baseline coordinates to the best baseline endpoint.

| start | optimizer | nfev | $\chi^2$ | dof | ladder | accepted | $y_V$ | distance |
|---|---|---|---|---|---|---|---|---|
| $y_V=0$ from `basin_probes/even_shared_widths` | `ftol` (status 2) | 2049 | 224.6477 | 42 | certified | yes | 0.001325 | 9.23 |
| $y_V=0$ from `basin_probes/babar` | `ftol` (status 2) | 2349 | 224.6477 | 42 | certified | yes | 0.001326 | 9.24 |
| $y_V=0$ from `continuation_4128128/stored_negative_c4` | `ftol` (status 2) | 1826 | **224.6474** | 42 | certified | yes (selected) | 0.001328 | 9.52 |
| $y_V=0$ from `basin_probes/odd` | `ftol` (status 2) | 1984 | 224.6480 | 42 | certified | yes | 0.001324 | 9.17 |
| $y_V=+0.022$ from best | evaluation cap (status 0) | 5000 | 224.6668 | 42 | not evaluated | no | 0.001351 | 12.96 |
| $y_V=-0.022$ from best | evaluation cap (status 0) | 5000 | 224.6485 | 42 | not evaluated | no | 0.001327 | 9.47 |
| `zero_control` ($y_V=0$ fixed) | `ftol` (status 2) | 27 | 225.7266 | 43 | certified | yes | 0 | 0.31 |

All seven runs reach the same region: the four accepted runs agree to $6\times10^{-4}$ in $\chi^2$, and both offset starts, though stopped by the evaluation cap, end within $0.02$ of it with the same $y_V$. No bound is active and no near-singular direction is found. The final optimality is $1.1$–$2.6$, so the optimizer stops on `ftol` rather than on a zero gradient, as in the earlier fits. The large scaled distances come from flat directions of the baseline minimum: the largest shift of any baseline coordinate is $0.92$ baseline local errors ($\eta'$ total-width pull) in every run.

**Control.** `zero_control` gives $\chi^2=225.7266$, equal to the baseline $225.727$ within $4\times10^{-4}$, with the largest coordinate shift $0.003$ local errors. Its scaled distance $0.307$ exceeds the control tolerance $0.01$, so the ledger flag `reproduces_baseline_endpoint` is false. The task 3.2 control (`joint_vector_isospin_breaking_vpp_fit.json`) has the same numbers ($\chi^2=225.72664$, distance $0.307$, flag false), as expected because both likelihoods equal the baseline bitwise at zero coefficient. The distance is set by the optimizer stopping on `ftol` along flat directions, not by the option.

## Selected endpoint and adequacy

Selected: the run from `continuation_4128128/stored_negative_c4`.

- $\chi^2_{(4,128,128)}=224.647$ for $75-33=42$ nominal dof; total nominal $p=8.4\times10^{-27}$. Baseline: $225.727$ for 43 dof, $p=1.3\times10^{-26}$.
- $\Delta\chi^2=1.08$ for one added coordinate.
- $y_V=(1.33\pm1.22)\times10^{-3}$ GeV$^{-1}$ (linearized local error from the $(4,128,128)$ Jacobian, not a profile interval), compatible with zero at about $1.1\sigma$. This is about $0.8\%$ of $|z_V|$ and about $6\%$ of the stage-3 scale $0.022$ GeV$^{-1}$, far inside the bounds $\pm0.5$ GeV$^{-1}$.

| block (rows) | $\chi^2$ | nominal $p$ (row-count dof) | baseline $\chi^2$ |
|---|---|---|---|
| odd (12) | 34.71 | $5.2\times10^{-4}$ | 35.18 |
| masses and leptonic widths (7) | 14.56 | $0.042$ | 12.94 |
| shared widths (7) | 104.81 | $1.1\times10^{-19}$ | 106.84 |
| $\phi\to K\bar K$ (2) | 17.99 | $1.2\times10^{-4}$ | 17.88 |
| BaBar (45) | 52.53 | $0.20$ | 52.83 |
| normalization constraints (2) | 0.06 | $0.97$ | 0.06 |

**Criterion** (total nominal $p>0.01$ and no block with nominal $p<10^{-3}$): **not met.** The total $p$ is $8.4\times10^{-27}$ and three blocks fail: shared widths (the worst), $\phi\to K\bar K$ and odd. The failing blocks are the same as in the baseline.

## Interpretation

- At fixed other coordinates the shared-width and $\phi\to K\bar K$ rows do not depend on $y_V$ ([option record](joint-vector-double-trace-singlet-source-likelihood-option.md), "Rows affected"), because $y_V$ changes only photon residues and no vector mass or $VPP$ coupling. These two blocks, which carry most of the baseline $\chi^2$, can move only through readjusting the other coordinates, and they move by less than 2 units.
- The odd block falls by $0.5$ and the BaBar block by $0.3$, while the masses-and-leptonic-widths block rises by $1.6$: the gain in the $\omega,\phi$ photon couplings of the radiative and $KK\pi$ rows is paid in the leptonic widths.
- The double-trace singlet-source structure therefore does not remove the inadequacy. This agrees with the [4.1 result](joint-vector-double-trace-source-mass-and-kkpi.md) that the term cannot move the $\omega,\phi$ masses, which was the reason it was listed as candidate cause 3 in the [candidate-causes record](joint-vector-fit-inadequacy-candidate-causes.md).

## Assumptions and limits

- Nominal $p$-values; the per-block values use the row count as dof (assumption of alp-uzh.6).
- $\chi^2$ values are upper bounds on the minimum within the `ftol` noise (about $10^{-3}$ here, from the spread of the four accepted runs).
- The $y_V$ error is linearized, not a profile interval.
- Model assumptions of the option record (order-$e^2$ mass shifts omitted, exact strong isospin, excited poles unchanged, baseline source-branch check kept).
- The extension was fitted alone, not combined with the isospin-breaking $VPP$ coordinates of [task 3.2](joint-vector-isospin-breaking-vpp-fit-outcome.md).

## Consequence for the plan

Task 4.2 asks for the implementation, the refit, the criterion check and the ledger entry; all are done, and the criterion fails. The next section of the change is section 5 (isospin breaking in the even sector).
