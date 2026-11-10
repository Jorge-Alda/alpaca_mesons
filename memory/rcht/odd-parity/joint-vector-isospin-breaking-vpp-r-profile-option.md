---
title: Tolerance options and the profile option in r of the isospin-breaking VPP joint fit driver (dry run only)
status: numerical
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh]
provenance: numerical-computation
---

# Tolerance options and profile option in $r$ of the isospin-breaking $VPP$ fit driver

**No fit and no profile point was run.** This record adds two capabilities to the driver of [joint-vector-isospin-breaking-vpp-fit-plan.md](joint-vector-isospin-breaking-vpp-fit-plan.md), after the restart option of [joint-vector-isospin-breaking-vpp-restart-option.md](joint-vector-isospin-breaking-vpp-restart-option.md):

- command-line optimizer tolerances;
- a profile in one of the coefficients $r_1, r_2, r_3$.

It reports one dry-run evaluation.

Driver: [`scripts/fit_joint_vector_isospin_breaking_vpp.py`](../../../scripts/fit_joint_vector_isospin_breaking_vpp.py); tests: [`scripts/test_fit_joint_vector_isospin_breaking_vpp.py`](../../../scripts/test_fit_joint_vector_isospin_breaking_vpp.py).

```bash
.venv/bin/python scripts/fit_joint_vector_isospin_breaking_vpp.py --dry-run --restart-from data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json --profile-r 1=-1,-0.5,0,0.5,1
# profile (not run): drop --dry-run; optional --workers N, --output PATH, --overwrite, --max-nfev, --ftol, --xtol, --gtol
```

## Conventions and units

These are the same as in the fit plan: GeV and nb, with $r_i$ in GeV$^{-2}$. Residuals and $\chi^2$ are dimensionless. Conventions follow Alda:2026nrd [@Alda:2026nrd].

## Inputs

- The likelihood option with coordinates $r_1, r_2, r_3$, bounded to $[-5,5]$ GeV$^{-2}$: [joint-vector-isospin-breaking-vpp-likelihood-option.md](joint-vector-isospin-breaking-vpp-likelihood-option.md). It is evaluated through the bitwise-equal class of [joint-vector-isospin-breaking-vpp-calculation-reuse.md](joint-vector-isospin-breaking-vpp-calculation-reuse.md).
- The restart endpoint, selected by the rule of [joint-vector-isospin-breaking-vpp-restart-option.md](joint-vector-isospin-breaking-vpp-restart-option.md).
- The data are the PDG rows [@ParticleDataGroup:2024cfk] and BaBar $KK\pi$ [@BaBar:2007ceh], entering through the baseline likelihood.
- **Assumption, for the time estimate only:** 9 minutes of wall time per profile point. This value was supplied with the task. The driver does not measure it, and it was not checked against [joint-vector-isospin-breaking-vpp-evaluation-cost.md](joint-vector-isospin-breaking-vpp-evaluation-cost.md).

## Options

1. **`--ftol`, `--xtol`, `--gtol`.** Each takes a finite positive float. The defaults are $10^{-8}$, $10^{-9}$ and $10^{-7}$, the values of `OPTIMIZER_SETTINGS`. They replace the corresponding `least_squares` arguments and are written to the ledger's `optimizer_settings`, keeping the key order. With the defaults the block is identical, as a JSON string, to the one stored in `data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json`. The options apply to the two-start fit, the restart and the profile.
2. **`--profile-r INDEX=V1,V2,...`.** INDEX is 1, 2 or 3, for $r_1$, $r_2$, $r_3$. The values are the grid in GeV$^{-2}$. The option can be repeated, and each flag is a separate 1-D scan. The driver refuses:
   - an INDEX other than 1, 2 or 3, a missing `=`, or a value that is not a number;
   - a value outside $[-5,5]$ GeV$^{-2}$ or not finite;
   - a repeated value within a grid, or a repeated INDEX across flags;
   - `--profile-r` without `--restart-from`, because the profile starts from the restart endpoint.

   For every grid point:
   - the chosen $r_i$ is held fixed at the grid value, through the same fixed-coordinate mechanism as the `zero_control` start;
   - the 32 baseline coordinates and the other two $r$ are refit at $(2,16,16)$ from the restart endpoint, with the fixed value substituted;
   - the endpoint is re-evaluated at $(4,128,128)$.

   Each point records:
   - the fixed value;
   - the optimizer state, status, message, `nfev`, `njev` and optimality;
   - $\chi^2$ at $(2,16,16)$ and at $(4,128,128)$, and the $(4,128,128)$ $\chi^2$ per block;
   - the 35 endpoint coordinates;
   - `profile: true` and the optimizer settings.

   The full refit record and the high-order evaluation are kept as well. The ledger has `study: isospin_breaking_vpp_r_profile`, a `profile_definition` and a `profile` list. Its default path is `data/vector_sector/joint_vector_isospin_breaking_vpp_r_profile.json`. The driver prints the number of points before any refit.
3. **`--profile-grid`.** This makes the points the Cartesian product of the `--profile-r` grids, so several $r$ are fixed at once. Without it, a point that fixes more than one $r$ is refused. `--profile-grid` without `--profile-r` is refused.
4. **`--overwrite`.** This is used only for the profile. Without it, the driver refuses to start if the output ledger exists. Finished points stored in `OUTPUT.parts/` are reused only when their start and optimizer settings are identical. Interrupted points resume from their progress file, with the same checkpoint rule as the restart.
5. **`--dry-run` with `--profile-r`.** The dry run does not call `least_squares`, and it writes nothing. It does the following:
   - checks the model and the source-hash gate;
   - builds the restart start and checks that every grid point lies inside the bounds;
   - evaluates the residual and Jacobian once at the first point of each 1-D scan (with `--profile-grid`, at the first product point);
   - recomputes the optimizer scales;
   - reports the number of points, whether the output exists, and the time estimate under the assumption above, both serial and with `--workers`.

## Choices (labelled, not results)

- **No basin switch, no `zero_control` run, no quadrature ladder.** A profile point is a refit with a constraint, not a new candidate endpoint. So the acceptance rule of the plan, which needs a certified quadrature ladder, is not applied, and the ledger drops `acceptance_rule`, `basin_switch_rule` and `control_rule`. The adequacy verdict at $(4,128,128)$ is still recorded for each point.
- **Degrees of freedom.** These are the retained rows minus the free coordinates: $75-34=41$ for a 1-D scan.
- **Optimizer scales.** These are as in restart mode. For the fixed $r$ the scale is unused.
- **Start.** Every point starts from the same restart endpoint, not from the endpoint of the neighbouring grid point.

## Dry-run result (numerical; one evaluation, not a fit)

The command above was run on 2026-11-10 at orders $(2,16,16)$ with four-bin blocks. It reported:

- `optimizer_called: false`;
- the source-hash gate passed, with 42 reference files checked;
- the restart source is run 2 of `scan_point`, as in the restart record;
- 5 points, 34 free coordinates and 41 dof;
- the added optimizer scales $9.763\times10^{-4}$, $5.825\times10^{-2}$ and $1.822\times10^{-3}$ GeV$^{-2}$;
- an output that does not yet exist;
- a time estimate of $0.75$ h serial, or $0.45$ h with 2 workers.

The first point fixes $r_1=-1$ GeV$^{-2}$, with $r_2, r_3$ and the 32 baseline coordinates at the restart endpoint, where $r_1=-0.31953$ GeV$^{-2}$. The $\chi^2$ at $(2,16,16)$ per block is:

| block (rows) | $\chi^2$ |
|---|---|
| odd (12) | 5538747.1 |
| masses and leptonic widths (7) | 0.792 |
| shared widths (7) | 12840.8 |
| $\phi\to K\bar K$ (2) | 441.07 |
| BaBar (45) | 252.80 |
| normalization constraints (2) | 0.162 |
| total | 5552282.7 |

The residual and Jacobian are finite. A shift of $r_1$ by $-0.68$ GeV$^{-2}$ from the restart endpoint raises $\chi^2$ from $98.0$ to $5.6\times10^6$, mostly in the odd block. This agrees with the small $r_1$ optimizer scale. A grid of spacing $0.5$ GeV$^{-2}$ in $r_1$ is therefore much wider than the region where $\chi^2$ changes by a few units. The grid choice is left open.

The dry runs without `--profile-r` were also run on the same date: plain `--dry-run`, and `--dry-run --restart-from data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json`. Both passed the gate and printed the default `optimizer_settings`, identical to the source ledger. The restart evaluation gave $\chi^2=98.0028256648142$ and `optimizer_scales_equal_source_ledger: true`. No file was added or changed in the repository by the three dry runs (checked with `git status`).

## Verification

`scripts.test_fit_joint_vector_isospin_breaking_vpp` has 51 tests, run on 2026-11-10, all passing. They use the cheap orders $(1,3,4)$, except where the driver's own settings apply, and the optimizer is always mocked. The new tests cover:

- default tolerances giving an `optimizer_settings` identical to `OPTIMIZER_SETTINGS`, including key order;
- parsing of custom tolerances;
- rejection of zero, negative, nan, inf and non-numeric tolerances;
- tolerances reaching the `least_squares` call;
- a mocked restart recording the tolerances in the ledger;
- profile parsing of one and of several scans;
- rejection of a bad INDEX or format, values outside $[-5,5]$, nan, inf, repeated values, a duplicate INDEX, a profile without `--restart-from`, and `--profile-grid` without `--profile-r`;
- the 1-D point list, the Cartesian product and the first points;
- the substitution of the fixed value leaving every other coordinate unchanged;
- refusal of a point that fixes several $r$ without `--profile-grid`;
- refusal to overwrite an existing output;
- the time-estimate arithmetic;
- a profile point with a mocked optimizer: 34 free coordinates, the recorded keys, high-order dof 41, and block $\chi^2$ summing to the total;
- `main` with `--profile-r` on a mocked point function: the ledger schema, the printed point count, refusal of a rerun without `--overwrite`, and reuse of finished parts with `--overwrite`;
- the profile dry run, 1-D and Cartesian, never calling `least_squares`, the point functions, the process pool or the ledger writer, with one evaluation per scan;
- `main --dry-run --profile-r` writing nothing.

Not checked:

- any profile refit;
- the 9-minute assumption;
- whether the grid values are suitable.
