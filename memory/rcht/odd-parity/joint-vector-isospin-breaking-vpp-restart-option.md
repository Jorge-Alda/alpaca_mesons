---
title: Restart and evaluation-limit options of the isospin-breaking VPP joint fit driver (dry run only)
status: numerical
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh]
provenance: numerical-computation
---

# Restart and evaluation-limit options of the isospin-breaking $VPP$ fit driver

**No fit was run.** This record adds two command-line options to the driver of [joint-vector-isospin-breaking-vpp-fit-plan.md](joint-vector-isospin-breaking-vpp-fit-plan.md). It reports one dry-run evaluation at the restart point.

Driver: [`scripts/fit_joint_vector_isospin_breaking_vpp.py`](../../../scripts/fit_joint_vector_isospin_breaking_vpp.py); tests: [`scripts/test_fit_joint_vector_isospin_breaking_vpp.py`](../../../scripts/test_fit_joint_vector_isospin_breaking_vpp.py).

```bash
.venv/bin/python scripts/fit_joint_vector_isospin_breaking_vpp.py --dry-run --restart-from data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json
# fit (not run): drop --dry-run; optional --max-nfev N and --output PATH
```

## Conventions and units

These are the same as in the fit plan: GeV and nb, with $r_i$ in GeV$^{-2}$. Residuals and $\chi^2$ are dimensionless. Conventions follow Alda:2026nrd [@Alda:2026nrd].

## Inputs

- Source ledger `data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json`, produced by the plan's driver. Entry `scan_point` has four unsuccessful runs that stopped at `max_nfev=5000`. Its starts, gates and bounds are given in the [fit plan](joint-vector-isospin-breaking-vpp-fit-plan.md).
- Likelihood option: [joint-vector-isospin-breaking-vpp-likelihood-option.md](joint-vector-isospin-breaking-vpp-likelihood-option.md). Since this change, the driver evaluates it through the bitwise-equal class of [joint-vector-isospin-breaking-vpp-calculation-reuse.md](joint-vector-isospin-breaking-vpp-calculation-reuse.md).
- The data are the PDG rows [@ParticleDataGroup:2024cfk] and BaBar $KK\pi$ [@BaBar:2007ceh], entering through the baseline likelihood.

## Options

1. `--restart-from PATH`. The driver reads the ledger PATH. It selects the run of entry `scan_point` with the lowest $(4,128,128)$ $\chi^2$, among finished runs with a valid high-order evaluation, whether accepted or not. That run's 35-coordinate endpoint becomes the single start, with $r$ free and the bounds of the plan. The start definition records:
   - the ledger path and its SHA-256;
   - the entry, the run index and the run's start endpoint;
   - the run's high-order and fit-order $\chi^2$, its acceptance flag, the optimizer message and the cumulative `nfev`.

   The driver refuses to start when:
   - a coordinate lies outside the bounds;
   - the residual at the fit order is nonfinite or equals the penalty value;
   - the ledger's coordinate names differ from the 35 names of the model.
2. `--max-nfev N`. The default is 5000, unchanged. The value replaces `max_nfev` in the optimizer settings, which are written to the ledger as `optimizer_settings`. This applies to both the original two-start fit and the restart.
3. `--dry-run` works with both options. With `--restart-from`, it does the following:
   - checks the model and the source-hash gate;
   - builds the restart start;
   - evaluates the residual and Jacobian once at the restart point and prints the $\chi^2$ per block;
   - recomputes the optimizer scales.

   It never calls `least_squares`.

## Choices (labelled, not results)

- **Single start, no basin switch, no `zero_control` run.** A restart continues one run. The basin-switch rule of the plan (refit from the other accepted endpoints when the scaled distance exceeds 1) belongs to the original multistart. The selected run lies at scaled distance $268.9$ from the best baseline endpoint, so applying that rule would not be a restart.
- **Optimizer scales.** These follow the definition of the plan: the 32 scales of the endpoint comparison, and reciprocal Jacobian column norms of $r_1,r_2,r_3$ at the `zero_control` start. The dry run checks that they are bitwise equal to the scales stored in the source ledger.
- **Output.** The default output in restart mode is `data/vector_sector/joint_vector_isospin_breaking_vpp_restart_fit.json`. An output equal to the source ledger is rejected. Checkpoints use the same per-run progress files as the plan.

## Dry-run result (numerical; one evaluation, not a fit)

The command above was run on 2026-11-10 at orders $(2,16,16)$ with four-bin blocks. The source-hash gate passed, with 42 reference files checked.

Restart source: run index 2 of `scan_point`, which started from `continuation_4128128/stored_negative_c4`. It was not accepted, and it stopped with "The maximum number of function evaluations is exceeded." after 5000 evaluations. Its $(4,128,128)$ $\chi^2$ is $98.1188$. The ledger SHA-256 is `f78624e6…906a`.

Restart $r=(-0.31953,\,1.58668,\,0.35146)$ GeV$^{-2}$. The table gives the $\chi^2$ at $(2,16,16)$, with 75 rows, 35 free coordinates and 40 dof:

| block (rows) | $\chi^2$ |
|---|---|
| odd (12) | 18.901 |
| masses and leptonic widths (7) | 0.792 |
| shared widths (7) | 10.537 |
| $\phi\to K\bar K$ (2) | 2.520 |
| BaBar (45) | 65.091 |
| normalization constraints (2) | 0.162 |
| total | 98.0028256648142 |

The total equals the `fit_order_chi2` stored for run 2. The residual and Jacobian are finite. The added optimizer scales are $9.763\times10^{-4}$, $5.825\times10^{-2}$ and $1.822\times10^{-3}$ GeV$^{-2}$; `optimizer_scales_equal_source_ledger` is true. The printed output was the same, apart from the elapsed time, before and after the driver switched to the reuse class.

## Verification

`scripts.test_fit_joint_vector_isospin_breaking_vpp` has 31 tests, run on 2026-11-10, all passing. They use the cheap orders $(1,3,4)$, except where the driver's own settings apply, and the optimizer is always mocked. The 12 new tests cover:

- argument defaults and parsing of `--restart-from` and `--max-nfev`;
- rejection of `--max-nfev 0` or non-integers, and of an output equal to the source ledger;
- run selection by lowest high-order $\chi^2$, skipping error runs and invalid high-order runs;
- selection of run 2 in the real ledger;
- the restart start being the stored endpoint, with $r$ free and the source path and run index recorded;
- refusal of endpoints outside the bounds (an $r$ coordinate, $G_V$, an even-sector coordinate);
- refusal of nonfinite and penalty residuals;
- refusal of different coordinate names;
- `dry_run` with a restart never calling `least_squares` (driver, sensitivity and multistart modules), the run functions or the process pool;
- `main --dry-run --restart-from` writing nothing;
- `main --restart-from --max-nfev 1500`, with a mocked run, writing the start and `optimizer_settings.max_nfev=1500` to the ledger.

Not checked: the restart fit itself.
