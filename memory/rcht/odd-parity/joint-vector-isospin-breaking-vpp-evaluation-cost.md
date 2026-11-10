---
title: Per-evaluation cost of the isospin-breaking VPP likelihood option against the baseline
status: numerical
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd]
provenance: numerical-computation
---

# Per-evaluation cost of the isospin-breaking $VPP$ option

This record contains timing and call counts only. No fit was run and the likelihood numerics are unchanged. It explains the run time of the `scan_point` refits of [joint-vector-isospin-breaking-vpp-fit-plan.md](joint-vector-isospin-breaking-vpp-fit-plan.md). That entry took $32438$ s in total, of which $30509$ s were the four refits. It also estimates the cost of the restart of [joint-vector-isospin-breaking-vpp-restart-option.md](joint-vector-isospin-breaking-vpp-restart-option.md).

Script: [`scripts/benchmark_isospin_breaking_vpp_cost.py`](../../../scripts/benchmark_isospin_breaking_vpp_cost.py); tests: [`scripts/test_benchmark_isospin_breaking_vpp_cost.py`](../../../scripts/test_benchmark_isospin_breaking_vpp_cost.py); output: [`data/vector_sector/isospin_breaking_vpp_cost_benchmark.json`](../../../data/vector_sector/isospin_breaking_vpp_cost_benchmark.json).

```bash
OMP_NUM_THREADS=1 .venv/bin/python scripts/benchmark_isospin_breaking_vpp_cost.py --output data/vector_sector/isospin_breaking_vpp_cost_benchmark.json
```

## Setup

- **Point.** The endpoint of run 2 of `scan_point`, the restart point, with $r=(-0.3195,1.5867,0.3515)$ GeV$^{-2}$.
- **Quadrature.** $(2,16,16)$ with four-bin blocks, giving 6 blocks per evaluation. `OMP_NUM_THREADS=1`. The machine and software versions are stored in the output file.
- **Configurations.**
  - Baseline: `GlobalLikelihood` with the first 32 coordinates, i.e. the optimized baseline path of [kkpi-global-fit-channel-jacobian.md](kkpi-global-fit-channel-jacobian.md).
  - Option at $r=0$ and option at the restart $r$: `IsospinBreakingVppLikelihood` with the default switches (the class of [joint-vector-isospin-breaking-vpp-likelihood-option.md](joint-vector-isospin-breaking-vpp-likelihood-option.md)). This is the class used for the original fit.
  - Reuse at $r=0$ and reuse at the restart $r$: the driver's `FitLikelihood`, built on [joint-vector-isospin-breaking-vpp-calculation-reuse.md](joint-vector-isospin-breaking-vpp-calculation-reuse.md).
- **Timing method.** This follows [kkpi-global-fit-channel-jacobian.md](kkpi-global-fit-channel-jacobian.md). There is one warm-up residual and Jacobian per configuration. Then come 5 rounds in which the configurations alternate, with residual rounds and Jacobian rounds as separate loops. The table gives medians of the wall time.
- **Counts.** One instrumented residual and one instrumented Jacobian per configuration, made in a separate pass. Wrappers count each call and accumulate its wall time; these accumulated times include the wrapper overhead. The test checks that the wrappers leave the residual and Jacobian bitwise unchanged.

## Timing (median of 5, seconds; ratio to baseline)

| configuration | residual | ratio | Jacobian | ratio |
|---|---|---|---|---|
| baseline (32 coordinates) | 0.0429 | 1.000 | 0.767 | 1.000 |
| option, $r=0$ | 0.0432 | 1.007 | 1.030 | 1.343 |
| option, restart $r$ | 0.0725 | 1.688 | 1.518 | 1.979 |
| reuse, $r=0$ | 0.0432 | 1.006 | 0.906 | 1.181 |
| reuse, restart $r$ | 0.0454 | 1.058 | 0.939 | 1.224 |

The spread of the 5 samples is below 1.5% in every case. A first benchmark run, made before the reuse configurations were added, gave the same values within 0.5%.

## Call counts per evaluation

The table counts calls per Jacobian; the counts for one residual are given after it. "Zero-excitation" calls are the second running-propagator calls with $\beta=(0,0)$ made by the option.

| per Jacobian | baseline | option $r=0$ | option restart | reuse restart |
|---|---|---|---|---|
| full residuals inside the inherited Jacobian | 11 | 11 | 11 | 11 |
| full residuals added by the option's Jacobian | 0 | 4 | 4 | 3 |
| channel-restricted BaBar evaluations | 6 | 6 | 6 | 6 |
| BaBar derivative pass (odd coordinates) | 1 | 1 | 1 | 1 |
| block integrations (`_babar_blocks_full`) | 114 | 138 | 138 | 132 |
| `old_blocks` | 14 | 18 | 18 | 17 |
| running propagator, $K^*$ / $\rho$ | 360 / 360 | 456 / 456 | 456 / 456 | 432 / 432 |
| zero-excitation running propagator, $K^*$ / $\rho$ | 0 / 0 | 72 / 36 | 456 / 228 | 432 / 216 (reused, not evaluated) |
| `reference_propagator_pair` | 1440 | 1932 | 2508 | 1728 |
| reference propagator (pole function) | 2160 | 2952 | 4104 | 2592 |
| `vertex_factors` / `s_channel_factors` | 0 / 0 | 24 / 18 | 182 / 138 | 13 / 8 |
| `s_channel_light_correction` | 0 | 18 | 108 | 102 |

For one residual at the restart point, the option makes 24 $K^*$ and 24 $\rho$ running-propagator calls, the same as the baseline. It adds 24 zero-excitation $K^*$ and 12 zero-excitation $\rho$ calls ($s$ channel), as well as 8 `vertex_factors`, 6 `s_channel_factors` and 6 light corrections. The 4 extra full residuals of the option's Jacobian are the duplicated reference residual and the three forward differences in $r_1,r_2,r_3$. At $r=0$ the option uses the baseline path, so its only zero-excitation calls come from the three residuals with $r$ shifted by one step.

## Sources of the extra cost (option at the restart $r$)

The extra time per Jacobian over the baseline is $1.518-0.767=0.751$ s. It is split below using the accumulated wrapper times; the split is approximate because those times include the wrapper overhead.

| source | s per Jacobian | share of the extra |
|---|---|---|
| zero-excitation $K^*$ running propagator in `mode_scaled` (456 calls) | 0.450 | 60% |
| zero-excitation $s$-channel $\rho$ running propagator (228 calls) | 0.025 | 3% |
| three forward-difference full residuals for $r_1,r_2,r_3$ (at baseline residual cost) | 0.129 | 17% |
| duplicate `self.residual(values)` in the option's `jacobian` (at baseline residual cost) | 0.043 | 6% |
| `vertex_factors` and `s_channel_factors` recomputed for every block (320 calls) | 0.089 | 12% |
| $\omega$ light correction in `light_scaled` (108 calls) | 0.006 | 1% |
| not attributed (identity lookups, cache resets, wrapper overhead) | 0.010 | 1% |

Each zero-excitation $K^*$ call costs about 1 ms, about as much as the full call it follows, because it recomputes the reference propagator pair and the on-shell pole. This one cost explains why the option's residual is 1.69 times the baseline at nonzero $r$ and 1.00 at $r=0$.

## Does the per-call cost explain the $\sim9$ h `scan_point` entry?

Assumption (labelled): one optimizer run costs $\mathrm{nfev}\times t_\mathrm{residual}+\mathrm{njev}\times t_\mathrm{Jacobian}$, with the costs measured at the restart point. The trf linear algebra, the callback and the checkpoint writes are neglected.

With the option's costs ($0.0725$ s and $1.518$ s), the four recorded runs ($\mathrm{nfev}=5000$; $\mathrm{njev}=4812, 4857, 4774, 4800$) are predicted to take $7667, 7736, 7610, 7649$ s. They took $8092, 7594, 7399, 7424$ s, which is within $-5\%$ to $+3\%$ of the prediction. So yes: the per-call cost of the option, with about 1 Jacobian per evaluation, accounts for the $30509$ s of refits. The remaining $\approx1930$ s of the entry are the $(4,128,128)$ re-evaluation and local errors of the four runs, which were not timed here. The Jacobians take about 95% of the refit time.

## Wall-time estimate for the restart run

Under the same assumption, with `max_nfev` $=5000$ and $\mathrm{njev}=4774$ (the ratio of run 2):

- with the option class as used for the original fit: $7610$ s, about 2.1 h;
- with the driver's current reuse class ($0.0454$ s per residual, $0.939$ s per Jacobian): $4708$ s, about 1.3 h.

These are upper estimates for a run that reaches the limit. A run that converges earlier stops sooner. The $(4,128,128)$ evaluation and, if the run succeeds, the quadrature ladder come on top.

## Verification

`scripts.test_benchmark_isospin_breaking_vpp_cost` has 6 tests at orders $(1,3,4)$, run on 2026-11-10, all passing. They check:

- the restart run index is 2;
- full-residual counts per Jacobian are 11+0 (baseline), 11+4 (option) and 11+3 (reuse), with 6 channel-restricted evaluations and 1 derivative pass;
- zero-excitation $K^*$ calls match the $K^*$ calls one to one at nonzero $r$ and are absent at $r=0$;
- the wrappers preserve the residual and Jacobian bitwise and restore all patched functions;
- the reuse class is bitwise equal to the option;
- the arithmetic of the wall-time formula.

The benchmark output also records `instrument_preserves_outputs` and `reuse_bitwise_equal_to_option`, both true at $(2,16,16)$.
