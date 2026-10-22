---
title: Gauss-Legendre quadrature for photon pseudoscalar-pair widths
status: numerical
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: []
provenance: numerical-computation
---

# Implemented Gauss-Legendre quadrature for $a\to\gamma PP$

## Scope and assumptions

The Gauss-Legendre method is the official implementation of the odd vector-exchange width, exposed as `photon_two_pseudoscalar_width`. The former nested-adaptive routine has been removed from the production module; the benchmark script retains a standalone adaptive reference for numerical comparisons. The official method preserves the same coherent amplitude, finite constant pole widths, Dalitz boundaries, normalization, and identical-particle factor. No literature source or physical input was changed.

## Inputs and provenance

- The amplitude, physical Dalitz domain, and width normalization are those of [alp-photon-pp-width.md](alp-photon-pp-width.md) and the implementation in `alpaca_chiral/photon_two_pseudoscalar.py`.
- The benchmark uses deterministic diagonal and Hermitian $ds$ source matrices from `scripts/benchmark_photon_pp_scalar.py` and fitted hadronic inputs read through the existing amplitude constructor.
- No external references were consulted for this numerical comparison.

## Method

For equal-mass daughters, the new routine evaluates both Dalitz variables with vectorized Gauss-Legendre rules. A sine-squared map in $s$ smooths the kinematic endpoints. For unequal-mass daughters, it retains adaptive outer integration and applies converged Gauss-Legendre quadrature to the inner variable. The $t$ and crossed $u$ integrations are subdivided at real pole centers and at offsets of $\pm2m_V\Gamma_V$ and $\pm8m_V\Gamma_V$. The quadrature order starts at 16 and doubles through 32, 64, and 128; successive results are compared at relative tolerance `epsrel/4`. Unresolved inner integrations fall back to adaptive quadrature. For amplitudes with only an $s$-channel term, the inner integral is a quadratic polynomial and is evaluated exactly with a two-node Gauss-Legendre rule.

The public routine accepts the same five positional inputs as the former width and the keyword-only controls `order` (default 16) and `epsrel` (default $10^{-5}$). Exact-zero/closed channels return before quadrature construction.

## Verification

Command: `.venv/bin/python -m unittest scripts.benchmark_photon_pp_quadrature`

All four tests pass. Nine nonzero cases compare against the standalone nested-adaptive reference tightened to $2\times10^{-7}$ in each variable. The suite also checks order 16 versus 32, exact zeros, meson exchange for unequal masses, and the forced adaptive fallback at order 1. Before promotion, the largest relative difference between the new method and the former production width over the nine benchmark cases was $3.17\times10^{-9}$.

## Benchmark

Command: `.venv/bin/python -m scripts.benchmark_photon_pp_quadrature --repeats 7`

Widths include amplitude preparation; timings are seven-run medians on the current environment. The measured adaptive-to-new times and speedups were:

| $m_a$ (GeV) | Final state | Adaptive (ms) | Gauss-Legendre (ms) | Speedup |
|---:|:---|---:|---:|---:|
| 0.28 | $\pi^+\pi^-$ | 7.45 | 5.92 | 1.26 |
| 0.60 | $\pi^+\pi^-$ | 8.22 | 4.75 | 1.73 |
| 1.05 | $\pi^+\pi^-$ | 9.53 | 8.11 | 1.17 |
| 1.50 | $\pi^+\pi^-$ | 11.30 | 9.18 | 1.23 |
| 1.05 | $K^+K^-$ | 9.71 | 10.44 | 0.93 |
| 1.50 | $K^+K^-$ | 35.07 | 12.05 | 2.91 |
| 1.50 | $K_LK_S$ | 31.63 | 11.88 | 2.66 |
| 1.50 | $K^+\pi^-$ | 26.91 | 28.13 | 0.96 |
| 2.00 | $K^+\pi^-$ | 40.45 | 31.04 | 1.30 |

The gain depends on channel and mass: seven cases were faster, while the near-threshold $K^+K^-$ case and the $1.50$ GeV unequal-mass case were respectively 7% and 4% slower in this repeated run. These environment-specific timings do not support a universal speedup claim. Reproduce them with the benchmark script; the machine-readable run is generated on standard output.
