---
title: Gauss-Legendre quadrature candidate for three-pseudoscalar ALP decays
status: numerical
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: []
provenance: numerical-computation
---

# Implemented Gauss-Legendre quadrature for $a\to3P$

## Scope and inputs

The official `three_pseudoscalar_width` uses the Gauss-Legendre method. The former nested-adaptive production function has been removed; its independent implementation remains in `scripts/benchmark_three_body_quadrature.py` as a comparison reference. The official function preserves all labelled permutations, finite vector widths, coherent terms, physical boundaries, and factorial symmetry factors. No literature inputs or physical amplitude terms were changed.

The amplitude and masses are prepared by `three_pseudoscalar_amplitude` and `_pseudoscalars()` in `alpaca_chiral/three_pseudoscalar.py`. Its callable retains its scalar behavior and now exposes array evaluation and invariant pole metadata for the fixed rule. Deterministic diagonal and Hermitian $ds$ couplings are defined in `scripts/benchmark_three_body_quadrature.py`.

## Numerical method

The full two-dimensional Dalitz region is integrated; no symmetry reduction is assumed. The inner and outer dimensions use vectorized Gauss--Legendre rules, with a sine-squared outer map. Intervals are split at accessible real pole centers and at offsets $\pm2m_V\Gamma_V$ and $\pm8m_V\Gamma_V$. Orders 16, 32, 64, and 128 are compared. If inner or outer order doubling does not meet the requested relative tolerance, the corresponding integration falls back to adaptive quadrature. The exact coherent squared amplitude is evaluated before integration; no narrow-width or propagator approximation is made.

The public function takes the same six positional inputs as before, with keyword-only `order=16` and `epsrel=10^{-5}`.

## Verification

Command: `.venv/bin/python -m unittest scripts.benchmark_three_body_quadrature scripts.benchmark_three_body_scalar`

All eleven tests pass, including the seven previous scalar-regression tests. The new tests cover eight physical $3P$ cases within a fourteen-case combined $3P/VPP$ set, an independently implemented nested adaptive reference at $2\times10^{-7}$ with pole splits, order 16 versus 32, threshold/zero cases, permutations, and forced fallback. The largest difference from the original default adaptive width is $6.20\times10^{-7}$ over the combined set; the largest difference from the tightened reference is $5.35\times10^{-9}$.

## Benchmark

Command: `.venv/bin/python -m scripts.benchmark_three_body_quadrature --repeats 7`

Seven-run medians include amplitude preparation. The eight $3P$ cases span speedups of $0.91$ to $11.65$. The largest gain is $11.65$ for $\pi^0K^+K^-$ at $m_a=1.7$ GeV. The $\pi^+\pi^-K_S$ cases at $m_a=1.3$ and $1.5$ GeV have speedups $1.11$ and $0.91$, respectively. Other benchmark cases include $\pi^0\pi^+\pi^-$ near threshold and at $0.6$ and $1.5$ GeV, $\pi^0\pi^0\eta$, and $\pi^0K_SK_L$.

The case-level medians, widths, and relative differences are emitted as JSON by the benchmark script. Results depend on the channel, mass, and runtime environment. The Gauss-Legendre calculation is now the canonical public width; the adaptive implementation remains only in the benchmark reference.
