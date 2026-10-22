---
title: Gauss-Legendre quadrature candidate for vector-pseudoscalar-pseudoscalar ALP decays
status: numerical
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: []
provenance: numerical-computation
---

# Implemented Gauss-Legendre quadrature for $a\to VPP$

## Scope and inputs

The official `vector_two_pseudoscalar_width` uses the Gauss-Legendre method. The former nested-adaptive production function has been removed; its independent implementation remains in `scripts/benchmark_three_body_quadrature.py` as a comparison reference. The official function retains the full coherent odd/even amplitude, polarization sum, finite-width poles, physical Dalitz boundaries, and identical-pseudoscalar factor. No physical amplitude or propagator term was changed.

The amplitude, pole groups, and fixed-$s$ coefficients are prepared by `vector_two_pseudoscalar_amplitude` in `alpaca_chiral/vector_two_pseudoscalar.py`. The benchmark includes diagonal sources and a Hermitian physical $ds$ source, along with fitted hadronic inputs read by the existing module. No external literature inputs were consulted.

## Numerical method

The complete two-dimensional Dalitz domain is integrated with vectorized Gauss--Legendre rules and a sine-squared outer map. Accessible real pole centers and offsets $\pm2m_V\Gamma_V$ and $\pm8m_V\Gamma_V$ split the intervals; zero-width poles are not broadened or regularized. Orders 16, 32, 64, and 128 are compared, with adaptive fallback when order doubling does not meet tolerance. The full coherent polarization sum, including odd-even interference and its $t$-dependent Gram determinant, is evaluated at every node.

The public function takes the same six positional inputs as before, with keyword-only `order=16` and `epsrel=10^{-5}`.

## Verification

Command: `.venv/bin/python -m unittest scripts.benchmark_three_body_quadrature scripts.benchmark_three_body_scalar`

All eleven tests pass, including the seven previous scalar-regression tests. The new suite covers six physical $VPP$ cases in the fourteen-case combined set, a tightened independent adaptive reference at $2\times10^{-7}$ with pole splits, order 16 versus 32, exact zero and closed channels, exchange of pseudoscalars, and forced adaptive fallback. Across all fourteen cases, the maximum difference from the original default adaptive widths is $6.20\times10^{-7}$; the maximum difference from the tightened reference is $5.35\times10^{-9}$.

## Benchmark

Command: `.venv/bin/python -m scripts.benchmark_three_body_quadrature --repeats 7`

Seven-run medians include amplitude preparation. The six $VPP$ cases span speedups of $0.67$ to $4.02$. The largest gain is $4.02$ for $\omega\pi^0K_S$ at $m_a=1.8$ GeV. The $\omega\pi^+\pi^-$ case at $m_a=1.8$ GeV is slower, with a speedup of $0.67$. Other cases cover $\omega\pi^+\pi^-$ near threshold, $\phi K^+K^-$, and the unequal-mass $\omega K^+\pi^-$ channel.

The benchmark script emits case-level medians, widths, and relative differences as JSON. The measured speedup is case-dependent. The Gauss-Legendre calculation is now the canonical public width; the adaptive implementation remains only in the benchmark reference.
