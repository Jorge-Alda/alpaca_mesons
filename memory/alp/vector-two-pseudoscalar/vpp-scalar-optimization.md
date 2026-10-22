---
title: Scalar evaluation optimization for vector-pseudoscalar-pseudoscalar decays
status: numerical
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: []
provenance: numerical-computation
---

# Scalar evaluation optimization for $a\to VPP$

## Scope and assumptions

The scalar-prepared amplitude and width are complete and approved for production use. They are now the canonical `vector_two_pseudoscalar_amplitude` and `vector_two_pseudoscalar_width` functions. The former slow production evaluators have been removed; independent slow references remain only in the benchmark script. The production implementation retains the coherent odd and even coefficients, finite-width poles, Dalitz boundaries, adaptive tolerances, pole subdivisions, and identical-pseudoscalar factor. No phase-space or propagator approximation is introduced.

## Inputs and provenance

- The reference amplitude terms and polarization sum are defined by `vector_two_pseudoscalar_amplitude` in `alpaca_chiral/vector_two_pseudoscalar.py`.
- The deterministic cases and comparison code are in `scripts/benchmark_three_body_scalar.py`.
- No external literature inputs were consulted for the optimization; vector, pseudoscalar, and fit inputs are supplied by the existing model code.

## Scalar polarization sum

At fixed $s=(p_1+p_2)^2$, let $t=(q_V+p_2)^2$ and let $t_-$ and $t_+$ be its physical limits. The Gram determinant in the odd polarization contraction is

$$\det G=\frac{s(t-t_-)(t_+-t)}{4}.$$

The optimized evaluator uses this full $t$-dependent expression. In particular, the angular zeros at both Dalitz boundaries are retained. The remaining transverse contractions are evaluated directly as scalars. It prepares squared masses and grouped complex poles once, computes the $s$-channel coefficients once per outer integration point, and reuses them across the inner quadrature.

`test_symbolic_gram` verifies the determinant identity by comparing it with the determinant of the reference $3\times3$ Gram matrix. Pointwise tests compare the full polarization sum including odd-even interference.

## Verification and benchmark

Run:

```sh
.venv/bin/python -m unittest scripts.benchmark_three_body_scalar
.venv/bin/python -m scripts.benchmark_three_body_scalar
```

The combined suite checks pointwise amplitudes and squared amplitudes, complete widths, identical-mass synthetic coefficients, thresholds, forbidden channels, and zero couplings. In the post-promotion 11-case benchmark, the five $VPP$ cases had maximum relative width difference $1.12\times10^{-15}$. Full-width speedups ranged from $4.52$ to $10.39$, and repeated integrand-evaluation speedups ranged from $5.87$ to $8.27$. Timings are three-run medians and include amplitude preparation in complete widths; they are environment- and case-dependent.

## Status

The optimization is complete and approved for production use. The canonical public functions use the scalar-prepared implementation; the old production evaluators and separate `_optimized` names have been removed.
