---
title: Scalar evaluation optimization for three-pseudoscalar decays
status: numerical
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: []
provenance: numerical-computation
---

# Scalar evaluation optimization for $a\to3P$

## Scope and assumptions

The scalar-prepared amplitude and width are complete and approved for production use. They are now the canonical `three_pseudoscalar_amplitude` and `three_pseudoscalar_width` functions. The former slow production evaluators have been removed; independent slow references remain only in the benchmark script. The production implementation retains all labelled permutations, flavour coefficients, finite vector widths, matching contributions, Dalitz boundaries, adaptive tolerances, and identical-particle factors. No zero-width approximation is introduced.

## Inputs and provenance

- Flavour traces, poles, and the reference amplitude are prepared by `three_pseudoscalar_amplitude` in `alpaca_chiral/three_pseudoscalar.py`.
- The benchmark inputs and independent pointwise/width comparisons are in `scripts/benchmark_three_body_scalar.py`.
- No external literature inputs were consulted for the optimization. The model inputs remain those read by the existing module.

## Algebraic reduction

For a chosen spectator $k$ and pair $ij$, write $s_{ij}=(p_i+p_j)^2$ and $s_{ik},s_{jk}$ for the crossed invariants. Use the on-shell identity

$$s_{ij}+s_{ik}+s_{jk}=m_a^2+m_i^2+m_j^2+m_k^2.$$

The matching polynomial in the reference expression reduces to

$$L_{k;ij}=\frac{s_{ij}(s_{ik}-s_{jk})+(m_j^2-m_i^2)(m_a^2-m_k^2)}{4}.$$

The pole numerator is $N_{k;ij}=s_{ij}L_{k;ij}$. Therefore the pole and matching terms combine exactly:

$$\frac{N_{k;ij}}{M_V^2-s_{ij}-iM_V\Gamma_V}+L_{k;ij}
=L_{k;ij}\frac{M_V^2-iM_V\Gamma_V}{M_V^2-iM_V\Gamma_V-s_{ij}}.$$

This identity is valid with the retained constant finite width. `test_symbolic_three_pseudoscalar` proves the polynomial identities with SymPy. The implementation also prepares the chiral polynomial, squared masses, and grouped vector residues once per amplitude rather than rebuilding an invariant dictionary and mass combinations at every Dalitz point.

## Verification and benchmark

Run:

```sh
.venv/bin/python -m unittest scripts.benchmark_three_body_scalar
.venv/bin/python -m scripts.benchmark_three_body_scalar
```

The suite checks symbolic identities, pointwise agreement at interior, near-boundary, and pole-region points, integrated widths, thresholds, and zero couplings. In the post-promotion 11-case benchmark, the six $3P$ cases had maximum relative width difference $6.67\times10^{-16}$. Full-width speedups ranged from $5.72$ to $9.26$, and repeated integrand-evaluation speedups ranged from $6.33$ to $12.58$. Timings are three-run medians and include amplitude preparation in complete widths; they are environment- and case-dependent.

The $\pi^0K_SK_L$ benchmark returns a reference residual of approximately $1.26\times10^{-43}$ GeV. It is reproduced by the optimized path and is recorded as a numerical cancellation residual, not a physical nonzero prediction.

## Status

The optimization is complete and approved for production use. The canonical public functions use the scalar-prepared implementation; the old production evaluators and separate `_optimized` names have been removed.
