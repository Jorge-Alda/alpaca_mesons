---
title: Scalar evaluation optimization for photon pseudoscalar-pair decays
status: numerical
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: []
provenance: numerical-computation
---

# Completed scalar evaluation optimization for $a\to\gamma PP$

## Approval and completion

The optimized scalar evaluation is approved for use as the standard implementation of $a\to\gamma PP$ in `alpaca_chiral`. The former pointwise evaluator and width routine have been removed; the public `photon_two_pseudoscalar_amplitude` and `photon_two_pseudoscalar_width` now use the optimized algebra. This completion applies to $a\to\gamma PP$ only.

## Assumptions and scope

This comparison changes only evaluation order and scalar arithmetic for the existing odd vector-exchange amplitude. It retains the same flavour preparation, pole masses and constant widths, Dalitz domain, adaptive quadrature tolerances, pole subdivisions, and identical-particle factors. It does not establish an improvement for other three-body channels.

## Inputs and provenance

- The production amplitude preparation and optimized public functions are in `alpaca_chiral/photon_two_pseudoscalar.py`.
- Deterministic diagonal and Hermitian $ds$ coupling matrices are defined in `scripts/benchmark_photon_pp_scalar.py`.
- No external literature inputs were consulted for this numerical comparison; fitted model inputs are read by the existing amplitude constructor.

## Method

The amplitude constructor prepares complex poles $m_V^2-i m_V\Gamma_V$, groups terms by $s,t,u$, and stores mass combinations used in the Dalitz limits. At fixed $s$, `at_fixed_s(s)` evaluates the $s$-channel contribution once and returns a scalar function of $t$. The width computes the physical $t$ boundaries once per outer integration point and uses scalar arithmetic for the photon polarization sum. It retains the established adaptive tolerances and real pole subdivisions.

The executable verification and timing procedure is:

```sh
.venv/bin/python -m unittest scripts.benchmark_photon_pp_scalar
.venv/bin/python -m scripts.benchmark_photon_pp_scalar
```

The script compares the production amplitude and polarization sums with a slow test-only scalar reference at fixed interior, endpoint, and pole-region Dalitz points; compares integrated widths; checks exact-zero cases; and measures medians over three runs. Width timings include amplitude construction. The reference evaluator and width quadrature use the production prepared flavour terms, so these checks validate scalar evaluation and integration, not an independent reconstruction of the flavour algebra.

## Results

All three tests pass. In the post-consolidation run across five selected channel and mass points, the largest relative width difference was $1.11\times10^{-15}$. Measured full-width speedups ranged from $2.24$ to $7.92$, while repeated fixed-$s$ integrand evaluations sped up by $4.05$ to $6.03$. These timings are environment- and channel-dependent and are not a universal performance guarantee. The benchmark set includes $\pi^+\pi^-$ at $m_a=0.6$ and $1.5$ GeV, $K^+K^-$, $K_LK_S$, and $K^+\pi^-$ with a complex $ds$ source at $m_a=1.5$ GeV.

## Verification and approval basis

Command: `.venv/bin/python -m scripts.benchmark_photon_pp_scalar`.

Result: pointwise comparisons passed at relative tolerance $2\times10^{-13}$ for amplitudes and $2\times10^{-12}$ for polarization sums; width comparisons passed at relative tolerance $2\times10^{-10}$. The optimized implementation is now the sole production path for this decay class. The former scalar evaluator remains only as an independent validation calculation in the benchmark script.
