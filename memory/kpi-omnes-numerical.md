---
title: Numerical elastic K pi Omnes functions
status: mixed
created: 2026-09-23
updated: 2026-09-23
bibtex_keys:
  - Omnes:1958hv
  - Watson:1952ji
  - Pelaez:2016tgi
  - Pelaez:2020gnd
  - Jamin:2001zq
provenance: numerical-computation
---

# Numerical Elastic $K\pi$ Omnes Functions

## Implementation

[`scripts/implementation/omnes.py`](../scripts/implementation/omnes.py) evaluates the normalized once-subtracted scalar function

$$
\Omega(s;s_0)=\exp\left\{
\frac{s-s_0}{\pi}\int_{s_{\mathrm{th}}}^{\infty}
\frac{\delta(s')\,ds'}{(s'-s_0)(s'-s-i0)}
\right\},
\qquad \Omega(s_0;s_0)=1,
$$

using a caller-supplied phase in radians [@Omnes:1958hv]. On the upper rim of the elastic cut, the principal-value real part is combined with the exact imaginary exponent $i\delta(s)$, so the phase of $\Omega$ equals the elastic scattering phase [@Watson:1952ji].

The phase input has a finite endpoint. The caller must explicitly choose either a constant continuation, matched continuously at that endpoint, or a truncated phase set to zero above it. Both are assumptions. For a constant continuation, the asymptotic part of the dispersion integral is evaluated analytically; the finite remainder is integrated with SciPy Cauchy principal-value quadrature. Sensitivity to this continuation should be varied as a systematic uncertainty.

`PairOmnesFSI` dispatches separate scalar functions by physical $K\pi$ pair and can be passed directly as the `fsi` argument of the $a\to K\bar K\pi$ width. This remains a diagonal independent-band approximation, not a charge-basis or inelastic Muskhelishvili--Omnes matrix [@Jamin:2001zq].

## Literature phase implementation

[`scripts/implementation/kpi_phases.py`](../scripts/implementation/kpi_phases.py) implements the central elastic $I=1/2$ and $I=3/2$ scalar conformal fits of Pelaez and Rodas, with selectable constrained or unconstrained parameters [@Pelaez:2016tgi]. It uses the isospin-limit masses of that analysis by default. For the $I=1/2$ Adler zero, the implementation uses the $4\Delta_{K\pi}^2$ expression in the later review, which reproduces the quoted $s_A\simeq0.236\,\mathrm{GeV}^2$; the 2016 printed expression omits the factor four and is numerically inconsistent with its quoted value [@Pelaez:2020gnd]. The $I=1/2$ phase is restricted to the $K\eta$ threshold, where the source switches to an inelastic representation; the $I=3/2$ phase is restricted to $1.8\,\mathrm{GeV}$.

Freezing the $I=1/2$ phase above the $K\eta$ threshold produces a mathematically defined scalar Omnes function but is an explicit approximation. It does not implement the $K_0^*(1430)$ region as a controlled coupled-channel calculation. Applying one isospin-symmetric phase to charge channels with physical split thresholds is likewise an approximation; a complete treatment requires the charge-basis Omnes matrix and consistent scattering input.

## Verification

Run

```bash
.venv/bin/python -m unittest scripts.implementation.test_omnes scripts.implementation.test_kpi_phases
```

The tests verify subtraction normalization, the exact constant-phase solution, reality below threshold, the phase on the upper cut, explicit tail behavior, invalid-input rejection, pair dispatch, phase-domain restrictions, phase signs, and composition of the literature phase with the Omnes integrator.
