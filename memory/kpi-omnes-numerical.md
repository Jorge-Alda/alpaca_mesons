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

The phase input has a finite endpoint. The caller must explicitly choose either a constant continuation, matched continuously at that endpoint, or a truncated phase set to zero above it. Both are assumptions. For a constant continuation, the asymptotic part of the dispersion integral is evaluated analytically. In the finite remainder, the substitution $s'=s_{\mathrm{th}}+u^2$ regularizes the threshold and an analytical principal-value subtraction removes the on-cut pole before SciPy quadrature. Sensitivity to the high-energy continuation should be varied as a systematic uncertainty.

`TabulatedOmnes` precomputes the direct function over a requested physical invariant interval. It uses a signed square-root coordinate at the two-body threshold, interpolates $\log|\Omega|$ with separate monotone cubic pieces below and above threshold, and retains the exact input phase on the cut. Two omitted points per interval are compared with direct dispersion-integral evaluations until the requested relative tolerance is met. The default usage tolerance is $10^{-6}$.

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

## Performance diagnosis and optimized width

The current direct implementation is unsuitable for repeated evaluation inside the nested Dalitz integral. On the present machine, 100 evaluations of the constrained Pelaez--Rodas $I=1/2$ Omnes function took approximately $0.72\,\mathrm{s}$, whereas 100,000 evaluations of the contact amplitude took approximately $0.063\,\mathrm{s}$. A contact-width integration at $m_a=1.5\,\mathrm{GeV}$ requested the FSI callback 16,758 times. Consequently, inserting the direct Omnes evaluator replaces each inexpensive amplitude evaluation by two adaptive principal-value dispersion integrals. The nesting is

$$
\text{outer Dalitz quadrature}
\longrightarrow
\text{inner Dalitz quadrature}
\longrightarrow
\Omega(s),\Omega(t)
\longrightarrow
\text{adaptive dispersion quadrature}.
$$

There is also avoidable repetition: for fixed outer $s$, the same $s$-band Omnes function is recomputed at every inner-$t$ node. `PairOmnesFSI` dispatches functions but does not cache their values.

The implemented optimization has three layers:

1. `TabulatedOmnes` precomputes each distinct $\Omega(s)$ over the full physical invariant range and validates the interpolation to the selected tolerance.
2. Exploit the separable band amplitude $\mathcal M(s,t)=A(s)+B(t)$. For fixed $s$,

$$
\int_{t_-}^{t_+}|A(s)+B(t)|^2dt
=(t_+-t_-)|A(s)|^2
+\int_{t_-}^{t_+}|B(t)|^2dt
+2\operatorname{Re}\left[A(s)\left(\int_{t_-}^{t_+}B(t)dt\right)^*\right].
$$

`s_wave_decay_width_fast` adaptively interpolates $B(t)$ with piecewise cubics. It integrates $B(t)$ and the resulting degree-six polynomial for $|B(t)|^2$ analytically interval by interval. The two-dimensional adaptive Dalitz integral then becomes one outer integral with constant-time evaluations of the two primitive differences.
3. Add exact memoization as a fallback for direct evaluation and hoist $A(s)$ outside the inner integral. This removes repeated $s$-band dispersion integrals but does not eliminate the many distinct $t$ evaluations, so it is a secondary optimization rather than the final solution.

For the documented $m_a=1.2\,\mathrm{GeV}$, $K^+K^-\pi^0$ example, a $10^{-6}$ Omnes table used 75 nodes and took approximately $2.9\,\mathrm{s}$ to construct. The optimized width then took approximately $3.2\,\mathrm{s}$. Its relative difference from the original nested width evaluated with the same tabulated Omnes input was $8\times10^{-9}$ when the internal width integration used one fifth of the requested $10^{-6}$ error budget. In the same benchmark, the original nested width with tabulated input took approximately $110\,\mathrm{s}$. Timings are machine-dependent numerical provenance; accuracy must be checked for each new mass range and phase model.
