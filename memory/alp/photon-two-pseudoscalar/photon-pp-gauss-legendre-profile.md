---
title: Gauss-Legendre profile for charged kaon and pion radiative decays
status: numerical
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: []
provenance: numerical-computation
---

# Assumptions and inputs

This profile measures the production `alpaca_chiral` Gauss–Legendre width for
$a\to\gamma\pi^+\pi^-$ and $a\to\gamma K^+K^-$. Both use $m_a=1.5$ GeV,
$f_a=100$ GeV, the same deterministic ALP coupling object
`coupling_inputs()[0]` from `scripts/benchmark_photon_pp_quadrature.py`, and
the production defaults `order=16` and `epsrel=10^{-5}`. The amplitude and
width conventions are given in the [amplitude](alp-photon-pp-amplitude.md)
and [width](alp-photon-pp-width.md) records; the numerical algorithm is
documented in [photon-pp-gauss-legendre-quadrature.md](photon-pp-gauss-legendre-quadrature.md).
No external literature input was consulted.

The coupling object is the fixture already used by the production quadrature
benchmark. It is not the direct no-mixing source vector used in the separate
[nested-adaptive diagnostic](photon-pp-kaon-pion-profile.md); therefore these
width values are not a cross-implementation agreement check. The timing ratio
here compares only pion and kaon channels within this Gauss–Legendre API.

# Method

The executable diagnostic is
[`scripts/profile_photon_pp_gauss_legendre.py`](../../../scripts/profile_photon_pp_gauss_legendre.py).
It warms each case, measures full width calls over seven repetitions, records
the cProfile hotspots, counts vectorized grid points and adaptive fallback
calls, and lists the prepared pole terms. Reproduce with:

```text
.venv/bin/python -m scripts.profile_photon_pp_gauss_legendre --repeats 7
```

# Results

On 2026-10-22, the median times were $10.80$ ms for $\gamma\pi^+\pi^-$ and
$16.84$ ms for $\gamma K^+K^-$, a ratio of $1.56$. The width values were
$8.1151\times10^{-12}$ GeV and $1.2170\times10^{-10}$ GeV, respectively.
Both channels converged using quadrature orders 16 and 32, with no adaptive
fallback. The pion calculation evaluated 70,400 vectorized Dalitz points; the
kaon calculation evaluated 98,560, a factor of $1.40$. The profiled Python
call totals were similar (3,623 for pions and 3,643 for kaons), since most
pointwise work is performed by NumPy array operations.

The pion amplitude has three broad $\rho$ pole terms. The kaon amplitude has
five terms: the $s$-channel $\rho$, $\omega$, and $\phi$, plus the crossed
$K^*$ terms. The $\phi$ pole, with width $0.004249$ GeV, lies in the physical
$s$ range; the $\rho$ and $\omega$ centers lie below the $K^+K^-$ threshold.
The crossed $K^*$ widths are about $0.0514$ GeV. Explicit resonance-centered
subdivision yields five $s$ panels for pions and seven for kaons; both use
eleven $t$ panels per $s$ node in the measured grids. The longer
kaon time is consistent with more grid points and more pole terms; the
profile does not isolate these effects independently.

# Interpretation and limits

For these equal-mass channels, the Gauss–Legendre calculation reduces the
measured kaon/pion timing contrast compared with the earlier nested-adaptive
profile, but the absolute times across the two profiles must not be compared
as a benchmark because their ALP coupling inputs differ. This result applies
to the stated mass, source fixture, fit, and machine only.
