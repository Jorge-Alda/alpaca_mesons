---
title: Adaptive profile comparison for charged kaon and pion radiative decays
status: numerical
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: []
provenance: numerical-computation
---

# Assumptions and inputs

This diagnostic compares the nested adaptive width implementation in
[`scripts/implementation/alp_photon_pp_width.py`](../../../scripts/implementation/alp_photon_pp_width.py)
for $a\to\gamma\pi^+\pi^-$ and $a\to\gamma K^+K^-$. It uses $m_a=1.5$ GeV,
$f_a=100$ GeV, the dimensionless diagonal axial source
$c_q^A=\operatorname{diag}(0.8,0,0)$, zero pseudoscalar source, the same
preliminary fit inputs, `epsrel=10^{-5}`, `epsabs=0`, and `limit=150` in both
channels. The amplitude and width conventions are those in the
[amplitude](alp-photon-pp-amplitude.md) and
[width](alp-photon-pp-width.md) derivations. No external literature input was
consulted.

# Method

The executable profile is
[`scripts/profile_photon_pp_channels.py`](../../../scripts/profile_photon_pp_channels.py).
It warms each channel, measures complete width calls, and separately runs
`cProfile` over one call. It reports the prepared pole terms and the adaptive
quadrature call counts. Reproduce with:

```text
.venv/bin/python -m scripts.profile_photon_pp_channels --repeats 3
```

# Results

Across two three-run invocations on 2026-10-22, the median complete-width
times ranged from $0.0407$--$0.0455$ s for $\gamma\pi^+\pi^-$ and
$0.2605$--$0.3051$ s for $\gamma K^+K^-$. The measured kaon-to-pion time ratio
was $5.7$--$7.5$; the individual runs show noticeable timing spread. The width
values were repeat-stable at $8.8797\times10^{-10}$ GeV and
$1.9315\times10^{-10}$ GeV, respectively.

In the profiled call, `integrate_t` ran 126 times for the pion channel and
378 times for the kaon channel. The polarization-summed scalar was evaluated
8,862 and 24,948 times, respectively. Total profiled Python call counts were
156,311 and 486,909. The prepared pion amplitude has one broad $\rho$ pole in
each of $s,t,u$ (widths approximately $0.148$ GeV). The kaon amplitude has
three $s$-channel poles and crossed $K^*$ poles: the $\phi$ pole has width
$0.00425$ GeV, the $\omega$ pole has width $0.00866$ GeV, and the $K^*$ poles
have width $0.0514$ GeV. Only poles inside the physical Dalitz region can
drive refinements; at this ALP mass the narrow $\phi$ lies in the physical
$s$ range. The kaon evaluation also coherently sums five prepared pole terms,
compared with three for the pion case.

# Interpretation and limits

The measured slowdown comes from both more adaptive evaluations and higher
cost per evaluation. The pole structure is a plausible source of the extra
quadrature work: the kaon channel exposes narrow resonance bands, while the
pion channel has broad $\rho$ poles. The prepared kaon amplitude also has more
terms, increasing scalar evaluation cost. This profile does not isolate the
individual contribution of each pole, so it does not establish that one pole
alone causes the slowdown. Results are specific to this fit, source choice,
ALP mass, machine, and adaptive implementation. In particular, the repository
also contains a separate Gauss–Legendre implementation for the
`alpaca_chiral` API; these timings do not characterize that API.
