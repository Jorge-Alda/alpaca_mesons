---
title: The BaBar 3pi spectrum has no one-multiplet-free mass window
status: numerical
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [BaBar:2004ytv, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# No Mass Window Isolates the One-Multiplet Baseline in the BaBar $3\pi$ Data

## Question

The [BaBar blocker record](babar-3pi-spectrum-fit-blocker.md) rejects the
$e^+e^-\to\pi^+\pi^-\pi^0$ table on general grounds: the $\omega$ and $\phi$ peaks
are absent and our baseline is a single multiplet. A natural objection is that
one might still use a restricted mass window far from $\omega'(1350)$ and
$\omega''(1660)$, with the extra states carried as nuisance parameters and their
error budget assessed. This record answers that objection quantitatively.

**Answer: no such window exists.** The heavier states own at least 73% of the
coherent sum at every one of the 78 bins.

## Method

Inside the coherent sum of BaBar Eq. (14) the magnitude of a resonance term is

$$
w_V(m)=\frac{\Gamma_V m_V^{3/2}\sqrt{{\cal B}(V\to e^+e^-){\cal B}(V\to3\pi)}}{|D_V(m)|},
\qquad D_V(m)=m_V^2-m^2-im\Gamma_V,
$$

which is independent of the common factor $F_{\rho\pi}(m)$ and of the phases.
Comparing $w_V$ therefore answers separability without needing $F_{\rho\pi}$ or
the interference phases, none of which are tabulated. This is **not** a
prediction of $\sigma$ and **not** a fit; it is a comparison of term magnitudes.

Reproduce with [`scripts/diagnose_three_pion_resonance_separation.py`](../../../scripts/diagnose_three_pion_resonance_separation.py).

## Inputs and provenance

- Bin centres, $\sigma$, and the statistical/systematic split from
  [`babar_3pi_born_spectrum.csv`](../../../data/odd_parity/babar_3pi_born_spectrum.csv),
  checksum-validated by `scripts/extract_odd_parity_tables.py`, 78 uniform 25 MeV
  bins spanning $\sqrt{s'}=1.0625$--$2.9875$ GeV.
- $m_V,\Gamma_V$ and ${\cal B}(V\to e^+e^-){\cal B}(V\to3\pi)$ for
  $\omega',\omega''$ from BaBar Summary Eqs. (19)--(20) [@BaBar:2004ytv].
- $\omega$ and $\phi$ masses and widths are **inputs** to our model, taken from
  the PDG listing [@ParticleDataGroup:2024cfk]; their tabulated branching
  products are used only for the weight scale.

## Result 1 — the heavier states are less than one width away

| state | $m_V$ (MeV) | $\Gamma_V$ (MeV) | $m-m_\phi$ | in units of $\Gamma_V$ |
|---|---|---|---|---|
| $\phi$ | 1019.5 | 4.25 | — | — |
| $\omega'$ | 1350 | 450 | 330.5 MeV | **0.73** |
| $\omega''$ | 1660 | 230 | 640.5 MeV | 2.78 |

## Result 2 — their Lorentzians never decay away

The minimum of $|D_V(m)|^{-2}$ over $m>0$, relative to the peak value, is
$(\Gamma_V/m_V)^2$:

| state | $\Gamma_V/m_V$ | floor, as % of peak |
|---|---|---|
| $\omega$ | 0.0111 | 0.01% |
| $\phi$ | 0.0042 | 0.002% |
| $\omega'$ | **0.3333** | **11.1%** |
| $\omega''$ | 0.1386 | 1.9% |

$\omega'$ is broad enough that no low-mass window excludes it. Its
$52\times$ smaller branching product is cancelled by its tail being only
$\sim14\times$ more suppressed than the $\phi$ tail at the same mass.

## Result 3 — weight fractions inside the coherent sum

| $\sqrt{s'}$ (GeV) | $\sigma$ (nb) | $\omega$ | $\phi$ | $\omega'$ | $\omega''$ |
|---|---|---|---|---|---|
| 1.0625 (first bin) | 1.61 | 6.3% | 21.1% | **50.1%** | 22.5% |
| 1.1625 | 4.17 | 4.6% | 6.3% | 62.2% | 26.9% |
| 1.3125 | 6.75 | 2.6% | 2.5% | 63.7% | 31.1% |
| 1.6375 | 5.33 | 1.2% | 0.8% | 27.3% | 70.7% |
| 1.9875 | 0.30 | 2.0% | 1.3% | 37.4% | 59.3% |
| 2.9875 (last bin) | 0.20 | 3.1% | 1.9% | 47.0% | 48.0% |

Minimum combined $\omega'+\omega''$ weight anywhere: **72.6%**, at the first bin.
Maximum $\phi$-tail weight anywhere: **21.1%**. Bins below a 10% heavy-state
threshold: **none**. The excess over the $\phi$ tail accumulates fastest in the
first bins — by $m\le1.3375$ GeV, 24% of it is already below that point — so
restricting to the lowest bins does not help either.

## Result 4 — the error-budget trap

If $\omega',\omega''$ are fixed at their fitted values and their quoted
uncertainties propagated, the induced shift of their combined weight fraction is
**small**: 1.28% worst case in the first bin, 0.18--0.40% over the rest,
against 10--53% statistical errors. Widths dominate the sensitivity (6.2% and
12.5% shifts in $w$) over masses (2.9% and 1.4%).

This benign-looking budget is **not** a rescue, and the reason matters: it is
benign precisely *because* the heavy states own 73--98% of the weight. A small
fractional perturbation of the dominant term is still a large absolute
systematic. Such a fit would report a tight, well-conditioned constraint that is
almost entirely a re-observation of the parameters supplied for $\omega'$ and
$\omega''$, propagated into the $c_i,d_i$. This is the same failure mode already
observed with $\mathcal B(\omega\to\pi^0\gamma)$, where the single connected
datum adds rank one yet has a normalized residual of $\approx0.045$, i.e. it is
already saturated. See
[the saturation probe record](radiative-channel-saturation-probe.md).

## Conclusion and its limit

Restricting the mass window does not rescue this dataset. Both the heavy-state
weight floor and the propagating-error argument must be stated together; either
alone would mislead.

**Limit of this calculation.** It compares term magnitudes *within BaBar's own
four-resonance model*. It shows that if that decomposition is adopted, the
heavier states own the range. It is not an independent proof that they are
physically required, although BaBar's finding that four isoscalar states
describe the spectrum well only up to 1.8 GeV is corroborating evidence
[@BaBar:2004ytv]. A genuine test of the one-multiplet hypothesis needs an
amplitude, not a weight comparison, and that is the dispersive treatment scoped
in [the treatment-plan record](dispersive-three-pion-treatment-plan.md).