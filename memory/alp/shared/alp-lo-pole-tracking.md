---
title: Pole tracking and avoided crossings along a parameter path
status: numerical
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: numerical-computation
---

# Pole tracking and avoided crossings along a parameter path

## Assumptions and inputs

- LO neutral blocks of the [neutral-reduction record](alp-lo-neutral-block-reduction.md) with the conventions of Alda:2026nrd [@Alda:2026nrd]; spectrum at each point from the [pole solver](alp-lo-pole-solver.md) and the [stability check](alp-lo-stability-conditions.md); an unstable point stops the scan.
- Paths: $m_{a,0}^2$ (GeV$^2$) at fixed $f_\pi/f_a$, or $\epsilon=f_\pi/f_a$ at fixed $f_\pi$ and $m_{a,0}^2$.
- Assumption (heuristic): labels by continuity are a convention that makes the plot readable. Near an avoided crossing the label of the "ALP-like" state is not an invariant property.
- Inputs are illustrative, not fits.

## Method

1. Branches are labelled at the first point by increasing pole.
2. Between consecutive points the $K$-weighted overlaps $O_{bn}=|v_b^{\rm prev\,\mathsf T}K_{\rm mid}v^{(n)}_{\rm cur}|$ with $K_{\rm mid}=(K_{\rm prev}+K_{\rm cur})/2$ are maximized by a one-to-one assignment; signs make the assigned overlaps positive. For an exactly degenerate group (such as the $K^0$ doublet for flavour-diagonal sources) the basis is first rotated to match the previous vectors, because the solver's basis is arbitrary there.
3. The step is bisected (down to a minimal step) if an assigned overlap is below a threshold, if two coupled branches exchange order, or if the largest pole change exceeds a stated factor of the smallest gap between neighbouring ALP-carrying poles. An exchange of order is accepted without refinement only when the two states are decoupled (their matrix element of $s\,\mathrm dK-\mathrm dM^2$ is small); an exchange between coupled states that persists at the minimal step is recorded as unresolved.
4. The ALP-like branch at each point is the branch with the largest $|v_a|$.
5. An avoided crossing is an interior local minimum of the gap between two ALP-carrying neighbouring branches, refined by a bounded one-dimensional minimization.

## Results (illustrative inputs)

- Isolated two-level check: one avoided crossing, with minimal gap equal to $2|M^2_{P'a,\rm can}|$ of the [two-level record](alp-lo-near-degeneracy-two-level.md), located where $M^2_{aa,\rm can}=m_\pi^2$.
- Scan in $m_{a,0}^2$, flavour-diagonal sources, $\epsilon=1/20$: three avoided crossings with ALP-like switches, with minimal gaps $3.4\times10^{-5}$, $7.5\times10^{-2}$ and $5.4\times10^{-2}$ GeV$^2$, and no unresolved crossing and no order exchange between ALP-carrying branches. Crossings with the exact $K^0$ pole and other decoupled states are decoupled crossings.
- Scan in $m_{a,0}^2$, neutral sources, $\epsilon=1/20$: also an avoided crossing between two mesons ($|v_a|\simeq0.12$ on both) induced by their common ALP mixing, which does not switch the ALP-like label.
- At $\epsilon=0.005$ and $m_{a,0}^2=0.3$ GeV$^2$ the ALP-like pole lies within $10^{-3}$ GeV$^2$ of $m_{a,0}^2$ and $|v_a|$ within $10^{-2}$ of $1$ (consistent with the [weak-mixing record](alp-lo-weak-mixing-expansion.md)).
- At large mixing ($\epsilon=23/50$) the ALP-like label switches between two branches with $|v_a|\simeq0.62$ and $0.65$ at the closest approach: the ALP-like label is ambiguous there (both poles carry comparable ALP weight), so statements such as "the ALP mass" are not meaningful without specifying the pole.

## Script and verification

[`scripts/track_alp_meson_poles.py`](../../../scripts/track_alp_meson_poles.py); invoke:

```bash
.venv/bin/python scripts/track_alp_meson_poles.py
```

All checks pass (no unresolved crossings; positive gaps; every ALP-like switch has an avoided crossing). Units: poles and gaps in GeV$^2$; $\epsilon$ dimensionless.
