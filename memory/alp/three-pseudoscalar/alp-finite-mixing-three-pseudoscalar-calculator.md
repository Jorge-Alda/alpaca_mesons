---
title: Finite-mixing LO three-pseudoscalar amplitude and width calculator
status: mixed
created: 2026-11-13
updated: 2026-11-13
bibtex_keys: [Alda:2026nrd]
provenance: numerical-computation
---

# Finite-mixing LO three-pseudoscalar amplitude and width calculator

## Inputs and provenance

- Analytical amplitude and external-pole construction: [finite-mixing three-pseudoscalar derivation](alp-finite-mixing-three-pseudoscalar-amplitudes.md), with LO quadratic blocks and source derivatives from the [joint-vertex derivation](../shared/alp-lo-joint-vertex-expansion.md), [exact pole equations](../shared/alp-lo-exact-pole-equations.md), and [pole solver](../shared/alp-lo-pole-solver.md), in the conventions of Alda:2026nrd [@Alda:2026nrd].
- The implementation is [`scripts/calculate_alp_three_pseudoscalar.py`](../../../scripts/calculate_alp_three_pseudoscalar.py). It imports `calculate_masses` from [`scripts/calculate_alp_meson_masses.py`](../../../scripts/calculate_alp_meson_masses.py). Inputs are the same $f_a,m_{a,0},\widehat S_n,\widehat P_n,\widehat\Theta_n,f_\pi,B_0,m_0^2$ as that mass calculator, plus an ordered triple of meson labels. The amplitude callable takes $(s_{12},s_{23})$ in GeV$^2$; the width is in GeV.
- The aligned vacuum and electric-charge-conserving source assumptions of the mass calculator apply. Lists are derivatives with respect to $x=a/f_a$ at $x=0$, not Taylor coefficients. This generic nonlinear-source amplitude uses $\widehat S_0,\ldots,\widehat S_4$, $\widehat P_0,\ldots,\widehat P_3$, and $\widehat\Theta_0,\ldots,\widehat\Theta_3$; absent derivatives must be supplied as zero entries. The first two orders used for the masses are passed to the mass calculator.
- Neutral labels use the maximum-overlap convention of that mass calculator, including `K0_1` and `K0_2` for its real neutral-kaon pole states. They are not the strangeness eigenstates $K^0,\bar K^0$ or the CP mixtures $K_S,K_L$.

## Method

1. Assemble the ten-real-field LO mass-squared matrix in the neutral and charged generator basis using $\chi_0=2B_0\widehat S_0$ and the linear invariant sources $\widehat P_1,\widehat\Theta_1$. With charge-conserving inputs, the charged generators have exactly zero entries connecting them to the ALP and neutral generators. The code extracts and diagonalizes the six-dimensional neutral block for external neutral poles, then diagonalizes the full matrix once to form the internal spectral propagator; the charged eigenstates are needed for charged-meson exchange. The charged and neutral sectors are block diagonal by construction, not separated by an approximation. External labels and masses match the imported mass calculator.
2. For each external leg, use its complete meson matrix and ALP component from the pole vector. The parent is the physical ALP pole; a daughter may itself have an ALP component when it is a neutral mixed pole. Charged daughters have no such component.
3. Evaluate the quartic contact and the three cubic-exchange channels from the finite-mixing pole-aligned derivation. The internal propagator is the full spectral inverse of the LO mass-squared matrix, not a selected meson or ALP pole. The implementation combines the cubic vertices with the spectral residues once during amplitude preparation.
4. Use $s_{13}=m_a^2+m_1^2+m_2^2+m_3^2-s_{12}-s_{23}$ and the on-shell dot-product identities to express the kinetic contact as a polynomial in the Dalitz variables. The flavour contractions, external vectors, and exchange residues are computed once, outside repeated amplitude evaluations.
5. Integrate $|\mathcal M|^2$ over the physical Dalitz domain and divide by $256\pi^3m_a^3$ and by the factorials of repeated identical labels. Return exactly zero for nonzero total electric charge, for a closed phase space, or when all prepared amplitude coefficients vanish by the tested traces.

## Efficiency assessment

The weak-mixing formula in the analytical record is useful as a leading-order limit, but it omits the relative $O((f_\pi/f_a)^2)$ terms: parent normalization, mixed daughter components, nonlinear contact terms, and all exchange diagrams. It cannot replace the exact amplitude for this request. Its evaluation strategy does carry over: reduce dot products with the on-shell identities and prepare matrix traces and pole residues outside the phase-space loop. No additional weak-mixing expansion is made.

The implementation uses the shared vectorized Gauss-Legendre integrator [`alpaca_chiral/_three_body_quadrature.py`](../../../alpaca_chiral/_three_body_quadrature.py), with adaptive fallback and pole subdivisions, following the [optimized decay-integrand guidelines](../shared/optimized-decay-integrand-guidelines.md). Flavor contractions and spectral residues are prepared once. Exactly equal floating-point pole values are combined with no tolerance, so no near-degenerate poles are merged. The amplitude has zero-width LO internal poles. If an exchanged invariant can reach a pole with nonzero total residue, the stable-particle tree width is singular and requires a width prescription beyond this calculation; the implementation raises `ValueError` instead of silently regularizing it. A roundoff-split exact degeneracy may be conservatively reported as an accessible pole.

## Invocation and verification status

```bash
PYTHONPATH=scripts .venv/bin/python -c 'import numpy as np; from calculate_alp_three_pseudoscalar import decay_width; z=np.zeros((3,3)); print(decay_width(1., .6, [np.diag([.003,.005,.1]),z,z,z,z], [z,z,z,z], [0.,0.,0.,0.], ("pi0","pi+","pi-"), f_pi=.092, B_0=2.5, m_0_squared=.75))'
```

A newly added test file received only a partial run: two checks reported success, while completion of the remaining check was not established. Therefore the implementation's full numerical verification remains pending. The calculation is tree-level LO and does not include loop effects, resonance exchange, or a prescription for unstable internal poles.
