---
title: Numerical implementation of only-up radial widths
status: numerical
created: 2026-11-19
updated: 2026-11-19
bibtex_keys: [DiLuzio:2026poh, Alda:2026nrd]
provenance: numerical-computation
---

# Numerical radial widths

The function `sigma_two_pseudoscalar_width` in [only_up_decays.py](../../../scripts/implementation/only_up_decays.py) implements the Cartesian vertex and two-body width of the [radial decay record](murayama-sigma-two-pseudoscalar-decays.md), using its assumptions and GeV units. Inputs are $\kappa_u$, $\lambda$ and two physical state labels; the remaining inputs and $f_a$ come from [only_up_inputs.py](../../../scripts/implementation/only_up_inputs.py).

It computes $m_\sigma^2=B_0m_uf_\pi^2/f_a^2+2\lambda f_a^2$. The existing finite-mixing mass calculator numerically diagonalizes the leading-order quadratic matrix; the three-pseudoscalar preparation routine constructs pole vectors with the same matrix and labelling convention. The width uses its returned masses and vectors, evaluates the Cartesian vertex and applies the two-body phase space and identical-particle factor. This diagonalization is a project calculation, not a generic diagonalized-mass formula from DiLuzio:2026poh. The up-quark projector gives zero coupling to a neutral-kaon pair, including any $K_S,K_L$ combination; all other pairs outside the neutral sector and the charged conjugate pairs violate charge or strangeness conservation. These channels return zero. The approximation is leading chiral order and tree level, with stable particles and no rescattering, as in the input record. See the [pole solver](alp-lo-pole-solver.md) for the finite-mixing conventions.

Verification: `.venv/bin/python -m unittest scripts.implementation.test_sigma_width` checks the charged analytic width, threshold, neutral ordering, an open identical neutral channel, and all pairs of recognized physical labels.

## Plot mass labels

The `sigma_PP` plot in [murayama.py](../../../scripts/implementation/murayama.py) takes the radial mass array from the same $\kappa_u$ scan and $\lambda$ used for the widths. The scan uses 500 logarithmically spaced points from $10^{-3}$ to $1$. The upper axis labels round masses in GeV using the values $1$, $2$ and $5$ times powers of ten. For each value, crossings are found by linear interpolation of the supplied mass array against $\log_{10}\kappa_u$. Crossings on both sides of a minimum are considered separately; labels closer than $0.16$ in $\log_{10}\kappa_u$ are suppressed to limit overlap. This interpolation approximates the mass curve between scan points; no inverse mass function is assumed. The array must correspond to the coupling supplied to `sigma_PP`.

Verification: `.venv/bin/python -m unittest scripts.implementation.test_murayama_sigma_plot` checks round labels and their positions on both branches of a nonmonotonic mass curve, as well as their spacing.
