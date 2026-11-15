---
title: Conditional low-energy form factors for matching tests
status: derived
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Inputs and conventions

Both forms below are **conditional symbolic benchmarks** for subsequent matching tests, not calibrated predictions. Natural units, $s_i=q_i^2$ and the unnormalized scalar form factor of [the tree-level record](alp-doubly-virtual-tff.md) are used; $F$ has units $\mathrm{GeV}^{-1}$. All couplings, source traces, masses, spectral densities and cutoff are independent symbolic inputs; no vector-fit coefficients are used. The tree form assumes stable diagonal neutral vectors and the truncation in that record. The spectral form assumes the finite-cutoff double-spectral representation of [the isovector construction](alp-isovector-tff-spectral-construction.md), with its phase and normalization conventions. It is a component, not a complete low-energy form factor.

# Conditional tree benchmark

With $W_i(s)=R_i/(M_i^2-s)$ and the explicit $V_{i\gamma}$ and $P_{ij}$ given in the linked tree record, retain the source-resolved expression
$$
\frac{\widehat c_\gamma}{4\pi^2f_a}+2\sqrt2\sum_i[W_i(s_1)V_{i\gamma}(s_1,s_2)+W_i(s_2)V_{i\gamma}(s_2,s_1)]-2\sum_{i,j}W_i(s_1)W_j(s_2)P_{ij}(s_1,s_2).
$$
This is the actual tree expression, including both photon placements and ordered vector pairs. It includes the contact coefficient, so no additional meson-pole WZW contribution is inserted. The executable uses two symbolic vector species as a declared finite example, with no numerical fit input; the displayed sum admits the full neutral basis. Stable poles do not describe finite-width cuts.

# Conditional dispersive benchmark

For $s_i$ off the spectral cuts, take the symbolic finite-cutoff isovector contribution
$$
\frac1{\pi^2}\int_{4M_\pi^2}^{\Lambda^2}\!dx\int_{4M_\pi^2}^{\Lambda^2}\!dy\,\rho_a^{I=1}(x,y)\left[\frac1{(x-s_1)(y-s_2)}+\frac1{(x-s_2)(y-s_1)}\right].
$$
Its zero-virtuality value is $2\pi^{-2}\int dx\,dy\,\rho_a^{I=1}(x,y)/(xy)$. The density is an unspecified source-dependent input and the cutoff is a model choice. The representation omits single-spectral and subtraction pieces, other channels and high-energy strength. It must not be added to the tree benchmark without a separate overlap prescription. For timelike arguments use the boundary values prescribed in the spectral record.

# Method and verification

The [script](../../../scripts/derive_alp_low_energy_benchmarks.py) builds the ordered tree vertices and the two photon placements from the linked inputs, and constructs the symbolic double integral. Run `.venv/bin/python scripts/derive_alp_low_energy_benchmarks.py` and append `--test` for four unit checks of Bose symmetry, the spectral normalization at zero and the contact limit. The tests do not establish a matching prescription or a numerical spectral density.
