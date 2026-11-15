---
title: Conditional matching-scale and effective-pole variations
status: derived
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Hoferichter:2018kwz]
provenance: analytical-derivation
---

# Inputs and assumptions

Use the [spectral matching prescription](alp-spectral-high-energy-matching.md), its ALP source trace, common asymptotic distribution amplitude and leading-order fixed-scale coefficient. Natural units are used; $Q^2=-q_1^2=-q_2^2>0$, and $s_m$ and $M_{\rm eff}^2$ have units $\mathrm{GeV}^2$. The [tree and dispersive forms](alp-low-energy-matching-benchmarks.md) are conditional independent benchmarks, not fit results. All variations below hold their low-energy input, source matrices, renormalization scheme and effective-pole parameters fixed. Changes of these inputs and possible spectral overlap are separate uncertainties.

# Matching-scale variation

The existing continuum formula gives
$$
F^{\rm asym}(-Q^2,-Q^2;s_m)=\frac{4f_\pi^2}{f_a}\left\langle Q^2\mathbb M_a^{-1}(\bm{\mathcal C})\right\rangle_\mu\frac{(Q^2)^2}{(s_m+Q^2)^3}.
$$
For two arbitrary positive endpoints $s_m,s_m'$, the exact finite change is
$$
\frac{4f_\pi^2}{f_a}\left\langle Q^2\mathbb M_a^{-1}(\bm{\mathcal C})\right\rangle_\mu (Q^2)^2\left[\frac1{(s_m'+Q^2)^3}-\frac1{(s_m+Q^2)^3}\right].
$$
Its derivative is
$$
\frac{\partial F^{\rm asym}}{\partial s_m}=-\frac{12f_\pi^2}{f_a}\left\langle Q^2\mathbb M_a^{-1}(\bm{\mathcal C})\right\rangle_\mu\frac{(Q^2)^2}{(s_m+Q^2)^4}.
$$
The variation vanishes at $Q^2=0$ and changes only subleading terms at large $Q^2$. It is the variation of the continuum alone: if the low-energy density, overlap subtraction or effective-pole parameters also depend on $s_m$, their variations must be added. The derivative is not a physical uncertainty until an admissible interval and correlated inputs are supplied.

# Analytic effective-pole comparison

Following the one- and two-pole variants of Messerli:2025rnv, compare the pole term in the spectral record with two stable poles of masses $M_1,M_2$ and residues $g_1,g_2$. Impose the same residual normalization and singly virtual $1/Q^2$ coefficient:
$$
g_1+g_2=R_{\rm eff},\qquad g_1M_1^2+g_2M_2^2=R_{\rm eff}M_{\rm eff}^2.
$$
For distinct masses the script solves these equations. Their exact equal-spacelike difference is
$$
\frac{g_1M_1^4}{(M_1^2+Q^2)^2}+\frac{g_2M_2^4}{(M_2^2+Q^2)^2}-\frac{R_{\rm eff}M_{\rm eff}^4}{(M_{\rm eff}^2+Q^2)^2},
\qquad
g_1=R_{\rm eff}\frac{M_{\rm eff}^2-M_2^2}{M_1^2-M_2^2},\quad
g_2=R_{\rm eff}\frac{M_1^2-M_{\rm eff}^2}{M_1^2-M_2^2}.
$$
The difference vanishes at the origin and is subleading to $1/Q^2$; its value is a conditional model spread, not a fitted uncertainty. The two masses and residual constraints are required inputs; no positivity of the residues is assumed. Both models introduce pole strength, and physical overlap must be resolved separately. The script applies the same construction independently to the explicit tree and symbolic dispersive benchmarks; their identical low-energy pieces cancel from these differences. The symbolic dispersive input is not a complete TFF.

# Verification

Run `.venv/bin/python scripts/derive_alp_matching_variation.py` and append `--test`. SymPy derives the scale difference and pole residues and checks both benchmarks, endpoint constraints and high-energy limits. The scale variation and pole-model spread are reported separately from hadronic inputs and scheme dependence; no numerical interval, dilepton-loop variation, or fitted uncertainty is claimed. No Euclidean-only diagnostic is used in this comparison.
