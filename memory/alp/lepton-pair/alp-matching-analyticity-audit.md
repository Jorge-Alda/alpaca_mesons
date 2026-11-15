---
title: Analyticity of low and high energy matching examples
status: derived
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv]
provenance: analytical-derivation
---

# Inputs and assumptions

Use natural units and the conditional [low-energy benchmarks](alp-low-energy-matching-benchmarks.md), the [light-cone result](alp-light-cone-high-energy-tff.md), and their stated conventions. The following one-variable examples use $Q^2=-s>0$ and do not specify the doubly virtual form factor. The constants $F_0$ and $C$ represent a chosen normalization and a leading $1/Q^2$ coefficient on that ray; they are independent symbolic inputs, not numerical predictions. No spectral density is fitted.

# Sharp switch

Take $F_0 M^2/(Q^2+M^2)$ below the matching scale and $C/Q^2$ above it. Their difference is generically nonzero at the switch, as calculated by the script. A condition on the inputs can enforce continuity at one point. Even then, if a single holomorphic form factor equalled the low-energy expression on an open spacelike interval, the identity theorem would fix it to that expression throughout their connected analytic domain. It cannot equal a different high-energy expression on another open interval. A sharp switch can therefore be a Euclidean integration prescription, but does not by itself supply an analytic continuation or a physical discontinuity.

# Analytic example with added spectral strength

For distinct positive $M^2,N^2$, solve for two pole residues using zero correction at $Q^2=0$ and leading coefficient $C$ at large $Q^2$. The script produces
$$
\frac{C Q^2}{(Q^2+M^2)(Q^2+N^2)}
=\frac{C M^2}{M^2-N^2}\frac{1}{Q^2+M^2}
-\frac{C N^2}{M^2-N^2}\frac{1}{Q^2+N^2}.
$$
It preserves the normalization and introduces no discontinuity on a low-energy cut that is disjoint from both pole locations. It adds two delta-function spectral contributions at the poles, so it does not preserve the complete discontinuity of the original form. If a pole lies on an existing cut, the additional distribution must be counted in that cut. This example does not identify or subtract any overlap with the existing high-energy spectral density, and it is not a physical completion by itself. The unequal residues have opposite signs for real $C$.

For an analytic spectral construction, one may instead keep the given low-energy discontinuity and append a separately specified high-energy spectral contribution with disjoint support, imposing the normalization and asymptotic moment constraints on its density and subtraction terms. Its existence and overlap with a given hard expression require those inputs; the finite-cutoff double-spectral benchmark alone also omits single-spectral and subtraction pieces. No claim of absence of double counting follows merely from an unspecified density.

# Verification

Run `.venv/bin/python scripts/audit_alp_matching_analyticity.py` and append `--test` for four checks of the example's normalization, asymptotic coefficient, nonzero pole residues and generic discontinuity at the sharp join. The identity-theorem argument is analytic; the symbolic tests do not prove a physical dispersive matching prescription.
