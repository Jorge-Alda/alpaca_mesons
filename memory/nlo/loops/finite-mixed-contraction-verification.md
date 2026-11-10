---
title: Verification of finite mixed tadpole contractions
status: derived
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, Gasser:1983yg]
provenance: analytical-derivation
---

# Verification of finite mixed tadpole contractions

## Inputs and assumptions

Use the source and field conventions, the $Q=sK-M^2+\Sigma$ sign, and the loop content of [mixed ALP–meson self-energy](mixed-alp-meson-self-energy.md), the expansion in [vertices and topologies](vertices-and-topologies.md), and the finite propagator of [loop functions and subtraction](loop-functions-and-subtraction.md). Masses and decay constants are in GeV and $s$ in GeV$^2$. Work to one loop and first order in $1/f_a$, with arbitrary noncommuting flavour matrices, no CP-violating sources and no constant $\theta$ term.

## Derivation and result

1. The script expands the leading Lagrangian using the independent building-block expansion in `scripts/derive_nlo_loop_vertices.py`. It replaces $\Phi$ by an arbitrary external matrix times a background field plus two independent arbitrary fluctuation matrices. It takes the derivative field to contain the external field: coincident field–derivative contractions vanish in dimensional regularization.
2. It extracts the coefficient linear in the background and in each fluctuation, separately for $c_q^A$ and $c_q^P$. It compares the coefficient to the sum of both orders of the fluctuation matrices in the formula in [mixed ALP–meson self-energy](mixed-alp-meson-self-energy.md). The axial coefficient multiplies $\partial_\mu a\,\partial^\mu P_i$ and therefore $s$ in the stated inverse-propagator convention.
3. Both differences vanish by exact symbolic trace algebra. The two fluctuation matrices are arbitrary, so this identity applies to every pair of nonet generators and every external generator. It applies to diagonal propagators by identifying the two fluctuation matrices: the coefficient of the squared fluctuation is half the coefficient of the two distinct fluctuations. Linearity then covers any symmetric finite propagator, including neutral-meson mixing and all mass eigenvalues. No commutation of sources or masses is assumed.

The absent $c_q^V$, $c_q^S$ and $c_G$ vertices are treated in [vertices and topologies](vertices-and-topologies.md); this comparison concerns the nonzero axial and pseudoscalar tadpoles. The result verifies the contraction of the stated vertex expansion, conditional on the inputs and subtraction; it is not an independent calculation of the loop integral.

## Verification

Run `.venv/bin/python scripts/verify_mixed_finite_contractions.py`. Its unit test compares the two extracted symbolic coefficients with the formula, retaining both fluctuation orders, and prints both sides using SymPy's LaTeX printer.
