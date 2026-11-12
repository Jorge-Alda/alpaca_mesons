# Proposal

## Why

Beads issue alp-sjp.6 asks how finite ALP-meson mixing changes the physical meson and ALP masses. The LO spectrum is solved numerically (`lo-finite-mixing-spectrum`), and the weak-mixing record gives implicit or first-order pole formulas. There is no closed analytical statement of the poles that holds at finite mixing, no expression in terms of the redefinition-invariant couplings, and no exact sum rules to test numerical spectra against. This change supplies that part of alp-sjp.6, at leading order and analytically.

## What Changes

Scope is leading order (the $O(p^2)$ chiral Lagrangian of Alda:2026nrd with the full $K$, $M^2$ blocks of the [finite-mixing record](../../../memory/alp/shared/finite-alp-meson-mixing.md)), tree level, neutral and charged mesons, and the ALP. Conventions are those of Alda:2026nrd. The deliverable is analytical expressions with symbolic proofs; scripts verify them (exact rationals and the existing solver are used only as checks, no numerical inputs, no scans). No expansion in $f_\pi/f_a$ is part of the change.

- **Exact relations**: all-orders statements that need no expansion: zero shift of the charged mesons (decoupling theorem), exact $K^0$ pole at $m_{K^0,0}^2$ (one combination), trace sum rule $\sum_n m_n^2=\operatorname{tr}(K^{-1}M^2)$ and product rule $\prod_n m_n^2=\det M^2$ (with $\det K=1$), each written with the explicit analytic $K^{-1}$ and in the sources $c_q^A$, $c_q^P$, $c_G$, $\chi_0$, $m_0^2$.
- **Exact pole equations**: the neutral meson and ALP poles at finite mixing, written in the invariants $c_q^{AP}$, $c_{GA}$, $\chi_0$, $m_0^2$, $m_{a,0}^2$ and $\operatorname{tr}\widehat S_2$. In the canonical frame the poles are the eigenvalues of one real symmetric matrix, with interlacing bounds against the unmixed masses. Invariance of the spectrum under constant axial redefinitions of alp-sjp.1 is proved symbolically.
- **Explicit forms**: for diagonal $\chi_0$ the exact secular polynomial, and the exact poles of $\pi^0$, $\eta$, $\eta'$ and the $K^0$ doublet.
- **Closed approximations**: an isolated meson pole and an isolated ALP pole, each in closed form, with an exact error bound for the meson and a fixed-point iteration for the ALP. They hold when the poles are far from each other, which is stated as an assumption.
- **Regimes**: in the small-mixing regime no shift at order $f_\pi/f_a$ is kept. In the large-mixing regime the finite shifts are the exact poles above or their closed approximations.
- **Limits and domain**: the $f_a\to\infty$ limit returns the unmixed masses, the level-repulsion sign rule, and the failure condition near degeneracy with the exact two-level formula as replacement.
- Memory records (one concept each, front matter and indexes per AGENTS.md); scripts written by the `coder` agent.
- **Out of scope** (alp-sjp.6 therefore stays open): decay constants, residues, mixing-angle observables, vector/scalar/axial audit, NLO and loop effects, uncertainties and numerical inputs (alp-sjp.2 NLO, alp-sjp.4, alp-sjp.7), charge-violating sources.

## Capabilities

### New Capabilities
- `lo-finite-mixing-mass-shifts`: exact LO pole equations of the pseudoscalar mesons and the ALP under finite mixing, exact sum rules, interlacing bounds, closed approximations for isolated poles, and invariance statements.

### Modified Capabilities
<!-- None. `lo-finite-mixing-spectrum` is consumed unchanged. -->

## Impact

- New records under `memory/alp/shared/` (exact relations, exact pole equations, arrow-matrix spectrum, isolated meson pole, isolated ALP pole, and the earlier second-order record kept as an illustration); `INDEX.md` there and in the parents updated.
- New scripts under `scripts/` (coder agent), reusing `derive_alp_lo_neutral_blocks.py`, `derive_alp_single_pole_alignment.py` and `solve_alp_meson_spectrum.py`.
- Beads: note on alp-sjp.6 listing what is done and what remains; do not close it. No commit unless asked.
