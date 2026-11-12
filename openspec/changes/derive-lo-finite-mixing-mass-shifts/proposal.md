# Proposal

## Why

Beads issue alp-sjp.6 asks how finite ALP-meson mixing changes the physical meson and ALP masses. The LO spectrum is solved numerically (`lo-finite-mixing-spectrum`), and the weak-mixing record gives implicit or first-order pole formulas. There is no closed analytical statement of the mass shifts of the individual mesons and of the ALP, no expression in terms of the redefinition-invariant couplings, and no exact sum rules to test numerical spectra against. This change supplies that part of alp-sjp.6, at leading order and analytically.

## What Changes

Scope is leading order (the $O(p^2)$ chiral Lagrangian of Alda:2026nrd with the full $K$, $M^2$ blocks of the [finite-mixing record](../../../memory/alp/shared/finite-alp-meson-mixing.md)), tree level, neutral and charged mesons, and the ALP. Conventions are those of Alda:2026nrd. The deliverable is analytical expressions with symbolic proofs; scripts verify them (exact rationals and the existing solver are used only as checks, no numerical inputs, no scans).

- **Exact relations**: all-orders statements that need no expansion: zero shift of the charged mesons (decoupling theorem), exact $K^0$ pole at $m_{K^0,0}^2$ (one combination), trace sum rule $\sum_n m_n^2=\operatorname{tr}(K^{-1}M^2)$ and product rule $\prod_n m_n^2=\det M^2$ (with $\det K=1$), each written with the explicit analytic $K^{-1}$ and in the sources $c_q^A$, $c_q^P$, $c_G$, $\chi_0$, $m_0^2$.
- **Neutral meson shifts**: explicit closed forms in $\epsilon=f_\pi/f_a$ for $\pi^0$, $\eta$, $\eta'$ and the mixed $K^0$ combination: leading $O(\epsilon^2)$ shift $\delta m_j^2=D_{ja}(m_{j,0}^2)^2/(m_{j,0}^2-m_{a,0}^2)$ and the next order, written through the redefinition-invariant combinations $c_q^{AP}=c_q^A-c_q^P$ and $c_{GA}=c_G+\operatorname{tr}c_q^A$ instead of the individual couplings. Explicit $\eta$-$\eta'$ rotation in the isospin limit; general $m_q$ through the eigenvectors of $M^2_{PP}$. Prove invariance of each shift under the constant chiral redefinitions of alp-sjp.1.
- **ALP mass shift**: $m_a^2-m_{a,0}^2$ at $O(\epsilon^2)$ and $O(\epsilon^4)$, with its decomposition into the kinetic, anomaly and meson-exchange parts, and its form in the pole-aligned frame (exact implicit equation from the aligned $D'_{aa}$).
- **Pole alignment**: use the single-pole-aligned coordinates to express each shift as the zero of one scalar function, to organize the expansion, and to separate direct level repulsion from the kinetic and anomaly renormalization of the bare ALP term.
- **Limits and domain**: no-ALP limit ($\epsilon\to0$, masses of the unmixed chiral Lagrangian), weak-mixing agreement with Alda:2026nrd and `alp-lo-weak-mixing-expansion.md`, level-repulsion sign rule, and the failure condition near degeneracy with the exact two-level formula as replacement.
- Memory records (one concept each, front matter and indexes per AGENTS.md); scripts written by the `coder` agent.
- **Out of scope** (alp-sjp.6 therefore stays open): decay constants, residues, mixing-angle observables, vector/scalar/axial audit, NLO and loop effects, uncertainties and numerical inputs (alp-sjp.2 NLO, alp-sjp.4, alp-sjp.7), charge-violating sources.

## Capabilities

### New Capabilities
- `lo-finite-mixing-mass-shifts`: analytical LO shifts of the pseudoscalar meson and ALP masses under finite mixing, exact sum rules, selection rules for zero shifts, and invariance statements.

### Modified Capabilities
<!-- None. `lo-finite-mixing-spectrum` is consumed unchanged. -->

## Impact

- New records under `memory/alp/shared/` (exact relations, meson shifts, ALP shift); `INDEX.md` there and in the parents updated.
- New scripts under `scripts/` (coder agent), reusing `derive_alp_weak_mixing` verification code, `derive_alp_lo_neutral_blocks.py`, `derive_alp_single_pole_alignment.py` and `solve_alp_meson_spectrum.py`.
- Beads: note on alp-sjp.6 listing what is done and what remains; do not close it. No commit unless asked.
