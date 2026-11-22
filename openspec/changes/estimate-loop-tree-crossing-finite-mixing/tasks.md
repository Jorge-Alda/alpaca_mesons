# Tasks

Conventions for every task: notation and invariants as in the existing records (the multiplet $\Phi_I$, $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$, $\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1$, $\mathbb M_{m_a^2}^{-1}$, matrix notation, no flavour indices); no new symbol without the user's approval; no assumed values or representatives for ALP coefficients; no parameter scan. Scripts are written by the `coder` agent, records are atomic with front matter and index entries.

## 1. Choice of observable

- [x] 1.1 Write the record comparing the meson pole mass, the ALP-meson mixing entry, the ALP-like pole mass and an on-shell ALP-sourced amplitude against the tests (a)-(d) of design D1. Verify with an exact-rational script that, for each candidate, the leading powers of $\epsilon$ of the tree and one-loop terms and its behaviour under a constant axial redefinition match what the record states.
- [ ] 1.2 State the selected observable, or that none passes, with the reasons. Verify that the record names, for the selected observable, the free parameters that could absorb its loop shift and shows that the chosen coefficient is not among them.

## 2. Tree quantity at finite mixing

- [ ] 2.1 Write the exact tree quantity of the selected observable from the finite-mixing pole solution in $\mathbb M_{m_a^2}^{-1}$ and the invariants. Verify against the existing exact solver with exact rational sources at several finite $\epsilon$, and record the difference from its leading expansion.

## 3. One-loop quantity

- [ ] 3.1 Derive the one-loop correction of the selected observable by loop function and by number of ALP lines, reusing the vertices of `alp-loop-power-counting.md` and `alp-invariant-loop-crossing-meson-poles.md`, with the power of $\epsilon$ of each class. Verify by comparing each coefficient with the exact one-loop expression built from the explicit Lagrangian to fourth order.
- [ ] 3.2 Separate the coefficient of the non-analytic logarithm from the scale-dependent remainder (design D3). Verify that the logarithm coefficient is independent of $\mu$ and that the remainder carries the dependence, with a symbolic check.
- [ ] 3.3 State which parts are absorbed by the sources or by NLO operators. Verify that the statement matches the structure of the local operators listed in `alp-dressed-mass-*` and `memory/nlo/` records, and name any part that cannot be classified.

## 4. Invariance

- [ ] 4.1 Check, in exact arithmetic, that the loop-to-tree ratio and the logarithm coefficient are unchanged under a constant axial redefinition, using the frame-change tools of the existing scripts. Verify that the test fails if a non-invariant quantity (for example the off-shell tadpole entry) is substituted.

## 5. Crossing formula

- [ ] 5.1 Derive the closed formula for the $f_\pi/f_a$ at which the ratio equals $R_{\max}$, in the invariants, the unmixed masses, $m_a$, $f_\pi$, $B_0$ and $m_0^2$, for a general $R_{\max}$. Verify against the exact ratio at exact rational sources at the formula's value, to the stated error.
- [ ] 5.2 Write the record with the formula, its criterion value, the dependence on $\mu$, and the loop classes that contribute. Verify that every symbol is defined in an existing record or flagged for approval.

## 6. Validity and limits

- [ ] 6.1 State the validity conditions, the $f_a\to\infty$ limit and the near-degenerate behaviour with the exact two-level replacement. Verify each statement against the exact solver for one exact rational example per case.
- [ ] 6.2 Add pointers from `alp-loop-power-counting.md` and `alp-loop-crossing-physical-inputs.md` to the new records, marking their representative-based numbers as superseded. Verify that the indexes in `memory/alp/shared/INDEX.md` and the parents list every new record and that all links resolve.

## 7. Integration

- [ ] 7.1 Run the new and existing test suites (`derive_alp_invariant_loop_crossing`, `derive_alp_loop_power_counting`, `derive_alp_dressed_frame_unitarity`) and verify that all pass. Record in beads the status of the loop issues under alp-gii that this change touches, without closing them.

## Workflow follow-up

- Review the records with the user before any drafting; the `drafter` agent is used only on request.
- Archive the change after the review requirements are satisfied.
