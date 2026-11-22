# Tasks

Conventions for every task: notation and invariants as in the existing records (the multiplet $\Phi_I$, $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$, $\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1$, $\mathbb M_{m_a^2}^{-1}$, matrix notation, no flavour indices); no new symbol without the user's approval (no primed polynomials); no assumed values, representatives or relations among ALP coefficients or quark masses; no parameter scan; no one-loop evaluation of any process. Scripts are written by the `coder` agent, records are atomic with front matter and index entries.

## 1. Choice of observable class

- [x] 1.1 Write the record comparing the meson pole mass, the ALP-meson mixing entry, the ALP-like pole mass and an on-shell ALP-sourced amplitude against the tests (a)-(d) of design D1. Verify with an exact-rational script that, for each candidate, the leading powers of $\epsilon$ of the tree and one-loop terms and its behaviour under a constant axial redefinition match what the record states.
- [x] 1.2 State the observable class used for the counting (an on-shell amplitude with the ALP at first order in $1/f_a$ at tree level) with the reasons, from the task 1.1 record, without evaluating any loop. Verify that the record states why the pole-mass candidates are not used, and add to the earlier selection record a pointer stating that its Dalitz and absorption analysis is not used by the estimate.

## 2. Tree quantity at finite mixing

- [x] 2.1 Write the exact tree quantity of the chosen class from the finite-mixing pole solution in $\mathbb M_{m_a^2}^{-1}$ and the invariants, without assuming relations among them. Verify against the existing exact solver with exact rational sources at several finite $\epsilon$, and record the difference from its leading expansion.

## 3. Loop size from power counting

- [x] 3.1 Derive, in matrix notation, the power counting of the dressed-frame vertices: the power of $\epsilon$ per closed ALP line, the loop factor and the dimensionless invariant contraction for each loop class, reusing `alp-loop-power-counting.md` and `alp-invariant-loop-crossing-meson-poles.md`. Verify the vertex scalings against the explicit fourth-order Lagrangian of the existing scripts, without evaluating any loop integral.
- [x] 3.2 State the naturalness assumption and the loop factor used (design D3). Verify that every symbol is defined in an existing record or flagged for approval.
- [x] 3.3 List the cancellations that make the estimate overestimate the neglected term, among classes at equal power of $\epsilon$ and from relations among the invariants. Verify each condition symbolically, keeping physical parameters unevaluated.

## 4. Invariance

- [ ] 4.1 Check, in exact arithmetic, that the invariant contractions and the loop-to-tree ratio are unchanged under a constant axial redefinition, using the frame-change tools of the existing scripts. Verify that the test fails if a non-invariant quantity (for example the off-shell tadpole entry) is substituted.

## 5. Crossing formula

- [ ] 5.1 Derive the closed formula for the $f_\pi/f_a$ at which the ratio of the leading ALP-line class to the tree equals $R_{\max}$, in the invariants, the unmixed masses, $m_a$, $f_\pi$, $B_0$, $m_0^2$ and the loop factor, for a general $R_{\max}$. Verify that the formula reproduces the stated powers of $\epsilon$ and the invariant contraction at exact rational sources.
- [ ] 5.2 Write the record with the formula, its criterion value, the naturalness assumption and the loop classes that contribute. Verify that every symbol is defined in an existing record or flagged for approval.

## 6. Validity and limits

- [ ] 6.1 State the validity conditions, the $f_a\to\infty$ limit and the near-degenerate behaviour with the exact two-level replacement. Verify each statement against the exact solver for one exact rational example per case.
- [ ] 6.2 Add pointers from `alp-loop-power-counting.md` and `alp-loop-crossing-physical-inputs.md` to the new records, marking their representative-based numbers as superseded. Verify that the indexes in `memory/alp/shared/INDEX.md` and the parents list every new record and that all links resolve.

## 7. Integration

- [ ] 7.1 Run the new and existing test suites (`derive_alp_invariant_loop_crossing`, `derive_alp_loop_power_counting`, `derive_alp_dressed_frame_unitarity`) and verify that all pass. Record in beads the status of the loop issues under alp-gii that this change touches, without closing them.

## Workflow follow-up

- Review the records with the user before any drafting; the `drafter` agent is used only on request.
- Archive the change after the review requirements are satisfied.
