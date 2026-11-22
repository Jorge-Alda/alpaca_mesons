# Design

## Context

See proposal.md for motivation. Existing material this change builds on and does not redo:

- exact LO finite-mixing poles and pole vectors (`alp-lo-pole-solver.md`, `alp-lo-exact-pole-equations.md`, `alp-lo-arrow-matrix-spectrum.md`);
- the dressed frame, in which the ALP enters only through $\widehat M(x)$, $\widehat\Theta(x)$ with $K=\mathbf 1$ and the invariants $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$ (`alp-dressed-frame-derivative-coupling-removal.md`, the dressed-mass and dressed-topological records);
- the loop classes of the meson pole masses written with these invariants, and the failure of the tadpole-ratio crossing under constant axial redefinitions (`alp-invariant-loop-crossing-meson-poles.md`);
- the ALP-line loop power counting (`alp-loop-power-counting.md`) and the one-loop LSZ structure (`loop-corrections-to-alp-meson-lsz-mixing.md`);
- the observable-candidate comparison of task 1.1 (`alp-loop-observable-candidates-comparison.md`), which stands.

Facts established in earlier work that constrain the design (each is recorded in the notes, not new): at leading order the meson pole-mass shift from tree-level mixing and the ALP-tadpole class are both $\propto\epsilon^2$, so their ratio is $\epsilon$-independent; the exact poles depend nonlinearly on $\epsilon$; ALP lines carry $1/f_a$ per line; meson-only loops give the usual $\epsilon$-independent chiral expansion.

## Goals / Non-Goals

**Goals:**
- State, with written arguments, the observable class for which the loop-to-tree ratio can depend on $\epsilon$.
- Derive the closed formula for the crossing in the invariants, with the tree part exact in $\epsilon$ and the loop part given by vertex power counting.

**Non-Goals:**
- Any condition defined by tree-level quantities alone (mixing angle, partial-wave unitarity).
- The evaluation of any loop integral or amplitude, and the analysis of what counterterms or $L_i$ absorb.
- Numerical scans, fits, a choice of ALP source, or values for the invariants or the quark masses; assumptions such as $m_u=m_d$ or $\langle\lambda_3\widehat P_1\rangle=0$.
- NLO constants and resonances; the near-degenerate regime beyond stating where the formula fails.

## Decisions

**D1. Observable class chosen by an explicit rule, not by convenience.** The candidates are compared against four tests: (a) invariant under constant chiral redefinitions; (b) the ALP enters at tree level at first order in $1/f_a$, so that ALP-line loops carry extra powers of $\epsilon$ relative to the tree; (c) the correction is not absorbed into a free parameter; (d) no fitted input is needed. The task 1.1 record shows that the meson pole mass fails (b) and the ALP-like pole mass fails (c). The class used for the counting is an on-shell amplitude with the ALP at first order in $1/f_a$ at tree level. Alternative considered: a specific process (such as $A\to3\pi^0$); rejected because it turns the estimate into a one-loop process calculation.

**D2. Tree quantity exact, loop side by vertex power counting.** The tree amplitude uses the exact pole solution (pole masses and pole vectors as functions of $\epsilon$), via $\mathbb M_{m_a^2}^{-1}$ and the invariants. The loop side is the size of each loop class, obtained from the scaling of the dressed-frame vertices: $\epsilon$ per ALP line, the loop factor $m^2/(16\pi^2f_\pi^2)$, and a dimensionless matrix contraction of the invariants. No loop integral is evaluated. Alternative considered: computing the coefficients of the loop functions of a chosen amplitude; rejected as defeating the purpose.

**D3. Naturalness instead of an absorption analysis.** The neglected next-order term is taken to have the size of the loop. This avoids deciding which parts counterterms absorb, which depends on the operator basis beyond what the notes fix. The price is an explicit assumption that no cancellation occurs, with the cancellations listed (among classes of equal power of $\epsilon$, or from relations among the invariants). Alternative considered: using only the part not absorbable by local operators; rejected because the task 1.2 analysis showed it requires evaluating loop functions and an order-by-order reading of the derivative expansion.

**D4. Frame invariance by explicit check.** The matrix contractions and the ratio are checked under the constant axial redefinition used in the existing scripts, in exact arithmetic. No off-shell entry is used.

**D5. Analytic result first, scripts verify.** The derivation is written in sympy with the closed formula as the output. Exact rational sources are used only to test vertex scalings and identities. No parameter scan is run or presented.

**D6. Criterion value.** The crossing is defined as ratio $=R_{\max}$ with $R_{\max}$ of order one, and the formula is written for general $R_{\max}$ so that the dependence on this choice is explicit.

## Risks / Trade-offs

- [Naturalness overestimates when classes cancel] → state the cancellation conditions and that the estimate does not apply there.
- [Several loop classes at the same power of $\epsilon$] → report each class and the sum, and the condition under which the sum vanishes.
- [Near an ALP-meson degeneracy the closed formula fails] → state the condition and give the exact two-level replacement.
- [Pole vectors in a frame with exact mixing can be singular for an unmixed pole] → treat the unmixed pole separately.
- [The estimate is of the size of the correction, not its value] → state that it is an estimate of where the tree-level approximation stops being controlled, not a prediction of the correction.

## Open Questions

- None that change the specs or tasks. CP-violating classes ($\widehat P_2$, $\widehat\Theta_2$) are included as separate matrix contractions.
