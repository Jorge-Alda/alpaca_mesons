# Design

## Context

See proposal.md for motivation. Existing material this change builds on and does not redo:

- exact LO finite-mixing poles and pole vectors (`alp-lo-pole-solver.md`, `alp-lo-exact-pole-equations.md`, `alp-lo-arrow-matrix-spectrum.md`);
- the dressed frame, in which the ALP enters only through $\widehat M(x)$, $\widehat\Theta(x)$ with $K=\mathbf 1$ and the invariants $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$ (`alp-dressed-frame-derivative-coupling-removal.md`, the dressed-mass and dressed-topological records);
- the loop classes of the meson pole masses written with these invariants, and the failure of the tadpole-ratio crossing under constant axial redefinitions (`alp-invariant-loop-crossing-meson-poles.md`);
- the ALP-line loop power counting (`alp-loop-power-counting.md`) and the one-loop LSZ structure (`loop-corrections-to-alp-meson-lsz-mixing.md`).

Facts established in earlier work that constrain the design (each is recorded in the notes, not new): at leading order the meson pole-mass shift from tree-level mixing and the ALP-tadpole class are both $\propto\epsilon^2$, so their ratio is $\epsilon$-independent; the exact poles depend nonlinearly on $\epsilon$; ALP-line loops carry $1/f_a$ per ALP line; meson-only loops give the usual $\epsilon$-independent chiral expansion.

## Goals / Non-Goals

**Goals:**
- Select, with written arguments, an observable for which the loop-to-tree ratio can depend on $\epsilon$ and is not absorbed into a free parameter.
- Derive the closed formula for the crossing in the invariants, exact in the tree part.

**Non-Goals:**
- Any condition defined by tree-level quantities alone (mixing angle, partial-wave unitarity).
- Numerical scans, fits, a choice of ALP source, or values for the invariants.
- NLO constants and resonances; the near-degenerate regime beyond stating where the formula fails.

## Decisions

**D1. Observable chosen by an explicit decision rule, not by convenience.** The candidates (meson pole mass; ALP-meson mixing entry; ALP-like pole mass; on-shell ALP-sourced amplitude) are compared against four tests: (a) invariant under constant chiral redefinitions; (b) the ALP enters at tree level at first order in $1/f_a$, so that ALP-line loops carry extra powers of $\epsilon$ relative to the tree; (c) the loop correction is not absorbed into a free parameter; (d) no fitted input is needed. Working hypotheses to be confirmed or refuted in task 1, not assumed:
- the meson pole mass fails (b) at leading order, since tree-level mixing and the tadpole scale alike;
- the ALP-like pole mass fails (c), because $m_{a,0}^2$ is a free parameter that absorbs its loop shift;
- an on-shell ALP-sourced amplitude, normalized to its exact tree value, passes (a), (b) and (d), and passes (c) only through the part of the loop that is not absorbed (see D3).
The outcome of the comparison, including a negative one (no candidate passes), is a deliverable. Alternative considered: choosing the candidate that gives the smallest crossing; rejected because it selects the answer before the argument.

**D2. Tree quantity exact, loop quantity at one loop.** The tree amplitude uses the exact pole solution (pole masses and pole vectors as functions of $\epsilon$), via $\mathbb M_{m_a^2}^{-1}$ and the invariants. The one-loop correction is organized as in the existing loop record: coefficients of independent loop functions, tagged by the number of ALP lines. Alternative considered: working at leading order in $\epsilon$ throughout; rejected because the question concerns $f_\pi/f_a\gtrsim1$ and the previous review showed that exact pole dependence changes the ratio.

**D3. Use the part of the loop that cannot be absorbed.** A tadpole proportional to $A_0$ in a quantity that has the form of an LO operator can be absorbed into the sources (renormalization of $\widehat S_n$, $\widehat P_n$) or into NLO constants. The scheme-independent content is the coefficient of the non-analytic logarithm, $\ln(m^2/\mu^2)$, whose coefficient equals the running of the corresponding counterterm. The estimate is built from this coefficient, and the analytic remainder (the part that depends on $\mu$ and the subtraction) is reported separately. Alternatives considered: using the full $A_0$ (rejected: depends on $\mu$ and on $\overline{\rm MS}$ conventions); comparing with an experimental precision (rejected: needs a numerical input that does not belong to an analytic derivation).

**D4. Frame invariance by explicit check.** Every ratio is checked under the constant axial redefinition used in the existing scripts, with an exact-rational test, before being used. The off-shell entry is never used on its own; on-shell combinations are used.

**D5. Analytic result first, scripts verify.** The derivation is written in sympy, with the closed formula as the output. Exact rational sources are used only to test identities. No parameter scan is run or presented.

**D6. Criterion value.** The crossing is defined as ratio $=1$ for the selected observable, stated as an order-one estimate, with the formula written for a general criterion $R_{\max}$ so that the dependence on this choice is explicit.

## Risks / Trade-offs

- [No candidate satisfies (a)-(d)] → report this as the result with the reasons, and state which quantity would have to be defined (for example, an amplitude with a fixed UV completion that fixes the invariants).
- [The non-absorbable coefficient depends on how the invariants are related by a UV model] → state the relations assumed, and keep the formula in the independent invariants.
- [Several loop classes contribute at the same power of $\epsilon$ and may cancel] → report each class and the sum, and check the cancellation condition explicitly.
- [Near an ALP-meson degeneracy the closed formula fails] → state the condition and give the exact two-level replacement.
- [Pole vectors in a frame with exact mixing can be singular for an unmixed pole] → treat the unmixed pole separately, as in the existing class-A analysis.
- [The loop part of the one-loop amplitude involves derivative vertices of charged mesons] → use the existing vertices from the loop record, not new ones.

## Open Questions

- Whether CP-violating classes ($\widehat P_2$, $\widehat\Theta_2$) should be included in the crossing or reported separately (does not change the approach).
