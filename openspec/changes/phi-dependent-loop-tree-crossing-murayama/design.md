# Design

## Context

See proposal.md for motivation. The ALP-only crossing (`estimate-loop-tree-crossing-finite-mixing`, records `alp-loop-crossing-exact-tree-condition.md` and `alp-loop-crossing-validity.md`) gives an implicit equation $R(\epsilon)=R_{\max}$ with the exact tree denominator. Its loop numerator is truncated at the leading mixing-dependent order; that truncation is the only approximation this change removes.

Reused without re-derivation (the ALP-only results and their hypotheses are recorded in the notes):

- exact LO finite-mixing poles, projectors and norms $N(s)$ (`alp-lo-pole-solver.md`, `alp-lo-exact-pole-equations.md`, `alp-lo-finite-mixing-operator-inverse.md`);
- the dressed-frame invariants $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$ and the operator $\mathbb M_s$ (`alp-dressed-vertex-power-counting-invariants.md`, `alp-dressed-mass-*`, `alp-dressed-topological-*`);
- the exact three-meson tree with contact and exchange terms (`alp-finite-mixing-three-meson-tree.md`) and the observable-class decision D1 of the ALP-only change (`alp-loop-observable-candidates-comparison.md`);
- the implicit crossing form and the fixed independent Dalitz invariants (`alp-loop-crossing-exact-tree-condition.md`);
- the scalar-moment reduction of loop momenta (`alp-loop-momentum-reduction.md`), reused for the estimate's weights;
- the constant axial redefinition check (`alp-loop-estimate-frame-invariance.md`).

New in this change, with the radial mode:

- the radial-phase sources and their Taylor coefficients (`alp-murayama-radial-phase-mass-coefficients.md`, `radial-decay-polar-cartesian-general-sources.md`);
- the radial mass relation $m_\sigma^2 = 2\lambda f_a^2 + \tfrac{B_0 f_\pi^2}{f_a^2}(\langle\widehat S_{10}\rangle-2\langle\widehat S_{20}\rangle)$ (`murayama-sigma-two-pseudoscalar-decays.md`);
- the exact tree with radial exchange (`alp-three-pseudoscalar-radial-exchange.md`), whose polar $1/m_\sigma^2$ decoupling is not assumed here, since $m_\sigma$ is not large in the regime of interest;
- the loop-size estimate for the radial classes, which is not in the notes. No loop is computed in this change; the estimate is a prescription built from the classes and their powers of $\epsilon$.

Checks already recorded in notes that motivate the regime (to be re-derived symbolically, not assumed): in the up-only model with $f_a=\sqrt2\,m_u/\kappa_u$, the $\lambda$ term in $m_\sigma^2$ relative to the first term is $8\lambda m_u^3/(B_0 f_\pi^2\kappa_u^4)$, and the first term is $B_0 f_\pi^2\kappa_u^2/(2m_u)$, so $m_\sigma$ is set by $\kappa_u$ for $\kappa_u\gtrsim 0.3$. Numerically $m_\sigma$ and $m_a$ from the existing finite-mixing calculator agree to about 1% at $\kappa_u=0.3$ and $1$.

## Goals / Non-Goals

**Goals:**
- Derive the implicit crossing in $\varphi$-valued invariants, then specialize to up-only.
- Replace the truncated numerator by a loop-size estimate that covers all classes, with their $\epsilon$ powers stated.
- State the naturalness assumption of the ALP-only change for the neglected term, including its $\mu$-dependent and counterterm parts, without analysing absorption.

**Non-Goals:**
- Meson-only loops, which are excluded (D6), and $\widehat\Theta_1\neq0$, which is excluded (D7).
- Computing any loop contribution: no loop integral, scalar loop function, loop coefficient, $\mu$ derivative of a loop, or loop amplitude in any frame.
- Any expansion in $\epsilon$ as the main result. Leading-order inversions may be quoted only as limits.
- Any numerical value for $\lambda$, $\kappa_u$, $R_{\max}$ or $\mu$. No parameter scan is part of this change.
- A Dalitz integral. The fixed independent invariants of the ALP-only derivation are kept.
- Computing the near-degeneracy of the ALP-like state with $\eta'$, or the validity of leading-order chiral perturbation theory at $m_a\approx1$–$2.3$ GeV. Both are stated as caveats only.
- Re-deriving the ALP-only results listed above.

## Decisions

**D1. General $\varphi$-valued result first.** The sources are written as functions of the complex radial-phase field $\varphi=(f_a+\sigma)e^{ia/f_a}/\sqrt2$ (polar coordinates, as in the radial-exchange record). The DiLuzio $n_q$ and up-only forms are obtained by substitution. Alternative rejected: starting from the $n_q$ coefficients, which would need the general case recovered afterwards.

**D2. Polar fields for the crossing.** The on-shell amplitude is field-independent (`alp-polar-cartesian-amplitude-equivalence.md`), but the split into a "no $\sigma$" part and a $\sigma$ part is not. Polar fields keep $\sigma$ as an internal line with no quadratic mixing with $a$. A Cartesian cross-check is kept as a verification step only.

**D3. Loop-size estimate of the loops with an ALP-like or radial line, no $\epsilon$ truncation, no loop computed.** The numerator is a size estimate, not a computed amplitude. Each loop that contains an internal ALP-like or radial line enters with its $\epsilon$ power, its reference scale and its invariant contraction stated. The weights come from the scalar-moment magnitude prescription of `alp-loop-momentum-reduction.md`, and the cancellation conditions of `alp-loop-estimate-cancellations.md` are its stated limits. Meson-only loops are excluded (D6). The full numerator is kept, so the two-line terms and the $\epsilon$-dependence of $N(m_a^2)$ that the ALP-only record drops by its leading-order truncation are present in the estimate.

**D6. Meson-only loops excluded.** They carry no power of $f_\pi/f_a$ beyond the ordinary chiral correction, so they are not part of the mixing-dependent crossing. This is the choice of the ALP-only crossing record.

**D7. $\widehat\Theta_1=0$.** The topological invariant is set to zero, as in the Murayama model and its generalizations. With this choice $\mathcal C=B_0\widehat P_1$, and the radial-exchange tree is stated in this case only.

**D4. Naturalness instead of an absorption analysis.** As in D3 of the ALP-only change, the neglected next-order term, including its renormalization-scale-dependent and counterterm parts, is taken to have the size of the loop-size estimate. The estimate contains no $\mu$. The cancellations that would break this are listed. What counterterms or $L_i$ absorb is not analysed, since that depends on the operator basis, which the notes do not fix.

**D5. Threshold symbolic.** $R_{\max}$ stays a symbol. The implicit equation is written as $R(\text{invariants})=R_{\max}$.

## Risks / Trade-offs

- **[The estimate may overstate the true loop size, or miss cancellations] → State the cancellation conditions of `alp-loop-estimate-cancellations.md` and the naturalness assumption as its limits; the estimate is not a computed loop.**
- **[Fixed Dalitz point may not represent the channel] → Keep the existing choice; state that integrated observables are not covered by this change.**
- **[$\sigma$ light relative to the meson scale, so tree and loop poles are close to the Dalitz edge] → Keep exact propagators; no decoupling argument is used.**
- **[Leading-order chiral perturbation theory outside its range for $m_a\approx1$–$2.3$ GeV] → Stated as a caveat; no value is derived there.**
- **[$a$–$\eta'$ near-degeneracy] → Stated as a caveat; not computed.**

## Migration Plan

Not applicable. The change adds records and scripts only. Records are written under `memory/alp/shared/` and `memory/alp/three-pseudoscalar/`, and symbolic scripts under `scripts/` by the coder agent.

## Open Questions

None. The counterterm treatment is settled by D4, as in the ALP-only change.
