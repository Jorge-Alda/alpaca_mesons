# Tasks

## 1. Phi-valued invariants and radial mass (reference only, not rederived)

These results are already in the notes. Each task checks that the cited record contains the statement and links to it. No derivation is repeated.

- [x] 1.1 Confirm that the dressed sources as functions of the radial-phase field $\varphi$, and their Taylor coefficients, are stated in `memory/alp/radial-mode/couplings-and-decays/alp-murayama-radial-phase-mass-coefficients.md` and `memory/alp/radial-mode/couplings-and-decays/finite-mixing/radial-decay-polar-cartesian-general-sources.md`, with the $n_q$ case as a special case. Verify: the cited statements are present and linked; no SymPy rederivation.
- [x] 1.2 Add an index entry for the $\varphi$-valued invariants in `memory/alp/shared/INDEX.md`, pointing to the records of 1.1 and 1.2. Verify: the entry links to each record and adds no new equation.

## 2. Vertices and tree with the radial mode (reference only, not rederived)

- [x] 2.1 Confirm that the vertices with up to four pseudoscalar, ALP-like and radial legs are listed in `memory/alp/radial-mode/couplings-and-decays/finite-mixing/murayama-sigma-two-pseudoscalar-decays.md` and `memory/alp/three-pseudoscalar/finite-mixing/alp-three-pseudoscalar-radial-exchange.md`. Verify: the cited vertices are present and linked.
- [x] 2.2 Confirm that the exact three-meson tree with radial exchange is in `memory/alp/three-pseudoscalar/finite-mixing/alp-three-pseudoscalar-radial-exchange.md`, with the numerical check against `three_pseudoscalar(..., radial_mass_squared)` recorded there. Verify: the cited record contains the amplitude and its check. If a $\varphi$-valued form is not stated verbatim there, record that gap and do not rederive it.

## 3. Loop-size estimate with the radial mode

**MUST NOT compute loop contributions.** Do not evaluate any loop integral, scalar loop function, loop coefficient, $\mu$ derivative of a loop, or a loop amplitude in either frame. Only classify the loop contributions, state their $\epsilon$ powers and invariant contractions, and propose a defensible size estimate using the existing naturalness prescription. Any script in this section may check bookkeeping only (counts, powers, symmetry factors, the algebra of the estimate's prescription), never loop values.

- [x] 3.1 List the loop classes of the estimate in a record, using the classes of `memory/alp/loop-corrections/finite-mixing/three-meson-amplitudes/alp-three-meson-one-loop-estimate.md` and `memory/alp/decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-dressed-vertex-power-counting-invariants.md`, extended by the radial lines of section 2. For each class, state the $\epsilon$ power, the reference scale and the invariant contraction left symbolic. Verify: the record contains the table with no loop value; a bookkeeping script checks only the counts and powers.
- [x] 3.2 Apply the existing scalar-moment magnitude prescription of `memory/alp/loop-corrections/finite-mixing/three-meson-amplitudes/alp-loop-momentum-reduction.md` as the estimate's weight, and the cancellation conditions of `memory/alp/loop-corrections/finite-mixing/power-counting-and-validity/crossing-and-validity/alp-loop-estimate-cancellations.md` as the stated limit of the estimate. Verify: the record states the weights and the cancellation conditions with no moment evaluated.
- [x] 3.3 State the naturalness assumption of `memory/alp/loop-corrections/finite-mixing/power-counting-and-validity/loop-counting/alp-loop-naturalness-assumption.md` for the neglected next-order term, including its renormalization-scale-dependent and counterterm parts, with no absorption analysis, and the cancellation conditions of 3.2 as its limits. Verify: the record contains the assumption and the list, and contains no $\mu$ and no $\mu$ derivative.
- [ ] 3.4 Check, symbolically on the estimate's invariant contractions only, that the estimate is unchanged under a constant axial redefinition. Verify: the check acts on the contractions and does not involve any loop value.

- [x] 3.5 Estimate only the factors that multiply the flavour traces of each loop (loop factor, powers of $\epsilon$ and $N(m_a^2)$, scalar-moment scale), leaving the traces exact and symbolic, with no loop algebra. Verify: `memory/alp/radial-mode/couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md` gives the estimate formula, and every prefactor used is the scalar-moment size of that record.

## 4. Implicit crossing in general form

- [ ] 4.1 Write the implicit crossing $R = R_{\max}$ with the loop-size estimate of section 3 (loops with an ALP-like or radial line only, $\widehat\Theta_1=0$) in the numerator and the exact tree in the denominator, in $\varphi$-valued invariants, with $R_{\max}$ symbolic and fixed independent Dalitz invariants. Verify: script substitutes the ALP-only limit (no radial line) and reproduces the ALP-only implicit equation of `memory/alp/loop-corrections/finite-mixing/power-counting-and-validity/crossing-and-validity/alp-loop-crossing-exact-tree-condition.md`.
- [ ] 4.2 Check frame invariance of the estimated $R$ under a constant axial redefinition. Verify: script confirms the unchanged $R$ in both frames, on the estimate's contractions.
- [x] 4.3 Write the record of the general crossing, with its validity conditions and the two caveats of the spec, stated but not computed. Verify: the record states both caveats and contains no numerical value for them.

## 5. Up-only specialization

- [ ] 5.1 Specialize the general crossing to up-only with $f_a=\sqrt2\,m_u/\kappa_u$, keeping $\kappa_u$ and $\lambda$ analytical. Verify: script substitutes the up-only tree and the estimate's contractions in the stated invariants.
- [x] 5.2 Derive $m_\sigma$ and $m_a$ analytically in $\kappa_u$ and $\lambda$, and state the ratio of the $\lambda$ term to the first term. Verify: script checks the ratio $8\lambda m_u^3/(B_0f_\pi^2\kappa_u^4)$ symbolically.
- [ ] 5.3 Write the up-only crossing record in `memory/alp/three-pseudoscalar/` with front matter, inputs linked to their records, and an explicit list of the ALP-only results reused. Verify: each reused result is linked to its record.

## 6. Cross-checks

- [x] 6.1 Cartesian cross-check of the magnitude prescription of 3.2 with exact projectors. Verify: the check acts on the estimate's weights only, and the record states that it depends on the field choice.
- [x] 6.2 Confirm the index entries of every new record in `memory/alp/shared/INDEX.md` and `memory/alp/three-pseudoscalar/INDEX.md`, and update `memory/plan.md`. Verify: the `memory/README.md` crawl reaches every new record.

## Workflow follow-up

- Archive this change after the tasks above are complete.
