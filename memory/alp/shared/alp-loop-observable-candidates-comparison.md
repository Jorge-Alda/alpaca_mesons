---
title: Candidate observables for the loop-to-tree comparison at finite f_pi/f_a, tested against invariance, order in 1/f_a, absorption and inputs
status: mixed
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd, Gasser:1983yg, Kaiser:1998ds]
provenance: analytical-derivation
---

# Candidate observables for the loop-to-tree comparison at finite $f_\pi/f_a$

This record compares four candidate observables against the four tests of decision D1 of the OpenSpec change `estimate-loop-tree-crossing-finite-mixing`. The powers of $\epsilon$ and the behaviour under a constant axial redefinition are verified by an exact script. The statements on absorption (test (c)) are partly verified and partly conjectures, labelled below. **This record does not select an observable**; the selection is a separate step (task 1.2).

## Tests (design D1)

- (a) The quantity is invariant under constant axial redefinitions $g_R=g_L^\dagger=e^{ix\kappa}$.
- (b) The ALP enters the tree term at first order in $1/f_a$. Then loops with ALP lines carry extra powers of $\epsilon$ relative to the tree, and the loop-to-tree ratio can depend on $\epsilon$.
- (c) The loop correction is not absorbed into a free parameter.
- (d) No fitted input is needed.

## Conventions and units

- Units GeV. Notation of [@Alda:2026nrd]: $U=\exp(i\Phi/f_\pi)$, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, $N_f=3$, $x=a/f_a$, $\epsilon=f_\pi/f_a$. Fields $\Phi_I=(P_i,a)$ with inverse propagator $D(s)=sK-M^2$ ([pole solver](alp-lo-pole-solver.md)).
- Dressed frame ([derivative-coupling removal](alp-dressed-frame-derivative-coupling-removal.md)): $K=\mathbf 1$, and the ALP enters only through $\widehat M(x)=\sum_n\frac{x^n}{n!}(\widehat S_n+i\widehat P_n)$ and $\widehat\Theta(x)=\sum_n\frac{x^n}{n!}\widehat\Theta_n$.
- $\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1$ (written $C$ in the [meson-pole loop record](alp-invariant-loop-crossing-meson-poles.md)), and $Y=\mathbb M^{-1}_{m_a^2}(\mathcal C)$, the solution of $\mathbb M_{m_a^2}(Y)=\mathcal C$.
- $A$ denotes the ALP-like pole and $m_A^2$ its mass. $m_{n,0}^2$ denotes the unmixed mass of meson $n$ and $m_{a,0}^2$ the bare ALP mass.
- Loop functions $A_0$ and $\bar B$ as in the [meson-pole loop record](alp-invariant-loop-crossing-meson-poles.md). The symbol $\bar B$ is still awaiting approval there. One-loop pole shifts are split into coefficients of independent loop functions. A "class" is one such coefficient.
- "Leading power $p$": all coefficients below $\epsilon^p$ vanish identically in the invariants, and the coefficient of $\epsilon^p$ is a nonzero polynomial in them.

## Assumptions

1. LO $U(3)$ Lagrangian in the dressed frame, one loop, as in the [meson-pole loop record](alp-invariant-loop-crossing-meson-poles.md) (assumptions 1–3 there). Flavour-diagonal invariants ($\widehat S_0=m_q$, $\widehat S_1=0$, $\widehat P_0=0$). The invariants are kept as symbols, and no value or relation among them is assumed.
2. $m_{a,0}^2$ is held fixed as $\epsilon$ varies, so $m_A^2$ is $\epsilon$-independent at leading order. If instead $m_{a,0}^2$ is not tuned and $m_A^2\propto\epsilon^2$, the powers below change. That regime is **not verified** here.
3. Poles are non-degenerate, and $m_A^2$ is away from the unmixed neutral masses.
4. Meson-sector numbers in the script are exact rationals used only for checks. They have $m_u=m_d$ (so that $\pi^0$ is exactly along $\lambda_3$), and $m_s$, $m_0^2$ are chosen so that the unmixed $\eta,\eta'$ masses are rational.
5. "CP-conserving" means that all CP-odd invariants vanish: $\widehat S_n$ with $n$ odd, and $\widehat P_n$, $\widehat\Theta_n$ with $n$ even.

## Inputs

- Classes A, B, C of the meson pole masses and their invariance at one loop: [meson-pole loop record](alp-invariant-loop-crossing-meson-poles.md).
- Vertex scaling $V\propto f_\pi^{2-n_P}f_a^{-n_a}$ and the tadpole and bubble combinatorics: [power-counting record](alp-loop-power-counting.md), sections L1–L3.
- Tree-level on-shell equality between frames, and the dressed frame: [derivative-coupling removal](alp-dressed-frame-derivative-coupling-removal.md).
- Pole vectors and residues: [pole solver](alp-lo-pole-solver.md).
- Beyond-LO derivative operators with singlet sources: [singlet-derivative record](alp-singlet-derivative-terms-energy-growth.md), from [@Kaiser:1998ds].
- Polynomial (subtraction-dependent) parts of the loops are absorbed by next-order local counterterms [@Gasser:1983yg].

## Candidates

1. **Meson pole mass $m_n^2$** ($\pi^0,\eta,\eta'$, charged and $K$ poles). The tree term is the exact pole. Its ALP-induced part is $m_n^2-m_{n,0}^2$.
2. **ALP–meson mixing entry** $D_{na}(s)$ (off shell), and the $(n,a)$ entry of the one-loop self-energy. The script uses $n=\pi^0$.
3. **ALP-like pole mass $m_A^2$**, with one-loop shift $v^{(A)\mathsf T}[\delta M^2-m_A^2\delta K]v^{(A)}$.
4. **On-shell ALP-sourced amplitude.** The case checked explicitly is $\mathcal M(A\to\pi^0\pi^0\pi^0)$, with CP-conserving sources and $m_u=m_d$ in the meson sector. The ALP sources may break isospin, and the amplitude is linear in that breaking. With CP-conserving sources all cubic vertices vanish ([power-counting record](alp-loop-power-counting.md), L1), so the tree term is the contact vertex contracted with the pole vectors. The general power counting for an amplitude with one external $A$ follows from the vertex scaling of L1 ([power-counting record](alp-loop-power-counting.md)); only $A\to3\pi^0$ is checked explicitly.

## Result: leading powers of $\epsilon$ (exact, verified)

The columns of ratio powers give the power of $\epsilon$ of loop/tree. For candidate 1 there are two possible tree references: (i) the ALP-induced tree part, and (ii) the full tree $m_n^2$.

| candidate | tree | loop class | loop power | ratio power |
|---|---|---|---|---|
| 1, $\pi^0,\eta,\eta'$ | full $\epsilon^0$; ALP part $\epsilon^2$ | meson lines, $\epsilon$-independent part | 0 | (i) $-2$, (ii) 0 |
| | | meson lines, $\epsilon$-dependent part (class C) | 2 | (i) 0, (ii) 2 |
| | | ALP line, $A_0(m_A^2)$ (class A) | 2 | (i) 0, (ii) 2 |
| | | ALP pair, $\bar B(m_n^2;m_A,m_A)$ (class B) | 4 (needs CP-odd invariants, zero for CP-conserving sources) | (i) 2, (ii) 4 |
| | | ALP–meson bubble | 6 | (i) 4, (ii) 6 |
| | | meson–meson bubble | 8 | — |
| | | loop-induced one-point terms | zero through $\epsilon^4$ | — |
| 1, charged poles | ALP part exactly zero | ALP line, $A_0(m_A^2)$ | 2 | (ii) 2 |
| 2, $D_{\pi^0a}$ | 1 | meson lines | 1 | 0 |
| | | ALP line, $A_0(m_A^2)$ | 3 | 2 |
| | | bubbles | 5 | 4 |
| 3, $m_A^2$ | ALP part $\epsilon^2$ (plus $m_{a,0}^2$) | meson lines, $A_0(m_P^2)$ | 2 | 0 |
| | | ALP line, $A_0(m_A^2)$ | 4 | 2 |
| | | bubbles $(A,P)$ / $(A,A)$ / $(P,P')$ | 4 / 6 / 6, all zero for CP-conserving sources | 2 / 4 / 4 |
| 4, $A\to3\pi^0$ | 1 (the $\epsilon^0,\epsilon^2,\epsilon^4$ terms vanish) | meson lines, $A_0(m_P^2)$ and $\bar B(s;m_P,m_P)$ | 1 | 0 |
| | | ALP line, $A_0(m_A^2)$ | 3 | 2 |
| | | $\bar B(s;m_A,m_P)$ | 3 | 2 |
| | | $\bar B(s;m_A,m_A)$ | 5 | 4 |
| | | external legs, $A_0(m_P^2)$ / $A_0(m_A^2)$ | 1 / 3 | 0 / 2 |

Further exact results:

- **Candidate 1, form of class A.** At $O(\epsilon^2)$ the coefficient of $A_0(m_A^2)$ for $\pi^0$, $\eta$ and $\eta'$ is $B_0\big\langle\lambda_n^2\,\big[-\frac{\epsilon^2}{2f_\pi^2}\big(\widehat S_2+4Y\widehat P_1-4Y^2m_q\big)\big]\big\rangle$, with the same diagonal matrix in brackets for the three poles. With $m_u=m_d$ this is three conditions on two unknowns, so the agreement is a nontrivial check. The bracket has the form of a shift of $m_q$ in the LO meson mass term. (The script evaluates $Y$ at $m_{a,0}^2$, which equals $m_A^2$ at this order.) For the charged pairs the remainder beyond this form is exactly $\frac{\epsilon^2}{2f_\pi^2}m_A^2\langle[Y,\lambda_n]^2\rangle$, which agrees with class A of the [meson-pole loop record](alp-invariant-loop-crossing-meson-poles.md).
- **Candidate 1, class B.** It does not vanish when only $\widehat P_2=\widehat\Theta_2=0$. The vertex $V_{\pi^0AA}$ then starts at $\epsilon^4$ through $\operatorname{tr}\widehat S_3$, which gives $b_{\pi^0,aa}=O(\epsilon^8)$. It vanishes only when all CP-odd invariants vanish. This refines, beyond the order given there, the statement of the meson-pole loop record that class B needs $\widehat P_2$, $\widehat\Theta_2$.
- **Candidate 4, tree.** The script verifies, through $\epsilon^5$,
$$
\mathcal M(A\to\pi^0\pi^0\pi^0)=\frac{\epsilon B_0}{f_\pi^2}\big\langle(2Ym_q-\widehat P_1)\lambda_3\big\rangle+O(\epsilon^3),
$$
  with the sign convention $\sum_kv^{(A)}_k\lambda_k=2\epsilon Y$ of the meson-pole loop record. For $m_u=m_d$ the script reduces the coefficient to $\frac{B_0}{f_\pi^2}\,m_A^2\langle\lambda_3\widehat P_1\rangle/(m_{\pi^0}^2-m_A^2)$, in which $\widehat\Theta_1$ drops out. The coefficient is generically nonzero, so the leading power is exactly 1.

## Result: behaviour under the constant axial redefinition (exact, verified)

- **Candidate 1.** Every on-shell coefficient (pole mass, each loop function, each bubble) is equal in the dressed frame and in the $\kappa$ frame, exactly through $\epsilon^4$, for $\pi^0,\eta,\eta'$ and the charged poles. Both frames are computed directly with symbolic $\epsilon$.
- **Candidate 2.** At tree level $D^{(\kappa)}_{\pi^0a}(s)-D_{\pi^0a}(s)=-\epsilon\langle\lambda_3\kappa\rangle(s-m_{\pi^0,0}^2)$, with $\kappa$ symbolic. This vanishes at the unmixed mass but not at the pole, where it starts at $\epsilon^3$. At one loop the frame difference of the $(\pi^0,a)$ entry at $s=m_{\pi^0}^2$ is nonzero, starting at $\epsilon^1$ (meson lines) and $\epsilon^3$ ($A_0(m_A^2)$). The entry is therefore frame dependent even on shell.
- **Candidate 3.** The on-shell coefficients are equal in both frames, exactly.
- **Candidate 4.** In the neutral sector the $\kappa$ map is exactly linear: a shift of each neutral field proportional to $\epsilon\,\kappa\,a$, with no residual nonlinear terms ($\kappa$ symbolic). The tree amplitude is equal in both frames, exactly. Because the map is linear, all neutral vertices are congruent, so the one-loop part with neutral internal lines is equal on shell. The part with charged internal lines was **not computed**. Its invariance rests on the field-redefinition assumption of section L6 of the [power-counting record](alp-loop-power-counting.md), which is checked there only for two-point functions.

## Result: absorption into free parameters (test (c))

- **Candidate 3 (exact).** $\partial m_A^2/\partial m_{a,0}^2=(v^{(A)}_a)^2=1+O(\epsilon^2)\neq0$. Hence any one-loop shift $\delta m_A^2$ is removed by $\delta m_{a,0}^2=-\delta m_A^2/(v^{(A)}_a)^2$.
- **Pure-ALP tadpole on any vertex (exact).** For the dressed Lagrangian with general Taylor coefficients, contracting one pair of pure-ALP fields with $\langle aa\rangle=-A_0(m_A^2)$ gives $-\frac{A_0(m_A^2)}{2f_a^2}\,\partial_x^2\mathcal L$. Verified for all terms with up to four fields and source Taylor order up to 6. This is the same as the shifts
$$
\widehat M_n\to\widehat M_n-\frac{A_0(m_A^2)}{2f_a^2}\widehat M_{n+2},\qquad
\widehat\Theta\to\widehat\Theta-\frac{A_0(m_A^2)}{2f_a^2}\widehat\Theta'',
$$
  plus a term containing only the ALP, $\frac{f_\pi^2m_0^2}{3f_a^2}A_0(m_A^2)\,(\widehat\Theta')^2$.
- **Candidate 1 (verified form, conjectured interpretation).** For neutral poles, class A has the form of a shift of $m_q$ (see above), so it can be absorbed into the quark masses. The charged remainder $\propto m_A^2\langle[Y,\lambda_n]^2\rangle$ does not have that form. *Conjecture:* classes B and C are not of the form of an LO parameter shift. This has not been checked.
- **Candidate 4 (partly conjecture).**
  - The pure-ALP part of the $A_0(m_A^2)$ class has the form of a shift of $\widehat P_1$ by $-\frac{A_0(m_A^2)}{2f_a^2}\widehat P_3$ (verified identity above). $\widehat P_1$ is a free coefficient, so this part is absorbed, including its dependence on $\mu$. Whether the finite $m_A^2\ln m_A^2$ part counts as absorbed depends on whether $\widehat P_1$ is fixed independently, for example by matching at a stated scale. For a single ALP with fixed mass, that part cannot be separated numerically from a change of the renormalized $\widehat P_1$ at scale $\mu$.
  - *Conjecture:* the parts that come from the meson components of $A$ and from the ALP components of the meson lines (coefficients with $Y$ and $1/(m_A^2-m_P^2)$), and the $A_0(m_P^2)$ logarithms, which are not analytic in $m_q$, are not of that form.
  - *Conjecture:* in the decay region the invariant masses satisfy $s_{ij}<(m_A+m_P)^2$, so $\bar B(s_{ij};m_A,m_P)$ and $\bar B(s_{ij};m_A,m_A)$ are real and analytic in $s_{ij}$. Their Taylor polynomials could then be reproduced by beyond-LO derivative operators with $x$-dependent coefficients ([singlet-derivative record](alp-singlet-derivative-terms-energy-growth.md)). Behaviour that is not analytic in $s$ appears only in crossed channels above the $AP$ threshold.
  - Remark: the tree term of any single ALP observable is linear in free invariants. "Absorbed into a free parameter" therefore has content only if it is read as "absorbed by a local counterterm polynomial in the Lagrangian parameters". Fixing this reading of test (c) is left to task 1.2.
  - Follow-up (task 1.2): the reading of test (c), the absorbed parts and the selected coefficients are in the [selection record](alp-loop-observable-selection-a-to-3pi0.md). There, the bubble conjecture above is resolved: local operators with at most four derivatives reproduce only a constant and a multiple of $\sum s_{ij}^2$, not the full $s_{ij}$ dependence of the bubbles.

## Comparison against tests (a)–(d)

| candidate | (a) | (b) | (c) | (d) |
|---|---|---|---|---|
| 1, $m_n^2$ | pass (exact) | fails: the ALP enters the tree at $\epsilon^2$; ratio to the ALP part is $\epsilon^0$ for classes A, C; only the CP-odd class B gives $\epsilon^2$ | partial: neutral class A has the form of an $m_q$ shift; charged remainder, B, C not (conjecture) | needs only the LO meson parameters and invariants, but the $m_q$-shift form of class A makes the ratio depend on how $m_q$ is fixed |
| 2, $D_{na}$ | fails, off shell at tree and at the pole at one loop (exact) | formally passes (tree $\epsilon$, ALP line $\epsilon^3$) | not meaningful, since a redefinition changes the entry | pass |
| 3, $m_A^2$ | pass (exact) | fails: ALP-induced tree $\epsilon^2$, meson lines $\epsilon^2$ | fails (exact: absorbed by $m_{a,0}^2$) | fails if $m_A$ is derived (needs the free $m_{a,0}^2$); if $m_A$ is an input, the quantity predicts nothing |
| 4, $A\to3\pi^0$ | pass at tree (exact) and for neutral loops; charged loops rest on the L6 assumption | pass: tree $\epsilon$, ALP-line classes give ratio $\epsilon^2$ ($A_0(m_A^2)$, $\bar B(s;m_A,m_P)$) and $\epsilon^4$ ($\bar B(s;m_A,m_A)$); meson-only classes $\epsilon^0$ | partial: the pure-ALP tadpole is a $\widehat P_1$ shift (exact); mixing-induced pieces and meson logarithms are not (conjecture); bubbles analytic in the decay region (conjecture) | pass: pole masses, $f_\pi$, LO meson parameters and invariants |

Remark on candidate 1 against the full tree. Every ALP class has ratio $\epsilon^2$ to the full $m_n^2$. That ratio compares ALP loops with the ALP-independent part of the meson mass, not with the ALP effect at tree level. In class A the prefactor is $\epsilon^2/f_\pi^2=1/f_a^2$ and $Y$ does not contain $f_\pi$, so this ratio is controlled by $f_a$ and not by the mixing.

## Comparison with the working hypotheses of D1

1. *Meson pole mass fails (b) at leading order.* Confirmed for CP-conserving sources at fixed $m_{a,0}^2$. Refinement: with CP-violating sources, class B has ratio $\epsilon^2$ to the ALP-induced tree term. Classes A and C dominate at leading order.
2. *ALP-like pole mass fails (c).* Confirmed exactly. It also fails (b).
3. *On-shell amplitude passes (a), (b), (d), and (c) only through the part that is not absorbed.* (a), (b) and (d) are confirmed within the limits stated: charged internal lines were not checked for (a), and $A\to3\pi^0$ is the only explicit case. **Partial contradiction with D3:** the pure-ALP part of the $A_0(m_A^2)$ class is exactly a shift of $\widehat P_1$. The coefficient of its logarithm is the running of $\widehat P_1$. It is an independent content only if $\widehat P_1$ is fixed by an independent condition, not as scheme-independent content in general. The candidate's content under (c) then rests on the mixing-induced pieces and the meson logarithms, which are conjectured, not verified, to be non-absorbable. If these are judged insufficient, the negative outcome allowed by D1 applies.
4. *Not anticipated by D1.* The one-loop mixing entry is frame dependent even at the pole. The ALP–meson bubble in $m_n^2$ starts at $\epsilon^6$ and is not absent. Class B also receives an $\operatorname{tr}\widehat S_3$ contribution at $\epsilon^8$.

## Not verified

- The regime where $m_{a,0}^2$ is not tuned ($m_A^2\propto\epsilon^2$).
- Charged internal lines in the $A\to3\pi^0$ loops, and in its external-leg corrections.
- The one-loop axial check of the full amplitude.
- The non-absorbability conjectures for classes B and C of candidate 1 and for the mixing-induced pieces of candidate 4.
- Only one meson-sector rational test set was used. The invariants are symbolic throughout.

## Verification

[`verify_alp_loop_observable_candidates.py`](../../../scripts/verify_alp_loop_observable_candidates.py) imports [`derive_alp_invariant_loop_crossing.py`](../../../scripts/derive_alp_invariant_loop_crossing.py) for the ten-field Lagrangian in both frames. It builds exact power series in symbolic $\epsilon$ for the poles and pole vectors from Rayleigh–Schrödinger perturbation theory with a rational reduced resolvent, keeping all ALP invariants as sympy symbols. For the amplitude it uses the closed-form neutral dressed Lagrangian to sixth order in the fields. The 22 tests in [`test_verify_alp_loop_observable_candidates.py`](../../../scripts/test_verify_alp_loop_observable_candidates.py) assert every power, identity and frame statement above in exact arithmetic. The existing suites `test_derive_alp_invariant_loop_crossing.py` (20 tests) and `test_derive_alp_loop_power_counting.py` (20 tests) still pass.

```sh
.venv/bin/python scripts/verify_alp_loop_observable_candidates.py
.venv/bin/python -m unittest discover -s scripts -p test_verify_alp_loop_observable_candidates.py
```
