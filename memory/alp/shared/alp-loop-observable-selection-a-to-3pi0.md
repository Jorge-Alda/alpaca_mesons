---
title: Selection of the observable for the loop-to-tree comparison at finite f_pi/f_a - coefficients of the bubbles with an ALP-like line in A -> 3 pi0
status: mixed
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd, Gasser:1983yg, Kaiser:1998ds, Herrera-Siklody:1996tqr]
provenance: analytical-derivation
---

# Selection of the observable: bubbles with an ALP-like line in $A\to3\pi^0$

This record is task 1.2 of the OpenSpec change `estimate-loop-tree-crossing-finite-mixing`. It selects the observable used for the loop-to-tree comparison, names the free parameters that could absorb its one-loop correction, and shows which part of that correction they cannot absorb. The comparison of the four candidates is in the [candidate record](alp-loop-observable-candidates-comparison.md) (task 1.1) and is not repeated.

## Observable class used for the estimate (revised scope, task 1.2)

The estimate of the OpenSpec change `estimate-loop-tree-crossing-finite-mixing` compares the tree amplitude with the power-counting size of the one-loop term, without evaluating any loop. The observable class is an **on-shell ALP-sourced amplitude in which the ALP enters at tree level at first order in $1/f_a$**. Its tree term is exact in $\epsilon$ (task 2.1), and the one-loop size is taken from the vertex scalings (task 3.1).

Reasons, from the [candidate record](alp-loop-observable-candidates-comparison.md):

- **Meson pole mass** fails test (b). Its tree shift and the ALP-tadpole class both scale as $\epsilon^2$, so their ratio does not depend on $\epsilon$ at leading order. It cannot locate a value of $f_\pi/f_a$.
- **ALP-meson mixing entry** fails test (a). It changes under a constant axial redefinition, even at the pole at one loop.
- **ALP-like pole mass** fails test (c). $m_{a,0}^2$ absorbs any shift of the pole mass exactly.
- **On-shell amplitude** passes (a) and (d). It has an ALP at first order in $1/f_a$ at tree level, so ALP-line loops carry extra powers of $\epsilon$. Its absorption is discussed in the Dalitz analysis below.

**The Dalitz and absorption analysis below is not used by the estimate.** It was a one-loop evaluation of specific bubble coefficients and is kept only as a record of the task 1.2 investigation. The estimate does not depend on the bubble selection, on the decision of the order-by-order reading (assumption 4), or on the assumptions $m_u=m_d$ and $\langle\lambda_3\widehat P_1\rangle=0$.

## Earlier selection (superseded for the estimate)

**Selection.** Candidate 4, the on-shell amplitude $\mathcal M(A\to\pi^0\pi^0\pi^0)$, restricted to the coefficients of the bubble functions that have an ALP-like internal line: $\bar B(s_{ij};m_A,m_P)$ with $P=\pi^0,\eta,\eta'$, and $\bar B(s_{ij};m_A,m_A)$. The rest of the one-loop amplitude is either absorbed (all $s_{ij}$-independent parts), or is the usual chiral expansion with ratio $\epsilon^0$ to the tree, which needs fitted constants (meson bubbles).

## Conventions and units

- Units GeV. Notation of [@Alda:2026nrd] and of the [candidate record](alp-loop-observable-candidates-comparison.md): $U=\exp(i\Phi/f_\pi)$, $N_f=3$, $x=a/f_a$, $\epsilon=f_\pi/f_a$, fields $\Phi_I=(P_i,a)$, dressed frame with $\widehat M(x)=\sum_n\frac{x^n}{n!}(\widehat S_n+i\widehat P_n)$, $\widehat\Theta(x)=\sum_n\frac{x^n}{n!}\widehat\Theta_n$, $\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1$, $Y=\mathbb M^{-1}_{m_a^2}(\mathcal C)$. $A$ is the ALP-like pole.
- Kinematics: $A(p)\to\pi^0(p_1)\pi^0(p_2)\pi^0(p_3)$, $s_{ij}=(p_i+p_j)^2$, on shell $s_{12}+s_{13}+s_{23}=m_A^2+3m_{\pi^0}^2$. Symmetric point $s_0=(m_A^2+3m_{\pi^0}^2)/3$.
- Loop functions as in the [meson-pole loop record](alp-invariant-loop-crossing-meson-poles.md): $A_0(m^2)=\frac{m^2}{16\pi^2}[1-\ln(m^2/\mu^2)]$ and $\bar B(s;m_1,m_2)=-\int_0^1dz\,\ln\frac{zm_1^2+(1-z)m_2^2-z(1-z)s-i0}{\mu^2}$ ($\overline{\rm MS}$). In this convention $\bar B$ contains the whole $\mu$ dependence of the finite part, and $d\bar B/d\ln\mu^2=1$. The symbol $\bar B$ is still awaiting approval there.
- **Local abbreviations flagged for approval:** $e_2'$, $e_3'$ denote the elementary symmetric polynomials of degree 2 and 3 of the three deviations $s_{ij}-s_0$ (the degree-1 one vanishes on shell). They are used only in this record.

## Assumptions

1. Those of the [candidate record](alp-loop-observable-candidates-comparison.md): LO $U(3)$ in the dressed frame, one loop, flavour-diagonal invariants kept as symbols, $m_{a,0}^2$ held fixed as $\epsilon$ varies, non-degenerate poles, CP-conserving sources.
2. The meson sector of the exact checks has $m_u=m_d$ (rational test inputs only). The ALP sources may break isospin.
3. Reading of test (c) of design D1, fixed here. A one-loop contribution is *absorbed* if, for all external kinematics, it equals the tree-level contribution of local operators from a stated set, with coefficients that do not depend on the kinematics. Since the tree of a single ALP amplitude is linear in free invariants, any $s_{ij}$-independent shift of one amplitude is absorbed. Non-absorbable content must therefore show up in the dependence on the kinematics (or in relations among several observables, which are not used here).
4. Order-by-order reading (assumption, to be confirmed by the user): operators with six or more derivatives are higher order in the chiral expansion and are not counted as absorbing a one-loop effect at the order considered.

## Inputs

- Tree amplitude, its leading power $\epsilon^1$, and the leading powers of the one-loop classes ($A_0(m_A^2)$, $\bar B(s;m_A,m_P)$ at $\epsilon^3$; $\bar B(s;m_A,m_A)$ at $\epsilon^5$; meson classes at $\epsilon^1$): [candidate record](alp-loop-observable-candidates-comparison.md).
- Pure-ALP tadpole as the source shift $\widehat M\to\widehat M-\frac{A_0}{2f_a^2}\widehat M''$, $\widehat\Theta\to\widehat\Theta-\frac{A_0}{2f_a^2}\widehat\Theta''$: [candidate record](alp-loop-observable-candidates-comparison.md).
- Absence of cubic vertices for CP-conserving sources: [power-counting record](alp-loop-power-counting.md), L1.
- Linearity of the constant axial map in the neutral sector: [candidate record](alp-loop-observable-candidates-comparison.md).
- Retained NLO operators $L_4,L_5,L_6,L_7,L_8,\Lambda_1,\Lambda_2,L_{18},L_{25}$ with the ALP entering through $\chi$ and $X_{TC}$: [NLO Lagrangian](../../nlo/lagrangian.md). Singlet-source operators $V_5\,D_\mu\theta D^\mu\theta$, $O_{46}$, $O_{47}$, $O_{50}$, $O_{52}$: [singlet-source derivative basis](../../nlo/singlet-source-derivative-basis.md), from [@Kaiser:1998ds; @Herrera-Siklody:1996tqr].
- Polynomial (subtraction-dependent) parts of loops are absorbed by local counterterms of the next order [@Gasser:1983yg]; renormalization of the NLO constants: [one-loop effects record](../../nlo/one-loop-effects-and-renormalized-lecs.md).

## Free parameters that could absorb the one-loop correction of $A\to3\pi^0$

| set | parameters | what they give for $A\to3\pi^0$ on shell |
|---|---|---|
| LO | all Taylor coefficients $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$ (the tree is $\propto\langle(2Ym_q-\widehat P_1)\lambda_3\rangle$ at $O(\epsilon)$), $m_{a,0}^2$, the pure-ALP potential, $m_q$, $m_0^2$, $B_0$, $f_\pi$ | an $s_{ij}$-independent constant (no derivative couplings in the neutral sector, no cubic vertices) |
| NLO in [memory/nlo](../../nlo/INDEX.md) | $L_4,L_5,L_6,L_7,L_8,\Lambda_1,\Lambda_2,L_{18},L_{25}$; $V_5$; $O_{46},O_{47},O_{50},O_{52}$ (at most two derivatives each) | an $s_{ij}$-independent constant |
| $O(p^4)$, not in memory/nlo | $L_1,L_2,L_3$ (fixed QCD constants, entering through the meson components of $A$); four-derivative singlet-source operators of the full $U(3)$ basis (unknown coefficients) | at most a constant plus a multiple of $\sum_{ij}s_{ij}^2$ |
| $O(p^6)$ and higher | e.g. six-derivative operators with $\chi(x)$, coefficients free for each ALP | can also give $e_3'$; excluded by assumption 4 |

Why the entries hold:

- **Neutral vertices carry no momenta.** For commuting diagonal $\Phi$, $\frac{f_\pi^2}{4}\langle\partial_\mu U\partial_\nu U^\dagger\rangle=\frac14\langle\partial_\mu\Phi\partial_\nu\Phi\rangle$ exactly (S1). The mass term and $X_{TC}=-\operatorname{tr}\Phi/f_\pi-2\widehat\Theta(x)$ have no derivatives, and the ALP kinetic term is $\frac12(\partial a)^2$. So every vertex among $\pi^0,\eta,\eta',a$ in the dressed frame is momentum independent (exact, all orders in the fields). In the charged sector this is not true.
- **Contact terms.** No cubic vertex exists for CP-conserving sources: at LO this is verified, and at NLO it follows from CP, because three CP-odd states cannot couple. The tree amplitude from any set of local operators is therefore a contact term times external-leg factors. If the operators have at most two derivatives, the amplitude is linear in the $p_I\cdot p_J$. Symmetrized over the three $\pi^0$ and taken on shell, it is a constant (S2a). With at most four derivatives it lies in the span of $1$ and $\sum_{ij}s_{ij}^2$ (S2b). $\sum s_{ij}^2$ is not constant on shell, and $\sum s_{ij}^3$ is not in that span, because $e_2'$ and $e_3'$ vary independently on the Dalitz plane (S2c).
- **Singlet-source four-derivative operators** such as $D_\mu\theta\langle u^\mu u^\nu u_\nu\rangle$. At leading order in $\epsilon$ and with $m_u=m_d$ they do not contribute, since $\langle\lambda_3^3\rangle=0$ (argued, not computed). In any case they lie in the span of $1$ and $\sum s^2$ by S2b.

## One-loop $A\to3\pi^0$: absorbed and non-absorbed parts

| part | dependence on $s_{ij}$ | status |
|---|---|---|
| tadpoles with neutral lines (contractions $aa$, $aP$, $PP$, including the $A$-pole part of the neutral propagators), coefficients of $A_0(m_r^2)$ | none (momentum-independent vertices) | absorbed (constant). The pure-ALP part equals the shift $\widehat P_1\to\widehat P_1-\frac{A_0}{2f_a^2}\widehat P_3$ (and the corresponding shifts of all $\widehat M_n$, $\widehat\Theta_n$). The coefficient of its logarithm is the running of free Taylor coefficients. It has no separate content unless a UV completion fixes $\widehat M(x)$, $\widehat\Theta(x)$ at a stated scale. The other neutral tadpoles also have LO form, with flavour-dependent coefficients (argued). |
| tadpoles with charged or strange lines (two-derivative vertices from $\langle[\Phi,\partial\Phi]^2\rangle$) | linear in $p_I\cdot p_J$, so constant on shell (S2a) | absorbed (constant); not of LO form (argued, not computed) |
| external-leg corrections | none | absorbed (constant) |
| UV divergence of every neutral bubble, and its $\mu$ dependence | $d\bar B/d\ln\mu^2=1$ (S3c) times a momentum-independent coefficient: constant | absorbed |
| meson bubbles $\bar B(s_{ij};m_P,m_{P'})$, neutral lines | not polynomial | not absorbed by LO + memory/nlo; ratio $\epsilon^0$ at leading order, the usual chiral expansion |
| charged and kaon bubbles ($\pi^+\pi^-$, $K\bar K$) | momentum-dependent vertices give UV divergences $\propto\sum s_{ij}^2$ | need $L_1,L_2,L_3$ (fitted, outside memory/nlo); ratio $\epsilon^0$; argued, **not computed** |
| **bubbles with an ALP-like line: coefficients of $\bar B(s_{ij};m_A,m_P)$, $P=\pi^0,\eta,\eta'$ ($\epsilon^3$) and $\bar B(s_{ij};m_A,m_A)$ ($\epsilon^5$)** | not polynomial; contains $e_3'$ | **not absorbed by any operator with at most four derivatives (selected)** |

Argument for the selected coefficients:

1. **Their functions contain an $e_3'$ term.** For a function $g$, the coefficient of $e_3'$ in $\sum_{ij}g(s_{ij})$, expanded about $s_0$, is exactly $g'''(s_0)/2$ (S3b). For $g=\bar B(\,\cdot\,;m_1,m_2)$, $g'''(s_0)>0$ whenever $s_0<(m_1+m_2)^2$. This uses $D-z(1-z)s=(zm_1-(1-z)m_2)^2+z(1-z)[(m_1+m_2)^2-s]$ with $D=zm_1^2+(1-z)m_2^2$. It holds for $(m_A,m_P)$ with $m_P\ge m_{\pi^0}$, and for $(m_A,m_A)$ with the decay open. The Taylor coefficients of $\bar B$ at $s^2$ and $s^3$ are positive (S3a).
2. **No operator with at most four derivatives gives $e_3'$** (S2b, S2c). LO and memory/nlo give only a constant. $L_{1,2,3}$ and four-derivative singlet-source operators give only $\sum s_{ij}^2$.
3. **Loop functions with different thresholds are independent functions of $s_{ij}$.** So the coefficient of $\bar B(s_{ij};m_A,m_P)$ is identifiable separately from the meson bubbles, which also contain $e_3'$ terms.
4. **Each coefficient is a product of two tree-level quartic vertices** contracted with the exact pole vectors. It contains no $\mu$ and no fitted number. A shift of a free parameter changes it only through the tree vertices, never independently of them. Hence none of the parameters in the table above is this coefficient, and no counterterm up to four derivatives can change it.
5. **Frame invariance.** The four coefficients are equal in the dressed frame and in the frame of a constant axial redefinition with symbolic $\kappa$, exactly through $\epsilon^5$ (S5). All their lines are neutral.

## Structure of the selected coefficients (exact, verified with symbolic invariants)

- **Leading powers** are confirmed: $\epsilon^3$ for $(A,\pi^0)$, $(A,\eta)$, $(A,\eta')$, and $\epsilon^5$ for $(A,A)$. The tree is $O(\epsilon)$, so the ratios are $\epsilon^2$ and $\epsilon^4$ (S4).
- **Invariants entering.** At leading order they depend on $\widehat P_1$, the diagonal entries of $\widehat S_2$ and $\widehat\Theta_1$; $(A,A)$ also depends on $\widehat P_3$.
- **Factor structure at the test point.**
  - $(A,\pi^0)$ is proportional to $\langle\lambda_3\widehat P_1\rangle$, the same isovector combination as the leading tree.
  - $(A,\eta)$ and $(A,\eta')$ factorize into two pieces. One is an isoscalar combination of $\widehat P_1$ and $\widehat\Theta_1$. The other contains $\langle\lambda_3\widehat P_1\rangle$ and the isovector part of $\widehat S_2$ (the $u$ entry minus the $d$ entry).
  - The numerical prefactors depend on the rational test inputs and carry no physical meaning.
- **The loop does not follow the tree.** With $\langle\lambda_3\widehat P_1\rangle=0$, the leading tree vanishes, and with $m_u=m_d$ the tree vanishes through $\epsilon^5$. The $(A,\pi^0)$ coefficient then vanishes through $\epsilon^5$. But $(A,\eta)$ and $(A,\eta')$ survive at $\epsilon^3$, proportional to the isovector part of $\widehat S_2$ times the isoscalar combination. $(A,A)$ survives at $\epsilon^5$. The same holds when, in addition, the isovector part of $\widehat P_3$ vanishes.
  - Consequence: the loop-to-tree ratio is $\epsilon^2$ times ratios of independent invariants, and it is not bounded by a function of the masses alone. The crossing formula of later tasks must keep these ratios as invariants. No relation among them is assumed.
- **Isospin limit.** All four coefficients vanish through $\epsilon^5$ when every isovector component of every invariant vanishes.

## Consequences for the OpenSpec design (not applied here)

- **D3** builds the estimate from the coefficient of $\ln(m^2/\mu^2)$. For the selected observable that coefficient belongs to the absorbed tadpoles: it is the running of free source coefficients. The content that cannot be absorbed is in the coefficients of the bubbles with an ALP-like line. These coefficients are finite and $\mu$-independent, and are identified through their $s_{ij}$ dependence. D3 should say this. The absorbed constants should be listed with their $\mu$ dependence, and the meson-bubble background (with $L_{1,2,3}$) should be reported separately.
- **D1, test (d).** $L_{1,2,3}$ enter only the excluded meson-bubble background.
- **Size measure.** The tree is constant while the selected loop part depends on $s_{ij}$. Comparing them needs a measure that does not change when a constant or a multiple of $\sum s_{ij}^2$ is added. Possible choices:
  - half the spread of the real part over the Dalitz region of the selected bubbles, after removing the span of $1$ and $\sum s^2$;
  - the coefficient of $e_3'$ at $s_0$.

  This is a choice for D6 and needs a user decision.

## Not verified

- Charged and kaon internal lines in the one-loop amplitude: the $\sum s^2$ divergence needing $L_{1,2,3}$, and frame invariance (L6 assumption). These lines have no $A$ threshold, so they do not enter the selected coefficients. A charged intermediate state $A\pi^\pm$ is forbidden by charge.
- The claim that the neutral tadpoles have LO form with flavour-dependent coefficients, and that the charged tadpoles are not of LO form (argued; for $A\to3\pi^0$ both are constants anyway).
- The zero of the singlet-source four-derivative operators at $m_u=m_d$ (argued from $\langle\lambda_3^3\rangle=0$). It is not needed, because of S2b.
- The regime with $m_A^2\propto\epsilon^2$ (untuned $m_{a,0}^2$), and $m_u\neq m_d$ in the meson sector.

## Verification

[`verify_alp_loop_observable_selection.py`](../../../scripts/verify_alp_loop_observable_selection.py) imports [`verify_alp_loop_observable_candidates.py`](../../../scripts/verify_alp_loop_observable_candidates.py) without modifying it. It runs checks S1–S5 in exact arithmetic with symbolic invariants and masses:

- S1: free neutral kinetic term;
- S2: on-shell contact structure;
- S3: $\bar B$ on the Dalitz plane and its $\mu$ dependence;
- S4: coefficients of the bubbles with an ALP-like line and the cases where the tree vanishes;
- S5: equality in the dressed and $\kappa$ frames.

The 17 tests in [`test_verify_alp_loop_observable_selection.py`](../../../scripts/test_verify_alp_loop_observable_selection.py) assert every statement above that is marked as verified, and run in about 75 s. The 22 tests of the task 1.1 suite still pass.

```sh
.venv/bin/python scripts/verify_alp_loop_observable_selection.py
.venv/bin/python -m unittest discover -s scripts -p test_verify_alp_loop_observable_selection.py
```
