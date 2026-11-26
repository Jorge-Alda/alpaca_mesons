# Design

## Context

See proposal.md for motivation. The estimate of `a→π⁺π⁻π⁰` (`memory/alp/shared/murayama-radial-up-only-estimate-result.md`) defines the vertex rules, the legs $(\lambda_x,\xi_x)$, the internal-line rules and the size formula
$$
M_1^{\rm est}=\sum_\Gamma S_\Gamma\,\frac{(m_\Gamma^2)^{2-r_\Gamma}}{16\pi^2}\,T_\Gamma .
$$
Its potential-part vertex is already written for generic $(r,s)$ and generic $\widehat S_{r,s}$, $\widehat P_{r,s}$. The up-only sources enter only afterwards. Related records:

- generator-sum form of the diagram products, tested with noncommuting sources: `murayama-loop-generator-trace-products.md`;
- external-leg residues: `murayama-radial-loop-external-residues.md`;
- tree $\sigma\to P_iP_j$ for general sources: `radial-decay-polar-cartesian-general-sources.md`;
- scalar-moment reduction and naturalness: `alp-loop-momentum-reduction.md`, `alp-loop-naturalness-assumption.md`, `alp-loop-estimate-cancellations.md`;
- axial-redefinition check: `alp-loop-estimate-frame-invariance.md`.

Observed: `σ→π⁺π⁻` has three external legs. At one loop the 1PI graphs have vertex valences $(5)$ for a tadpole, $(3,4)$ for a bubble and $(3,3,3)$ for a triangle. A bubble $(2,5)$ is the self-energy of an external leg, and a box would need a vertex without an external leg and with only two lines. Counting is to be confirmed by a script, not assumed.

## Goals / Non-Goals

**Goals:**
- A closed symbolic expression for $M_1^{\rm est}$ in generic $\widehat S_{r,s}$, $\widehat P_{r,s}$, in which all integrals are replaced by $(m^2)^{2-r}/(16\pi^2)$.
- Remove the generator sums analytically, so that the result is a list of single and double traces.

**Non-Goals:**
- Computing any loop integral, scalar function, coefficient, $\mu$ derivative or amplitude.
- Absorptive parts, thresholds, renormalization and counterterm analysis.
- Numerical scans, up-only numbers, integrated widths, the draft in `latex/` and the diagram atlas.
- Reducible graphs (tree with loop on an internal propagator); there is no internal propagator in a three-point tree with a contact vertex.

## Decisions

**D1. Same prescription as for $a\to\pi^+\pi^-\pi^0$.** The size of each scalar moment with $r$ denominators is $(m^2)^{2-r}/(16\pi^2)$, as requested. Alternative rejected: any evaluation of a Passarino–Veltman function, even in a limit.

**D2. Traces are exact; loop-momentum products are not.** Products of the loop momentum with any momentum are set to $m_{\rm term}^2$, the loop scale of the term (D3). Products of external momenta are exact: $p_\sigma^2=m_\sigma^2$, $p_{\pi^\pm}^2=m_{\pi^\pm}^2$, $p_{\pi^+}\!\cdot p_{\pi^-}=(m_\sigma^2-2m_{\pi^\pm}^2)/2$.

**D3. Loop scale, set per term (assumption).** Each internal line has the scale of its poles. Radial line: $m_\sigma$. ALP-like line: $m_a$. A meson line, with residue $(\mathbf 1_9\oplus1)-v_av_a^{\rm T}$, is split exactly into (i) the bare $U(3)$ generator sum, whose poles are the nine meson-like states, with scale $m_P$; and (ii)+(iii)+(iv) $=(0\oplus1)-v_av_a^{\rm T}$ (ALP-field entry $(N-1)/N$, mixed entries, meson block of the ALP-like pole), which is the remainder of the ALP-like propagator and has the scale $m_a$. For every term of a graph (one choice of piece on each meson line) the scale is the largest scale among its lines, $m_{\rm term}=\max(\text{line scales})$; bare meson lines therefore also enter the maximum in graphs with a radial or ALP-like line. $m_P$ is a single symbolic mass, the largest pole mass of the nine meson-like states; no numerical value is used. A mass per generator is not possible, because the completeness sum (D5) runs over all nine generators on every bare line. No lower limit on $m_{\rm term}$ other than the line scales is imposed. Consequences: with a common scale $m_a$, the ALP-like line plus pieces (ii)+(iii)+(iv) of a meson line is the bare ALP-field line, so the internal ALP–meson mass mixing is neglected inside the loops (it is kept in the vertices and on the external legs); the factor $m_a^2-m_\sigma^2$ in $\delta Z_{\pi^\pm}$ does not appear. The routing prescription is D12. Alternative rejected: one scale per graph.

**D12. Routing: largest magnitude per graph, termwise (assumption).** For a graph $\Gamma$ and each of its $2r_\Gamma$ routings $\rho$ (loop momentum on one internal line, either orientation), write the contribution as $\sum_tC_t^\rho\,t$, with $t$ the independent products of traces of the canonical generic-source basis of the graph scripts and $C_t^\rho=\sum_{\rm terms}S_\Gamma(m_{\rm term}^2)^{2-r_\Gamma}c_{t,{\rm term}}^\rho/(16\pi^2)$ (scales of D3 included, pieces summed before maximizing). For a self-energy the same is done on $\delta Z^\rho$, with the derivative rule applied in each routing. For each $t$ the routing $\rho_t^*$ with the largest $|C_t^\rho|$ is chosen, and the graph enters $M_1^{\rm est}$ as $\sum_tC_t^{\rho_t^*}t$ (sign kept). This is an upper bound per graph on the routing ambiguity of the size rule, not on the true loop: $|\sum_tC_t^\rho t|\le\sum_t|C_t^{\rho_t^*}|\,|t|$ for every routing $\rho$. Justification: the routing that maximizes the magnitude of the whole $T_\Gamma$ is not defined without numerical values of the source traces; the termwise maximum is defined symbolically and is at least as large. Ties of equal magnitude and opposite sign take the sign of the sum of the coefficients over all routings; if this sum is zero, the first routing in the canonical order is used and the term is marked (a choice at a degenerate mass relation). If $\rho_t^*$ depends on the masses, the branches are recorded as inequalities; at large $\epsilon$ the choice is made by the leading power in $\epsilon$, then by the size of the leading coefficient. The routing spread is kept as a diagnostic. Alternative rejected: a single routing per graph.

**D4. Meson-only chiral background excluded as an exact split.** The residue of each internal meson line is written exactly as $(\mathbf 1_9\oplus1)-v_av_a^{\rm T}$: bare $U(3)$ generator sum, ALP-field entry, minus the ALP-like pole. Only the term in which every internal line takes the bare generator sum and the only power of $\epsilon$ is the single $\sigma$ end is excluded. With the scales of D3 it has the scale $m_P$. Criterion: the removed term contains no letter of the ALP sector ($\mathbb M^{-1}_{m_a^2}(\mathcal C)$, $N$, $\widehat P_{0,1}$, $\widehat S_{r,s\ge1}$, $\widehat P_{r,s\ge1}$), and its ratio to the tree depends on $\epsilon$ only through the kinematics $p_\sigma^2=m_\sigma^2$ (checked by script for graphs 3, 15 and π2; at large $\epsilon$ their powers relative to the tree are $\epsilon^0$, $\epsilon^2$, $\epsilon^0$, below $R\sim\epsilon^4$). It is the chiral background of the model without the ALP, recorded but not estimated, and matters only if it alone approaches $R_{\max}$. Everything else in graphs with only meson lines is kept. Alternative rejected: excluding whole meson-only graphs.

**D10. Valid for any $\epsilon$, organised for $\epsilon\gg1$.** No step assumes $\epsilon\ll1$. $N(m_a^2)$ and $\mathbb M^{-1}_{m_a^2}(\mathcal C)$ stay exact and symbolic, and no graph is kept or dropped by a series in $\epsilon$. The large-$\epsilon$ behaviour is taken at fixed sources and fixed quartic coupling $\lambda>0$, as in the up-only model, with $m_P$ fixed: the ALP mass entry is $\epsilon^2A$, $A=-2B_0\langle\widehat S_{0,2}\rangle>0$, so $m_a^2\simeq\epsilon^2A$, $N\to1$, $\xi_a\to1$ and the meson part of $v_a$ falls as $\epsilon^{-1}$; $m_\sigma^2=2\lambda f_\pi^2/\epsilon^2+\epsilon^2\kappa_\sigma$ with $\kappa_\sigma=B_0(\langle\widehat S_{1,0}\rangle-2\langle\widehat S_{2,0}\rangle)>0$. The leading powers are labels of the exact estimate. Alternative rejected: a small-$\epsilon$ truncation.

**D11. External-leg mixing (assumptions).** The source terms are assumed CP-invariant ($\sigma$ CP-even, $a$ and $\Phi$ CP-odd). This gives $\langle\widehat S_{1,1}\rangle=0$ and removes $\sigma$ mixing with the ALP-like state and with any neutral pseudoscalar. $\pi^+\to K^+$ vanishes if the sources and the quark mass matrix are flavour-diagonal. Both are stated assumptions; the trace algebra is still done for generic noncommuting sources. The one-particle-reducible external-leg mixing graphs are therefore not included.

**D5. Generator sums by completeness.** With $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$ and $\lambda_0=\sqrt{2/3}\,\mathbf 1$, $\sum_i\lambda_i\otimes\lambda_i$ gives $\sum_i\langle A\lambda_iB\lambda_i\rangle=2\langle A\rangle\langle B\rangle$ and $\sum_i\langle A\lambda_i\rangle\langle\lambda_iB\rangle=2\langle AB\rangle$. The vertex products become traces, which is what keeps the result in terms of the invariants. The completeness identities are applied before any expansion, so no component notation appears.

**D9. Compact trace notation.** Every trace is written in matrix form only: single or double traces of products of the sources, the external matrices and the pole matrices; traces of commutators $[A,B]$ or anticommutators $\{A,B\}$; or sums over permutations or cyclic permutations of a fixed word. Matrix indices, explicit generator components $(\lambda_i)_{ab}$ and the unit-matrix basis $E_{ij}$ are not allowed in any record or script output.

**D6. Internal lines reused.** The ALP-like line uses the exact pole vector $v_a$ with $\mathbb M^{-1}_{m_a^2}(\mathcal C)$ and $N(m_a^2)$ left symbolic; the meson line is the generator sum minus the ALP pole; the radial line has no meson component, assuming $\widehat P_{1,0}=0$ and $\langle\widehat S_{1,1}\rangle=0$ as in the tree record. $\sigma$ self-couplings come from $V_\varphi$ and $\langle\widehat S_{r,0}\rangle$.

**D7. Charged pions do not mix with the ALP.** The ALP is neutral and electric charge is conserved, so $\pi^\pm$ legs carry $(\lambda_{\pi^\pm},0)$ exactly. This is a fact, not an assumption.

**D8. Check, not scan.** The estimate is checked (a) against the tree limit of the vertex rules, (b) by exact symbolic verification of the completeness identities and trace equalities for generic noncommuting symbols (no random or numerical matrices; all tests are analytical), (c) under a constant axial redefinition. The two completeness identities are checked once with the explicit $U(3)$ generators and symbolic $3\times3$ matrices inside the test; these matrices are never printed or written to a record. Neither check evaluates a loop.

## Risks / Trade-offs

- [Same-order graphs cancel, so the sum overestimates] → State the cancellation conditions of `alp-loop-estimate-cancellations.md` as limits.
- [$s=m_\sigma^2$ lies above the $\pi\pi$ and possibly $aa$, $a\pi$ thresholds, where the estimate has no absorptive part and can miss a rescattering enhancement] → State as a limit; give magnitude only.
- [The all-generator terms of meson-only loops are dropped, though this chiral correction could exceed $R_{\max}$ on its own] → State the background and its size scale $\max(m_P^2,m_\sigma^2)/(16\pi^2f_\pi^2)$ as a limit; do not compute it.
- [The termwise routing maximum (D12) can overstate graphs whose coefficients vanish in some routings] → Stated as an upper bound per graph on the routing ambiguity, not on the loop.
- [The size rule cannot resolve on-shell cancellations between the kinetic $\sigma aa$ vertex and the $\widehat S_{1,2}$ source vertex in graphs 2, 13 and σ4, which give the whole up-only coefficient] → Stated as a limit.
- [At large $\epsilon$, $m_a$ and $m_\sigma$ grow beyond the range of leading-order chiral perturbation theory, and $R\ll1$ bounds $\epsilon$ from above] → Stated as limits.
- [$\delta Z=-\Pi'(m^2)$ with the scalar-moment size and the vanishing one-loop $\sigma$ tadpole] → Assumptions, stated in the limits record.
- [Generic $r,s$ gives many vertices] → Truncate by explicit selection rules on the number of $\sigma$ and ALP-like ends, stated in the record, and verify the count with a script.
- [Neglected terms may carry counterterms] → Naturalness assumption, with no absorption analysis.
