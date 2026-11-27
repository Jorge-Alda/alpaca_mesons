# Tasks

All scripts are written by the `coder` agent only. No task may evaluate a loop integral, scalar loop function, loop coefficient, $\mu$ derivative or loop amplitude. All tests are exact symbolic identities (no random or numerical matrices). Traces in compact matrix form only. The calculation is done from scratch in the Cartesian fields; the only polar comparison used in the derivation is the tree equivalence (task 2.2). Records go in `memory/alp/radial-mode/finite-mixing/two-pion-loop/cartesian-estimate/` with an `INDEX.md`, linked from the parent index.

## 1. Frame, legs and tree

- [x] 1.1 Define $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$ by the expansion of $\widehat M(s,p)$, record the Cartesian Lagrangian, legs $s$, $p$, mesons and external $\sigma$, $\pi^+$, $\pi^-$, and restate each assumption (no tadpoles, no quadratic mixing, $\widehat\Theta_1=0$, CP, flavour-diagonality, $f_\pi$ independent of the radial fields) as a condition on the Cartesian invariants, derived from the quadratic Lagrangian and linked to the polar assumption. Verify: script derives each condition; exact tree $M_0$ in the invariants.

## 2. Vertices

- [x] 2.1 Derive the vertex rules for $n=3,4,5$ from the Cartesian Lagrangian in generic $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$, including $V_\varphi(s,p)$ and the meson kinetic vertices, by constructions 1 and 2. Verify: the constructions agree as exact identities of cyclic words for every leg pattern; relabelling symmetry; dimensional scaling $c^{4-n}$; up-only limit.
- [x] 2.2 Tree equivalence, the only comparison with the polar results: check that the on-shell $V(s;\pi^+,\pi^-)$ equals the polar tree of the existing notes, using the relation between the Cartesian and polar coefficients, stationarity and the pole equation. Verify: exact symbolic script and test; no off-shell, quartic or quintic comparison.
- [x] 2.3 Establish the large-$\epsilon$ powers of the vertices at fixed sources and fixed $\lambda$ and state the regime (D10). Verify: symbolic script; no truncation in $\epsilon$.

## 3. Graphs and traces

- [x] 3.1 Enumerate from the Cartesian vertices the one-loop 1PI graphs and external-leg corrections, with $S_\Gamma$, $r_\Gamma$ and the D4 excluded terms; check the absence of a box and of a 1PI $(2,5)$ graph. Verify: script counts per class. The comparison with the polar graph list is made only afterwards (task 4.4).
- [x] 3.2 Remove the generator sums by completeness and write $T_\Gamma$ in compact trace form, with the D3 split of the meson line and the D12 routing. Verify: exact symbolic tests (routing-independent graphs, $|C_t^{\rho_t^*}|\ge|C_t^\rho|$, $\delta Z_{\pi^\pm}$ without $m_a^2-m_\sigma^2$, dimensional rescaling with $\max$).
- [x] 3.3 Write the trace records per graph class. Verify: no loop value, no component notation, no $E_{ij}$.

## 4. Estimate and comparison

- [x] 4.1 Assemble $M_1^{\rm est}$ and $R$ with the D3 scales; check dimensions. Verify: script.
- [x] 4.2 Check invariance under a constant axial redefinition of the sources. Verify: script.
- [x] 4.3 Compute the large-$\epsilon$ leading powers and coefficients of each graph, $M_1^{\rm est}$ and $R$ (generic and up-only), the D4 criterion, and whether same-scaling graphs have leading coefficients that could cancel. Verify: symbolic script.
- [x] 4.4 Diagnostic comparison with the polar estimate, graph class by graph class, including the up-only limit. State that agreement is not expected, that the results are orders of magnitude that can miss cancellations differently in each parametrization, and that neither is authoritative. Verify: record values equal script output.
- [x] 4.5 Record the limits (restated polar limits plus the order-of-magnitude caveat). Verify: no loop value in the record.

## 5. Indexes

- [x] 5.1 Create `cartesian-estimate/INDEX.md`, update `memory/alp/radial-mode/finite-mixing/two-pion-loop/INDEX.md`, parents if affected, and `memory/plan.md`. Verify: the `memory/README.md` crawl reaches every new record.
- [ ] 5.2 Propose a commit message and archive the change after review.
