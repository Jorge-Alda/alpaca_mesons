# Tasks

All scripts are written by the `coder` agent only. No task may evaluate a loop integral, scalar loop function, loop coefficient, $\mu$ derivative or loop amplitude. Scripts may check bookkeeping and trace algebra only. All tests are analytical (exact symbolic); no random or numerical matrices.

## 1. External legs and tree

- [x] 1.1 Record the legs $\sigma$, $\pi^+$, $\pi^-$, the kinematics ($p_\sigma^2=m_\sigma^2$, external dot products exact) and the exact tree amplitude $M_0=V(\sigma;\pi^+,\pi^-)$ for generic sources, linking `radial-decay-polar-cartesian-general-sources.md` and `murayama-radial-up-only-estimate-result.md`. Verify: a script reproduces the tree from the vertex rule and agrees with `murayama-sigma-pp` for the up-only sources.
- [x] 1.2 State the assumptions (no tadpoles, no quadratic radial mixing, $\widehat\Theta_1=0$). Verify: each is linked to the record that states it.

## 2. Graphs and vertices

- [x] 2.1 Enumerate the one-loop 1PI graphs with legs $\sigma\pi^+\pi^-$ (valences $(5)$, $(3,4)$, $(3,3,3)$) containing at least one ALP-like or radial internal line (a meson line may accompany it), with $S_\Gamma$ and $r_\Gamma$, and the external-leg corrections; list the excluded meson-only chiral-background terms (D4, exact split) separately. Verify: script counts the graphs per class and shows no box and no $(2,5)$ 1PI graph.
- [x] 2.2 Write the vertex rules for generic $\widehat S_{r,s}$, $\widehat P_{r,s}$ up to five legs, including the kinetic parts, the $\sigma$ self-couplings and the selection rules on $\sigma$ and ALP-like ends that truncate the sum. Verify: the cubic and quartic rules agree with `murayama-loop-generator-trace-products.md` for the up-only sources.

- [x] 2.3 Establish the $\epsilon$-scaling at large $\epsilon$ of $N(m_a^2)$, $\xi_a$, $\mathbb M^{-1}_{m_a^2}(\mathcal C)$ and the meson-line residue at fixed sources and fixed $\lambda$, and state the graph selection and D4 under D10 (no truncation in $\epsilon$). Verify: symbolic script; the selection record states the regime and no graph is dropped by a small-$\epsilon$ series.

## 3. Traces

- [x] 3.1 Remove the generator sums by completeness for each graph and write $T_\Gamma$ in the compact trace form of D9 (single and double traces, commutators, anticommutators, sums over permutations or cyclic permutations). Verify: the symbolic script confirms the completeness identities and the trace equalities exactly for generic noncommuting symbols; no numerical matrices and no tolerance.
- [x] 3.2 Write the records of the traces per graph class in `memory/alp/shared/`, with front matter and linked inputs. Verify: no loop value appears, and no component notation or $E_{ij}$ appears.

## 4. Estimate

- [x] 4.1 Assemble $M_1^{\rm est}=\sum_\Gamma\sum_{\rm terms}S_\Gamma\,(m_{\rm term}^2)^{2-r_\Gamma}/(16\pi^2)\,T_{\rm term}$ and $R=|M_1^{\rm est}|/|M_0|$, with $m_{\rm term}$ as in D3. Verify: dimensional analysis by rescaling all dimensionful inputs, as for $a\to\pi^+\pi^-\pi^0$.
- [x] 4.2 Check invariance under a constant axial redefinition of the sources on the traces. Verify: script.
- [x] 4.3 Record the limits (including the CP and flavour-diagonality assumptions of D11): absorptive parts and thresholds, cancellations, the naturalness assumption without counterterm analysis. Verify: the record contains no value of a loop.

- [x] 4.4 Implement the per-term loop scales of D3 (bare generator part at $m_P$, rest of the meson line at $m_a$, maximum over the lines of each term) and the termwise routing maximum of D12 in the graph, trace and estimate scripts, and regenerate the trace, estimate and graph-list records. Verify: exact symbolic tests that (a) graphs independent of the routing enter with their single-routing value; (b) $|C_t^{\rho_t^*}|\ge|C_t^\rho|$ for every routing and every $t$; (c) the commutator factor $m_a^2-m_\sigma^2$ is absent from $\delta Z_{\pi^\pm}$; (d) with a common routing and common scale $m_a$, ALP-like line plus pieces (ii)+(iii)+(iv) equals the bare ALP-field line; (e) dimensional rescaling with $\max$ is exact.
- [x] 4.5 Compute the large-$\epsilon$ leading powers (fixed sources, fixed $\lambda$) of every graph, of $M_1^{\rm est}$ and of $R$ with the scales of D3 and the routing rule of D12, the leading coefficient (generic and up-only sources), and check the restated D4 criterion. Verify: symbolic script; the D4 terms contain no ALP-sector letter and depend on $\epsilon$ only through $m_\sigma^2$.
- [x] 4.6 Derive the up-only limit $M_1^{\rm est}/M_0\to-21B_0m_u\epsilon^4/(64\pi^2f_\pi^2)$ graph by graph (graphs 2, 13, σ4) in its own record, with the up-only substitution of the sources and the region boundary $\kappa_\sigma=A$. Verify: symbolic script and tests; the record values equal the script output.

## 5. Indexes

- [x] 5.1 Update `memory/alp/shared/INDEX.md` (and parents if affected) and `memory/plan.md`. Verify: the `memory/README.md` crawl reaches every new record.
- [ ] 5.2 Propose a commit message and archive the change after review.
