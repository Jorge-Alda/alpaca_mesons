---
title: Compact flavour traces of the kept one-loop 1PI graphs of sigma to pi+ pi-
status: derived
created: 2026-11-26
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Compact traces of the kept 1PI graphs of $\sigma\to\pi^+\pi^-$

Tasks 3.1, 3.2 (1PI part) and 4.4 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama`. No loop integral, loop function or loop amplitude is evaluated. The selection holds for any $\epsilon$; no step assumes $\epsilon\ll1$.

## Inputs

The 26 kept graphs and their numbering come from the [graph-list record](murayama-sigma-pipi-loop-graph-list.md), the selection rules (D3, D4, D11, D12) from the [selection record](murayama-sigma-pipi-loop-graph-selection.md), the vertex rules from the [vertex-rule record](murayama-sigma-pipi-loop-vertex-rules.md), the pole vector $v_a$ and the meson-line residue from the [large-$\epsilon$ record](murayama-sigma-pipi-loop-large-epsilon-scaling.md), and the legs, kinematics, charge rule and assumptions from the [tree and legs record](murayama-sigma-pipi-loop-tree-and-legs.md). Units are GeV. The quantity shown for each graph is its contribution to $M_1^{\rm est}$,
$$
X_\Gamma^{\max}=\sum_tC_t^{\rho_t^*}\,t,\qquad \sum_tC_t^\rho\,t=\sum_{\rm terms}S_\Gamma\,\frac{(m_{\rm term}^2)^{2-r_\Gamma}}{16\pi^2}\,T_{\rm term}^\rho ,
$$
of mass dimension 1, with $T_{\rm term}^\rho$ the product of vertices of one term in routing $\rho$ (mass dimension $2r_\Gamma-3$).

## Method

For each kept graph and each routing of the loop momentum, the product of the vertex rules is built with the internal legs of the lines. A meson line carries the residue $(\mathbf 1_9\oplus1)-v_av_a^{\rm T}$, split exactly into three pieces: the bare $U(3)$ generator sum, the ALP-field entry, and minus the ALP-like pole $v_a$ at both ends. A "meson" line is therefore not a pure meson line: which pieces contribute to each graph (bare generators, ALP-field entry $1-1/N(m_a^2)$, mixed entries, meson block of the pole) is listed in the [meson-line piece record](murayama-sigma-pipi-loop-meson-line-pieces.md). The generator sums are removed by $U(3)$ completeness, which follows from $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$ and $\lambda_0=\sqrt{2/3}\,\mathbf 1$:
$$
\sum_{i=0}^8\langle A\lambda_iB\lambda_i\rangle=2\langle A\rangle\langle B\rangle,\qquad \sum_{i=0}^8\langle A\lambda_i\rangle\langle\lambda_iB\rangle=2\langle AB\rangle .
$$
The rewrite is applied to the compact vertex expressions before any expansion. When the summed generator stands inside a sum over permutations, the trace is first brought to the form $\langle Z\,D\rangle$, with $Z$ the generator, using identities valid for any matrix $Z$: $\langle X\sum_{\rm perm}(Z\,y_2\cdots y_m)\rangle=\langle Z\sum_{\rm perm}(y_2\cdots y_m\,X)\rangle$, and for commutators $\langle[Z,y]R\rangle=\langle Z[y,R]\rangle$, $\langle\{Z,y\}R\rangle=\langle Z\{y,R\}\rangle$. Then the charge rule and the assumptions ($\widehat P_{0,0}=\widehat P_{1,0}=0$, $\langle\widehat S_{0,1}\rangle=\langle\widehat S_{1,1}\rangle=0$) are applied.

**D4 (exact split, restated criterion).** In four meson-only graphs (the closed meson line on the five-leg vertex, the bubble with two meson lines and $\sigma$ alone at one vertex, and the closed meson line on each $\pi^\pm$ leg) the term in which every meson line takes the bare generator sum is removed, in every routing. This term has the scale $m_P$. The script checks the restated criterion of D4 on it, routing by routing: it contains no letter of the ALP sector ($\mathbb M^{-1}_{m_a^2}(\mathcal C)$, $N(m_a^2)$, $\widehat P_{0,1}$, $\widehat S_{r,s\ge1}$, $\widehat P_{r,s\ge1}$), and its ratio to $M_0$ (for a $\pi^\pm$ leg, its $\delta Z$) depends on $\epsilon$ only through the kinematics $p_\sigma^2=m_\sigma^2$; see the [selection record](murayama-sigma-pipi-loop-graph-selection.md), rule 3. All other pieces are kept. The test checks that the kept value equals the full value minus the removed term, exactly, in every routing.

**Loop scale per term (D3, assumption).** Every term of a graph, i.e. every choice of one piece on each meson line, has its own scale $m_{\rm term}$, the largest of the scales of its lines: $m_\sigma$ for a radial line, $m_a$ for an ALP-like line and for the ALP-field entry and the ALP-like pole of a meson line (together the pieces (ii)+(iii)+(iv) $=(0\oplus1)-v_av_a^{\rm T}$, the remainder of the ALP-like propagator), and $m_P$ for the bare generator sum. $m_P$ is a single symbolic mass, the largest pole mass of the nine meson-like states; no number is used for it. Bare meson lines enter the maximum also in graphs with a radial or an ALP-like line, so that the scales that occur are $m_\sigma$, $m_a$, $m_P$, $\max(m_P,m_a)$, $\max(m_P,m_\sigma)$, $\max(m_a,m_\sigma)$ and $\max(m_P,m_a,m_\sigma)$ (listed per graph below). In each term the products with the loop momentum are set to $k\cdot k\to m_{\rm term}^2$ and $k\cdot p\to m_{\rm term}^2$ after expanding, with the external products exact.

**Routing (D12).** For each of the $2r_\Gamma$ routings $\rho$ the size-rule value of the graph is written as $\sum_tC_t^\rho\,t$, with $t$ the independent products of traces of words (the canonical basis for generic sources). The graph enters $M_1^{\rm est}$ as $\sum_tC_t^{\rho_t^*}t$, with $\rho_t^*$ the routing of largest $|C_t^\rho|$; the sign is kept. This is an upper bound per graph on the routing ambiguity of the size rule, not on the true loop: $|\sum_tC_t^\rho t|\le\sum_t|C_t^{\rho_t^*}|\,|t|$ for every routing $\rho$. Ties of equal magnitude and opposite sign take the sign of the sum of the coefficients over all routings; if this sum vanishes, the first routing in the canonical order is used and the product is marked (a choice at a degenerate mass relation). The routings are numbered $\rho=1,\dots,2r_\Gamma$: $\rho=2j+1$ ($\rho=2j+2$) has the loop momentum on line $j$ along (against) the loop, where line $j$ joins vertex $j$ to vertex $j+1$ in the vertex list, starting from $0$.

**Mass-dependent choices.** In 80 products of the whole estimate (all in graphs 13, 14, 15, 20, 21, 25, 26) the largest routing depends on the masses. In each of them the routing coefficients have a single loop scale $s$ ($m_a$ or $\max(m_P,m_a)$) and are, up to a common factor, affine in one ratio $w=(m_\sigma^2+\gamma m_{\pi^\pm}^2)/s^2$. The comparison is done exactly on $w>0$ (which holds because $m_\sigma\ge2m_{\pi^\pm}$, the decay being open), and the coefficient is written as a Piecewise expression in homogeneous inequalities between squared masses. Three patterns occur:
- bubbles 13, 14, 15, $\gamma=-1$, $s=m_a$ or $\max(m_P,m_a)$: routing 3 for $m_\sigma^2-m_{\pi^\pm}^2\le s^2$ and routing 4 for $m_\sigma^2-m_{\pi^\pm}^2>s^2$ (59 products; at the boundary the two routings are equal and opposite and their sum vanishes, so routing 3 is taken and the product is marked: 59 marks, all at this boundary point);
- bubbles 13, 14, 15, $\gamma=-1$: $s=m_a$: routing 3 for $m_\sigma^2-m_{\pi^\pm}^2\le2s^2$ or $\ge5s^2$, routing 4 in between (9 products, in which routing 2 gives exactly zero);
- triangles 20, 21, 25, 26, $\gamma=-2$, $s=\max(m_P,m_a)$: routing 1 for $m_\sigma^2-2m_{\pi^\pm}^2\le2s^2$ and routing 2 above (12 products).
No product needs the tie rule on an open region.

**Display.** For a graph whose maximum is attained by one routing for every product, that routing is shown in compact form. Otherwise the compact form of the routing that attains the maximum for most products is shown, followed by the correction $\sum_t(C_t^{\rho_t^*}-C_t^{\rho})\,t$ in products of traces of words, where pairs of words that differ by $\lambda_{\pi^+}\lambda_{\pi^-}\leftrightarrow\lambda_{\pi^-}\lambda_{\pi^+}$ with equal coefficients are written with the anticommutator. Piecewise coefficients are printed as a common factor times the cases (the "otherwise" branch is the complement of the conditions above it).

**Routing spread (diagnostic).** For the products with numerical ratios, the record gives the range of $C_t^\rho/C_t^{\rho_t^*}$ over the routings. The routing ambiguity of the size rule is covered by the maximum.

Notation: $\sum_{\rm perm}(y_1\cdots y_m)$ is the sum over all $m!$ orderings of the product $y_1\cdots y_m$ (identical factors are counted as distinct), $[\cdot,\cdot]$ and $\{\cdot,\cdot\}$ are the commutator and anticommutator, and $\langle\cdot\rangle$ is the flavour trace. No loop value appears. Everything is exact in $N(m_a^2)$, $\mathbb M^{-1}_{m_a^2}(\mathcal C)$, $\epsilon$ and the source traces; there is no expansion in $\epsilon$. The large-$\epsilon$ powers of the graphs are in the [large-$\epsilon$ graph record](murayama-sigma-pipi-loop-large-epsilon-graphs.md).

## Verification

Construction 1 expands every vertex into words and applies completeness to the words. Construction 2 works on the compact expressions, as above. For every kept graph and every routing the two agree as exact identities of cyclic words. A separate test confirms the two completeness identities once, with the explicit $U(3)$ generators and two symbolic $3\times3$ matrices; it also checks $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$, $\lambda_0=\sqrt{2/3}\,\mathbf 1$, and that the octet sum alone fails. These matrices are never printed. Further tests (graph and trace scripts): graphs that do not depend on the routing enter with their single-routing value; $|C_t^{\rho_t^*}|\ge|C_t^\rho|$ for every routing and product, exactly for numerical ratios and branch by branch for the mass-dependent ones (checked on the original coefficients with $m_\sigma^2=ws^2-\gamma m_{\pi^\pm}^2$); the three branch patterns; with a common routing and the common scale $m_a$ (and also with the per-term scales), an ALP-like line plus pieces (ii)+(iii)+(iv) of a meson line equals the bare ALP-field line, for every graph with an ALP-like line; the D4 removal is exact; the displayed form (base routing plus correction) reproduces the maximum; and the closed meson line on the five-leg vertex plus the closed ALP-like line equals the bare ALP-field line alone. All tests use exact rationals and generic noncommuting symbols, and the output contains no generator, no matrix index and no element of the unit-matrix basis.

## Routing-dependent 1PI graphs (diagnostic)

Graphs 13, 14 and 15 (the bubbles with $\pi^+\pi^-$ at one vertex and two ALP-like, ALP-like and meson, or two meson lines): 20 products with numerical ratios $C_t^\rho/C_t^{\rho^*}$ between $-1/3$ and $1$ (some routings give exactly zero), and 12, 28, 28 products with a mass-dependent largest routing. Graphs 20, 21, 25 and 26 (the triangles with a meson line between $\pi^+$ and $\pi^-$ and two lines of ALP-like or meson type at $\sigma$): 3 products each with a mass-dependent largest routing (routing 1 or 2). No single routing attains the maximum for all products of these seven graphs. The other 19 graphs do not depend on the routing.

## Result

**Graph 1.** Vertices $(\pi^+,\pi^-,\sigma)$; lines radial (closed on the vertex); $S_\Gamma=\frac{1}{2}$, $r_\Gamma=1$, $m_{\rm term}\in\{m_\sigma\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(- \frac{3 B_{0} \epsilon^{3} m_\sigma^{2}}{32 \pi^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{3,0}\rangle
\end{aligned}
$$

**Graph 2.** Vertices $(\pi^+,\pi^-,\sigma)$; lines ALP-like (closed on the vertex); $S_\Gamma=\frac{1}{2}$, $r_\Gamma=1$, $m_{\rm term}\in\{m_a\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(- \frac{B_{0} \epsilon^{3} m_a^{2}}{32 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{1,2}\rangle \\
&+\left(- \frac{B_{0} \epsilon^{3} m_a^{2}}{48 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat P_{1,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} m_a^{2}}{192 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat S_{1,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle
\end{aligned}
$$

**Graph 3.** Vertices $(\pi^+,\pi^-,\sigma)$; lines meson (closed on the vertex); $S_\Gamma=\frac{1}{2}$, $r_\Gamma=1$, $m_{\rm term}\in\{m_a\}$; D4 term removed; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(- \frac{B_{0} \epsilon^{3} m_a^{2} \left(N(m_a^2) - 1\right)}{32 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{1,2}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} m_a^{2}}{48 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat P_{1,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0} \epsilon^{3} m_a^{2}}{192 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat S_{1,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle
\end{aligned}
$$

**Graph 4.** Vertices $(\pi^+,\sigma)$, $(\pi^-)$; lines meson, radial; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(\frac{B_{0}^{2} \epsilon^{3}}{16 \pi^{2} f_\pi^{3}}\right)\langle \{\widehat S_{1,0},\lambda_{\pi^-}\}\,\{\widehat S_{2,0},\lambda_{\pi^+}\}\rangle
\end{aligned}
$$

**Graph 5.** Vertices $(\pi^+,\sigma)$, $(\pi^-)$; lines ALP-like, meson; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(\frac{B_{0}^{2} \epsilon^{3}}{32 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^-}\}\,\{\widehat S_{1,1},\lambda_{\pi^+}\}\rangle
\end{aligned}
$$

**Graph 6.** Vertices $(\pi^+,\sigma)$, $(\pi^-)$; lines meson, meson; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(\frac{B_{0}^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)}{32 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^-}\}\,\{\widehat S_{1,1},\lambda_{\pi^+}\}\rangle
\end{aligned}
$$

**Graph 7.** Vertices $(\pi^+)$, $(\pi^-,\sigma)$; lines meson, radial; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(\frac{B_{0}^{2} \epsilon^{3}}{16 \pi^{2} f_\pi^{3}}\right)\langle \{\widehat S_{1,0},\lambda_{\pi^+}\}\,\{\widehat S_{2,0},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

**Graph 8.** Vertices $(\pi^+)$, $(\pi^-,\sigma)$; lines ALP-like, meson; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(\frac{B_{0}^{2} \epsilon^{3}}{32 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{1,1},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

**Graph 9.** Vertices $(\pi^+)$, $(\pi^-,\sigma)$; lines meson, meson; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(\frac{B_{0}^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)}{32 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{1,1},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

**Graph 10.** Vertices $(\pi^+,\pi^-)$, $(\sigma)$; lines radial, radial; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{m_\sigma\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(- \frac{3 B_{0}^{2} \epsilon^{5}}{32 \pi^{2} f_\pi^{3}}\right)\langle \widehat S_{1,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{2,0}\rangle \\
&+\left(\frac{3 B_{0}^{2} \epsilon^{5}}{16 \pi^{2} f_\pi^{3}}\right)\langle \widehat S_{2,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{2,0}\rangle \\
&+\left(- \frac{3 B_{0}^{2} \epsilon^{5}}{16 \pi^{2} f_\pi^{3}}\right)\langle \widehat S_{3,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{2,0}\rangle \\
&+\left(\frac{3 B_{0} \epsilon^{3} m_\sigma^{2}}{32 \pi^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{2,0}\rangle
\end{aligned}
$$

**Graph 11.** Vertices $(\pi^+,\pi^-)$, $(\sigma)$; lines ALP-like, radial; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(- \frac{B_{0}^{2} \epsilon^{5}}{16 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat S_{2,1}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{1,1}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5}}{8 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{2,0}\rangle
\end{aligned}
$$

**Graph 12.** Vertices $(\pi^+,\pi^-)$, $(\sigma)$; lines meson, radial; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(- \frac{B_{0}^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat S_{2,1}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{1,1}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5}}{8 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{2,0}\rangle
\end{aligned}
$$

**Graph 13.** Vertices $(\pi^+,\pi^-)$, $(\sigma)$; lines ALP-like, ALP-like; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{m_a\}$; routing dependent in 32 of 92 products of traces of words; 20 with numerical ratios: $C_t^\rho/C_t^{\rho^*}$ from $- \frac{1}{3}$ to $1$; 12 with a mass-dependent largest routing (Piecewise); 9 decided by the first routing at a boundary point (sum over routings zero).
$$
\begin{aligned}
X_\Gamma^{\max}\big|_{\rho=3}={}&\left(- \frac{B_{0}^{2} \epsilon^{5}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{1,2}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,2}\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{1,2}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5}}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{1,2}\rangle\langle \widehat P_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5}}{96 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{1,2}\rangle\langle \widehat S_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{3 B_{0} \epsilon^{3} m_a^{2}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,2}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{\epsilon^{3} m_a^{2} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_{0} \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5}}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \widehat P_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \widehat S_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5}}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle\langle \widehat P_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle\langle \widehat S_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0} \epsilon^{3} m_a^{2}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat P_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} m_a^{2}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle
\end{aligned}
$$
The routing $\rho=3$ attains the maximum for most products; the other products are corrected by
$$
\begin{aligned}
X_\Gamma^{\max}-X_\Gamma^{\max}\big|_{\rho=3}={}&\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \widehat S_{1,2}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{6} & \text{otherwise} \end{cases}\right)\right)\langle \widehat S_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{\epsilon^{3} m_a^{2} \left(2 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \frac{\left(m_{\pi^\pm}\right)^{2}}{2} + \frac{m_\sigma^{2}}{2} \\\frac{1}{6} & \text{for}\: m_a^{2} > - \frac{\left(m_{\pi^\pm}\right)^{2}}{5} + \frac{m_\sigma^{2}}{5} \\0 & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{6} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{6} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{3} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{3} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{\epsilon^{3} m_a^{2} \left(2 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \frac{\left(m_{\pi^\pm}\right)^{2}}{2} + \frac{m_\sigma^{2}}{2} \\- \frac{1}{3} & \text{for}\: m_a^{2} > - \frac{\left(m_{\pi^\pm}\right)^{2}}{5} + \frac{m_\sigma^{2}}{5} \\0 & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle
\end{aligned}
$$

**Graph 14.** Vertices $(\pi^+,\pi^-)$, $(\sigma)$; lines ALP-like, meson; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{m_a,\max\left(m_P, m_a\right)\}$; routing dependent in 48 of 192 products of traces of words; 20 with numerical ratios: $C_t^\rho/C_t^{\rho^*}$ from $- \frac{1}{3}$ to $1$; 28 with a mass-dependent largest routing (Piecewise); 25 decided by the first routing at a boundary point (sum over routings zero).
$$
\begin{aligned}
X_\Gamma^{\max}\big|_{\rho=3}={}&\left(- \frac{B_{0}^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{1,2}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,2}\rangle \\
&+\left(- \frac{B_{0} \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{1,2}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5} \left(N(m_a^2) - 2\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{1,2}\rangle\langle \widehat P_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{1,2}\rangle\langle \widehat S_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\,[[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}],\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\,[[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}],\lambda_{\pi^+}]\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{3}}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{0,0}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{3}}{48 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\,\sum_{\rm perm}\big(\widehat P_{0,1}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{3 B_{0} \epsilon^{3} m_a^{2} \left(N(m_a^2) - 1\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,2}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5}}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5} \left(N(m_a^2) - 2\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{B_{0} \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle [[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}],\lambda_{\pi^-}]\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{B_{0} \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle [[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}],\lambda_{\pi^+}]\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{\epsilon^{3} m_a^{2} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(- \frac{B_{0} \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5}}{6 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \widehat P_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5}}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \widehat S_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5} \left(N(m_a^2) - 4\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle\langle \widehat P_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5}}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle\langle \widehat S_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0} \epsilon^{3} m_a^{2} \left(N(m_a^2) - 2\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat P_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{3}}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat P_{1,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{0,0}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{3}}{48 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat P_{1,1}\,\sum_{\rm perm}\big(\widehat P_{0,1}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0} \epsilon^{3} m_a^{2}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle
\end{aligned}
$$
The routing $\rho=3$ attains the maximum for most products; the other products are corrected by
$$
\begin{aligned}
X_\Gamma^{\max}-X_\Gamma^{\max}\big|_{\rho=3}={}&\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{6} & \text{otherwise} \end{cases}\right)\right)\langle \widehat S_{1,2}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{3} & \text{otherwise} \end{cases}\right)\right)\langle \widehat S_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{\epsilon^{3} m_a^{2} \left(2 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \frac{\left(m_{\pi^\pm}\right)^{2}}{2} + \frac{m_\sigma^{2}}{2} \\- \frac{1}{3} & \text{for}\: m_a^{2} > - \frac{\left(m_{\pi^\pm}\right)^{2}}{5} + \frac{m_\sigma^{2}}{5} \\0 & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{3} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{3} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{24} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{24} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{24} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat P_{1,1}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{24} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{1,0}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{2}{3} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{2}{3} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\,\lambda_{\pi^-}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\rangle \\
&+\left(\frac{\epsilon^{3} m_a^{2} \left(2 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \frac{\left(m_{\pi^\pm}\right)^{2}}{2} + \frac{m_\sigma^{2}}{2} \\\frac{2}{3} & \text{for}\: m_a^{2} > - \frac{\left(m_{\pi^\pm}\right)^{2}}{5} + \frac{m_\sigma^{2}}{5} \\0 & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\widehat P_{1,1}\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\widehat S_{1,0}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\,\widehat P_{1,1}\,\lambda_{\pi^+}\rangle
\end{aligned}
$$

**Graph 15.** Vertices $(\pi^+,\pi^-)$, $(\sigma)$; lines meson, meson; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{m_a,\max\left(m_P, m_a\right)\}$; D4 term removed; routing dependent in 48 of 192 products of traces of words; 20 with numerical ratios: $C_t^\rho/C_t^{\rho^*}$ from $- \frac{1}{3}$ to $1$; 28 with a mass-dependent largest routing (Piecewise); 25 decided by the first routing at a boundary point (sum over routings zero).
$$
\begin{aligned}
X_\Gamma^{\max}\big|_{\rho=3}={}&\left(- \frac{B_{0}^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)^{2}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{1,2}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,2}\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{1,2}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{1,2}\rangle\langle \widehat P_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5}}{96 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{1,2}\rangle\langle \widehat S_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0} \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\,[[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}],\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_{0} \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\,[[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}],\lambda_{\pi^+}]\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{3}}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{0,0}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{3}}{48 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\,\sum_{\rm perm}\big(\widehat P_{0,1}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{3 B_{0} \epsilon^{3} m_a^{2} \left(N(m_a^2) - 1\right)^{2}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,2}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle [[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}],\lambda_{\pi^-}]\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle [[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}],\lambda_{\pi^+}]\,\widehat P_{1,1}\rangle \\
&+\left(\frac{\epsilon^{3} m_a^{2} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_{0} \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5}}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \widehat P_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{5}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \widehat S_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5} \left(N(m_a^2) - 2\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle\langle \widehat P_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{5}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle\langle \widehat S_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} m_a^{2} \left(N(m_a^2) - 1\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat P_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{3}}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat P_{1,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{0,0}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)}{48 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat P_{1,1}\,\sum_{\rm perm}\big(\widehat P_{0,1}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} m_a^{2}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle
\end{aligned}
$$
The routing $\rho=3$ attains the maximum for most products; the other products are corrected by
$$
\begin{aligned}
X_\Gamma^{\max}-X_\Gamma^{\max}\big|_{\rho=3}={}&\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \widehat S_{1,2}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{6} & \text{otherwise} \end{cases}\right)\right)\langle \widehat S_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{\epsilon^{3} m_a^{2} \left(2 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \frac{\left(m_{\pi^\pm}\right)^{2}}{2} + \frac{m_\sigma^{2}}{2} \\\frac{1}{6} & \text{for}\: m_a^{2} > - \frac{\left(m_{\pi^\pm}\right)^{2}}{5} + \frac{m_\sigma^{2}}{5} \\0 & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{6} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{6} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{24} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{24} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{24} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat P_{1,1}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{24} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{1,0}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{3} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{5} m_a^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{3} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\,\lambda_{\pi^-}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\rangle \\
&+\left(\frac{\epsilon^{3} m_a^{2} \left(2 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: m_a^{2} \geq - \frac{\left(m_{\pi^\pm}\right)^{2}}{2} + \frac{m_\sigma^{2}}{2} \\- \frac{1}{3} & \text{for}\: m_a^{2} > - \frac{\left(m_{\pi^\pm}\right)^{2}}{5} + \frac{m_\sigma^{2}}{5} \\0 & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\widehat P_{1,1}\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\- \frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\widehat S_{1,0}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_{0} \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{\pi^{2} N(m_a^2) f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} \\\frac{1}{12} & \text{otherwise} \end{cases}\right)\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^-}\,\widehat P_{1,1}\,\lambda_{\pi^+}\rangle
\end{aligned}
$$

**Graph 16.** Vertices $(\pi^+)$, $(\pi^-)$, $(\sigma)$; lines $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: radial; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(\frac{3 B_{0}^{3} \epsilon^{5}}{32 \pi^{2} f_\pi^{3} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \widehat S_{1,0}\rangle\langle \{\widehat S_{1,0},\lambda_{\pi^+}\}\,\{\widehat S_{1,0},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{3 B_{0}^{3} \epsilon^{5}}{16 \pi^{2} f_\pi^{3} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \widehat S_{2,0}\rangle\langle \{\widehat S_{1,0},\lambda_{\pi^+}\}\,\{\widehat S_{1,0},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{3 B_{0}^{3} \epsilon^{5}}{16 \pi^{2} f_\pi^{3} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \widehat S_{3,0}\rangle\langle \{\widehat S_{1,0},\lambda_{\pi^+}\}\,\{\widehat S_{1,0},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{3 B_{0}^{2} \epsilon^{3} m_\sigma^{2}}{32 \pi^{2} f_\pi^{3} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \{\widehat S_{1,0},\lambda_{\pi^+}\}\,\{\widehat S_{1,0},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

**Graph 17.** Vertices $(\pi^+)$, $(\pi^-)$, $(\sigma)$; lines $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: radial; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(\frac{B_{0}^{3} \epsilon^{5}}{16 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \widehat S_{2,1}\rangle\langle \{\widehat S_{0,1},\lambda_{\pi^-}\}\,\{\widehat S_{1,0},\lambda_{\pi^+}\}\rangle \\
&+\left(\frac{B_{0}^{3} \epsilon^{5}}{8 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^-}\}\,\{\widehat S_{1,0},\lambda_{\pi^+}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{2,0}\rangle
\end{aligned}
$$

**Graph 18.** Vertices $(\pi^+)$, $(\pi^-)$, $(\sigma)$; lines $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: radial; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(\frac{B_{0}^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \widehat S_{2,1}\rangle\langle \{\widehat S_{0,1},\lambda_{\pi^-}\}\,\{\widehat S_{1,0},\lambda_{\pi^+}\}\rangle \\
&+\left(- \frac{B_{0}^{3} \epsilon^{5}}{8 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^-}\}\,\{\widehat S_{1,0},\lambda_{\pi^+}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{2,0}\rangle
\end{aligned}
$$

**Graph 19.** Vertices $(\pi^+)$, $(\pi^-)$, $(\sigma)$; lines $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: ALP-like; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(\frac{B_{0}^{3} \epsilon^{5}}{16 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \widehat S_{2,1}\rangle\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{1,0},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_{0}^{3} \epsilon^{5}}{8 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{1,0},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{2,0}\rangle
\end{aligned}
$$

**Graph 20.** Vertices $(\pi^+)$, $(\pi^-)$, $(\sigma)$; lines $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: ALP-like; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing dependent in 3 of 12 products of traces of words; 3 with a mass-dependent largest routing (Piecewise).
$$
\begin{aligned}
X_\Gamma^{\max}\big|_{\rho=1}={}&\left(\frac{B_{0}^{3} \epsilon^{5}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{3} \left(2 \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 8 \max\left(m_P, m_a\right)^{2}\right)}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_{0}^{3} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0}^{3} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle
\end{aligned}
$$
The routing $\rho=1$ attains the maximum for most products; the other products are corrected by
$$
\begin{aligned}
X_\Gamma^{\max}-X_\Gamma^{\max}\big|_{\rho=1}={}&\left(\frac{B_{0}^{2} \epsilon^{3}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + \frac{m_\sigma^{2}}{2} \\- \frac{3}{8} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,1}\,\widehat S_{0,1}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{3}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + \frac{m_\sigma^{2}}{2} \\- \frac{3}{4} & \text{otherwise} \end{cases}\right)\right)\langle \widehat S_{0,1}\,\lambda_{\pi^+}\,\widehat S_{0,1}\,\lambda_{\pi^-}\rangle
\end{aligned}
$$

**Graph 21.** Vertices $(\pi^+)$, $(\pi^-)$, $(\sigma)$; lines $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: ALP-like; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing dependent in 3 of 12 products of traces of words; 3 with a mass-dependent largest routing (Piecewise).
$$
\begin{aligned}
X_\Gamma^{\max}\big|_{\rho=1}={}&\left(\frac{B_{0}^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{3} \left(N(m_a^2) - 1\right) \left(2 \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 8 \max\left(m_P, m_a\right)^{2}\right)}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_{0}^{3} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0}^{3} \epsilon^{5} \left(N(m_a^2) - 2\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle
\end{aligned}
$$
The routing $\rho=1$ attains the maximum for most products; the other products are corrected by
$$
\begin{aligned}
X_\Gamma^{\max}-X_\Gamma^{\max}\big|_{\rho=1}={}&\left(- \frac{B_{0}^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + \frac{m_\sigma^{2}}{2} \\\frac{3}{8} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,1}\,\widehat S_{0,1}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + \frac{m_\sigma^{2}}{2} \\\frac{3}{4} & \text{otherwise} \end{cases}\right)\right)\langle \widehat S_{0,1}\,\lambda_{\pi^+}\,\widehat S_{0,1}\,\lambda_{\pi^-}\rangle
\end{aligned}
$$

**Graph 22.** Vertices $(\pi^+)$, $(\pi^-)$, $(\sigma)$; lines $\pi^+$–$\pi^-$: radial; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: meson; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(- \frac{B_{0}^{3} \epsilon^{3}}{32 \pi^{2} f_\pi^{3} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \{\{\widehat S_{1,0},\lambda_{\pi^-}\},\widehat S_{1,0}\}\,\{\widehat S_{1,0},\lambda_{\pi^+}\}\rangle
\end{aligned}
$$

**Graph 23.** Vertices $(\pi^+)$, $(\pi^-)$, $(\sigma)$; lines $\pi^+$–$\pi^-$: ALP-like; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: meson; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(- \frac{B_{0}^{3} \epsilon^{3}}{32 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\{\widehat S_{0,1},\lambda_{\pi^-}\},\widehat S_{1,0}\}\,\{\widehat S_{0,1},\lambda_{\pi^+}\}\rangle
\end{aligned}
$$

**Graph 24.** Vertices $(\pi^+)$, $(\pi^-)$, $(\sigma)$; lines $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: meson; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_\Gamma^{\max}={}&\left(\frac{B_{0}^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \widehat S_{2,1}\rangle\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{1,0},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_{0}^{3} \epsilon^{5}}{8 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{1,0},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{2,0}\rangle
\end{aligned}
$$

**Graph 25.** Vertices $(\pi^+)$, $(\pi^-)$, $(\sigma)$; lines $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: meson; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing dependent in 3 of 12 products of traces of words; 3 with a mass-dependent largest routing (Piecewise).
$$
\begin{aligned}
X_\Gamma^{\max}\big|_{\rho=1}={}&\left(\frac{B_{0}^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{3} \left(N(m_a^2) - 1\right) \left(2 \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 8 \max\left(m_P, m_a\right)^{2}\right)}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_{0}^{3} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0}^{3} \epsilon^{5} \left(N(m_a^2) - 2\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle
\end{aligned}
$$
The routing $\rho=1$ attains the maximum for most products; the other products are corrected by
$$
\begin{aligned}
X_\Gamma^{\max}-X_\Gamma^{\max}\big|_{\rho=1}={}&\left(- \frac{B_{0}^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + \frac{m_\sigma^{2}}{2} \\\frac{3}{8} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,1}\,\widehat S_{0,1}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + \frac{m_\sigma^{2}}{2} \\\frac{3}{4} & \text{otherwise} \end{cases}\right)\right)\langle \widehat S_{0,1}\,\lambda_{\pi^+}\,\widehat S_{0,1}\,\lambda_{\pi^-}\rangle
\end{aligned}
$$

**Graph 26.** Vertices $(\pi^+)$, $(\pi^-)$, $(\sigma)$; lines $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: meson; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing dependent in 3 of 20 products of traces of words; 3 with a mass-dependent largest routing (Piecewise).
$$
\begin{aligned}
X_\Gamma^{\max}\big|_{\rho=1}={}&\left(\frac{B_{0}^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)^{2}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_{0}^{3} \epsilon^{3} \left(N(m_a^2) - 1\right)}{32 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\{\widehat S_{0,1},\lambda_{\pi^-}\},\widehat S_{1,0}\}\,\{\widehat S_{0,1},\lambda_{\pi^+}\}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)^{2} \left(2 \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 8 \max\left(m_P, m_a\right)^{2}\right)}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_{0}^{3} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(- \frac{B_{0}^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle
\end{aligned}
$$
The routing $\rho=1$ attains the maximum for most products; the other products are corrected by
$$
\begin{aligned}
X_\Gamma^{\max}-X_\Gamma^{\max}\big|_{\rho=1}={}&\left(\frac{B_{0}^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + \frac{m_\sigma^{2}}{2} \\- \frac{3}{8} & \text{otherwise} \end{cases}\right)\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S_{0,1}\,\widehat S_{0,1}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)^{2}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}} \left(\begin{cases} 0 & \text{for}\: \max\left(m_P, m_a\right)^{2} \geq - \left(m_{\pi^\pm}\right)^{2} + \frac{m_\sigma^{2}}{2} \\- \frac{3}{4} & \text{otherwise} \end{cases}\right)\right)\langle \widehat S_{0,1}\,\lambda_{\pi^+}\,\widehat S_{0,1}\,\lambda_{\pi^-}\rangle
\end{aligned}
$$

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_traces.py
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_traces.py --markdown
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_traces test_derive_sigma_pipi_loop_graphs
```

[Script](../../../scripts/derive_sigma_pipi_loop_traces.py), [tests](../../../scripts/test_derive_sigma_pipi_loop_traces.py), [graph script](../../../scripts/derive_sigma_pipi_loop_graphs.py) and its [tests](../../../scripts/test_derive_sigma_pipi_loop_graphs.py), [shared algebra](../../../scripts/sigma_pipi_loop_algebra.py). The `--markdown` option prints the result blocks of this record and of the [external-leg record](murayama-sigma-pipi-loop-traces-external-legs.md).
