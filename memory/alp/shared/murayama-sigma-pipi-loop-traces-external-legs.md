---
title: Compact flavour traces and residue corrections of the external legs in the one-loop sigma to pi+ pi- estimate
status: derived
created: 2026-11-26
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Compact traces and residue corrections of the external legs

Tasks 3.1, 3.2 (external legs) and 4.4 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama`. No loop integral, loop function or loop amplitude is evaluated. The selection holds for any $\epsilon$; no step assumes $\epsilon\ll1$.

## Inputs and assumptions

The 6 kept $\sigma$ self-energies and the 5 kept self-energies of each of $\pi^\pm$ come from the [graph-list record](murayama-sigma-pipi-loop-graph-list.md), with the rules of the [selection record](murayama-sigma-pipi-loop-graph-selection.md). The vertices come from the [vertex-rule record](murayama-sigma-pipi-loop-vertex-rules.md), the pole vector $v_a$ and the meson-line residue from the [large-$\epsilon$ record](murayama-sigma-pipi-loop-large-epsilon-scaling.md), and the legs and assumptions from the [tree and legs record](murayama-sigma-pipi-loop-tree-and-legs.md). The two legs of a self-energy carry incoming momenta $p$ and $-p$, and $p^2$ is kept symbolic until the end. Self-energies whose $\delta Z$ is exactly zero (tadpoles with a vertex that does not depend on $p$) only shift the masses, and the physical pole masses are used. One-loop tadpoles (one-point functions): those of $\pi^0$, $\eta$, $\eta'$ and of the ALP field vanish under CP (D11); the $\sigma$ tadpole is set to zero by assumption (a vacuum condition on $f_a$, consistent with eliminating $\lambda$ through $m_\sigma^2$). The one-particle-reducible graphs in which an external leg mixes with another state vanish under the CP and flavour-diagonal assumptions of D11 and are not included.

## Method

For each kept graph and each routing of the loop momentum, the product of the vertex rules is built with the internal legs of the lines. A meson line carries the residue $(\mathbf 1_9\oplus1)-v_av_a^{\rm T}$, split exactly into three pieces: the bare $U(3)$ generator sum, the ALP-field entry, and minus the ALP-like pole $v_a$ at both ends. A "meson" line is therefore not a pure meson line: which pieces contribute to each graph (bare generators, ALP-field entry $1-1/N(m_a^2)$, mixed entries, meson block of the pole) is listed in the [meson-line piece record](murayama-sigma-pipi-loop-meson-line-pieces.md). The generator sums are removed by $U(3)$ completeness, which follows from $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$ and $\lambda_0=\sqrt{2/3}\,\mathbf 1$:
$$
\sum_{i=0}^8\langle A\lambda_iB\lambda_i\rangle=2\langle A\rangle\langle B\rangle,\qquad \sum_{i=0}^8\langle A\lambda_i\rangle\langle\lambda_iB\rangle=2\langle AB\rangle .
$$
The rewrite is applied to the compact vertex expressions before any expansion. When the summed generator stands inside a sum over permutations, the trace is first brought to the form $\langle Z\,D\rangle$, with $Z$ the generator, using identities valid for any matrix $Z$: $\langle X\sum_{\rm perm}(Z\,y_2\cdots y_m)\rangle=\langle Z\sum_{\rm perm}(y_2\cdots y_m\,X)\rangle$, and for commutators $\langle[Z,y]R\rangle=\langle Z[y,R]\rangle$, $\langle\{Z,y\}R\rangle=\langle Z\{y,R\}\rangle$. Then the charge rule and the assumptions ($\widehat P_{0,0}=\widehat P_{1,0}=0$, $\langle\widehat S_{0,1}\rangle=\langle\widehat S_{1,1}\rangle=0$) are applied.

**D4 (exact split, restated criterion).** In four meson-only graphs (the closed meson line on the five-leg vertex, the bubble with two meson lines and $\sigma$ alone at one vertex, and the closed meson line on each $\pi^\pm$ leg) the term in which every meson line takes the bare generator sum is removed, in every routing. This term has the scale $m_P$. The script checks the restated criterion of D4 on it, routing by routing: it contains no letter of the ALP sector ($\mathbb M^{-1}_{m_a^2}(\mathcal C)$, $N(m_a^2)$, $\widehat P_{0,1}$, $\widehat S_{r,s\ge1}$, $\widehat P_{r,s\ge1}$), and its ratio to $M_0$ (for a $\pi^\pm$ leg, its $\delta Z$) depends on $\epsilon$ only through the kinematics $p_\sigma^2=m_\sigma^2$; see the [selection record](murayama-sigma-pipi-loop-graph-selection.md), rule 3. All other pieces are kept. The test checks that the kept value equals the full value minus the removed term, exactly, in every routing.

**Loop scale per term (D3, assumption).** Every term of a graph, i.e. every choice of one piece on each meson line, has its own scale $m_{\rm term}$, the largest of the scales of its lines: $m_\sigma$ for a radial line, $m_a$ for an ALP-like line and for the ALP-field entry and the ALP-like pole of a meson line (together the pieces (ii)+(iii)+(iv) $=(0\oplus1)-v_av_a^{\rm T}$, the remainder of the ALP-like propagator), and $m_P$ for the bare generator sum. $m_P$ is a single symbolic mass, the largest pole mass of the nine meson-like states; no number is used for it. Bare meson lines enter the maximum also in graphs with a radial or an ALP-like line, so that the scales that occur are $m_\sigma$, $m_a$, $m_P$, $\max(m_P,m_a)$, $\max(m_P,m_\sigma)$, $\max(m_a,m_\sigma)$ and $\max(m_P,m_a,m_\sigma)$ (listed per graph below). In each term the products with the loop momentum are set to $k\cdot k\to m_{\rm term}^2$ and $k\cdot p\to m_{\rm term}^2$ after expanding, with the external products exact.

**Routing (D12).** For each of the $2r_\Gamma$ routings $\rho$ the size-rule value of the graph (for a self-energy, $\delta Z^\rho$, with the derivative rule below applied term by term in each routing) is written as $\sum_tC_t^\rho\,t$, with $t$ the independent products of traces of words (the canonical basis for generic sources). The graph enters $M_1^{\rm est}$ as $\sum_tC_t^{\rho_t^*}t$, with $\rho_t^*$ the routing of largest $|C_t^\rho|$; the sign is kept. This is an upper bound per graph on the routing ambiguity of the size rule, not on the true loop: $|\sum_tC_t^\rho t|\le\sum_t|C_t^{\rho_t^*}|\,|t|$ for every routing $\rho$. Ties of equal magnitude and opposite sign take the sign of the sum of the coefficients over all routings; if this sum vanishes, the first routing in the canonical order is used and the product is marked (a choice at a degenerate mass relation). The routings are numbered $\rho=1,\dots,2r_\Gamma$: $\rho=2j+1$ ($\rho=2j+2$) has the loop momentum on line $j$ along (against) the loop, where line $j$ joins vertex $j$ to vertex $j+1$ in the vertex list, starting from $0$.

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

## Residue prescription (assumption)

With the vertex convention $M=V$, the one-loop two-point function $\Pi(p^2)$ enters the inverse propagator as $p^2-m^2+\Pi(p^2)$. The residue is then $Z=1+\delta Z$ with $\delta Z=-\Pi'(m^2)$, and the external legs add
$$
M_0\Big(\tfrac12\delta Z_\sigma^{\max}+\tfrac12\delta Z_{\pi^+}^{\max}+\tfrac12\delta Z_{\pi^-}^{\max}\Big)
$$
to $M_1^{\rm est}$, where each $\delta Z^{\max}$ is the sum over the self-energies of that leg of their termwise routing maxima. For one term of a self-energy graph, with numerator $T_{\rm term}(p^2)$ and scale $m_{\rm term}$, the derivative is replaced by
$$
\delta Z_{\rm term}=-\frac{S_\Gamma}{16\pi^2}\Big[\theta(r_\Gamma\ge2)\,(m_{\rm term}^2)^{2-(r_\Gamma+1)}\,T_{\rm term}+(m_{\rm term}^2)^{2-r_\Gamma}\,\frac{\partial T_{\rm term}}{\partial p^2}\Big]_{p^2=m_{\rm ext}^2},
$$
and $\delta Z_\Gamma^\rho$ is the sum over the terms in routing $\rho$. The first term is the moment with one more propagator, from differentiating a denominator that depends on $p$. The second is the explicit $p^2$ dependence of the numerator; $k\cdot p\to m_{\rm term}^2$ does not contribute to it. For a tadpole no denominator depends on $p$, and only the second term is present. For the closed ALP-like line on the four-meson kinetic vertex it gives $\frac{\partial V(e,\bar e,a,a)}{\partial p^2}=\langle[\lambda_a,\lambda_e][\lambda_a,\lambda_{\bar e}]\rangle/(12f_\pi^2)$ with $\lambda_a=2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)/\sqrt{N(m_a^2)}$. This agrees with [the external-residue record](murayama-radial-loop-external-residues.md), which is checked in the tests. The overall sign $\delta Z=-\Pi'$ follows the vertex convention $M=V$; the size estimate does not fix the relative signs of the loops.

The meson-only $\sigma$ self-energy (graph σ6) has two $\sigma$ ends, so its all-generator-sum term is kept (D4 does not apply); with the scale $m_P$ this term gives $\delta Z\sim\epsilon^2$ at large $\epsilon$ ([large-$\epsilon$ graph record](murayama-sigma-pipi-loop-large-epsilon-graphs.md)). The closed meson line on each $\pi^\pm$ leg (graph π2) has its D4 term removed; what remains in $\delta Z$ is only the meson block of the ALP-like pole, $-4\epsilon^2\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\otimes\mathbb M^{-1}_{m_a^2}(\mathcal C)/N(m_a^2)$, because the ALP-field entry and the mixed entries give vertices without $p^2$ dependence ([meson-line piece record](murayama-sigma-pipi-loop-meson-line-pieces.md)). It has the scale $m_a$, the same as the closed ALP-like line (graph π1), and $\delta Z_{\pi1}+\delta Z_{\pi2}=0$ exactly (tested): no factor $m_a^2-m_\sigma^2$ appears. What remains of $\delta Z_{\pi^\pm}$ is the sum of the three bubbles π3, π4 and π5. This is the statement, at the level of the integrands, that with a common scale the ALP-like line plus pieces (ii)+(iii)+(iv) of a meson line is the bare ALP-field line, whose closed loop on the four-pion vertex has no $p^2$ dependence.

## Routing-dependent self-energies (diagnostic)

Graphs σ4, σ5 and σ6, the $\sigma$ self-energies with two ALP-like, ALP-like and meson, or two meson lines. Their kinetic $\sigma aa$ vertices depend on the momenta of the ALP-field ends. For four products (the term without traces, $\langle\widehat S_{1,2}\rangle$, $\langle\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle$ and $\langle\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle$) two of the four routings give zero and the other two give the maximum. Routing 2 attains the maximum for every product. No $\pi^\pm$ self-energy depends on the routing.

## Result

**Graph σ1** (self-energy of $\sigma$). Lines radial, radial; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{m_\sigma\}$; routing independent.
$$
\begin{aligned}
\delta Z_\Gamma^{\max}={}&\left(- \frac{9 \epsilon^{2} m_\sigma^{2}}{32 \pi^{2} f_\pi^{2}}\right) \\
&+\left(\frac{9 B_{0} \epsilon^{4}}{16 \pi^{2} f_\pi^{2}}\right)\langle \widehat S_{1,0}\rangle \\
&+\left(- \frac{9 B_{0}^{2} \epsilon^{6}}{32 \pi^{2} f_\pi^{2} m_\sigma^{2}}\right)\langle \widehat S_{1,0}\rangle\langle \widehat S_{1,0}\rangle \\
&+\left(\frac{9 B_{0}^{2} \epsilon^{6}}{8 \pi^{2} f_\pi^{2} m_\sigma^{2}}\right)\langle \widehat S_{1,0}\rangle\langle \widehat S_{2,0}\rangle \\
&+\left(- \frac{9 B_{0}^{2} \epsilon^{6}}{8 \pi^{2} f_\pi^{2} m_\sigma^{2}}\right)\langle \widehat S_{1,0}\rangle\langle \widehat S_{3,0}\rangle \\
&+\left(- \frac{9 B_{0} \epsilon^{4}}{8 \pi^{2} f_\pi^{2}}\right)\langle \widehat S_{2,0}\rangle \\
&+\left(- \frac{9 B_{0}^{2} \epsilon^{6}}{8 \pi^{2} f_\pi^{2} m_\sigma^{2}}\right)\langle \widehat S_{2,0}\rangle\langle \widehat S_{2,0}\rangle \\
&+\left(\frac{9 B_{0}^{2} \epsilon^{6}}{4 \pi^{2} f_\pi^{2} m_\sigma^{2}}\right)\langle \widehat S_{2,0}\rangle\langle \widehat S_{3,0}\rangle \\
&+\left(\frac{9 B_{0} \epsilon^{4}}{8 \pi^{2} f_\pi^{2}}\right)\langle \widehat S_{3,0}\rangle \\
&+\left(- \frac{9 B_{0}^{2} \epsilon^{6}}{8 \pi^{2} f_\pi^{2} m_\sigma^{2}}\right)\langle \widehat S_{3,0}\rangle\langle \widehat S_{3,0}\rangle
\end{aligned}
$$

**Graph σ2** (self-energy of $\sigma$). Lines ALP-like, radial; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
\delta Z_\Gamma^{\max}={}&\left(- \frac{B_{0}^{2} \epsilon^{6}}{4 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_a, m_\sigma\right)^{2}}\right)\langle \widehat S_{2,1}\rangle\langle \widehat S_{2,1}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{6}}{\pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_a, m_\sigma\right)^{2}}\right)\langle \widehat S_{2,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{2,0}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{6}}{\pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_a, m_\sigma\right)^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{2,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{2,0}\rangle
\end{aligned}
$$

**Graph σ3** (self-energy of $\sigma$). Lines meson, radial; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_\sigma\right),\max\left(m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
\delta Z_\Gamma^{\max}={}&\left(- \frac{B_{0}^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{4 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_a, m_\sigma\right)^{2}}\right)\langle \widehat S_{2,1}\rangle\langle \widehat S_{2,1}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{6}}{\pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_a, m_\sigma\right)^{2}}\right)\langle \widehat S_{2,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{2,0}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{6}}{\pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_a, m_\sigma\right)^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{2,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{2,0}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{4}}{2 \pi^{2} f_\pi^{2} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \widehat P_{2,0}\,\widehat P_{2,0}\rangle
\end{aligned}
$$

**Graph σ4** (self-energy of $\sigma$). Lines ALP-like, ALP-like; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{m_a\}$; routing dependent in 4 of 10 products of traces of words; 4 with numerical ratios: $C_t^\rho/C_t^{\rho^*}$ from $0$ to $1$; attained by routing 2 for every product.
$$
\begin{aligned}
\delta Z_\Gamma^{\max}={}&\left(- \frac{\epsilon^{2} m_a^{2}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2}}\right) \\
&+\left(- \frac{B_{0} \epsilon^{4}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2}}\right)\langle \widehat S_{1,2}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{6}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \widehat S_{1,2}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0} \epsilon^{4}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{6}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{B_{0} \epsilon^{4}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle
\end{aligned}
$$

**Graph σ5** (self-energy of $\sigma$). Lines ALP-like, meson; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{m_a,\max\left(m_P, m_a\right)\}$; routing dependent in 4 of 15 products of traces of words; 4 with numerical ratios: $C_t^\rho/C_t^{\rho^*}$ from $0$ to $1$; attained by routing 2 for every product.
$$
\begin{aligned}
\delta Z_\Gamma^{\max}={}&\left(- \frac{\epsilon^{2} m_a^{2} \left(N(m_a^2) - 1\right)}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2}}\right) \\
&+\left(- \frac{B_{0} \epsilon^{4} \left(N(m_a^2) - 1\right)}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2}}\right)\langle \widehat S_{1,2}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \widehat S_{1,2}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{6}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{6} \left(N(m_a^2) - 2\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{4}}{8 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\,\{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{4}}{4 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{2 B_{0} \epsilon^{4}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{6}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(- \frac{2 B_{0}^{2} \epsilon^{6}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{B_{0} \epsilon^{4} \left(N(m_a^2) - 2\right)}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{6} \left(N(m_a^2) - 4\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{4}}{8 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat P_{1,1}\,\widehat P_{1,1}\rangle
\end{aligned}
$$

**Graph σ6** (self-energy of $\sigma$). Lines meson, meson; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{m_P,m_a,\max\left(m_P, m_a\right)\}$; routing dependent in 4 of 17 products of traces of words; 4 with numerical ratios: $C_t^\rho/C_t^{\rho^*}$ from $0$ to $1$; attained by routing 2 for every product.
$$
\begin{aligned}
\delta Z_\Gamma^{\max}={}&\left(- \frac{\epsilon^{2} m_a^{2} \left(N(m_a^2) - 1\right)^{2}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2}}\right) \\
&+\left(- \frac{B_{0}^{2} \epsilon^{2}}{16 \pi^{2} f_\pi^{2} m_P^{2}}\right)\langle \widehat S_{1,0}\rangle\langle \widehat S_{1,0}\rangle \\
&+\left(- \frac{B_{0} \epsilon^{4} \left(N(m_a^2) - 1\right)^{2}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2}}\right)\langle \widehat S_{1,2}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)^{2}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \widehat S_{1,2}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{4}}{8 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\,\{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{4}}{4 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S_{1,0}\}\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0} \epsilon^{4}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{6}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0} \epsilon^{4} \left(N(m_a^2) - 1\right)}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(\frac{B_{0}^{2} \epsilon^{6} \left(N(m_a^2) - 2\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{B_{0}^{2} \epsilon^{4} \left(N(m_a^2) - 1\right)}{8 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat P_{1,1}\,\widehat P_{1,1}\rangle \\
&+\left(- \frac{3 B_{0}^{2} \epsilon^{2}}{16 \pi^{2} f_\pi^{2} m_P^{2}}\right)\langle \widehat S_{1,0}\,\widehat S_{1,0}\rangle
\end{aligned}
$$

**Graph π1** (self-energy of $\pi^+$). Lines ALP-like (closed on the vertex); $S_\Gamma=\frac{1}{2}$, $r_\Gamma=1$, $m_{\rm term}\in\{m_a\}$; routing independent.
$$
\begin{aligned}
\delta Z_\Gamma^{\max}={}&\left(- \frac{\epsilon^{2} m_a^{2}}{96 \pi^{2} N(m_a^2) f_\pi^{2}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle
\end{aligned}
$$

**Graph π2** (self-energy of $\pi^+$). Lines meson (closed on the vertex); $S_\Gamma=\frac{1}{2}$, $r_\Gamma=1$, $m_{\rm term}\in\{m_a\}$; D4 term removed; routing independent.
$$
\begin{aligned}
\delta Z_\Gamma^{\max}={}&\left(\frac{\epsilon^{2} m_a^{2}}{96 \pi^{2} N(m_a^2) f_\pi^{2}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle
\end{aligned}
$$

**Graph π3** (self-energy of $\pi^+$). Lines meson, radial; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
\delta Z_\Gamma^{\max}={}&\left(- \frac{B_{0}^{2} \epsilon^{2}}{32 \pi^{2} f_\pi^{2} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \{\widehat S_{1,0},\lambda_{\pi^+}\}\,\{\widehat S_{1,0},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

**Graph π4** (self-energy of $\pi^+$). Lines ALP-like, meson; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
\delta Z_\Gamma^{\max}={}&\left(- \frac{B_{0}^{2} \epsilon^{2}}{32 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

**Graph π5** (self-energy of $\pi^+$). Lines meson, meson; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
\delta Z_\Gamma^{\max}={}&\left(- \frac{B_{0}^{2} \epsilon^{2} \left(N(m_a^2) - 1\right)}{32 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S_{0,1},\lambda_{\pi^+}\}\,\{\widehat S_{0,1},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_traces.py --markdown
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_traces test_derive_sigma_pipi_loop_graphs
```

[Script](../../../scripts/derive_sigma_pipi_loop_traces.py), [tests](../../../scripts/test_derive_sigma_pipi_loop_traces.py), [graph script](../../../scripts/derive_sigma_pipi_loop_graphs.py) and its [tests](../../../scripts/test_derive_sigma_pipi_loop_graphs.py). The π graphs are those of $\pi^+$; the $\pi^-$ graphs follow from $\lambda_{\pi^+}\leftrightarrow\lambda_{\pi^-}$, and $\delta Z_{\pi^+}^{\max}=\delta Z_{\pi^-}^{\max}$ exactly.
