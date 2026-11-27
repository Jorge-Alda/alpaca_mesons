---
title: Cartesian frame - list of the one-loop sigma to pi+ pi- graphs and external self-energies with symmetry factors, loop scales and D4 terms
status: derived
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Cartesian frame: one-loop graphs of $\sigma\to\pi^+\pi^-$ and of the external self-energies

Task 3.1 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama-cartesian`. The graphs are enumerated by script from the Cartesian vertices; no polar graph is used. **No loop integral, scalar loop function, loop coefficient or loop amplitude is computed.** The estimate built on these graphs is an order of magnitude per graph class and can miss cancellations between graphs with the same $\epsilon$ scaling; neither it nor the polar estimate is authoritative. The comparison with the polar graph list is not made here (task 4.4).

## Inputs

- Rules R1–R10 (graph set, lines, four-piece split, D3 scales, size rule, completeness, D4, D11, tadpoles, routing, residues): [Cartesian selection record](cartesian-sigma-pipi-loop-graph-selection.md).
- Cartesian vertices with $\lambda$ eliminated through $m_\sigma^2$: [Cartesian vertex record](cartesian-sigma-pipi-loop-vertex-rules.md). Legs, pole vector and assumptions: [Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md).

## Method

For $V=1,2,3,4$ vertices on the loop, the script distributes the external legs over the vertices (each vertex gets at least one), assigns a line type ($s$, ALP-like, meson) to each line and keeps one representative per class under rotations and reflections of the cycle. $S_\Gamma=1/|{\rm automorphisms}|$ and $r_\Gamma=V$. For every graph and every routing it builds the product of the Cartesian vertices with the pieces of R2, removes the generator sums by completeness, applies the assumptions and the charge rule, and multiplies each term by its size factor (1PI) or takes its residue correction (self-energies). A graph "vanishes identically" if its value is zero in every routing before any size factor; the vertex that vanishes in each term is recorded. The D4 term is identified and tested with the criterion of R6.

## Result: topologies

- External legs $s$, $\pi^+$, $\pi^-$: one way to put them on one vertex, six ordered ways on two vertices, six on three, none on four. Hence the 1PI graphs are tadpoles (one five-leg vertex), bubbles (vertices with 3 and 4 legs) and triangles (three 3-leg vertices). **No box** exists (four vertices with an external leg each cannot be built from three legs) and **no 1PI bubble with a two-leg and a five-leg vertex** exists (a two-leg vertex is the quadratic Lagrangian, which is in the propagators); six such configurations would appear if two-leg vertices were allowed.
- Self-energies (legs in, out with momenta $p$, $-p$): tadpoles (one four-leg vertex) and bubbles (two 3-leg vertices).

## Result: counts per class

| class | graphs | kept | vanish identically | zero $\delta Z$ (mass shift only) | D4 term excluded |
|---|---|---|---|---|---|
| 1PI tadpoles | 3 | 3 | 0 | – | 1 ($C3$) |
| 1PI bubbles | 18 | 12 | 6 | – | 1 ($C15$) |
| 1PI triangles | 27 | 11 | 16 | – | 0 |
| $s$ self-energy tadpoles | 3 | 0 | 0 | 3 | 0 |
| $s$ self-energy bubbles | 6 | 6 | 0 | 0 | 0 |
| $\pi^\pm$ self-energy tadpoles | 3 | 2 | 0 | 1 | 1 ($C\pi3$) |
| $\pi^\pm$ self-energy bubbles | 6 | 3 | 3 | 0 | 0 |

In total 26 kept 1PI graphs, 6 kept $s$ self-energies and 5 kept self-energies for each of $\pi^+$ and $\pi^-$. The $\pi^-$ self-energies equal the $\pi^+$ ones graph by graph (tested), so one list $C\pi1$–$C\pi9$ serves both legs. Every graph that vanishes identically does so because one of its vertices has one charged pion and line ends without a generator (the pole vector $v_a$, $p$, $s$ or $\mathbb M^{-1}_{m_a^2}(\mathcal C)$); its trace has charge $\pm1$ and vanishes by the charge rule (the last column of the tables lists the vertex, $V(\text{legs})=0$).

## Result: graph tables

Line $j$ joins vertex $j$ to vertex $j+1$ in the order of the vertex column (from 0). "$m_{\rm term}$" lists the scales of the terms whose size-rule value is nonzero in some routing (D4 term removed). "Pieces per meson line" lists the pieces of R2 whose terms contribute in some routing. "Routing": whether the value depends on the routing of the loop momentum (details in the [1PI trace record](cartesian-sigma-pipi-loop-traces-1pi.md)); for vanishing graphs this column gives the vanishing vertex.

### One-particle-irreducible graphs

| graph | vertices (external legs) | lines | $S_\Gamma$ | $r_\Gamma$ | status | $m_{\rm term}$ | pieces per meson line | routing |
|---|---|---|---|---|---|---|---|---|
| $C1$ | $\{\pi^+,\pi^-,s\}$ | $s$ | $\frac{1}{2}$ | 1 | kept | $m_\sigma$ | – | independent |
| $C2$ | $\{\pi^+,\pi^-,s\}$ | ALP-like | $\frac{1}{2}$ | 1 | kept | $m_a$ | – | independent |
| $C3$ | $\{\pi^+,\pi^-,s\}$ | meson | $\frac{1}{2}$ | 1 | kept, D4 term excluded | $m_a$ | line 0: (ii),(iii),(iv) | independent |
| $C4$ | $\{\pi^+\}$, $\{\pi^-,s\}$ | $s$, $s$ | $\frac{1}{2}$ | 2 | vanishes identically | – | – | $V(\pi^+,s,s)=0$ |
| $C5$ | $\{\pi^+\}$, $\{\pi^-,s\}$ | ALP-like, $s$ | $1$ | 2 | vanishes identically | – | – | $V(\pi^+,s,v_a)=0$ |
| $C6$ | $\{\pi^+\}$, $\{\pi^-,s\}$ | ALP-like, ALP-like | $\frac{1}{2}$ | 2 | vanishes identically | – | – | $V(\pi^+,v_a,v_a)=0$ |
| $C7$ | $\{\pi^+\}$, $\{\pi^-,s\}$ | ALP-like, meson | $1$ | 2 | kept | $\max\left(m_P, m_a\right)$ | line 1: (i) | independent |
| $C8$ | $\{\pi^+\}$, $\{\pi^-,s\}$ | meson, $s$ | $1$ | 2 | kept | $\max\left(m_P, m_\sigma\right)$ | line 0: (i) | independent |
| $C9$ | $\{\pi^+\}$, $\{\pi^-,s\}$ | meson, meson | $\frac{1}{2}$ | 2 | kept | $\max\left(m_P, m_a\right)$ | line 0: (i),(ii); line 1: (i),(ii) | independent |
| $C10$ | $\{\pi^+,\pi^-\}$, $\{s\}$ | $s$, $s$ | $\frac{1}{2}$ | 2 | kept | $m_\sigma$ | – | independent |
| $C11$ | $\{\pi^+,\pi^-\}$, $\{s\}$ | ALP-like, $s$ | $1$ | 2 | kept | $\max\left(m_a, m_\sigma\right)$ | – | independent |
| $C12$ | $\{\pi^+,\pi^-\}$, $\{s\}$ | ALP-like, ALP-like | $\frac{1}{2}$ | 2 | kept | $m_a$ | – | depends on the routing |
| $C13$ | $\{\pi^+,\pi^-\}$, $\{s\}$ | ALP-like, meson | $1$ | 2 | kept | $m_a$, $\max\left(m_P, m_a\right)$ | line 1: (i),(ii),(iii),(iv) | depends on the routing |
| $C14$ | $\{\pi^+,\pi^-\}$, $\{s\}$ | meson, $s$ | $1$ | 2 | kept | $\max\left(m_a, m_\sigma\right)$ | line 0: (ii),(iii) | independent |
| $C15$ | $\{\pi^+,\pi^-\}$, $\{s\}$ | meson, meson | $\frac{1}{2}$ | 2 | kept, D4 term excluded | $m_a$, $\max\left(m_P, m_a\right)$ | line 0: (i),(ii),(iii),(iv); line 1: (i),(ii),(iii),(iv) | depends on the routing |
| $C16$ | $\{\pi^+,s\}$, $\{\pi^-\}$ | $s$, $s$ | $\frac{1}{2}$ | 2 | vanishes identically | – | – | $V(\pi^+,s,s,s)=0$ |
| $C17$ | $\{\pi^+,s\}$, $\{\pi^-\}$ | ALP-like, $s$ | $1$ | 2 | vanishes identically | – | – | $V(\pi^+,s,s,v_a)=0$ |
| $C18$ | $\{\pi^+,s\}$, $\{\pi^-\}$ | ALP-like, ALP-like | $\frac{1}{2}$ | 2 | vanishes identically | – | – | $V(\pi^+,s,v_a,v_a)=0$ |
| $C19$ | $\{\pi^+,s\}$, $\{\pi^-\}$ | ALP-like, meson | $1$ | 2 | kept | $\max\left(m_P, m_a\right)$ | line 1: (i) | independent |
| $C20$ | $\{\pi^+,s\}$, $\{\pi^-\}$ | meson, $s$ | $1$ | 2 | kept | $\max\left(m_P, m_\sigma\right)$ | line 0: (i) | independent |
| $C21$ | $\{\pi^+,s\}$, $\{\pi^-\}$ | meson, meson | $\frac{1}{2}$ | 2 | kept | $\max\left(m_P, m_a\right)$ | line 0: (i),(ii); line 1: (i),(ii) | independent |
| $C22$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | $s$, $s$, $s$ | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,s,s)=0$ |
| $C23$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | $s$, $s$, ALP-like | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,v_a,s)=0$ |
| $C24$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | $s$, $s$, meson | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,\mathbb M^{-1}_{m_a^2}(\mathcal C),s)=0$, $V(\pi^+,p,s)=0$, $V(\pi^-,s,s)=0$ |
| $C25$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | $s$, ALP-like, $s$ | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,s,s)=0$ |
| $C26$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | $s$, ALP-like, ALP-like | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,v_a,s)=0$ |
| $C27$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | $s$, ALP-like, meson | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,\mathbb M^{-1}_{m_a^2}(\mathcal C),s)=0$, $V(\pi^+,p,s)=0$, $V(\pi^-,s,v_a)=0$ |
| $C28$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | $s$, meson, $s$ | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,s,s)=0$ |
| $C29$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | $s$, meson, ALP-like | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,v_a,s)=0$ |
| $C30$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | $s$, meson, meson | $1$ | 3 | kept | $\max\left(m_P, m_\sigma\right)$ | line 1: (i); line 2: (i) | independent |
| $C31$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | ALP-like, $s$, $s$ | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,s,v_a)=0$ |
| $C32$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | ALP-like, $s$, ALP-like | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,v_a,v_a)=0$ |
| $C33$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | ALP-like, $s$, meson | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,\mathbb M^{-1}_{m_a^2}(\mathcal C),v_a)=0$, $V(\pi^+,p,v_a)=0$, $V(\pi^-,v_a,s)=0$ |
| $C34$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | ALP-like, ALP-like, $s$ | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,s,v_a)=0$ |
| $C35$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | ALP-like, ALP-like, ALP-like | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,v_a,v_a)=0$ |
| $C36$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | ALP-like, ALP-like, meson | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,\mathbb M^{-1}_{m_a^2}(\mathcal C),v_a)=0$, $V(\pi^+,p,v_a)=0$, $V(\pi^-,v_a,v_a)=0$ |
| $C37$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | ALP-like, meson, $s$ | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,s,v_a)=0$ |
| $C38$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | ALP-like, meson, ALP-like | $1$ | 3 | vanishes identically | – | – | $V(\pi^+,v_a,v_a)=0$ |
| $C39$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | ALP-like, meson, meson | $1$ | 3 | kept | $\max\left(m_P, m_a\right)$ | line 1: (i); line 2: (i) | independent |
| $C40$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | meson, $s$, $s$ | $1$ | 3 | kept | $\max\left(m_P, m_\sigma\right)$ | line 0: (i) | independent |
| $C41$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | meson, $s$, ALP-like | $1$ | 3 | kept | $\max\left(m_P, m_a, m_\sigma\right)$ | line 0: (i) | independent |
| $C42$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | meson, $s$, meson | $1$ | 3 | kept | $\max\left(m_P, m_a, m_\sigma\right)$ | line 0: (i); line 2: (ii),(iii) | independent |
| $C43$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | meson, ALP-like, $s$ | $1$ | 3 | kept | $\max\left(m_P, m_a, m_\sigma\right)$ | line 0: (i) | independent |
| $C44$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | meson, ALP-like, ALP-like | $1$ | 3 | kept | $\max\left(m_P, m_a\right)$ | line 0: (i) | independent |
| $C45$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | meson, ALP-like, meson | $1$ | 3 | kept | $\max\left(m_P, m_a\right)$ | line 0: (i); line 2: (ii),(iii) | independent |
| $C46$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | meson, meson, $s$ | $1$ | 3 | kept | $\max\left(m_P, m_a, m_\sigma\right)$ | line 0: (i); line 1: (ii),(iii) | independent |
| $C47$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | meson, meson, ALP-like | $1$ | 3 | kept | $\max\left(m_P, m_a\right)$ | line 0: (i); line 1: (ii),(iii) | independent |
| $C48$ | $\{\pi^+\}$, $\{\pi^-\}$, $\{s\}$ | meson, meson, meson | $1$ | 3 | kept | $\max\left(m_P, m_a\right)$ | line 0: (i),(ii); line 1: (i),(ii),(iii); line 2: (i),(ii),(iii) | independent |

### Self-energies of $s$

The two legs are $s$ (in, momentum $p$) and $s$ (out, momentum $-p$).

| graph | vertices (external legs) | lines | $S_\Gamma$ | $r_\Gamma$ | status | $m_{\rm term}$ | pieces per meson line | routing |
|---|---|---|---|---|---|---|---|---|
| $C\sigma1$ | $\{\text{in},\text{out}\}$ | $s$ | $\frac{1}{2}$ | 1 | mass shift only ($\delta Z=0$) | – | – | independent |
| $C\sigma2$ | $\{\text{in},\text{out}\}$ | ALP-like | $\frac{1}{2}$ | 1 | mass shift only ($\delta Z=0$) | – | – | independent |
| $C\sigma3$ | $\{\text{in},\text{out}\}$ | meson | $\frac{1}{2}$ | 1 | mass shift only ($\delta Z=0$) | – | – | independent |
| $C\sigma4$ | $\{\text{in}\}$, $\{\text{out}\}$ | $s$, $s$ | $\frac{1}{2}$ | 2 | kept | $m_\sigma$ | – | independent |
| $C\sigma5$ | $\{\text{in}\}$, $\{\text{out}\}$ | ALP-like, $s$ | $1$ | 2 | kept | $\max\left(m_a, m_\sigma\right)$ | – | independent |
| $C\sigma6$ | $\{\text{in}\}$, $\{\text{out}\}$ | ALP-like, ALP-like | $\frac{1}{2}$ | 2 | kept | $m_a$ | – | independent |
| $C\sigma7$ | $\{\text{in}\}$, $\{\text{out}\}$ | ALP-like, meson | $1$ | 2 | kept | $m_a$, $\max\left(m_P, m_a\right)$ | line 1: (i),(ii),(iii),(iv) | independent |
| $C\sigma8$ | $\{\text{in}\}$, $\{\text{out}\}$ | meson, $s$ | $1$ | 2 | kept | $\max\left(m_P, m_\sigma\right)$, $\max\left(m_a, m_\sigma\right)$ | line 0: (i),(ii),(iii),(iv) | independent |
| $C\sigma9$ | $\{\text{in}\}$, $\{\text{out}\}$ | meson, meson | $\frac{1}{2}$ | 2 | kept | $m_P$, $m_a$, $\max\left(m_P, m_a\right)$ | line 0: (i),(ii),(iii),(iv); line 1: (i),(ii),(iii),(iv) | independent |

### Self-energies of $\pi^\pm$

The two legs are $\pi^+$ (in, momentum $p$) and $\pi^-$ (out, momentum $-p$); exchanging the roles gives the same $\delta Z$, graph by graph.

| graph | vertices (external legs) | lines | $S_\Gamma$ | $r_\Gamma$ | status | $m_{\rm term}$ | pieces per meson line | routing |
|---|---|---|---|---|---|---|---|---|
| $C\pi1$ | $\{\text{in},\text{out}\}$ | $s$ | $\frac{1}{2}$ | 1 | mass shift only ($\delta Z=0$) | – | – | independent |
| $C\pi2$ | $\{\text{in},\text{out}\}$ | ALP-like | $\frac{1}{2}$ | 1 | kept | $m_a$ | – | independent |
| $C\pi3$ | $\{\text{in},\text{out}\}$ | meson | $\frac{1}{2}$ | 1 | kept, D4 term excluded | $m_a$ | line 0: (iv) | independent |
| $C\pi4$ | $\{\text{in}\}$, $\{\text{out}\}$ | $s$, $s$ | $\frac{1}{2}$ | 2 | vanishes identically | – | – | $V(\pi^+,s,s)=0$ |
| $C\pi5$ | $\{\text{in}\}$, $\{\text{out}\}$ | ALP-like, $s$ | $1$ | 2 | vanishes identically | – | – | $V(\pi^+,s,v_a)=0$ |
| $C\pi6$ | $\{\text{in}\}$, $\{\text{out}\}$ | ALP-like, ALP-like | $\frac{1}{2}$ | 2 | vanishes identically | – | – | $V(\pi^+,v_a,v_a)=0$ |
| $C\pi7$ | $\{\text{in}\}$, $\{\text{out}\}$ | ALP-like, meson | $1$ | 2 | kept | $\max\left(m_P, m_a\right)$ | line 1: (i) | independent |
| $C\pi8$ | $\{\text{in}\}$, $\{\text{out}\}$ | meson, $s$ | $1$ | 2 | kept | $\max\left(m_P, m_\sigma\right)$ | line 0: (i) | independent |
| $C\pi9$ | $\{\text{in}\}$, $\{\text{out}\}$ | meson, meson | $\frac{1}{2}$ | 2 | kept | $\max\left(m_P, m_a\right)$ | line 0: (i),(ii); line 1: (i),(ii) | independent |

## Result: excluded D4 terms (recorded, not estimated)

Three terms satisfy the criterion of R6 in every routing and are removed exactly: the closed meson line on the five-leg vertex ($C3$), the bubble with $\{\pi^+,\pi^-\}$ at one vertex and $s$ alone at the other with two meson lines ($C15$), and the closed meson line on the $\pi^\pm$ leg ($C\pi3$). They are the chiral background of the model without the ALP: no ALP-sector letter, no $N$, no $\lambda$, no $m_\varphi^2$, and their ratio to $M_0=-\frac{B_0\epsilon}{2f_\pi}\langle\widehat S^C_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle$ (for $C\pi3$, the $\delta Z$ itself) depends on $\epsilon$ only through $p_\sigma^2=m_\sigma^2$. With the size factor included they read:

$C3$ (all routings):
$$
\begin{aligned}
D4={}&\left(\frac{B_0 \epsilon m_P^{2}}{96 \pi^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \lambda_{\pi^+}\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_0 \epsilon m_P^{2}}{64 \pi^{2} f_\pi^{3}}\right)\langle \sum_{\rm perm}\big(\widehat S^C_{1,0}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle
\end{aligned}
$$

$C15$, routing $\rho=1$:
$$
\begin{aligned}
D4\big|_{\rho=1}={}&\left(- \frac{B_0^{2} \epsilon}{192 \pi^{2} f_\pi^{3}}\right)\langle \widehat S^C_{0,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0 \epsilon \left(2 m_P^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{96 \pi^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \lambda_{\pi^+}\,\lambda_{\pi^-}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon}{192 \pi^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \sum_{\rm perm}\big(\widehat S^C_{0,0}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon}{96 \pi^{2} f_\pi^{3}}\right)\langle \widehat S^C_{0,0}\,\widehat S^C_{1,0}\rangle\langle \lambda_{\pi^+}\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_0 \epsilon \left(2 m_P^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{64 \pi^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\rangle \\
&+\left(\frac{B_0 \epsilon \left(2 m_P^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{64 \pi^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\,\lambda_{\pi^-}\,\lambda_{\pi^+}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon}{64 \pi^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\,\sum_{\rm perm}\big(\widehat S^C_{0,0}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle
\end{aligned}
$$

In the other routings of $C15$ the factor $2m_P^2$ in the three kinematic coefficients is replaced by $0$ ($\rho=2$), $3m_P^2$ ($\rho=3$) and $-m_P^2$ ($\rho=4$); the other coefficients do not change (script output).

$C\pi3$ (all routings):
$$
\begin{aligned}
D4={}&\left(\frac{m_P^{2}}{32 \pi^{2} f_\pi^{2}}\right)\langle \lambda_{\pi^+}\,\lambda_{\pi^-}\rangle
\end{aligned}
$$


In five meson-only graphs the all-generator-sum term is exactly zero ($C9$, $C21$, $C48$, $C\sigma3$, $C\pi9$), because a vertex with three meson ends is proportional to $\widehat P^C_{0,0}=0$, or, for the tadpole $C\sigma3$, because the term has no $p^2$ dependence. In $C\sigma9$ (two $s$ ends, $\epsilon^2$) the criterion fails and the all-generator-sum term, with scale $m_P$, is kept.

## Result: pieces of the meson line and the quartic $V_\varphi$ vertices

- No meson line contributes through its $p$ entry (ii) alone. In $C3$ only (ii), (iii), (iv) remain after the D4 removal, and in $C\pi3$ only (iv) remains: the $p$ entry and the mixed entries reach the vertex $(\pi^+,\pi^-,p,p)$ or $(\pi^+,\pi^-,p,\mathbb M^{-1}_{m_a^2}(\mathcal C))$, which has no $p^2$ dependence.
- With the common scale $m_a$, $\delta Z_{C\pi2}+\delta Z_{C\pi3}=0$ exactly (tested): the ALP-like line plus pieces (ii)+(iii)+(iv) is the bare $p$ line, and the vertex $(\pi^+,\pi^-,p,p)$ comes only from $\widehat S^C_{0,2}$, without momentum. No factor $m_a^2-m_\sigma^2$ enters $\delta Z_{\pi^\pm}$.
- **Quartic $V_\varphi$ vertices ($\lambda$).** A four-leg vertex can contain a quartic $V_\varphi$ vertex only if its external legs have no meson matrix; this happens only in the tadpoles $C\sigma1$, $C\sigma2$, $C\sigma3$ of the $s$ self-energy (script, all graphs and line ends). Their $\lambda$-dependent part is nonzero and has no $p^2$ dependence (checked with $\lambda$ kept), so $\delta Z_\sigma$ gets nothing from them: they only shift $m_\sigma$. The cubic $V_\varphi$ vertices $V(s,s,s)=-6\lambda f_a$ and $V(s;x,y)=-2\lambda f_a\xi_x\xi_y$ do enter kept graphs: with $\lambda$ kept in the vertices, $T$ depends on $\lambda$ in $C10$, $C12$, $C13$, $C15$, $C40$, $C44$, $C45$, $C47$, $C48$, $C\sigma4$, $C\sigma6$, $C\sigma7$, $C\sigma9$ (script). There $\lambda$ is eliminated through $m_\sigma^2$, which is why $m_\sigma^2$ appears in their coefficients also outside the kinematics.
- **Routing.** Only $C12$, $C13$, $C15$ depend on the routing: they are the graphs with the four-meson kinetic vertex at $\{\pi^+,\pi^-\}$ and two ALP-like or meson lines. $C\pi2$ and $C\pi3$ contain the kinetic vertex but do not depend on the routing (the terms linear in $k\cdot p$ cancel). All other graphs have no momentum-dependent vertex and are routing independent (tested).

## Verification

Exact symbolic tests in [`test_derive_sigma_pipi_cartesian_loop_graphs.py`](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_graphs.py): the counts per class and the kept labels; no box and no two-leg vertex; $S_\Gamma=1/2$ exactly for tadpoles and bubbles with two lines of the same type; every vanishing graph has a vanishing charged vertex without a generator end; momentum conservation of the external products; the four-piece split equals the residue $(\mathbf 1_9\oplus1)-v_av_a^{\rm T}$ in every graph and routing; ALP-like line plus (ii)+(iii)+(iv) equals the bare $p$ line; $\delta Z_{C\pi2}+\delta Z_{C\pi3}=0$; routing independence without the kinetic vertex; the D4 terms (criterion with $\lambda$ kept, exact removal, failure for $C\sigma9$ with its $\epsilon^2$); the quartic $V_\varphi$ statement; $\pi^+$ and $\pi^-$ self-energies equal; exact rational coefficients; U(3) completeness once with explicit generators and symbolic $3\times3$ matrices (never printed). No loop value, no component notation and no element of the unit-matrix basis appear (checked by the trace tests on this record).

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_graphs.py
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_graphs.py --markdown
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_graphs
```

[Script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_graphs.py) (about 10 s), [tests](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_graphs.py). The tables and the D4 blocks above are the `--markdown` output.
