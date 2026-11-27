---
title: Cartesian frame - compact traces and residue corrections of the external legs in the one-loop sigma to pi+ pi- estimate
status: derived
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Cartesian frame: compact traces and residue corrections of the external legs

Task 3.2 (external legs) of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama-cartesian`. **No loop integral, scalar loop function, loop coefficient or loop amplitude is computed.** Each scalar moment with $r$ denominators is replaced by $(m_{\rm term}^2)^{2-r}/(16\pi^2)$ and each product with the loop momentum by $m_{\rm term}^2$ (D1, D2). The values below are an order of magnitude per graph class and can miss cancellations between graphs with the same $\epsilon$ scaling; neither this estimate nor the polar one is authoritative. The traces are derived from the Cartesian vertices; no polar vertex, graph or trace is used.

## Inputs

- Kept graphs, $S_\Gamma$, $r_\Gamma$, line types, D4 terms: [Cartesian graph-list record](cartesian-sigma-pipi-loop-graph-list.md).
- Rules R1–R10 (four-piece split, D3 scales, size rule, completeness, D4, D11, routing, residue prescription): [Cartesian selection record](cartesian-sigma-pipi-loop-graph-selection.md).
- Cartesian vertices: [Cartesian vertex record](cartesian-sigma-pipi-loop-vertex-rules.md). Legs, $v_a$, $N(m_a^2)$, $\mathbb M^{-1}_{m_a^2}(\mathcal C)$, assumptions: [Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md).

## Residue prescription and assumptions

Rule R10 of the [Cartesian selection record](cartesian-sigma-pipi-loop-graph-selection.md) (the polar residue prescription restated): the two legs of a self-energy carry incoming momenta $p$ and $-p$, $p^2$ is kept symbolic, and for each term
$$
\delta Z_{\rm term}=-\frac{S_\Gamma}{16\pi^2}\Big[\theta(r_\Gamma\ge2)\,(m_{\rm term}^2)^{1-r_\Gamma}\,T_{\rm term}+(m_{\rm term}^2)^{2-r_\Gamma}\,\frac{\partial T_{\rm term}}{\partial p^2}\Big]_{p^2=m_{\rm ext}^2},
$$
with $m_{\rm ext}=m_\sigma$ or $m_{\pi^\pm}$; $k\cdot p\to m_{\rm term}^2$ does not contribute to $\partial/\partial p^2$. The external legs add $M_0\big(\tfrac12\delta Z_\sigma^{\max}+\tfrac12\delta Z_{\pi^+}^{\max}+\tfrac12\delta Z_{\pi^-}^{\max}\big)$ to the estimate, each $\delta Z^{\max}$ being the sum over the kept self-energies of that leg of their termwise routing maxima. Self-energies with $\delta Z=0$ exactly only shift the masses; the physical masses are used (R8). One-loop tadpoles: those of $\pi^0$, $\eta$, $\eta'$ and $p$ vanish under CP; the one-loop $s$ tadpole is set to zero (assumption 9 of the [frame record](cartesian-sigma-pipi-loop-frame-and-tree.md)). The one-particle-reducible graphs in which an external leg mixes with another state are removed by the CP and flavour-diagonality assumptions (R7) and are not included.

## Method

For every kept graph and every routing $\rho$ the product of the Cartesian vertices is built with the line pieces of R2 and the loop momentum on one line (routing $2j+1$ along line $j$, $2j+2$ against it). Each term (one piece per meson line) gets its scale $m_{\rm term}$ (D3), with $k\cdot k,\ k\cdot p_j\to m_{\rm term}^2$. The generator sums are removed by completeness (D5) in two independent ways: construction 1 expands the vertices into words and applies completeness to words; construction 2 applies it to the compact vertex trees (sums over permutations, commutators, anticommutators) before any expansion, using $\langle X\sum_{\rm perm}(Z\,y_2\cdots y_m)\rangle=\langle Z\sum_{\rm perm}(y_2\cdots y_m\,X)\rangle$, $\langle[Z,y]R\rangle=\langle Z[y,R]\rangle$ and $\langle\{Z,y\}R\rangle=\langle Z\{y,R\}\rangle$ to bring the summed generator $Z$ to the front. The two agree as exact identities of cyclic words in every routing of every kept graph. Then the assumptions ($\widehat P^C_{0,0}=\widehat P^C_{1,0}=0$, $\langle\widehat S^C_{0,1}\rangle=\langle\widehat S^C_{1,1}\rangle=0$) and the charge rule are applied, and the D4 terms of the [graph-list record](cartesian-sigma-pipi-loop-graph-list.md) are removed.

The routing rule R9 is applied on the products of traces of words: for each product $t$ the routing of largest $|C_t^\rho|$ is chosen (termwise maximum, computed with the graph-agnostic helper `termwise_maximum` of `derive_sigma_pipi_loop_graphs.py`, which acts on a list of sums of products of traces and contains no graph, vertex or letter of the polar frame). The displayed compact form is the compact value of one routing, plus, where the choice depends on the masses, $\theta(\text{branch condition})$ times the compact difference between the routing of that branch and the displayed routing; the script checks that this reproduces the termwise maximum exactly, branch by branch.

## Result: structure

- **Routing.** No self-energy depends on the routing (tested). In the Cartesian frame the only momentum-dependent vertex is the four-meson kinetic vertex; it appears only in the closed loops $C\pi2$, $C\pi3$ on the $\pi^\pm$ leg, where the terms linear in $k\cdot p$ cancel between the two orientations.
- **$s$ leg.** The tadpoles $C\sigma1$–$C\sigma3$ have no $p^2$ dependence and give $\delta Z=0$; they contain the quartic $V_\varphi$ vertices, which therefore only shift $m_\sigma$ and do not enter $\delta Z_\sigma$. $\delta Z_\sigma^{\max}$ is the sum of the bubbles $C\sigma4$–$C\sigma9$, all through the derivative of one propagator ($r_\Gamma=2$, no $p^2$ in the numerators). The all-generator-sum term of $C\sigma9$ has two $s$ ends ($\epsilon^2$) and is kept (D4 does not apply).
- **$\pi^\pm$ legs.** $\delta Z_{\pi^+}=\delta Z_{\pi^-}$ graph by graph (tested). $C\pi1$ gives $\delta Z=0$. $\delta Z_{C\pi2}+\delta Z_{C\pi3}=0$ exactly: with the common scale $m_a$ the ALP-like line plus pieces (ii)+(iii)+(iv) of the meson line is the bare $p$ line, whose vertex $(\pi^+,\pi^-,p,p)$ has no momentum dependence; the D4 term of $C\pi3$ is removed. Hence $\delta Z_{\pi^\pm}^{\max}$ comes from the bubbles $C\pi7$, $C\pi8$, $C\pi9$. No factor $m_a^2-m_\sigma^2$ appears in any $\delta Z_{\pi^\pm}$ (tested, graph by graph, in every routing and in the maximum).

## Result: residue corrections per graph

Notation: $\sum_{\rm perm}(y_1\cdots y_m)$ is the sum over the $m!$ orderings of the product (identical factors counted as distinct), $[\cdot,\cdot]$ and $\{\cdot,\cdot\}$ are the commutator and anticommutator, $\langle\cdot\rangle$ is the flavour trace, $\theta(\cdot)$ is 1 if the inequality holds and 0 otherwise, and $N(m_a^2)=1+2\epsilon^2\langle(\mathbb M^{-1}_{m_a^2}(\mathcal C))^2\rangle$. Everything is exact in $N(m_a^2)$, $\mathbb M^{-1}_{m_a^2}(\mathcal C)$, $\epsilon$ and the source traces; there is no expansion in $\epsilon$. $\lambda$ is eliminated through $m_\sigma^2=2\lambda f_a^2+\epsilon^2B_0(\langle\widehat S^C_{1,0}\rangle-2\langle\widehat S^C_{2,0}\rangle)$, so $m_\sigma^2$ also appears outside the kinematics. Terms of different loop scales are written on separate lines.

### Self-energies of $s$

**$C\sigma4$.** Lines $s$, $s$; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{m_\sigma\}$; routing independent.
$$
\begin{aligned}
\delta Z_{C\sigma4}^{\max}={}&\left(- \frac{9 \epsilon^{2} m_\sigma^{2}}{32 \pi^{2} f_\pi^{2}}\right) \\
&+\left(\frac{9 B_0 \epsilon^{4}}{16 \pi^{2} f_\pi^{2}}\right)\langle \widehat S^C_{1,0}\rangle \\
&+\left(- \frac{9 B_0^{2} \epsilon^{6}}{32 \pi^{2} f_\pi^{2} m_\sigma^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{1,0}\rangle \\
&+\left(\frac{9 B_0^{2} \epsilon^{6}}{8 \pi^{2} f_\pi^{2} m_\sigma^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{2,0}\rangle \\
&+\left(- \frac{9 B_0^{2} \epsilon^{6}}{8 \pi^{2} f_\pi^{2} m_\sigma^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{3,0}\rangle \\
&+\left(- \frac{9 B_0 \epsilon^{4}}{8 \pi^{2} f_\pi^{2}}\right)\langle \widehat S^C_{2,0}\rangle \\
&+\left(- \frac{9 B_0^{2} \epsilon^{6}}{8 \pi^{2} f_\pi^{2} m_\sigma^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \widehat S^C_{2,0}\rangle \\
&+\left(\frac{9 B_0^{2} \epsilon^{6}}{4 \pi^{2} f_\pi^{2} m_\sigma^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \widehat S^C_{3,0}\rangle \\
&+\left(\frac{9 B_0 \epsilon^{4}}{8 \pi^{2} f_\pi^{2}}\right)\langle \widehat S^C_{3,0}\rangle \\
&+\left(- \frac{9 B_0^{2} \epsilon^{6}}{8 \pi^{2} f_\pi^{2} m_\sigma^{2}}\right)\langle \widehat S^C_{3,0}\rangle\langle \widehat S^C_{3,0}\rangle
\end{aligned}
$$

**$C\sigma5$.** Lines ALP-like, $s$; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
\delta Z_{C\sigma5}^{\max}={}&\left(- \frac{B_0^{2} \epsilon^{6}}{4 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_a, m_\sigma\right)^{2}}\right)\langle \widehat S^C_{2,1}\rangle\langle \widehat S^C_{2,1}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{\pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_a, m_\sigma\right)^{2}}\right)\langle \widehat S^C_{2,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{2,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{\pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_a, m_\sigma\right)^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{2,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{2,0}\rangle
\end{aligned}
$$

**$C\sigma6$.** Lines ALP-like, ALP-like; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{m_a\}$; routing independent.
$$
\begin{aligned}
\delta Z_{C\sigma6}^{\max}={}&\left(- \frac{\epsilon^{2} m_\sigma^{4}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right) \\
&+\left(\frac{B_0 \epsilon^{4} m_\sigma^{2}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{1,2}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{2,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0 \epsilon^{4} m_\sigma^{2}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \widehat S^C_{1,2}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \widehat S^C_{2,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0 \epsilon^{4} m_\sigma^{2}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{2,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \widehat S^C_{2,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0 \epsilon^{4} m_\sigma^{2}}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0 \epsilon^{4} m_\sigma^{2}}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle
\end{aligned}
$$

**$C\sigma7$.** Lines ALP-like, meson; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{m_a,\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
\delta Z_{C\sigma7}^{\max}={}&\left(- \frac{\epsilon^{2} m_\sigma^{4} \left(N(m_a^2) - 1\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right) \\
&+\left(\frac{B_0 \epsilon^{4} m_\sigma^{2} \left(N(m_a^2) - 1\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{1,2}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{2,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 2\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0 \epsilon^{4} m_\sigma^{2} \left(N(m_a^2) - 1\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \widehat S^C_{1,2}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \widehat S^C_{2,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 2\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0 \epsilon^{4} m_\sigma^{2} \left(N(m_a^2) - 1\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{2,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \widehat S^C_{2,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 2\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0 \epsilon^{4} m_\sigma^{2}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{2 B_0^{2} \epsilon^{6}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0 \epsilon^{4} m_\sigma^{2} \left(N(m_a^2) - 2\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 4\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{4}}{8 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,\{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{4}}{4 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{4}}{8 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat P^C_{1,1}\,\widehat P^C_{1,1}\rangle
\end{aligned}
$$

**$C\sigma8$.** Lines meson, $s$; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_\sigma\right),\max\left(m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
\delta Z_{C\sigma8}^{\max}={}&\left(- \frac{B_0^{2} \epsilon^{4}}{2 \pi^{2} f_\pi^{2} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \widehat P^C_{2,0}\,\widehat P^C_{2,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{4 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_a, m_\sigma\right)^{2}}\right)\langle \widehat S^C_{2,1}\rangle\langle \widehat S^C_{2,1}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{\pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_a, m_\sigma\right)^{2}}\right)\langle \widehat S^C_{2,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{2,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{\pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_a, m_\sigma\right)^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{2,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{2,0}\rangle
\end{aligned}
$$

**$C\sigma9$.** Lines meson, meson; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{m_P,m_a,\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
\delta Z_{C\sigma9}^{\max}={}&\left(- \frac{B_0^{2} \epsilon^{2}}{16 \pi^{2} f_\pi^{2} m_P^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{1,0}\rangle \\
&+\left(- \frac{3 B_0^{2} \epsilon^{2}}{16 \pi^{2} f_\pi^{2} m_P^{2}}\right)\langle \widehat S^C_{1,0}\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{\epsilon^{2} m_\sigma^{4} \left(N(m_a^2) - 1\right)^{2}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right) \\
&+\left(\frac{B_0 \epsilon^{4} m_\sigma^{2} \left(N(m_a^2) - 1\right)^{2}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)^{2}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)^{2}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{1,2}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)^{2}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{2,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0 \epsilon^{4} m_\sigma^{2} \left(N(m_a^2) - 1\right)^{2}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)^{2}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \widehat S^C_{1,2}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)^{2}}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \widehat S^C_{2,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0 \epsilon^{4} m_\sigma^{2} \left(N(m_a^2) - 1\right)^{2}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{2,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)^{2}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \widehat S^C_{2,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 1\right)}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0 \epsilon^{4} m_\sigma^{2}}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{6}}{2 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6}}{\pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0 \epsilon^{4} m_\sigma^{2} \left(N(m_a^2) - 1\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{6} \left(N(m_a^2) - 2\right)}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{2} m_a^{2}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{4}}{8 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,\{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{4}}{4 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{4} \left(N(m_a^2) - 1\right)}{8 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat P^C_{1,1}\,\widehat P^C_{1,1}\rangle
\end{aligned}
$$

### Self-energies of $\pi^\pm$ ($\pi^+$ shown; $\pi^-$ identical)

**$C\pi2$.** Lines ALP-like; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=1$, $m_{\rm term}\in\{m_a\}$; routing independent.
$$
\begin{aligned}
\delta Z_{C\pi2}^{\max}={}&\left(- \frac{\epsilon^{2} m_a^{2}}{96 \pi^{2} N(m_a^2) f_\pi^{2}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle
\end{aligned}
$$

**$C\pi3$.** Lines meson; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=1$, $m_{\rm term}\in\{m_a\}$; D4 term excluded; routing independent.
$$
\begin{aligned}
\delta Z_{C\pi3}^{\max}={}&\left(\frac{\epsilon^{2} m_a^{2}}{96 \pi^{2} N(m_a^2) f_\pi^{2}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle
\end{aligned}
$$

**$C\pi7$.** Lines ALP-like, meson; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
\delta Z_{C\pi7}^{\max}={}&\left(- \frac{B_0^{2} \epsilon^{2}}{32 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

**$C\pi8$.** Lines meson, $s$; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
\delta Z_{C\pi8}^{\max}={}&\left(- \frac{B_0^{2} \epsilon^{2}}{32 \pi^{2} f_\pi^{2} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \{\widehat S^C_{1,0},\lambda_{\pi^+}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

**$C\pi9$.** Lines meson, meson; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
\delta Z_{C\pi9}^{\max}={}&\left(- \frac{B_0^{2} \epsilon^{2} \left(N(m_a^2) - 1\right)}{32 \pi^{2} N(m_a^2) f_\pi^{2} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

## Verification

Exact symbolic tests in [`test_derive_sigma_pipi_cartesian_loop_traces.py`](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_traces.py) and [`test_derive_sigma_pipi_cartesian_loop_graphs.py`](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_graphs.py): construction 1 equals construction 2 in every routing; routing independence of every self-energy; $\delta Z_{\pi^+}=\delta Z_{\pi^-}$ per graph; $\delta Z_{C\pi2}+\delta Z_{C\pi3}=0$; no factor $m_a^2-m_\sigma^2$ in $\delta Z_{\pi^\pm}$; the quartic $V_\varphi$ vertices only in $C\sigma1$–$C\sigma3$, with no $p^2$ dependence; the D4 criterion and exact removal for $C\pi3$, and its failure for $C\sigma9$; dimensional rescaling ($\delta Z^{\max}$ scales as $c^0$); the sum of the $\pi^\pm$ maxima returned by `estimate_inputs()`; this record contains no loop value, no component notation and no element of the unit-matrix basis (test on the file).

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_traces.py --markdown
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_traces test_derive_sigma_pipi_cartesian_loop_graphs
```

[Trace script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_traces.py), [graph script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_graphs.py). The blocks above are the `--markdown` output.
