---
title: Cartesian frame - compact traces of the kept 1PI one-loop sigma to pi+ pi- graphs with per-term scales and termwise routing maximum
status: derived
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Cartesian frame: compact traces of the kept 1PI graphs of $\sigma\to\pi^+\pi^-$

Task 3.2 (1PI graphs) of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama-cartesian`. **No loop integral, scalar loop function, loop coefficient or loop amplitude is computed.** Each scalar moment with $r$ denominators is replaced by $(m_{\rm term}^2)^{2-r}/(16\pi^2)$ and each product with the loop momentum by $m_{\rm term}^2$ (D1, D2). The values below are an order of magnitude per graph class and can miss cancellations between graphs with the same $\epsilon$ scaling; neither this estimate nor the polar one is authoritative. The traces are derived from the Cartesian vertices; no polar vertex, graph or trace is used.

## Inputs

- Kept graphs, $S_\Gamma$, $r_\Gamma$, line types, D4 terms: [Cartesian graph-list record](cartesian-sigma-pipi-loop-graph-list.md).
- Rules R1–R10 (four-piece split, D3 scales, size rule, completeness, D4, D11, routing, residue prescription): [Cartesian selection record](cartesian-sigma-pipi-loop-graph-selection.md).
- Cartesian vertices: [Cartesian vertex record](cartesian-sigma-pipi-loop-vertex-rules.md). Legs, $v_a$, $N(m_a^2)$, $\mathbb M^{-1}_{m_a^2}(\mathcal C)$, assumptions: [Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md).

## Method

For every kept graph and every routing $\rho$ the product of the Cartesian vertices is built with the line pieces of R2 and the loop momentum on one line (routing $2j+1$ along line $j$, $2j+2$ against it). Each term (one piece per meson line) gets its scale $m_{\rm term}$ (D3), with $k\cdot k,\ k\cdot p_j\to m_{\rm term}^2$. The generator sums are removed by completeness (D5) in two independent ways: construction 1 expands the vertices into words and applies completeness to words; construction 2 applies it to the compact vertex trees (sums over permutations, commutators, anticommutators) before any expansion, using $\langle X\sum_{\rm perm}(Z\,y_2\cdots y_m)\rangle=\langle Z\sum_{\rm perm}(y_2\cdots y_m\,X)\rangle$, $\langle[Z,y]R\rangle=\langle Z[y,R]\rangle$ and $\langle\{Z,y\}R\rangle=\langle Z\{y,R\}\rangle$ to bring the summed generator $Z$ to the front. The two agree as exact identities of cyclic words in every routing of every kept graph. Then the assumptions ($\widehat P^C_{0,0}=\widehat P^C_{1,0}=0$, $\langle\widehat S^C_{0,1}\rangle=\langle\widehat S^C_{1,1}\rangle=0$) and the charge rule are applied, and the D4 terms of the [graph-list record](cartesian-sigma-pipi-loop-graph-list.md) are removed.

The routing rule R9 is applied on the products of traces of words: for each product $t$ the routing of largest $|C_t^\rho|$ is chosen (termwise maximum, computed with the graph-agnostic helper `termwise_maximum` of `derive_sigma_pipi_loop_graphs.py`, which acts on a list of sums of products of traces and contains no graph, vertex or letter of the polar frame). The displayed compact form is the compact value of one routing, plus, where the choice depends on the masses, $\theta(\text{branch condition})$ times the compact difference between the routing of that branch and the displayed routing; the script checks that this reproduces the termwise maximum exactly, branch by branch.

## Result: routing

- 23 of the 26 kept 1PI graphs have no momentum-dependent vertex and are routing independent (tested); they enter with their single value.
- $C12$, $C13$, $C15$ (bubbles with $\{\pi^+,\pi^-\}$ at the four-meson kinetic vertex and $s$ alone at the other vertex) depend on the routing for 18 of 138, 34 of 238 and 34 of 238 products. All 86 are mass dependent with one pattern (tested): up to a factor that does not depend on the routing, the coefficients in the routings $\rho=1,2,3,4$ are
$$
2m_{\rm term}^2-(m_\sigma^2-m_{\pi^\pm}^2),\qquad -(m_\sigma^2-m_{\pi^\pm}^2),\qquad 3m_{\rm term}^2-(m_\sigma^2-m_{\pi^\pm}^2),\qquad -m_{\rm term}^2-(m_\sigma^2-m_{\pi^\pm}^2),
$$
with $m_{\rm term}$ the single loop scale of the product ($m_a$, or $\max(m_P,m_a)$ in $C13$ and $C15$). With $m_\sigma\ge2m_{\pi^\pm}$, $m_\sigma^2-m_{\pi^\pm}^2>0$, and the largest routing is
$$
\rho_t^*=3\quad\text{for } m_\sigma^2-m_{\pi^\pm}^2\le m_{\rm term}^2,\qquad \rho_t^*=4\quad\text{for } m_\sigma^2-m_{\pi^\pm}^2>m_{\rm term}^2 .
$$
- At the boundary $m_\sigma^2-m_{\pi^\pm}^2=m_{\rm term}^2$, routings 3 and 4 are equal and opposite and the sum over the four routings vanishes, so routing 3 is taken and the product is marked (86 marks, all at this boundary point; no tie on an open region). No product has a mass-independent ratio between routings.
- **Routing spread (diagnostic).** Over the routings and both branches, $C_t^\rho/C_t^{\rho_t^*}$ ranges from $-1$ (at the boundary) to $1$; in every routing-dependent product the coefficient of routing 1 (routing 3) vanishes at $m_\sigma^2-m_{\pi^\pm}^2=2m_{\rm term}^2$ ($3m_{\rm term}^2$). The routing ambiguity of the size rule is covered by the termwise maximum, which is an upper bound per graph on that ambiguity, not on the loop.

## Result: termwise routing maximum per graph

Notation: $\sum_{\rm perm}(y_1\cdots y_m)$ is the sum over the $m!$ orderings of the product (identical factors counted as distinct), $[\cdot,\cdot]$ and $\{\cdot,\cdot\}$ are the commutator and anticommutator, $\langle\cdot\rangle$ is the flavour trace, $\theta(\cdot)$ is 1 if the inequality holds and 0 otherwise, and $N(m_a^2)=1+2\epsilon^2\langle(\mathbb M^{-1}_{m_a^2}(\mathcal C))^2\rangle$. Everything is exact in $N(m_a^2)$, $\mathbb M^{-1}_{m_a^2}(\mathcal C)$, $\epsilon$ and the source traces; there is no expansion in $\epsilon$. $\lambda$ is eliminated through $m_\sigma^2=2\lambda f_a^2+\epsilon^2B_0(\langle\widehat S^C_{1,0}\rangle-2\langle\widehat S^C_{2,0}\rangle)$, so $m_\sigma^2$ also appears outside the kinematics. Terms of different loop scales are written on separate lines.

**$C1$.** Lines $s$; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=1$, $m_{\rm term}\in\{m_\sigma\}$; routing independent.
$$
\begin{aligned}
X_{C1}^{\max}={}&\left(- \frac{3 B_0 \epsilon^{3} m_\sigma^{2}}{32 \pi^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{3,0}\rangle
\end{aligned}
$$

**$C2$.** Lines ALP-like; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=1$, $m_{\rm term}\in\{m_a\}$; routing independent.
$$
\begin{aligned}
X_{C2}^{\max}={}&\left(- \frac{B_0 \epsilon^{3} m_a^{2}}{32 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{1,2}\rangle \\
&+\left(- \frac{B_0 \epsilon^{3} m_a^{2}}{48 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat P^C_{1,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0 \epsilon^{3} m_a^{2}}{192 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle
\end{aligned}
$$

**$C3$.** Lines meson; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=1$, $m_{\rm term}\in\{m_a\}$; D4 term excluded; routing independent.
$$
\begin{aligned}
X_{C3}^{\max}={}&\left(- \frac{B_0 \epsilon^{3} m_a^{2} \left(N(m_a^2) - 1\right)}{32 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{1,2}\rangle \\
&+\left(\frac{B_0 \epsilon^{3} m_a^{2}}{48 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat P^C_{1,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0 \epsilon^{3} m_a^{2}}{192 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle
\end{aligned}
$$

**$C7$.** Lines ALP-like, meson; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_{C7}^{\max}={}&\left(\frac{B_0^{2} \epsilon^{3}}{32 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{1,1},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

**$C8$.** Lines meson, $s$; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_{C8}^{\max}={}&\left(\frac{B_0^{2} \epsilon^{3}}{16 \pi^{2} f_\pi^{3}}\right)\langle \{\widehat S^C_{1,0},\lambda_{\pi^+}\}\,\{\widehat S^C_{2,0},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

**$C9$.** Lines meson, meson; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_{C9}^{\max}={}&\left(\frac{B_0^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)}{32 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{1,1},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

**$C10$.** Lines $s$, $s$; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{m_\sigma\}$; routing independent.
$$
\begin{aligned}
X_{C10}^{\max}={}&\left(- \frac{3 B_0^{2} \epsilon^{5}}{32 \pi^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{2,0}\rangle \\
&+\left(\frac{3 B_0^{2} \epsilon^{5}}{16 \pi^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{2,0}\rangle \\
&+\left(- \frac{3 B_0^{2} \epsilon^{5}}{16 \pi^{2} f_\pi^{3}}\right)\langle \widehat S^C_{3,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{2,0}\rangle \\
&+\left(\frac{3 B_0 \epsilon^{3} m_\sigma^{2}}{32 \pi^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{2,0}\rangle
\end{aligned}
$$

**$C11$.** Lines ALP-like, $s$; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_{C11}^{\max}={}&\left(- \frac{B_0^{2} \epsilon^{5}}{16 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat S^C_{2,1}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{1,1}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{8 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{2,0}\rangle
\end{aligned}
$$

**$C12$.** Lines ALP-like, ALP-like; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{m_a\}$; 18 of 138 products depend on the routing (mass-dependent choice, routing 3 below the boundary).
$$
\begin{aligned}
X_{C12}^{\max}={}&\left(- \frac{B_0^{2} \epsilon^{5}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle \\
&+\left(\frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{96 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{192 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle \\
&+\left(\frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{96 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{96 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0 \epsilon^{3} m_\sigma^{2}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{\epsilon^{3} m_\sigma^{2} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{96 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0 \epsilon^{3} m_\sigma^{2}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0 \epsilon^{3} m_\sigma^{2}}{192 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\theta\big(- \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} > m_a^{2}\big)\Big[\left(- \frac{B_0 \epsilon^{5} m_a^{2}}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} m_a^{2}}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_0 \epsilon^{5} m_a^{2}}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{\epsilon^{3} m_a^{2} m_\sigma^{2}}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_0 \epsilon^{5} m_a^{2}}{6 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} m_a^{2}}{6 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle\Big]
\end{aligned}
$$

**$C13$.** Lines ALP-like, meson; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{m_a,\max\left(m_P, m_a\right)\}$; 34 of 238 products depend on the routing (mass-dependent choice, routing 3 below the boundary).
$$
\begin{aligned}
X_{C13}^{\max}={}&\left(- \frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 2\right)}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{96 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 2\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle \\
&+\left(\frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 2\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0 \epsilon^{3} m_\sigma^{2} \left(N(m_a^2) - 1\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{4 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 2\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{\epsilon^{3} m_\sigma^{2} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{6 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 4\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0 \epsilon^{3} m_\sigma^{2} \left(N(m_a^2) - 2\right)}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0 \epsilon^{3} m_\sigma^{2}}{96 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0 \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,[[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}],\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_0 \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,[[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}],\lambda_{\pi^+}]\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{3}}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{0,0}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{3}}{48 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,\sum_{\rm perm}\big(\widehat P^C_{0,1}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0 \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle [[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}],\lambda_{\pi^-}]\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0 \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle [[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}],\lambda_{\pi^+}]\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{3}}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat P^C_{1,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{0,0}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{3}}{48 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat P^C_{1,1}\,\sum_{\rm perm}\big(\widehat P^C_{0,1}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\theta\big(- \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} > m_a^{2}\big)\Big[\left(\frac{B_0 \epsilon^{5} m_a^{2}}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_0 \epsilon^{5} m_a^{2}}{6 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} m_a^{2}}{6 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{\epsilon^{3} m_a^{2} m_\sigma^{2}}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} m_a^{2}}{3 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0 \epsilon^{5} m_a^{2}}{3 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle\Big] \\
&+\theta\big(- \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} > \max\left(m_P, m_a\right)^{2}\big)\Big[\left(- \frac{B_0 \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{24 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,[[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}],\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0 \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{24 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,[[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}],\lambda_{\pi^+}]\rangle \\
&+\left(\frac{B_0 \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{24 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle [[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}],\lambda_{\pi^-}]\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0 \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{24 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle [[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}],\lambda_{\pi^+}]\,\widehat P^C_{1,1}\rangle\Big]
\end{aligned}
$$

**$C14$.** Lines meson, $s$; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_{C14}^{\max}={}&\left(- \frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat S^C_{2,1}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{1,1}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{8 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{1,1}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{2,0}\rangle
\end{aligned}
$$

**$C15$.** Lines meson, meson; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{m_a,\max\left(m_P, m_a\right)\}$; D4 term excluded; 34 of 238 products depend on the routing (mass-dependent choice, routing 3 below the boundary).
$$
\begin{aligned}
X_{C15}^{\max}={}&\left(- \frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)^{2}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle \\
&+\left(\frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{96 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{192 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)^{2}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle \\
&+\left(\frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{96 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)^{2}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{96 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0 \epsilon^{3} m_\sigma^{2} \left(N(m_a^2) - 1\right)^{2}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 1\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \{\lambda_{\pi^+},\lambda_{\pi^-}\}\,\widehat S^C_{0,2}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{\epsilon^{3} m_\sigma^{2} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{96 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0 \epsilon^{5} \left(3 m_a^{2} + \left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2}\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{5}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5} \left(N(m_a^2) - 2\right)}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{5}}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0 \epsilon^{3} m_\sigma^{2} \left(N(m_a^2) - 1\right)}{48 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat P^C_{0,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0 \epsilon^{3} m_\sigma^{2}}{192 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{0,0}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0 \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,[[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}],\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0 \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,[[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}],\lambda_{\pi^+}]\rangle \\
&+\left(\frac{B_0^{2} \epsilon^{3}}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{0,0}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{3}}{48 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,\sum_{\rm perm}\big(\widehat P^C_{0,1}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(\frac{B_0 \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle [[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}],\lambda_{\pi^-}]\,\widehat P^C_{1,1}\rangle \\
&+\left(\frac{B_0 \epsilon^{3} \left(\left(m_{\pi^\pm}\right)^{2} - m_\sigma^{2} + 3 \max\left(m_P, m_a\right)^{2}\right)}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle [[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}],\lambda_{\pi^+}]\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{3}}{96 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat P^C_{1,1}\,\sum_{\rm perm}\big(\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{0,0}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)}{48 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \widehat P^C_{1,1}\,\sum_{\rm perm}\big(\widehat P^C_{0,1}\,\lambda_{\pi^+}\,\lambda_{\pi^-}\big)\rangle \\
&+\theta\big(- \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} > m_a^{2}\big)\Big[\left(- \frac{B_0 \epsilon^{5} m_a^{2}}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,0}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} m_a^{2}}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{1,2}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_0 \epsilon^{5} m_a^{2}}{12 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle \widehat S^C_{2,0}\rangle\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{\epsilon^{3} m_a^{2} m_\sigma^{2}}{24 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_0 \epsilon^{5} m_a^{2}}{6 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0 \epsilon^{5} m_a^{2}}{6 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3}}\right)\langle [\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}]\,[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}]\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle\Big] \\
&+\theta\big(- \left(m_{\pi^\pm}\right)^{2} + m_\sigma^{2} > \max\left(m_P, m_a\right)^{2}\big)\Big[\left(\frac{B_0 \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{24 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,[[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}],\lambda_{\pi^-}]\rangle \\
&+\left(\frac{B_0 \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{24 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\mathbb M^{-1}_{m_a^2}(\mathcal C),\widehat S^C_{1,0}\}\,[[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}],\lambda_{\pi^+}]\rangle \\
&+\left(- \frac{B_0 \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{24 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle [[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^+}],\lambda_{\pi^-}]\,\widehat P^C_{1,1}\rangle \\
&+\left(- \frac{B_0 \epsilon^{3} \max\left(m_P, m_a\right)^{2}}{24 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle [[\mathbb M^{-1}_{m_a^2}(\mathcal C),\lambda_{\pi^-}],\lambda_{\pi^+}]\,\widehat P^C_{1,1}\rangle\Big]
\end{aligned}
$$

**$C19$.** Lines ALP-like, meson; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_{C19}^{\max}={}&\left(\frac{B_0^{2} \epsilon^{3}}{32 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^-}\}\,\{\widehat S^C_{1,1},\lambda_{\pi^+}\}\rangle
\end{aligned}
$$

**$C20$.** Lines meson, $s$; $S_\Gamma=1$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_{C20}^{\max}={}&\left(\frac{B_0^{2} \epsilon^{3}}{16 \pi^{2} f_\pi^{3}}\right)\langle \{\widehat S^C_{1,0},\lambda_{\pi^-}\}\,\{\widehat S^C_{2,0},\lambda_{\pi^+}\}\rangle
\end{aligned}
$$

**$C21$.** Lines meson, meson; $S_\Gamma=\frac{1}{2}$, $r_\Gamma=2$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_{C21}^{\max}={}&\left(\frac{B_0^{2} \epsilon^{3} \left(N(m_a^2) - 1\right)}{32 \pi^{2} N(m_a^2) f_\pi^{3}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^-}\}\,\{\widehat S^C_{1,1},\lambda_{\pi^+}\}\rangle
\end{aligned}
$$

**$C30$.** Lines $s$, meson, meson; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_{C30}^{\max}={}&\left(- \frac{B_0^{3} \epsilon^{3}}{32 \pi^{2} f_\pi^{3} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \{\{\widehat S^C_{1,0},\lambda_{\pi^-}\},\widehat S^C_{1,0}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^+}\}\rangle
\end{aligned}
$$

**$C39$.** Lines ALP-like, meson, meson; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_{C39}^{\max}={}&\left(- \frac{B_0^{3} \epsilon^{3}}{32 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\{\widehat S^C_{0,1},\lambda_{\pi^-}\},\widehat S^C_{1,0}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^+}\}\rangle
\end{aligned}
$$

**$C40$.** Lines meson, $s$, $s$; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_{C40}^{\max}={}&\left(\frac{3 B_0^{3} \epsilon^{5}}{32 \pi^{2} f_\pi^{3} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \{\widehat S^C_{1,0},\lambda_{\pi^+}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{3 B_0^{3} \epsilon^{5}}{16 \pi^{2} f_\pi^{3} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \{\widehat S^C_{1,0},\lambda_{\pi^+}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{3 B_0^{3} \epsilon^{5}}{16 \pi^{2} f_\pi^{3} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \widehat S^C_{3,0}\rangle\langle \{\widehat S^C_{1,0},\lambda_{\pi^+}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{3 B_0^{2} \epsilon^{3} m_\sigma^{2}}{32 \pi^{2} f_\pi^{3} \max\left(m_P, m_\sigma\right)^{2}}\right)\langle \{\widehat S^C_{1,0},\lambda_{\pi^+}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^-}\}\rangle
\end{aligned}
$$

**$C41$.** Lines meson, $s$, ALP-like; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_{C41}^{\max}={}&\left(\frac{B_0^{3} \epsilon^{5}}{16 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \widehat S^C_{2,1}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_0^{3} \epsilon^{5}}{8 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{2,0}\rangle
\end{aligned}
$$

**$C42$.** Lines meson, $s$, meson; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_{C42}^{\max}={}&\left(\frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \widehat S^C_{2,1}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_0^{3} \epsilon^{5}}{8 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{2,0}\rangle
\end{aligned}
$$

**$C43$.** Lines meson, ALP-like, $s$; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_{C43}^{\max}={}&\left(\frac{B_0^{3} \epsilon^{5}}{16 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \widehat S^C_{2,1}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^-}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^+}\}\rangle \\
&+\left(\frac{B_0^{3} \epsilon^{5}}{8 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^-}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^+}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{2,0}\rangle
\end{aligned}
$$

**$C44$.** Lines meson, ALP-like, ALP-like; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_{C44}^{\max}={}&\left(\frac{B_0^{3} \epsilon^{5}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_0^{3} \epsilon^{5}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_0^{3} \epsilon^{5}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{3} m_\sigma^{2}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_0^{3} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0^{3} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle
\end{aligned}
$$

**$C45$.** Lines meson, ALP-like, meson; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_{C45}^{\max}={}&\left(\frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{3} m_\sigma^{2} \left(N(m_a^2) - 1\right)}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_0^{3} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 2\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle
\end{aligned}
$$

**$C46$.** Lines meson, meson, $s$; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a, m_\sigma\right)\}$; routing independent.
$$
\begin{aligned}
X_{C46}^{\max}={}&\left(\frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \widehat S^C_{2,1}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^-}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^+}\}\rangle \\
&+\left(- \frac{B_0^{3} \epsilon^{5}}{8 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a, m_\sigma\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^-}\}\,\{\widehat S^C_{1,0},\lambda_{\pi^+}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{2,0}\rangle
\end{aligned}
$$

**$C47$.** Lines meson, meson, ALP-like; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_{C47}^{\max}={}&\left(\frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{3} m_\sigma^{2} \left(N(m_a^2) - 1\right)}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_0^{3} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(\frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 2\right)}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle
\end{aligned}
$$

**$C48$.** Lines meson, meson, meson; $S_\Gamma=1$, $r_\Gamma=3$, $m_{\rm term}\in\{\max\left(m_P, m_a\right)\}$; routing independent.
$$
\begin{aligned}
X_{C48}^{\max}={}&\left(\frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)^{2}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S^C_{1,0}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(\frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)^{2}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S^C_{1,2}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)^{2}}{16 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \widehat S^C_{2,0}\rangle\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_0^{3} \epsilon^{3} \left(N(m_a^2) - 1\right)}{32 \pi^{2} N(m_a^2) f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\{\widehat S^C_{0,1},\lambda_{\pi^-}\},\widehat S^C_{1,0}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^+}\}\rangle \\
&+\left(- \frac{B_0^{2} \epsilon^{3} m_\sigma^{2} \left(N(m_a^2) - 1\right)^{2}}{32 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle \\
&+\left(- \frac{B_0^{3} \epsilon^{5}}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat S^C_{1,0}\rangle \\
&+\left(- \frac{B_0^{3} \epsilon^{5} \left(N(m_a^2) - 1\right)}{8 \pi^{2} \left(N(m_a^2)\right)^{2} f_\pi^{3} \max\left(m_P, m_a\right)^{2}}\right)\langle \{\widehat S^C_{0,1},\lambda_{\pi^+}\}\,\{\widehat S^C_{0,1},\lambda_{\pi^-}\}\rangle\langle \mathbb M^{-1}_{m_a^2}(\mathcal C)\,\widehat P^C_{1,1}\rangle
\end{aligned}
$$

## Verification

Exact symbolic tests in [`test_derive_sigma_pipi_cartesian_loop_traces.py`](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_traces.py) and [`test_derive_sigma_pipi_cartesian_loop_graphs.py`](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_graphs.py):
- construction 1 equals construction 2 in every routing of every kept graph;
- graphs without the kinetic vertex are routing independent, and their maximum is their single value;
- the routing pattern above, and $|C_t^{\rho_t^*}|\ge|C_t^\rho|$ for every routing and product: the routing coefficients are first checked to have the four forms above exactly, then the inequality is proved on each branch with an exact univariate inequality;
- the displayed form equals the termwise maximum branch by branch;
- dimensional rescaling: with $B_0$, $f_\pi$, $m_\sigma$, $m_a$, $m_P$, $m_{\pi^\pm}$ and each source letter scaled by $c$ and $\epsilon$, $N$, $\lambda_{\pi^\pm}$, $\mathbb M^{-1}_{m_a^2}(\mathcal C)$ unscaled, every $X^{\max}$ scales as $c^1$ (with $\max(cx,cy)=c\max(x,y)$), and every branch condition is homogeneous of degree 2;
- the four-piece split, the ALP-like line plus (ii)+(iii)+(iv) against the bare $p$ line, and the D4 removal (graph tests);
- this record contains no loop value, no component notation and no element of the unit-matrix basis (test on the file).

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_traces.py
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_traces.py --markdown
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_traces
```

[Script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_traces.py) (about 20 s), [tests](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_traces.py). The blocks above are the `--markdown` output. The function `estimate_inputs()` returns $X^{\max}_\Gamma$ per graph, $\sum\delta Z^{\max}$ for $s$ and for $\pi^\pm$, and $M_0$, for the assembly of the estimate (task 4.1); `entries()` returns the values in every routing.
