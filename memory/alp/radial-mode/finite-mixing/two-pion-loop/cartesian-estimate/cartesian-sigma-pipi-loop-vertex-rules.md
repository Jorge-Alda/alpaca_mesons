---
title: Cartesian vertex rules up to five legs for generic sources in the one-loop sigma to pi+ pi- estimate
status: derived
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Cartesian vertex rules up to five legs for generic $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$

Task 2.1 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama-cartesian`. The Cartesian estimate is an independent order-of-magnitude estimate; neither it nor the polar estimate is authoritative. No loop integral, scalar loop function, loop coefficient or loop amplitude is evaluated. The vertices are derived from the Cartesian Lagrangian only; no polar vertex is used as input.

## Conventions and inputs

- Lagrangian, legs, vertex convention, units and assumptions: [Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md). The vertex convention is the one of the [polar vertex record](../one-loop-estimate/murayama-sigma-pipi-loop-vertex-rules.md) ($n$-th derivative of $\mathcal L$, $\partial_\mu\to-ik_\mu$, all momenta incoming); only this convention is shared.
- $m_\sigma^2=2\lambda f_a^2+\epsilon^2B_0(\langle\widehat S^C_{1,0}\rangle-2\langle\widehat S^C_{2,0}\rangle)$ from the Cartesian quadratic Lagrangian ([frame record](cartesian-sigma-pipi-loop-frame-and-tree.md), Step 2).
- Assumptions used: $f_\pi$ independent of $s$ and $p$; anomalous singlet term independent of $s$ and $p$; no derivative coupling of $s$, $p$ to the mesons (assumptions 3, 4, 8 of the frame record).

## Method

**Construction 1.** Direct noncommutative expansion of $\mathcal L$ with formal legs: $s=\sum\phi_j$ over the radial legs, $p=\sum_j\xi_j\phi_j$, $\Phi=\sum_j\phi_j\lambda_j$, $\partial_\mu\phi_j\to-ik_{j,\mu}\phi_j$. The coefficient of $\phi_1\cdots\phi_n$ is extracted as exact cyclic words. The source part is taken in the exponential form $\frac{f_\pi^2}{4}2B_0\langle\widehat MU^\dagger+U\widehat M^\dagger\rangle$ with $\widehat M=\sum(\widehat S^C_{r,s}+i\widehat P^C_{r,s})(s/f_a)^r(p/f_a)^s$. $V_\varphi$ is expanded as a polynomial in $s$ and $p$. The kinetic terms of $s$, $p$ and the mesons are expanded as well.

**Construction 2.** The closed rules below, written both as words and as compact trees.

## Result

For legs $x_1,\dots,x_n$, let $R$ be the radial legs, $A$ the legs that contribute their $p$ component $\xi_x$, and $P$ the legs that contribute their meson matrix $\lambda_x$, with $r=|R|$, $s=|A|$, $m=|P|$. Let $X^C_{r,s}=\widehat S^C_{r,s}$ for even $m$ and $\widehat P^C_{r,s}$ for odd $m$. A vertex is the sum of the parts that apply to its legs.

**Source part.**
$$
V^{(n)}_{\rm src}=B_0f_\pi^{2-n}\sum_{(R,A,P)}r!\,s!\,\epsilon^{r+s}\frac{(-1)^{m(m-1)/2}}{m!}\prod_{x\in A}\xi_x\sum_{\pi\in S_m}\big\langle X^C_{r,s}\lambda_{\pi_1}\cdots\lambda_{\pi_m}\big\rangle .
$$

**Radial potential $V_\varphi$.** Above the quadratic order, $-V_\varphi\supset-\lambda f_a(s^3+sp^2)-\frac{\lambda}{4}(s^2+p^2)^2$. $V_\varphi$ does not contain the mesons, so every pseudoscalar leg enters through its $p$ component:
$$
V(s,s,s)=-6\lambda f_a,\quad V(s;x,y)=-2\lambda f_a\,\xi_x\xi_y,\quad V(s,s,s,s)=-6\lambda,\quad V(s,s;x,y)=-2\lambda\,\xi_x\xi_y,\quad V(w,x,y,z)=-6\lambda\,\xi_w\xi_x\xi_y\xi_z .
$$
There is no quintic term, and the $V_\varphi$ vertices do not depend on the momenta. $\lambda$ is eliminated through $m_\sigma^2$ (decision D-A3), in the cubic and in the quartic vertices:
$$
\lambda f_a=\frac{\epsilon}{2f_\pi}\Big[m_\sigma^2-\epsilon^2B_0\big(\langle\widehat S^C_{1,0}\rangle-2\langle\widehat S^C_{2,0}\rangle\big)\Big],\qquad \lambda=\frac{\epsilon}{f_\pi}\,\lambda f_a .
$$
The traces $\langle\widehat S^C_{1,0}\rangle$, $\langle\widehat S^C_{2,0}\rangle$ then enter the trace basis.

**Four-meson kinetic vertex** (from $\frac{f_\pi^2}{4}\langle u_\mu u^\mu\rangle$; four legs with meson matrices):
$$
V^{(4)}_{\rm kin}=\frac{1}{48f_\pi^2}\sum_{\pi\in S_4}\big(-k_{\pi_2}\cdot k_{\pi_4}\big)\big\langle[\lambda_{\pi_1},\lambda_{\pi_2}][\lambda_{\pi_3},\lambda_{\pi_4}]\big\rangle .
$$
There is no cubic or quintic meson kinetic vertex. The kinetic term $\frac12(\partial s)^2+\frac12(\partial p)^2$ gives no vertex, and $f_\pi$ does not depend on $s$ or $p$, so $s$ and $p$ have no derivative coupling. The only momentum-dependent vertex is $V^{(4)}_{\rm kin}$.

The six-meson kinetic vertex is not needed: no one-loop graph of $\sigma\to\pi^+\pi^-$ has a vertex with more than five legs. For the generic source letters, $\widehat P^C_{0,0}=\widehat P^C_{1,0}=0$ are applied as matrix identities, and $\langle\widehat S^C_{0,1}\rangle=\langle\widehat S^C_{1,1}\rangle=0$ as trace identities.

**Up-only limit (D-A4).** With $\widehat S^C_{0,0}=\widehat S^C_{1,0}=m_uP_u$ and $\widehat P^C_{0,1}=-m_uP_u$, only $(r,s)\in\{(0,0),(1,0),(0,1)\}$ survive. So a source vertex has at most one end that is radial or $p$. For example,
$$
V^{\text{up-only}}(s;x,y)=-\frac{B_0m_u\epsilon}{2f_\pi}\big\langle P_u\{\lambda_x,\lambda_y\}\big\rangle-2\lambda f_a\,\xi_x\xi_y .
$$

## Verification

All checks are exact symbolic identities of cyclic words. No random or numerical matrix is used.
- **Constructions 1 and 2 agree**, both as words and as trees, for every pattern $n=3,4,5$ with $r=0,\dots,n$ radial legs and generic pseudoscalar legs, each with both $\lambda_x$ and $\xi_x$ (15 patterns, $\lambda$ kept).
- **Relabelling symmetry**: every vertex is invariant under the exchange of two legs of the same kind.
- **Dimensions**: under $X\to c^{d_X}X$, with $d=1$ for $B_0$, $f_\pi$, $m_\sigma$ and each source letter, $d=2$ for each $k_i\cdot k_j$, and $d=0$ for $\epsilon$, $\lambda_x$, $\xi_x$, every vertex with $\lambda$ eliminated scales as $c^{4-n}$.
- **Up-only limit**: for every pattern $n=3,4,5$, the closed rules after the substitution of D-A4 equal the direct expansion with $\widehat M=m_uP_u\varphi^*/f_a$.
- **Momentum dependence**: only the pattern of four pseudoscalar legs depends on the momenta, through $V^{(4)}_{\rm kin}$. The $V_\varphi$ vertices, including the quartic ones, have none, and the expansion of $\frac12(\partial s)^2+\frac12(\partial p)^2$ has no term with three or more legs.
- The $\lambda$ elimination reproduces $\lambda f_a$ when $m_\sigma^2$ is inserted.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_vertices.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_vertices
```

[Script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_vertices.py), [tests](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_vertices.py), [Cartesian algebra](../../../../../../scripts/sigma_pipi_cartesian_loop_algebra.py). The run takes about 50 s.
