---
title: Cartesian large epsilon at fixed sources - regime, pole data, meson-line residue pieces and leg factors
status: mixed
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Cartesian large $\epsilon$: regime, pole data, residue pieces and leg factors

Task 2.3 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama-cartesian`, decision D10 restated in the Cartesian invariants. The Cartesian estimate is an independent order-of-magnitude estimate; neither it nor the polar estimate is authoritative. No loop integral, scalar loop function, loop coefficient or loop amplitude is evaluated, and no number is used. The vertex powers and the tree are in the [Cartesian large-$\epsilon$ vertex record](cartesian-sigma-pipi-loop-large-epsilon-vertices-tree.md).

## Inputs

All exact, derived from the Cartesian quadratic Lagrangian in the [Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md), Steps 2 and 3:
$$
D_a(s)=s-\epsilon^2A+2\epsilon^2\big\langle\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\big\rangle,\qquad N(s)=1+2\epsilon^2\big\langle(\mathbb M_s^{-1}(\mathcal C))^2\big\rangle,\qquad m_\sigma^2=\frac{2\lambda f_\pi^2}{\epsilon^2}+\epsilon^2\kappa_\sigma,
$$
$$
v_a=\Big(\frac{2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)}{\sqrt N},\frac1{\sqrt N}\Big),\qquad \mathcal C=B_0\widehat P^C_{0,1},\qquad A=B_0\big(\langle\widehat S^C_{1,0}\rangle-2\langle\widehat S^C_{0,2}\rangle\big),\qquad \kappa_\sigma=B_0\big(\langle\widehat S^C_{1,0}\rangle-2\langle\widehat S^C_{2,0}\rangle\big).
$$
Here $s$ is the squared momentum of the line (not the radial field). $\mathbb M_s=L-s$ with $L(X)=B_0\{\widehat S^C_{0,0},X\}+\frac{m_1^2}{3}\langle X\rangle\mathbf 1$, and $c_k=\langle\mathcal C\,L^k(\mathcal C)\rangle$. The residue of an internal meson line is split as in the polar estimate (rules of the [polar selection record](../one-loop-estimate/murayama-sigma-pipi-loop-graph-selection.md), rule 3): $(\mathbf 1_9\oplus1)-v_av_a^{\rm T}$, written as four pieces.

## Regime (D10, assumptions)

1. **Held fixed** as $\epsilon=f_\pi/f_a$ grows: $f_\pi$, $B_0$, $m_1^2$, $m_P$, every source letter $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$, the matrices $\lambda_{\pi^\pm}$, and the quartic coupling $\lambda$ of $V_\varphi$. $f_a=f_\pi/\epsilon$ varies.
2. $A>0$ and $\kappa_\sigma>0$.
3. $\epsilon^2A$ lies above the spectrum of $L$, so $\mathbb M_s^{-1}=-\sum_kL^k/s^{k+1}$ at $s=m_a^2$.
4. Generic sources: every independent product of traces is nonzero, so the powers are upper bounds unless a trace combination vanishes.
5. Leading-order chiral perturbation theory is used at $m_a^2,m_\sigma^2\sim\epsilon^2$; this is a limit of the estimate, as in the polar estimate.

## Method

Each quantity is the exact expression above. A series in $u=1/\epsilon$ is used only to read the leading power and its coefficient. The sums over $c_k$ are kept up to $k=3$; the omitted terms are of higher order, because $c_k$ enters with $\epsilon^{-2k}$ relative to $c_0$. The root of $D_a$ is found order by order. All pieces of a quantity are combined before the leading term is read. For example, $1$ and $-1/N$ combine to $1-1/N=O(\epsilon^{-2})$. Up-only values: decision D-A4, $\widehat S^C_{0,0}=\widehat S^C_{1,0}=m_uP_u$, $\widehat P^C_{0,1}=-m_uP_u$, reduced with $P_u^2=P_u$, $\langle P_u\rangle=1$ and $\langle\mathbf 1\rangle=3$.

## Result: pole data (generic sources)

$$
m_a^2=\epsilon^2A+\frac{2c_0}{A}+\frac{1}{\epsilon^2}\Big(\frac{2c_1}{A^2}-\frac{4c_0^2}{A^3}\Big)+O(\epsilon^{-4}),\qquad c_0=B_0^2\langle\widehat P^C_{0,1}\widehat P^C_{0,1}\rangle .
$$

| quantity | leading power | leading coefficient |
|---|---|---|
| $m_a^2$ | $\epsilon^2$ | $A$ |
| $m_\sigma^2$ | $\epsilon^2$ | $\kappa_\sigma$ |
| $m_a^2-m_\sigma^2$, $A\neq\kappa_\sigma$ | $\epsilon^2$ | $A-\kappa_\sigma$ |
| $m_a^2-m_\sigma^2$, $A=\kappa_\sigma$ | $\epsilon^0$ | $2B_0^2\langle\widehat P^C_{0,1}\widehat P^C_{0,1}\rangle/A$ |
| $\mathbb M^{-1}_{m_a^2}(\mathcal C)$ | $\epsilon^{-2}$ | $-(B_0/A)\,\widehat P^C_{0,1}$ |
| $N(m_a^2)$, $1/N(m_a^2)$ | $\epsilon^0$ | $1$ |
| $N(m_a^2)-1$ | $\epsilon^{-2}$ | $2B_0^2\langle\widehat P^C_{0,1}\widehat P^C_{0,1}\rangle/A^2$ |
| $\xi_a=1/\sqrt{N(m_a^2)}$ ($p$ component of $v_a$) | $\epsilon^0$ | $1$ |
| meson matrix of $v_a$, $2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)/\sqrt N$ | $\epsilon^{-1}$ | $-(2B_0/A)\,\widehat P^C_{0,1}$ |

## Result: up-only pole data

$A=\kappa_\sigma=B_0m_u$, $c_0=B_0^2m_u^2$, $c_1=B_0^2m_u^2(2B_0m_u+m_1^2/3)$, and
$$
m_a^2-m_\sigma^2=2B_0m_u+\frac{1}{\epsilon^2}\Big(\frac{2m_1^2}{3}-2\lambda f_\pi^2\Big)+O(\epsilon^{-4}) .
$$
The first nonvanishing order is $\epsilon^0$, with coefficient $2B_0m_u>0$: the ALP-like pole lies above the radial pole for $\epsilon$ large enough. The poles are degenerate at order $\epsilon^2$. In the up-only limit, $\mathbb M^{-1}_{m_a^2}(\mathcal C)\to\epsilon^{-2}P_u$, $N-1\to2\epsilon^{-2}$, and the meson matrix of $v_a$ goes to $2\epsilon^{-1}P_u$.

## Result: pieces of the meson-line residue

| piece | $p$ ends | leading power | leading coefficient (generic) | up-only |
|---|---|---|---|---|
| (i) bare generator sum $\mathbf 1_9$ | 0 | $\epsilon^0$ | $1$ | $1$ |
| (ii) $p$ entry $1-1/N(m_a^2)$ | 2 | $\epsilon^{-2}$ | $2B_0^2\langle\widehat P^C_{0,1}\widehat P^C_{0,1}\rangle/A^2$ | $2$ |
| (iii) mixed entries $-2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)/N$ | 1 | $\epsilon^{-1}$ | $(2B_0/A)\,\widehat P^C_{0,1}$ | $-2P_u$ |
| (iv) meson block of the pole, $-4\epsilon^2\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\otimes\mathbb M^{-1}_{m_a^2}(\mathcal C)/N$ | 0 | $\epsilon^{-2}$ | $-(4B_0^2/A^2)\,\widehat P^C_{0,1}\otimes\widehat P^C_{0,1}$ | $-4P_u\otimes P_u$ |

## Result: leg factors with their placement

In the Cartesian frame a $p$ end sits either at a source vertex, which carries one power of $\epsilon$ per $p$ end (vertex rule $\epsilon^{r+s}$), or at a $V_\varphi$ vertex. The power of a $V_\varphi$ vertex ($\epsilon^{-1}$ cubic, $\epsilon^0$ quartic, with $\lambda$ fixed) does not depend on the number of its $p$ ends; it is listed in the [vertex record](cartesian-sigma-pipi-loop-large-epsilon-vertices-tree.md). There is no kinetic vertex with a $p$ end.

| leg factor | power | net power by placement of the $p$ ends |
|---|---|---|
| $\xi_a$ | $\epsilon^0$ | $\epsilon^1$ at a source vertex; $\epsilon^0$ at a $V_\varphi$ vertex |
| meson matrix of $v_a$ | $\epsilon^{-1}$ | $\epsilon^{-1}$ |
| (i) $\mathbf 1_9$ | $\epsilon^0$ | $\epsilon^0$ |
| (ii) $1-1/N$ | $\epsilon^{-2}$ | $\epsilon^0$ (both ends at source vertices); $\epsilon^{-1}$ (one at a $V_\varphi$ vertex); $\epsilon^{-2}$ (both at $V_\varphi$ vertices) |
| (iii) mixed | $\epsilon^{-1}$ | $\epsilon^0$ ($p$ end at a source vertex); $\epsilon^{-1}$ (at a $V_\varphi$ vertex) |
| (iv) meson block | $\epsilon^{-2}$ | $\epsilon^{-2}$ |

## Verification

Exact symbolic checks in the script and the tests:
- the series identity $(L-s)\big(-\frac1s\sum_{k\le K}(L/s)^k\big)=1-(L/s)^{K+1}$, and the self-adjointness of $L$, which gives $\langle L^j(\mathcal C)L^k(\mathcal C)\rangle=c_{j+k}$;
- $N$ from $\langle(\mathbb M^{-1})^2\rangle$ equals $D_a'(m_a^2)$ as series;
- the coefficients of the root;
- the up-only series root is a zero of the exact up-only $D_a(s)$ through relative order $\epsilon^{-6}$. The exact $D_a$ is obtained on the span of $P_u$ and $\mathbf 1$, which $L$ maps into itself. The test also checks that the leading term alone fails;
- the tables above and the sign $2B_0m_u>0$.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_large_epsilon_vertices.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_large_epsilon_vertices
```

[Script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_large_epsilon_vertices.py) (sections 1–3 and 5), [tests](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_large_epsilon_vertices.py).
