---
title: Large epsilon at fixed sources - assumptions, pole data, meson-line residue pieces and leg factors of the sigma to pi+ pi- loop estimate
status: mixed
created: 2026-11-26
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Large $\epsilon$: assumptions, pole data and leg factors

OpenSpec change `estimate-sigma-pipi-one-loop-murayama` (tasks 2.3 and 4.5, decisions D3, D10, D12). The large-$\epsilon$ behaviour, $\epsilon=f_\pi/f_a\to\infty$, is taken at fixed sources and fixed quartic coupling, as in the up-only model. This record gives the exact leading power of $\epsilon$ and the leading coefficient of the inputs of the estimate. Vertices and tree: [large-$\epsilon$ vertex record](murayama-sigma-pipi-loop-large-epsilon-vertices-tree.md). Graphs and $R$: [large-$\epsilon$ graph record](murayama-sigma-pipi-loop-large-epsilon-graphs.md). No loop integral is evaluated, and no number is used.

## Assumptions

1. **Held fixed** (independent of $\epsilon$): $f_\pi$, $B_0$, $m_1^2$, every source letter $\widehat S_{r,s}$, $\widehat P_{r,s}$, the matrices $\lambda_{\pi^\pm}$, and the quartic coupling $\lambda$ of $V_\varphi$. Then $f_a=f_\pi/\epsilon$ varies. The sources are generic: every independent product of traces is taken to be nonzero, so the powers below are upper bounds, valid unless a particular trace combination vanishes. $m_{\pi^\pm}$ is not an independent input: it comes from $\widehat S_{0,0}$ alone and does not depend on $\epsilon$. It is kept as a fixed symbol in the kinematics, and no relation between $m_{\pi^\pm}$ and $\widehat S_{0,0}$ is used. All leading coefficients below turn out to be free of $m_{\pi^\pm}$.
2. **ALP mass entry:** $\epsilon^2A$ with $A=-2B_0\langle\widehat S_{0,2}\rangle>0$ fixed ([large-$\epsilon$ record](murayama-sigma-pipi-loop-large-epsilon-scaling.md), Step 4). $c_k=\langle\mathcal C\,L^k(\mathcal C)\rangle$, $\mathcal C=B_0\widehat P_{0,1}$, $L(X)=B_0\{\widehat S_{0,0},X\}+\frac{m_1^2}{3}\langle X\rangle\mathbf 1$.
3. **Radial mass with $\lambda$ fixed:** from the [vertex-rule record](murayama-sigma-pipi-loop-vertex-rules.md),
$$
m_\sigma^2=\frac{2\lambda f_\pi^2}{\epsilon^2}+\epsilon^2\kappa_\sigma,\qquad \kappa_\sigma=B_0\big(\langle\widehat S_{1,0}\rangle-2\langle\widehat S_{2,0}\rangle\big)>0\ \text{(assumption)} .
$$
In $V(\sigma,\sigma,\sigma)$ and $V(\sigma^4)$ the coupling $\lambda$ enters, so $m_\sigma^2-\epsilon^2\kappa_\sigma=2\lambda f_\pi^2/\epsilon^2$ exactly. The script implements this by the exact rewriting $\langle\widehat S_{1,0}\rangle=\kappa_\sigma/B_0+2\langle\widehat S_{2,0}\rangle$ and $\langle\widehat S_{0,2}\rangle=-A/(2B_0)$ in every trace product before any series is taken. $\lambda$ is fixed and positive (bounded potential); the term $2\lambda f_a^2=2\lambda f_\pi^2/\epsilon^2$ falls, and the growth $m_\sigma^2\to\epsilon^2\kappa_\sigma$ comes from $\kappa_\sigma>0$.
4. **Kinematics:** $p_1\cdot p_2=(m_\sigma^2-2m_{\pi^\pm}^2)/2$ and $p_\sigma\cdot p_i=-m_\sigma^2/2$ ([tree and legs record](murayama-sigma-pipi-loop-tree-and-legs.md)), with $m_\sigma^2$ as above and $m_{\pi^\pm}^2$ kept exactly.
5. **Loop scale per term** as in the [selection record](murayama-sigma-pipi-loop-graph-selection.md) (D3): $m_{\rm term}$ is the largest of the scales of the lines of the term, $m_\sigma$ (radial), $m_a$ (ALP-like line, and pieces (ii)+(iii)+(iv) of a meson line) and $m_P$ (bare generator sum). $m_P$, the largest pole mass of the nine meson-like states, is held fixed ($\epsilon^0$, symbolic); since $m_a,m_\sigma\sim\epsilon$, $m_P$ drops out of every maximum that contains $m_a$ or $m_\sigma$ at large $\epsilon$, and only the terms with every line in its bare generator sum keep the scale $m_P$. $\max(m_\sigma,m_a)$ is $m_a$ for $A>\kappa_\sigma$ and $m_\sigma$ for $A<\kappa_\sigma$; for $A=\kappa_\sigma$, $m_a^2-m_\sigma^2\to2c_0/A>0$, which is the first branch. Loop-momentum products are set to $m_{\rm term}^2$ term by term.
6. **Routing (D12) at large $\epsilon$.** For each product of traces $t$ of the canonical generic-source basis (before the rewriting of item 3), the coefficient $C_t^\rho$ of each routing is expanded at large $\epsilon$. The routing with the largest leading power of $\epsilon$ is chosen; among those, the one with the largest magnitude of the leading coefficient (exact ratios, with $\kappa_\sigma=rA$); for equal leading terms, the next coefficients times the sign of the leading one decide; ties of equal magnitude and opposite sign take the sign of the sum of the coefficients over all routings. The choice is constant in regions of $r=\kappa_\sigma/A$ bounded by the roots of $L_i=\pm L_j$, of $L_i=0$ (leading coefficients) and by $r=1$; the script checks it at three points per region. The regions are $r<1$, $1<r<2$, $2<r<5$ and $r>5$; the boundary values (degenerate mass relations) are not resolved. No tie and no undecided comparison occurs.
7. The large-$s$ expansion $\mathbb M_s^{-1}=-\sum_k s^{-(k+1)}L^k$ holds, i.e. $\epsilon^2A$ lies above the spectrum of $L$ ([large-$\epsilon$ record](murayama-sigma-pipi-loop-large-epsilon-scaling.md), Step 4).
8. Leading-order chiral perturbation theory is used at $m_a^2,m_\sigma^2\sim\epsilon^2$; this is a limit of the estimate ([limits record](murayama-sigma-pipi-loop-limits.md)). Since $R\sim\epsilon^4$ ([graph record](murayama-sigma-pipi-loop-large-epsilon-graphs.md)), the condition $R\ll1$ fails before the condition of leading-order chiral perturbation theory.

## Method

Every quantity is a series in $1/\epsilon^2$: $m_a^2$ is the root of $D_a$, $\mathbb M^{-1}_{m_a^2}(\mathcal C)=-\sum_k (m_a^2)^{-(k+1)}L^k(\mathcal C)$ with the series for $m_a^2$ inserted, and $N(m_a^2)=1+2\epsilon^2\langle(\mathbb M^{-1}_{m_a^2}(\mathcal C))^2\rangle$. Each series is truncated by total order in $1/\epsilon^2$, not by $k$. The pieces of the meson-line residue are combined exactly before any power is read off: the ALP-field entry $1$ and the ALP-field component $-1/N$ of $-v_av_a^{\rm T}$ are each of order one and combine to $1-1/N=O(\epsilon^{-2})$. Per independent trace product, the leading coefficient is taken. If it cancels, the next order is used; a letter $\mathbb M^{-1}_{m_a^2}(\mathcal C)$ and a direct letter $\widehat P_{0,1}$ become the same trace at leading order, and the script handles this. If a cancellation persisted to the truncation order (relative $\epsilon^{-4}$), the script would stop with an error; this never happened.

## Result: pole data and inputs

| quantity | leading power | leading coefficient |
|---|---|---|
| $m_a^2$ | $\epsilon^2$ | $A$ |
| $m_\sigma^2$ | $\epsilon^2$ | $\kappa_\sigma$ |
| $m_a^2-m_\sigma^2$, $A\neq\kappa_\sigma$ | $\epsilon^2$ | $A-\kappa_\sigma$ |
| $m_a^2-m_\sigma^2$, $A=\kappa_\sigma$ | $\epsilon^0$ | $2B_0^2\langle\widehat P_{0,1}\widehat P_{0,1}\rangle/A=2c_0/A$ |
| $m_a^2-m_\sigma^2$, up-only sources | $\epsilon^0$ | $2B_0m_u+\big(\tfrac{2m_1^2}{3}-2\lambda f_\pi^2\big)\epsilon^{-2}+\dots$ |
| $\mathbb M^{-1}_{m_a^2}(\mathcal C)$ | $\epsilon^{-2}$ | $-(B_0/A)\,\widehat P_{0,1}$ |
| $N(m_a^2)$, $1/N(m_a^2)$ | $\epsilon^0$ | $1$ |
| $N(m_a^2)-1$ | $\epsilon^{-2}$ | $2B_0^2\langle\widehat P_{0,1}\widehat P_{0,1}\rangle/A^2$ |
| $\xi_a=1/\sqrt{N(m_a^2)}$ (ALP-field component of $v_a$) | $\epsilon^0$ | $1$ |
| meson matrix of $v_a$, $2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)/\sqrt{N}$ | $\epsilon^{-1}$ | $-(2B_0/A)\,\widehat P_{0,1}$ |

In the up-only model ($A=\kappa_\sigma=B_0m_u$) the radial and ALP-like poles are degenerate at order $\epsilon^2$; the loop scale $\max(m_\sigma,m_a)$ is then $m_a$, decided at order $\epsilon^0$, and $r=1$ lies at the boundary of the regions of item 6 (the up-only coefficient of $R$ is the same in all regions).

## Result: pieces of the meson-line residue

The residue $(\mathbf 1_9\oplus1)-v_av_a^{\rm T}$ is split exactly into four pieces ([meson-line piece record](murayama-sigma-pipi-loop-meson-line-pieces.md)):

| piece | ALP-field ends | leading power | leading coefficient |
|---|---|---|---|
| (i) bare generator sum $\mathbf 1_9$ | 0 | $\epsilon^0$ | $1$ |
| (ii) ALP-field entry $1-1/N(m_a^2)$ | 2 | $\epsilon^{-2}$ | $2B_0^2\langle\widehat P_{0,1}\widehat P_{0,1}\rangle/A^2$ |
| (iii) mixed entries $-2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)/N$ | 1 | $\epsilon^{-1}$ | $(2B_0/A)\,\widehat P_{0,1}$ |
| (iv) meson block of the pole, $-4\epsilon^2\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\otimes\mathbb M^{-1}_{m_a^2}(\mathcal C)/N$ | 0 | $\epsilon^{-2}$ | $-(4B_0^2/A^2)\,\widehat P_{0,1}\otimes\widehat P_{0,1}$ |

## Result: leg factors with their placement

Each ALP-field end at a source vertex carries one power of $\epsilon$; at a kinetic $\sigma aa$ or $\sigma\sigma aa$ vertex it carries none, and a kinetic vertex gets in addition $k\cdot k\to m_{\rm term}^2\sim\epsilon^2$ (or $m_P^2\sim\epsilon^0$ in a term with only bare generator sums) ([large-$\epsilon$ vertex record](murayama-sigma-pipi-loop-large-epsilon-vertices-tree.md)).

| leg factor | power | net power by placement |
|---|---|---|
| $\xi_a$ | $\epsilon^0$ | $\epsilon^1$ at a source vertex; $\epsilon^0$ at a kinetic vertex |
| meson matrix of $v_a$ | $\epsilon^{-1}$ | $\epsilon^{-1}$ |
| (i) $\mathbf 1_9$ | $\epsilon^0$ | $\epsilon^0$ |
| (ii) $1-1/N$ | $\epsilon^{-2}$ | $\epsilon^0$ (both ends at source vertices); $\epsilon^{-1}$ (one at a kinetic vertex); $\epsilon^{-2}$ (both at kinetic vertices) |
| (iii) mixed | $\epsilon^{-1}$ | $\epsilon^0$ (ALP-field end at a source vertex); $\epsilon^{-1}$ (at a kinetic vertex) |
| (iv) meson block | $\epsilon^{-2}$ | $\epsilon^{-2}$ |

The ALP-like end is mostly the ALP field: $\xi_a\to1$, and its meson part is suppressed by $\epsilon^{-1}$.

## Verification

Exact symbolic checks in the script: the defining relation $\mathbb M_s(\mathbb M_s^{-1}(\mathcal C))=\mathcal C$ order by order, the pole condition $D_a(m_a^2)=0$, agreement of the series root with its closed form to order $\epsilon^{-2}$, the exact combination $1-1/N$, and the leading powers and coefficients of the tables.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_large_epsilon.py --markdown
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_large_epsilon
```

[Script](../../../scripts/derive_sigma_pipi_loop_large_epsilon.py) (sections 5 to 11), [tests](../../../scripts/test_derive_sigma_pipi_loop_large_epsilon.py).
