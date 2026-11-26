---
title: Large epsilon at fixed sources - leading powers of every kept graph, of M1_est and of R in sigma to pi+ pi-
status: derived
created: 2026-11-26
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Large $\epsilon$: graphs, $M_1^{\rm est}$ and $R$

OpenSpec change `estimate-sigma-pipi-one-loop-murayama` (tasks 2.3 and 4.5, decisions D3, D10, D12). The powers are labels of the exact estimate; no graph and no term is dropped. No loop integral is evaluated.

## Inputs

- Assumptions of the large-$\epsilon$ behaviour (fixed sources and $\lambda$, $m_a^2\simeq\epsilon^2A$, $m_\sigma^2\simeq\epsilon^2\kappa_\sigma$, $m_P$ fixed, $N\to1$, $\xi_a\to1$), the per-term scales at large $\epsilon$, the large-$\epsilon$ routing rule and its regions of $\kappa_\sigma/A$, pole data, residue pieces and leg factors: [large-$\epsilon$ input record](murayama-sigma-pipi-loop-large-epsilon-inputs.md).
- Vertices and $M_0\sim\epsilon$: [large-$\epsilon$ vertex record](murayama-sigma-pipi-loop-large-epsilon-vertices-tree.md).
- Graphs, $S_\Gamma$, $r_\Gamma$, scales of the terms: [graph-list record](murayama-sigma-pipi-loop-graph-list.md). Per-routing values, $X_\Gamma^{\max}$, $\delta Z^{\max}$: [1PI trace record](murayama-sigma-pipi-loop-traces-1pi.md), [external-leg record](murayama-sigma-pipi-loop-traces-external-legs.md). $M_1^{\rm est}$ and $R$: [estimate record](murayama-sigma-pipi-loop-estimate.md).

## Method

For every kept graph and every routing the exact size-rule value with the per-term scales (1PI: $\sum_{\rm terms}S_\Gamma(m_{\rm term}^2)^{2-r_\Gamma}T_{\rm term}/(16\pi^2)$; self-energies: $\delta Z$ term by term) is taken from the graph scripts. For each product of traces $t$ of the canonical generic-source basis the routing is chosen by the large-$\epsilon$ rule (largest leading power of $\epsilon$ of $C_t^\rho$, then largest magnitude of the leading coefficient; [input record](murayama-sigma-pipi-loop-large-epsilon-inputs.md), item 6). The choice is constant in four regions of $r=\kappa_\sigma/A$: region 1 $\kappa_\sigma<A$, region 2 $A<\kappa_\sigma<2A$, region 3 $2A<\kappa_\sigma<5A$, region 4 $\kappa_\sigma>5A$. Ten graphs depend on the routing (13, 14, 15, 20, 21, 25, 26, σ4, σ5, σ6); no tie occurs. The resulting $X_\Gamma^{\max}$, $\delta Z^{\max}$ (and, with the same routing choices, the unsized $T_\Gamma$) are rewritten with the large-$\epsilon$ series of all inputs and the leading power is read off per independent trace product with the cancellation handling of the input record. Two source variants: generic fixed sources (primary) and $r\le1$ sources ($\widehat S_{r,s}=\widehat P_{r,s}=0$ for $r\ge2$).

## Result: 1PI graphs

Entries: leading power of $T_\Gamma$ / of $T_\Gamma/M_0$ / of the contribution $X_\Gamma^{\max}/M_0$ to $R$. "0" means exactly zero. Lines as in the [graph-list record](murayama-sigma-pipi-loop-graph-list.md); for triangles in the order $\pi^+$–$\pi^-$; $\pi^-$–$\sigma$; $\sigma$–$\pi^+$. The powers are the same in the four regions.

| no. | vertices | lines | generic | $r\le1$ |
|---|---|---|---|---|
| 1 | $(\pi^+,\pi^-,\sigma)$ | radial | $\epsilon^3/\epsilon^2/\epsilon^4$ | 0 |
| 2 | $(\pi^+,\pi^-,\sigma)$ | ALP-like | $\epsilon^3/\epsilon^2/\epsilon^4$ | same |
| 3 | $(\pi^+,\pi^-,\sigma)$ | meson (D4 term removed) | $\epsilon^1/\epsilon^0/\epsilon^2$ | same |
| 4, 7 | one pion with $\sigma$ | meson, radial | $\epsilon^3/\epsilon^2/\epsilon^2$ | 0 |
| 5, 8 | one pion with $\sigma$ | ALP-like, meson | $\epsilon^3/\epsilon^2/\epsilon^2$ | same |
| 6, 9 | one pion with $\sigma$ | meson, meson | $\epsilon^1/\epsilon^0/\epsilon^0$ | same |
| 10 | $(\pi^+,\pi^-)$, $(\sigma)$ | radial, radial | $\epsilon^5/\epsilon^4/\epsilon^4$ | 0 |
| 11 | $(\pi^+,\pi^-)$, $(\sigma)$ | ALP-like, radial | $\epsilon^5/\epsilon^4/\epsilon^4$ | 0 |
| 12 | $(\pi^+,\pi^-)$, $(\sigma)$ | meson, radial | $\epsilon^3/\epsilon^2/\epsilon^2$ | 0 |
| 13 | $(\pi^+,\pi^-)$, $(\sigma)$ | ALP-like, ALP-like | $\epsilon^5/\epsilon^4/\epsilon^4$ | same |
| 14 | $(\pi^+,\pi^-)$, $(\sigma)$ | ALP-like, meson | $\epsilon^3/\epsilon^2/\epsilon^2$ | same |
| 15 | $(\pi^+,\pi^-)$, $(\sigma)$ | meson, meson (D4 term removed) | $\epsilon^3/\epsilon^2/\epsilon^2$ | same |
| 16 | triangle | meson; radial; radial | $\epsilon^5/\epsilon^4/\epsilon^2$ | $\epsilon^1/\epsilon^0/\epsilon^{-2}$ |
| 17 | triangle | meson; ALP-like; radial | $\epsilon^5/\epsilon^4/\epsilon^2$ | 0 |
| 18 | triangle | meson; meson; radial | $\epsilon^3/\epsilon^2/\epsilon^0$ | 0 |
| 19 | triangle | meson; radial; ALP-like | $\epsilon^5/\epsilon^4/\epsilon^2$ | 0 |
| 20 | triangle | meson; ALP-like; ALP-like | $\epsilon^5/\epsilon^4/\epsilon^2$ | same |
| 21 | triangle | meson; meson; ALP-like | $\epsilon^3/\epsilon^2/\epsilon^0$ | same |
| 22 | triangle | radial; meson; meson | $\epsilon^3/\epsilon^2/\epsilon^0$ | same |
| 23 | triangle | ALP-like; meson; meson | $\epsilon^3/\epsilon^2/\epsilon^0$ | same |
| 24 | triangle | meson; radial; meson | $\epsilon^3/\epsilon^2/\epsilon^0$ | 0 |
| 25 | triangle | meson; ALP-like; meson | $\epsilon^3/\epsilon^2/\epsilon^0$ | same |
| 26 | triangle | meson; meson; meson | $\epsilon^1/\epsilon^0/\epsilon^{-2}$ | same |

Which pieces of each meson line contribute is in the [meson-line piece record](murayama-sigma-pipi-loop-meson-line-pieces.md).

## Result: self-energies

Entries: leading power of $T_\Gamma$ (the self-energy numerator) / of $\delta Z^{\max}$. Same in the four regions.

| no. | lines | generic | $r\le1$ |
|---|---|---|---|
| σ1 | radial, radial | $\epsilon^6/\epsilon^4$ | $\epsilon^{-2}/\epsilon^{-4}$ |
| σ2 | ALP-like, radial | $\epsilon^6/\epsilon^4$ | 0 |
| σ3 | meson, radial | $\epsilon^4/\epsilon^2$ | 0 |
| σ4 | ALP-like, ALP-like | $\epsilon^6/\epsilon^4$ | same |
| σ5 | ALP-like, meson | $\epsilon^4/\epsilon^2$ | same |
| σ6 | meson, meson | $\epsilon^2/\epsilon^2$ | same |
| π1 | ALP-like (closed) | $\epsilon^2/\epsilon^0$ | same |
| π2 | meson (closed; D4 term removed) | $\epsilon^0/\epsilon^0$ | same |
| π3 | meson, radial | $\epsilon^2/\epsilon^0$ | same |
| π4 | ALP-like, meson | $\epsilon^2/\epsilon^0$ | same |
| π5 | meson, meson | $\epsilon^0/\epsilon^{-2}$ | same |

In σ6 the all-generator-sum term (kept, two $\sigma$ ends) has the scale $m_P$, fixed, and sets the leading power $\epsilon^2$ of $\delta Z$; its leading coefficient is $-\frac{\epsilon^2}{16\pi^2f_\pi^2m_P^2}\big(3B_0^2\langle\widehat S_{1,0}\widehat S_{1,0}\rangle+4B_0^2\langle\widehat S_{2,0}\rangle^2+4B_0\kappa_\sigma\langle\widehat S_{2,0}\rangle+\kappa_\sigma^2\big)$. It stays below $\epsilon^4$. For $r\le1$ sources the $\sigma^3$ vertex is $-6\lambda f_\pi/\epsilon$ only, so σ1 drops by eight powers. $\delta Z_{\pi1}+\delta Z_{\pi2}=0$ exactly (both at the scale $m_a$).

**Regions.** Only the coefficients change between regions: σ2 ($1/A$ in region 1, $1/\kappa_\sigma$ in regions 2–4, from $\max(m_\sigma,m_a)$) at the leading power of $R$, and graph 15 (routing choice; e.g. factors $3A-\kappa_\sigma$ in region 1 and $A+\kappa_\sigma$ in region 2) and the other routing-dependent graphs below it.

## Result: graphs that set the leading power

Leading terms (generic sources) of the graphs at the maximal power, with $\{\cdot\}=\{\lambda_{\pi^+},\lambda_{\pi^-}\}$:

- graph 1: $X^{\max}\to-\frac{3B_0\kappa_\sigma}{32\pi^2f_\pi^3}\langle\widehat S_{3,0}\{\cdot\}\rangle\epsilon^5$;
- graph 2: $X^{\max}\to-\frac{AB_0}{32\pi^2f_\pi^3}\langle\widehat S_{1,2}\{\cdot\}\rangle\epsilon^5$;
- graph 10: $X^{\max}\to-\frac{3B_0^2}{16\pi^2f_\pi^3}\langle\widehat S_{3,0}\rangle\langle\widehat S_{2,0}\{\cdot\}\rangle\epsilon^5$;
- graph 11: $X^{\max}\to-\frac{B_0^2}{16\pi^2f_\pi^3}\langle\widehat S_{2,1}\rangle\langle\widehat S_{1,1}\{\cdot\}\rangle\epsilon^5$;
- graph 13: $X^{\max}\to-\frac{B_0}{16\pi^2f_\pi^3}\big(B_0\langle\widehat S_{1,2}\rangle+3A\big)\langle\widehat S_{0,2}\{\cdot\}\rangle\epsilon^5$;
- σ1: $\delta Z^{\max}\to-\frac{9B_0^2}{8\pi^2f_\pi^2\kappa_\sigma}\langle\widehat S_{3,0}\rangle^2\epsilon^4$;
- σ2: $\delta Z^{\max}\to-\frac{B_0^2}{4\pi^2f_\pi^2}\langle\widehat S_{2,1}\rangle^2\epsilon^4/A$ in region 1, $/\kappa_\sigma$ in regions 2–4;
- σ4: $\delta Z^{\max}\to-\frac{1}{2\pi^2f_\pi^2}\Big(A+B_0\langle\widehat S_{1,2}\rangle+\frac{B_0^2}{4A}\langle\widehat S_{1,2}\rangle^2\Big)\epsilon^4$; the terms $A$ and $B_0\langle\widehat S_{1,2}\rangle$ come from the four routing-dependent products, all attained by routing 2.

## Result: $M_1^{\rm est}$ and $R$

$M_1^{\rm est}\sim\epsilon^5$ and $R=|M_1^{\rm est}|/|M_0|\sim\epsilon^4$, for generic and for $r\le1$ sources, in every region; the attaining graphs are 1, 2, 10, 11, 13, σ1, σ2, σ4 (generic) and 2, 13, σ4 ($r\le1$). In region 1 ($\kappa_\sigma<A$) the signed ratio tends to
$$
\frac{M_1^{\rm est}}{M_0}\to\frac{\epsilon^4}{16\pi^2f_\pi^2}\Bigg[\frac{6A\langle\widehat S_{0,2}\{\cdot\}\rangle+A\langle\widehat S_{1,2}\{\cdot\}\rangle+2B_0\langle\widehat S_{0,2}\{\cdot\}\rangle\langle\widehat S_{1,2}\rangle+2B_0\langle\widehat S_{1,1}\{\cdot\}\rangle\langle\widehat S_{2,1}\rangle+6B_0\langle\widehat S_{2,0}\{\cdot\}\rangle\langle\widehat S_{3,0}\rangle+3\kappa_\sigma\langle\widehat S_{3,0}\{\cdot\}\rangle}{\langle\widehat S_{1,0}\{\cdot\}\rangle}
-4A-4B_0\langle\widehat S_{1,2}\rangle-\frac{9B_0^2\langle\widehat S_{3,0}\rangle^2}{\kappa_\sigma}-\frac{B_0^2\langle\widehat S_{1,2}\rangle^2}{A}-\frac{2B_0^2\langle\widehat S_{2,1}\rangle^2}{A}\Bigg] .
$$
In regions 2, 3 and 4 the last term has $1/\kappa_\sigma$ in place of $1/A$; nothing else changes between regions. For $r\le1$ sources the bracket loses its $\widehat S_{1,1}$, $\widehat S_{2,0}$, $\widehat S_{3,0}$ and $\widehat S_{2,1}$ terms. For the up-only sources it reduces to $M_1^{\rm est}/M_0\to-21B_0m_u\epsilon^4/(64\pi^2f_\pi^2)$, from graphs 2, 13 and σ4; the derivation is in the [up-only estimate record](murayama-sigma-pipi-loop-up-only-estimate.md).

## Result: D4 terms (restated criterion)

The removed all-generator-sum terms have the scale $m_P$. The script checks, in every routing, that they contain no letter of the ALP sector ($\mathbb M^{-1}_{m_a^2}(\mathcal C)$, $N$, $\widehat P_{0,1}$, $\widehat S_{r,s\ge1}$, $\widehat P_{r,s\ge1}$): their letters are $\widehat S_{1,0}$, $\lambda_{\pi^\pm}$ (graph 3), $\widehat S_{0,0}$, $\widehat S_{1,0}$, $\lambda_{\pi^\pm}$ (graph 15) and $\lambda_{\pi^\pm}$ (π2), and their masses are $m_P$ (and $m_\sigma$ from the kinematics for graph 15). Their ratio to $M_0$ ($\delta Z$ for π2) depends on $\epsilon$ only through $p_\sigma^2=m_\sigma^2$. The criterion holds for all three. At large $\epsilon$ their leading power relative to $M_0$ is $\epsilon^0$ (graph 3), $\epsilon^2$ (graph 15, from $p_{\pi^+}\cdot p_{\pi^-}$ and $p_\sigma\cdot p_{\pi^\pm}\sim m_\sigma^2$) and $\epsilon^0$ (π2). All stay below $\epsilon^4$, so the D4 terms do not affect the leading power of $R$.

## Notes

1. $\langle\widehat S_{2,0}\rangle$ appears alone in the leading coefficient of σ6: it is a genuine $\langle\widehat S_{1,0}\rangle$ of a completeness sum written with $\langle\widehat S_{1,0}\rangle=\kappa_\sigma/B_0+2\langle\widehat S_{2,0}\rangle$.
2. The boundaries of the regions ($\kappa_\sigma=A,2A,5A$) are degenerate mass relations and are not resolved; the up-only model sits at $\kappa_\sigma=A$, where the up-only coefficient is the same on both sides.
3. The limits of these powers (on-shell cancellations between the kinetic $\sigma aa$ vertex and the $\widehat S_{1,2}$ source vertex in graphs 2, 13 and σ4; validity of $R\ll1$ and of leading-order chiral perturbation theory) are stated in the [limits record](murayama-sigma-pipi-loop-limits.md).

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_large_epsilon.py --markdown
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_large_epsilon
```

[Script](../../../scripts/derive_sigma_pipi_loop_large_epsilon.py), [tests](../../../scripts/test_derive_sigma_pipi_loop_large_epsilon.py): regions, powers per graph and per source variant in every region, the leading coefficients of graph 13 and σ4, $R$ and the attaining graphs, the generic and up-only coefficients, and the D4 criterion and powers, all exact. The per-routing values come from [`derive_sigma_pipi_loop_graphs.py`](../../../scripts/derive_sigma_pipi_loop_graphs.py).
