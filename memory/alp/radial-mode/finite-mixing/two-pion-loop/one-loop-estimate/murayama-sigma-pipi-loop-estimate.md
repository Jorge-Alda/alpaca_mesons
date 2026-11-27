---
title: Size estimate of the one-loop sigma to pi+ pi- amplitude for generic sources and arbitrary epsilon
status: derived
created: 2026-11-26
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Size estimate $M_1^{\rm est}$ of the one-loop $\sigma\to\pi^+\pi^-$ amplitude

Tasks 4.1 and 4.4 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama`. **No loop integral, scalar loop function, loop coefficient, $\mu$ derivative or loop amplitude is computed.** Each loop enters only through the factor $(m_{\rm term}^2)^{2-r_\Gamma}/(16\pi^2)$ per term (D1, D3). $R$ is a size ratio, not a computed correction. The result is valid for any $\epsilon=f_\pi/f_a$: no graph and no term is dropped by a series in $\epsilon$ or $1/\epsilon$, and $N(m_a^2)$, $\mathbb M^{-1}_{m_a^2}(\mathcal C)$, $\epsilon$ and the source traces are kept exact.

## Inputs

- $M_0$, the legs and the assumptions: [tree and legs record](murayama-sigma-pipi-loop-tree-and-legs.md).
- The kept graphs, $S_\Gamma$, $r_\Gamma$ and the scales of their terms: [graph-list record](murayama-sigma-pipi-loop-graph-list.md). Selection rules (D3, D4, D11, D12): [selection record](murayama-sigma-pipi-loop-graph-selection.md). Pole data and meson-line residue: [large-$\epsilon$ record](../asymptotic-behavior/murayama-sigma-pipi-loop-large-epsilon-scaling.md).
- $X_\Gamma^{\max}$ (termwise routing maximum, per-term scales) and the routing spread: [1PI trace record](murayama-sigma-pipi-loop-traces-1pi.md). $\delta Z^{\max}$, the residue prescription and the routing spread of the self-energies: [external-leg record](murayama-sigma-pipi-loop-traces-external-legs.md).
- Size rule: [the estimate record for $a\to\pi^+\pi^-\pi^0$](../../../couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md).

## Result

$$
M_1^{\rm est}=\sum_{\Gamma=1}^{26}X_\Gamma^{\max}+M_0\Big(\tfrac12\delta Z_\sigma^{\max}+\tfrac12\delta Z_{\pi^+}^{\max}+\tfrac12\delta Z_{\pi^-}^{\max}\Big),\qquad
M_0=-\frac{B_0\epsilon}{2f_\pi}\big\langle\widehat S_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\big\rangle,\qquad
R=\frac{|M_1^{\rm est}|}{|M_0|},
$$
with
$$
X_\Gamma^{\max}=\sum_tC_t^{\rho_t^*}\,t,\qquad \sum_tC_t^\rho\,t=\sum_{\rm terms}S_\Gamma\,\frac{(m_{\rm term}^2)^{2-r_\Gamma}}{16\pi^2}\,T_{\rm term}^\rho,\qquad |C_t^{\rho_t^*}|=\max_\rho|C_t^\rho| .
$$
Here $t$ runs over the independent products of traces (canonical basis for generic sources), $\rho$ over the $2r_\Gamma$ routings of the loop momentum, and "terms" over the choices of one piece on each meson line; $m_{\rm term}$ is the largest scale among the lines of the term ($m_\sigma$ radial, $m_a$ ALP-like and pieces (ii)+(iii)+(iv) of a meson line, $m_P$ bare generator sum). $\delta Z^{\max}$ of a leg is the sum over its self-energies of the termwise routing maxima of $\delta Z^\rho$, with the derivative rule applied term by term at $m_{\rm term}$ in each routing. Each graph value is an upper bound per graph on the routing ambiguity of the size rule, not on the true loop: $|\sum_tC_t^\rho t|\le\sum_t|C_t^{\rho_t^*}|\,|t|$ for every routing $\rho$. Where the largest routing depends on the masses (80 products, graphs 13, 14, 15, 20, 21, 25, 26) the coefficient is a Piecewise expression in homogeneous inequalities between squared masses; the branches are listed in the [1PI trace record](murayama-sigma-pipi-loop-traces-1pi.md).

$\delta Z_\sigma^{\max}$ is the sum over σ1–σ6, and $\delta Z_{\pi^+}^{\max}=\delta Z_{\pi^-}^{\max}$ (exactly) the sum over π1–π5; $\delta Z_{\pi1}+\delta Z_{\pi2}=0$ exactly, so $\delta Z_{\pi^\pm}^{\max}$ comes from π3, π4 and π5. The four D4 terms (graphs 3, 15 and π2 for each charged leg; scale $m_P$) are removed exactly. Expanded into products of traces of words, $M_1^{\rm est}$ has 144 terms. Eleven dressed sources enter: $\widehat S_{0,1}$, $\widehat S_{0,2}$, $\widehat S_{1,0}$, $\widehat S_{1,1}$, $\widehat S_{1,2}$, $\widehat S_{2,0}$, $\widehat S_{2,1}$, $\widehat S_{3,0}$, $\widehat P_{0,1}$, $\widehat P_{1,1}$ and $\widehat P_{2,0}$. $\widehat S_{0,0}$ appears in graphs 13, 14 and 15 and cancels exactly in their sum. There is no up-only substitution.

### 1PI graphs

Numbering and lines as in the [graph-list record](murayama-sigma-pipi-loop-graph-list.md). "$m_{\rm term}$" lists the scales of the terms of the graph that contribute (a scale whose terms vanish exactly, such as $m_P$ in graphs whose all-generator-sum term is zero or removed, is not listed). "Routing" gives the number of products whose coefficient depends on the routing with numerical ratios or with a mass-dependent largest routing, and whether one routing attains the maximum for all products.

| no. | vertices | lines | $S_\Gamma$ | $r_\Gamma$ | $m_{\rm term}$ | routing |
|---|---|---|---|---|---|---|
| 1 | $(\pi^+,\pi^-,\sigma)$ | radial | $\frac{1}{2}$ | 1 | $m_\sigma$ | independent |
| 2 | $(\pi^+,\pi^-,\sigma)$ | ALP-like | $\frac{1}{2}$ | 1 | $m_a$ | independent |
| 3 | $(\pi^+,\pi^-,\sigma)$ | meson (D4 term removed) | $\frac{1}{2}$ | 1 | $m_a$ | independent |
| 4 | $(\pi^+,\sigma)$, $(\pi^-)$ | meson, radial | $1$ | 2 | $\max(m_P, m_\sigma)$ | independent |
| 5 | $(\pi^+,\sigma)$, $(\pi^-)$ | ALP-like, meson | $1$ | 2 | $\max(m_P, m_a)$ | independent |
| 6 | $(\pi^+,\sigma)$, $(\pi^-)$ | meson, meson | $\frac{1}{2}$ | 2 | $\max(m_P, m_a)$ | independent |
| 7 | $(\pi^+)$, $(\pi^-,\sigma)$ | meson, radial | $1$ | 2 | $\max(m_P, m_\sigma)$ | independent |
| 8 | $(\pi^+)$, $(\pi^-,\sigma)$ | ALP-like, meson | $1$ | 2 | $\max(m_P, m_a)$ | independent |
| 9 | $(\pi^+)$, $(\pi^-,\sigma)$ | meson, meson | $\frac{1}{2}$ | 2 | $\max(m_P, m_a)$ | independent |
| 10 | $(\pi^+,\pi^-)$, $(\sigma)$ | radial, radial | $\frac{1}{2}$ | 2 | $m_\sigma$ | independent |
| 11 | $(\pi^+,\pi^-)$, $(\sigma)$ | ALP-like, radial | $1$ | 2 | $\max(m_a, m_\sigma)$ | independent |
| 12 | $(\pi^+,\pi^-)$, $(\sigma)$ | meson, radial | $1$ | 2 | $\max(m_a, m_\sigma)$ | independent |
| 13 | $(\pi^+,\pi^-)$, $(\sigma)$ | ALP-like, ALP-like | $\frac{1}{2}$ | 2 | $m_a$ | 20 numerical; 12 mass-dependent; no single routing |
| 14 | $(\pi^+,\pi^-)$, $(\sigma)$ | ALP-like, meson | $1$ | 2 | $m_a$, $\max(m_P, m_a)$ | 20 numerical; 28 mass-dependent; no single routing |
| 15 | $(\pi^+,\pi^-)$, $(\sigma)$ | meson, meson (D4 term removed) | $\frac{1}{2}$ | 2 | $m_a$, $\max(m_P, m_a)$ | 20 numerical; 28 mass-dependent; no single routing |
| 16 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: radial | $1$ | 3 | $\max(m_P, m_\sigma)$ | independent |
| 17 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: radial | $1$ | 3 | $\max(m_P, m_a, m_\sigma)$ | independent |
| 18 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: radial | $1$ | 3 | $\max(m_P, m_a, m_\sigma)$ | independent |
| 19 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: ALP-like | $1$ | 3 | $\max(m_P, m_a, m_\sigma)$ | independent |
| 20 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: ALP-like | $1$ | 3 | $\max(m_P, m_a)$ | 3 mass-dependent; no single routing |
| 21 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: ALP-like | $1$ | 3 | $\max(m_P, m_a)$ | 3 mass-dependent; no single routing |
| 22 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: radial; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_P, m_\sigma)$ | independent |
| 23 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: ALP-like; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_P, m_a)$ | independent |
| 24 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_P, m_a, m_\sigma)$ | independent |
| 25 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_P, m_a)$ | 3 mass-dependent; no single routing |
| 26 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_P, m_a)$ | 3 mass-dependent; no single routing |

For the triangles the lines are listed in the order $\pi^+$–$\pi^-$; $\pi^-$–$\sigma$; $\sigma$–$\pi^+$. A "meson" line carries the full residue $(\mathbf 1_9\oplus1)-v_av_a^{\rm T}$, not a pure meson propagator; the pieces that contribute to each graph are in the [meson-line piece record](murayama-sigma-pipi-loop-meson-line-pieces.md).

### Self-energies of the external legs

The size factor is inside $\delta Z_\Gamma$ (residue prescription of the [external-leg record](murayama-sigma-pipi-loop-traces-external-legs.md)). The $\pi^-$ graphs equal the $\pi^+$ graphs with $\lambda_{\pi^+}\leftrightarrow\lambda_{\pi^-}$.

| no. | leg | lines | $S_\Gamma$ | $r_\Gamma$ | $m_{\rm term}$ | routing |
|---|---|---|---|---|---|---|
| σ1 | $\sigma$ | radial, radial | $\frac{1}{2}$ | 2 | $m_\sigma$ | independent |
| σ2 | $\sigma$ | ALP-like, radial | $1$ | 2 | $\max(m_a, m_\sigma)$ | independent |
| σ3 | $\sigma$ | meson, radial | $1$ | 2 | $\max(m_P, m_\sigma)$, $\max(m_a, m_\sigma)$ | independent |
| σ4 | $\sigma$ | ALP-like, ALP-like | $\frac{1}{2}$ | 2 | $m_a$ | 4 numerical; attained by routing 2 |
| σ5 | $\sigma$ | ALP-like, meson | $1$ | 2 | $m_a$, $\max(m_P, m_a)$ | 4 numerical; attained by routing 2 |
| σ6 | $\sigma$ | meson, meson | $\frac{1}{2}$ | 2 | $m_P$, $m_a$, $\max(m_P, m_a)$ | 4 numerical; attained by routing 2 |
| π1 | $\pi^\pm$ | ALP-like (closed on the vertex) | $\frac{1}{2}$ | 1 | $m_a$ | independent |
| π2 | $\pi^\pm$ | meson (closed on the vertex) (D4 term removed) | $\frac{1}{2}$ | 1 | $m_a$ | independent |
| π3 | $\pi^\pm$ | meson, radial | $1$ | 2 | $\max(m_P, m_\sigma)$ | independent |
| π4 | $\pi^\pm$ | ALP-like, meson | $1$ | 2 | $\max(m_P, m_a)$ | independent |
| π5 | $\pi^\pm$ | meson, meson | $\frac{1}{2}$ | 2 | $\max(m_P, m_a)$ | independent |

### Large-$\epsilon$ behaviour

The estimate is the full exact sum above and is valid for any $\epsilon$; no term is dropped. Its large-$\epsilon$ behaviour at fixed sources and fixed $\lambda$ is in the [large-$\epsilon$ graph record](../asymptotic-behavior/murayama-sigma-pipi-loop-large-epsilon-graphs.md): $R\sim\epsilon^4$, from graphs 1, 2, 10, 11, 13, σ1, σ2 and σ4, in every region of $\kappa_\sigma/A$. For the up-only sources $M_1^{\rm est}/M_0\to-21B_0m_u\epsilon^4/(64\pi^2f_\pi^2)$, from graphs 2, 13 and σ4 ([up-only estimate record](../asymptotic-behavior/murayama-sigma-pipi-loop-up-only-estimate.md)).

**Vanishing tree.** Where $\langle\widehat S_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle=0$, $M_0=0$ and the external-leg terms vanish identically. $R$ is then undefined, and only $M_1^{\rm est}$, which is nonzero, is reported. The script returns $R$ as undefined there and does not divide.

## Decisions and assumptions entering the result

1. **Routing (D12).** Termwise maximum over the $2r_\Gamma$ routings, as above. Ties of equal magnitude and opposite sign take the sign of the sum of the coefficients over all routings; if this sum vanishes, the first routing in the canonical order is used and the product is marked (a choice at a degenerate mass relation). In 80 of the products of traces of $M_1^{\rm est}$ (all in the 1PI graphs 13, 14, 15, 20, 21, 25, 26; none in a self-energy) the largest routing depends on the masses. For each of them the routing coefficients share one loop scale $s$ and are, up to a common positive factor, affine in one ratio $w$; $M_1^{\rm est}$ stores the coefficient as a Piecewise expression whose branches are the coefficients of the chosen routings and whose conditions are homogeneous inequalities between squared masses (no number is inserted for any mass):
   - 59 products (27 with $s=m_a$, 32 with $s=\max(m_P,m_a)$; graphs 13, 14, 15): routing 3 for $m_\sigma^2-m_{\pi^\pm}^2\le s^2$, routing 4 for $m_\sigma^2-m_{\pi^\pm}^2>s^2$;
   - 9 products ($s=m_a$; graphs 13, 14, 15): routing 3 for $m_\sigma^2-m_{\pi^\pm}^2\le2s^2$ or $\ge5s^2$, routing 4 for $2s^2<m_\sigma^2-m_{\pi^\pm}^2<5s^2$;
   - 12 products ($s=\max(m_P,m_a)$; graphs 20, 21, 25, 26, three each): routing 1 for $m_\sigma^2-2m_{\pi^\pm}^2\le2s^2$, routing 2 above.

   The comparison is exact on $w>0$, which follows from the physical constraint $m_\sigma\ge2m_{\pi^\pm}$ (open decay). At each boundary point the chosen branch is the one whose value coincides with the maximum there. The tie rule with a vanishing sum over routings is needed only at the boundary point $m_\sigma^2-m_{\pi^\pm}^2=s^2$ of the first pattern (59 products, marked as "first routing at a boundary point"); no product with numerical routing ratios needs it, and no tie occurs on an open region. The routing spread is kept as a diagnostic in the trace records.
2. **Loop scale per term (D3, assumption).** $m_P$ is the largest pole mass of the nine meson-like states, kept as a symbol. The mixed entries (iii) of a meson line, like pieces (ii) and (iv), have the scale $m_a$; with a common scale the ALP–meson mass mixing is neglected inside the loops, and it is kept in the vertices and on the external legs. No lower limit on $m_{\rm term}$ other than the line scales is imposed.
3. **Residue prescription.** $\delta Z=-\Pi'(m^2)$ with the size rule applied to the derivative, term by term, as stated in the [external-leg record](murayama-sigma-pipi-loop-traces-external-legs.md).
4. **$\sigma$ tadpole.** The one-loop $\sigma$ tadpole is set to zero by assumption (a vacuum condition on $f_a$, consistent with eliminating $\lambda$ through $m_\sigma^2$). The tadpoles of $\pi^0$, $\eta$, $\eta'$ and of the ALP field vanish exactly under CP.

The one-particle-reducible graphs in which an external leg mixes with another state are not included. They vanish under D11: CP-invariant source terms with $\sigma$ CP-even and $a$, $\Phi$ CP-odd (no $\sigma$–$a$ or $\sigma$–pseudoscalar mixing at any loop order), and flavour-diagonal sources and quark masses (no $\pi^\pm$–$K^\pm$ mixing). Both are assumptions, stated in the [tree and legs record](murayama-sigma-pipi-loop-tree-and-legs.md) (assumption 10) and in the [selection record](murayama-sigma-pipi-loop-graph-selection.md) (rule 4).

The limits (absorptive parts and thresholds, cancellations, naturalness without counterterm analysis, the scale assignment and the D11 assumptions) are in the [limits record](murayama-sigma-pipi-loop-limits.md).

## Verification

- Dimensions, with $X\to c^{d_X}X$: $d=1$ for $B_0$, $f_\pi$, $m_\sigma$, $m_a$, $m_{\pi^\pm}$, $m_P$ and for every source letter in a trace; $d=0$ for $\epsilon$, $N(m_a^2)$, $\lambda_{\pi^\pm}$ and $\mathbb M^{-1}_{m_a^2}(\mathcal C)$. Every maximum of masses (two or three arguments) is rescaled with the exact identity $\max(cx,cy,\dots)=c\max(x,y,\dots)$ for $c>0$: the script substitutes a mass symbol of dimension one for it during the rescaling. Every condition of a Piecewise coefficient is checked to be homogeneous, so it is unchanged by the rescaling (a non-homogeneous condition raises an error, tested). Under this rescaling $M_0\to cM_0$, each $X_\Gamma^{\max}\to cX_\Gamma^{\max}$, $\delta Z^{\max}\to\delta Z^{\max}$, $M_1^{\rm est}\to cM_1^{\rm est}$ and hence $R\to R$. All of these are exact symbolic identities, checked branch by branch.
- $M_1^{\rm est}$ equals the sum of the 26 values $X_\Gamma^{\max}$ and $M_0\,(\delta Z_\sigma^{\max}+\delta Z_{\pi^+}^{\max}+\delta Z_{\pi^-}^{\max})/2$; $\delta Z_{\pi^+}^{\max}=\delta Z_{\pi^-}^{\max}$ exactly. No floating-point number appears.
- The routing tests of task 4.4 (single-routing graphs unchanged, $|C_t^{\rho_t^*}|\ge|C_t^\rho|$, absence of $m_a^2-m_\sigma^2$ from $\delta Z_{\pi^\pm}$, field-basis identity) are described in the [1PI trace record](murayama-sigma-pipi-loop-traces-1pi.md) and run in the graph tests.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_estimate.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_estimate
```

[Script](../../../../../../scripts/derive_sigma_pipi_loop_estimate.py), [tests](../../../../../../scripts/test_derive_sigma_pipi_loop_estimate.py).
