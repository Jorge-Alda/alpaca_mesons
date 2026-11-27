---
title: Which pieces of the meson-line residue contribute to each kept sigma to pi+ pi- graph
status: derived
created: 2026-11-26
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Pieces of the meson-line residue in the kept $\sigma\to\pi^+\pi^-$ graphs

OpenSpec change `estimate-sigma-pipi-one-loop-murayama`. A meson line is not a pure meson line once the ALP mixes; this record lists which pieces of its residue contribute to each kept graph. The statements are exact, for generic sources, at any $\epsilon$. No loop is evaluated.

## Split of the residue

An internal meson line carries the residue $(\mathbf 1_9\oplus1)-v_av_a^{\rm T}$ on the space (nine meson matrices) $\oplus$ (ALP field), with $v_a=\big(2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)/\sqrt N,\,1/\sqrt N\big)$ and $N=N(m_a^2)$ ([large-$\epsilon$ record](../asymptotic-behavior/murayama-sigma-pipi-loop-large-epsilon-scaling.md)). It is split exactly into four pieces:

- (i) bare generator sum $\mathbf 1_9$ (both ends meson);
- (ii) ALP-field entry $1-1/N=(N-1)/N$ (both ends ALP field);
- (iii) mixed entries $-2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)/N$ (one end ALP field, one end meson);
- (iv) meson block of the ALP-like pole, $-4\epsilon^2\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\otimes\mathbb M^{-1}_{m_a^2}(\mathcal C)/N$ (both ends meson).

A vertex with three meson ends vanishes ($\widehat P_{0,0}=0$), and so does $\sigma$ with three meson ends ($\widehat P_{1,0}=0$) ([vertex-rule record](murayama-sigma-pipi-loop-vertex-rules.md)). The four pieces reproduce every graph value exactly, routing by routing (tested). Loop scales (D3, [selection record](murayama-sigma-pipi-loop-graph-selection.md)): piece (i) has the scale $m_P$ (largest pole mass of the nine meson-like states), pieces (ii), (iii) and (iv), which together are $(0\oplus1)-v_av_a^{\rm T}$, the remainder of the ALP-like propagator, have the scale $m_a$; a term has the largest scale among its lines. Their large-$\epsilon$ behaviour is in the [large-$\epsilon$ input record](../asymptotic-behavior/murayama-sigma-pipi-loop-large-epsilon-inputs.md).

## Result

Measure: the size-rule value with the per-term scales ($X_\Gamma$ for 1PI graphs, $\delta Z$ for self-energies); a piece "contributes" if its value is nonzero in some routing of the loop momentum (exact, after $N=1+2\epsilon^2\langle(\mathbb M^{-1}_{m_a^2}(\mathcal C))^2\rangle$ is substituted). "Line $j$" joins vertex $j$ to vertex $j+1$ in the vertex list of the [graph-list record](murayama-sigma-pipi-loop-graph-list.md), starting from $0$; for triangles line 0 is $\pi^+$–$\pi^-$, line 1 is $\pi^-$–$\sigma$, line 2 is $\sigma$–$\pi^+$.

| graph (meson line) | pieces that contribute |
|---|---|
| 4 (0), 5 (1), 7 (0), 8 (1), 16–20 (0), 21 (0), 22 (1, 2), 23 (1, 2), 24 (0), 25 (0), π3, π4 | (i) only |
| 6, 9, π5 (both lines); 26 (line 0) | (i) and (ii); the term with (i) on every meson line is exactly zero, so each contribution has at least one line in its ALP-field entry |
| 3 (closed line, D4 term removed) | (ii), (iii), (iv) |
| 12 (0), 18 (1), 21 (1), 24 (2), 25 (2) | (ii) and (iii), both at the leading power of the graph |
| 14 (1), σ3, σ5, σ6 | (i), (ii), (iii), (iv) |
| 15 (both lines, D4 term removed) | (i), (ii), (iii), (iv); at large $\epsilon$ (region $\kappa_\sigma<A$) the leading $\epsilon^3$ of the termwise routing maximum comes from (iii) on both lines and from (i) on one line with (iii) on the other |
| 26 (lines 1, 2) | (i), (ii), (iii) |
| π2 (closed line, D4 term removed) | (iv) only: the ALP-field entry and the mixed entries give no $p^2$ dependence, so no $\delta Z$; with the scale $m_a$ it cancels exactly against the closed ALP-like line π1 |

## Consequences for the wording

- No meson line contributes through its ALP-field entry (ii) alone.
- The description "bare generators, whose term with every line bare vanishes at a three-meson vertex, plus the ALP-field part weighted by $(N-1)/N$" is exact for graphs 4–9, 16, 17, 19, 20, 22, 23, π3, π4, π5 and line 0 of graph 26.
- It is not exact for graphs 3, 12, 14, 15, 18, 21, 24, 25, lines 1 and 2 of 26, σ3, σ5, σ6 and π2: there the mixed entries (iii) and/or the meson block (iv) contribute, often at the leading power.
- In the tables of the [graph-list record](murayama-sigma-pipi-loop-graph-list.md) and the [estimate record](murayama-sigma-pipi-loop-estimate.md) a line called "meson" means the full residue above; the pieces that contribute are those of this table.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_large_epsilon.py --markdown
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_large_epsilon
```

[Script](../../../../../../scripts/derive_sigma_pipi_loop_large_epsilon.py) (survival table; the large-$\epsilon$ powers of each piece are those of its termwise routing maximum in the region $\kappa_\sigma<A$), [graph values](../../../../../../scripts/derive_sigma_pipi_loop_graphs.py) with the optional `piece_function` and `choice_filter` arguments (the four-piece split carries the scale of each piece).
