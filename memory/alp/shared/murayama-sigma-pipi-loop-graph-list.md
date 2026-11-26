---
title: One-loop graphs of sigma to pi+ pi-: list, symmetry factors, loop scales and excluded terms
status: derived
created: 2026-11-26
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# One-loop graphs of $\sigma\to\pi^+\pi^-$

Tasks 2.1, 2.3 and 4.4 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama`. No loop is evaluated.

## Inputs and assumptions

- Vertices: [vertex-rule record](murayama-sigma-pipi-loop-vertex-rules.md). Legs, kinematics, charge rule and assumptions: [tree and legs record](murayama-sigma-pipi-loop-tree-and-legs.md). Selection rules, CP and flavour assumption (D11), per-term loop scale (D3) and routing rule (D12): [selection record](murayama-sigma-pipi-loop-graph-selection.md). Pole data and meson-line residue: [large-$\epsilon$ record](murayama-sigma-pipi-loop-large-epsilon-scaling.md).
- Internal lines, as in [the estimate record](murayama-radial-up-only-estimate-result.md): a radial line has $\sigma$ at both ends. An ALP-like line has $v_a=\big(2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)/\sqrt{N(m_a^2)},\,1/\sqrt{N(m_a^2)}\big)$ at both ends. A meson line has the residue $(\mathbf 1_9\oplus1)-v_av_a^{\rm T}$, split exactly into three terms: the $U(3)$ generator sum at both ends, the ALP-field entry at both ends, and minus $v_a$ at both ends. "Meson" in the lines column means this full residue, not a pure meson line; the pieces that contribute to each graph are in the [meson-line piece record](murayama-sigma-pipi-loop-meson-line-pieces.md). Leading powers at large $\epsilon$: [large-$\epsilon$ graph record](murayama-sigma-pipi-loop-large-epsilon-graphs.md).
- $S_\Gamma=1/|{\rm Aut}_\Gamma|$, where ${\rm Aut}_\Gamma$ is the set of rotations and reflections of the loop, acting on the line ends, that keep the external legs and the line types. $r_\Gamma$ is the number of internal lines.

## Topologies

A one-loop 1PI graph is a loop of $V$ vertices. Each vertex carries two line ends and at least one external leg, because there is no two-point vertex: quadratic terms are part of the propagators. With the three legs $\sigma$, $\pi^+$, $\pi^-$, the number of ways to distribute the legs over $V$ vertices is $1,6,6,0,0$ for $V=1,\dots,5$. So the 1PI graphs are tadpoles (one five-leg vertex), bubbles (a three-leg and a four-leg vertex) and triangles (three three-leg vertices). **No box exists**, since it would need four vertices with an external leg each. **No $(2,5)$ 1PI bubble exists**, since it would need a vertex with two lines and no external leg. The self-energy of an external leg is one-particle reducible; it enters through the residue $\delta Z=-\Pi'(m^2)$. A tree with a loop on an internal propagator does not exist, because the tree is a single contact vertex. The one-particle-reducible graphs in which an external leg mixes with another state vanish under D11 (CP invariance with $\sigma$ CP-even and $a$, $\Phi$ CP-odd; flavour-diagonal sources and quark masses for $\pi^\pm$), and are not listed.

## Counts

- 1PI: 48 graphs (3 tadpoles, 18 bubbles, 27 triangles). 26 kept (3 tadpoles, 12 bubbles, 11 triangles); 22 vanish identically (6 bubbles, 16 triangles).
- Self-energies: 9 graphs per leg (3 tadpoles, 6 bubbles). Kept: 6 for $\sigma$, 5 for each of $\pi^\pm$.
- Total kept: 42. No graph is dropped by a series in $\epsilon$ or $1/\epsilon$.
- An excluded all-generator-sum term (D4; scale $m_P$; no letter of the ALP sector, and its ratio to $M_0$ depends on $\epsilon$ only through $p_\sigma^2=m_\sigma^2$) occurs in 4 graphs (1PI nos. 3 and 15, and the closed meson line on each $\pi^\pm$ leg).

## 1PI graphs

The vertices column lists the external legs at each vertex. The lines column lists the internal lines around the loop. "$m_{\rm term}$" lists the loop scales of the terms that contribute (D3: for each choice of one piece on every meson line, the largest of $m_\sigma$ for a radial line, $m_a$ for an ALP-like line and for the ALP-field entry and ALP-like pole of a meson line, and $m_P$ for the bare generator sum); for graphs that vanish identically it lists the scales of all terms. The all-generator-sum term column applies to meson-only graphs: "excluded (D4)" means the term is removed and satisfies the restated D4 criterion; "exactly zero" means the term vanishes, and only the ALP-field and ALP-like-pole terms remain.

| no. | vertices | lines | $S_\Gamma$ | $r_\Gamma$ | $m_{\rm term}$ | all-generator-sum term | status |
|---|---|---|---|---|---|---|---|
| 1 | $(\pi^+,\pi^-,\sigma)$ | radial (closed on the vertex) | $\frac{1}{2}$ | 1 | $m_\sigma$ | — | kept |
| 2 | $(\pi^+,\pi^-,\sigma)$ | ALP-like (closed on the vertex) | $\frac{1}{2}$ | 1 | $m_a$ | — | kept |
| 3 | $(\pi^+,\pi^-,\sigma)$ | meson (closed on the vertex) | $\frac{1}{2}$ | 1 | $m_a$ | excluded (D4) | kept (D4 term removed) |
|  | $(\pi^+,\sigma)$, $(\pi^-)$ | radial, radial | $\frac{1}{2}$ | 2 | $m_\sigma$ | — | vanishes identically |
|  | $(\pi^+,\sigma)$, $(\pi^-)$ | ALP-like, radial | $1$ | 2 | $\max(m_a, m_\sigma)$ | — | vanishes identically |
| 4 | $(\pi^+,\sigma)$, $(\pi^-)$ | meson, radial | $1$ | 2 | $\max(m_P, m_\sigma)$ | — | kept |
|  | $(\pi^+,\sigma)$, $(\pi^-)$ | ALP-like, ALP-like | $\frac{1}{2}$ | 2 | $m_a$ | — | vanishes identically |
| 5 | $(\pi^+,\sigma)$, $(\pi^-)$ | ALP-like, meson | $1$ | 2 | $\max(m_P, m_a)$ | — | kept |
| 6 | $(\pi^+,\sigma)$, $(\pi^-)$ | meson, meson | $\frac{1}{2}$ | 2 | $\max(m_P, m_a)$ | exactly zero | kept |
|  | $(\pi^+)$, $(\pi^-,\sigma)$ | radial, radial | $\frac{1}{2}$ | 2 | $m_\sigma$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-,\sigma)$ | ALP-like, radial | $1$ | 2 | $\max(m_a, m_\sigma)$ | — | vanishes identically |
| 7 | $(\pi^+)$, $(\pi^-,\sigma)$ | meson, radial | $1$ | 2 | $\max(m_P, m_\sigma)$ | — | kept |
|  | $(\pi^+)$, $(\pi^-,\sigma)$ | ALP-like, ALP-like | $\frac{1}{2}$ | 2 | $m_a$ | — | vanishes identically |
| 8 | $(\pi^+)$, $(\pi^-,\sigma)$ | ALP-like, meson | $1$ | 2 | $\max(m_P, m_a)$ | — | kept |
| 9 | $(\pi^+)$, $(\pi^-,\sigma)$ | meson, meson | $\frac{1}{2}$ | 2 | $\max(m_P, m_a)$ | exactly zero | kept |
| 10 | $(\pi^+,\pi^-)$, $(\sigma)$ | radial, radial | $\frac{1}{2}$ | 2 | $m_\sigma$ | — | kept |
| 11 | $(\pi^+,\pi^-)$, $(\sigma)$ | ALP-like, radial | $1$ | 2 | $\max(m_a, m_\sigma)$ | — | kept |
| 12 | $(\pi^+,\pi^-)$, $(\sigma)$ | meson, radial | $1$ | 2 | $\max(m_a, m_\sigma)$ | — | kept |
| 13 | $(\pi^+,\pi^-)$, $(\sigma)$ | ALP-like, ALP-like | $\frac{1}{2}$ | 2 | $m_a$ | — | kept |
| 14 | $(\pi^+,\pi^-)$, $(\sigma)$ | ALP-like, meson | $1$ | 2 | $m_a$, $\max(m_P, m_a)$ | — | kept |
| 15 | $(\pi^+,\pi^-)$, $(\sigma)$ | meson, meson | $\frac{1}{2}$ | 2 | $m_a$, $\max(m_P, m_a)$ | excluded (D4) | kept (D4 term removed) |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: radial; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: radial | $1$ | 3 | $m_\sigma$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: radial; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: radial | $1$ | 3 | $\max(m_a, m_\sigma)$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: radial; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: radial | $1$ | 3 | $\max(m_P, m_\sigma)$, $\max(m_a, m_\sigma)$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: ALP-like; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: radial | $1$ | 3 | $\max(m_a, m_\sigma)$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: ALP-like; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: radial | $1$ | 3 | $\max(m_a, m_\sigma)$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: ALP-like; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: radial | $1$ | 3 | $\max(m_a, m_\sigma)$, $\max(m_P, m_a, m_\sigma)$ | — | vanishes identically |
| 16 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: radial | $1$ | 3 | $\max(m_P, m_\sigma)$ | — | kept |
| 17 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: radial | $1$ | 3 | $\max(m_P, m_a, m_\sigma)$ | — | kept |
| 18 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: radial | $1$ | 3 | $\max(m_P, m_a, m_\sigma)$ | — | kept |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: radial; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: ALP-like | $1$ | 3 | $\max(m_a, m_\sigma)$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: radial; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: ALP-like | $1$ | 3 | $\max(m_a, m_\sigma)$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: radial; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: ALP-like | $1$ | 3 | $\max(m_a, m_\sigma)$, $\max(m_P, m_a, m_\sigma)$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: ALP-like; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: ALP-like | $1$ | 3 | $\max(m_a, m_\sigma)$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: ALP-like; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: ALP-like | $1$ | 3 | $m_a$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: ALP-like; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: ALP-like | $1$ | 3 | $m_a$, $\max(m_P, m_a)$ | — | vanishes identically |
| 19 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: ALP-like | $1$ | 3 | $\max(m_P, m_a, m_\sigma)$ | — | kept |
| 20 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: ALP-like | $1$ | 3 | $\max(m_P, m_a)$ | — | kept |
| 21 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: ALP-like | $1$ | 3 | $\max(m_P, m_a)$ | — | kept |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: radial; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_P, m_\sigma)$, $\max(m_a, m_\sigma)$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: radial; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_a, m_\sigma)$, $\max(m_P, m_a, m_\sigma)$ | — | vanishes identically |
| 22 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: radial; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_P, m_\sigma)$ | — | kept |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: ALP-like; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_a, m_\sigma)$, $\max(m_P, m_a, m_\sigma)$ | — | vanishes identically |
|  | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: ALP-like; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: meson | $1$ | 3 | $m_a$, $\max(m_P, m_a)$ | — | vanishes identically |
| 23 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: ALP-like; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_P, m_a)$ | — | kept |
| 24 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: radial; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_P, m_a, m_\sigma)$ | — | kept |
| 25 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: ALP-like; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_P, m_a)$ | — | kept |
| 26 | $(\pi^+)$, $(\pi^-)$, $(\sigma)$ | $\pi^+$–$\pi^-$: meson; $\pi^-$–$\sigma$: meson; $\sigma$–$\pi^+$: meson | $1$ | 3 | $\max(m_P, m_a)$ | exactly zero | kept |

## Self-energies of the external legs

The external legs carry incoming momenta $p$ and $-p$. For $\pi^+$ the two legs are $\lambda_{\pi^+}$ and $\lambda_{\pi^-}$. The $\pi^-$ graphs are the same with the two exchanged; they have the same statuses. "Zero $\delta Z$ (exact)" means a tadpole whose vertex does not depend on $p$. It only shifts the mass, and the physical masses are used.

| leg | no. | lines | $S_\Gamma$ | $r_\Gamma$ | $m_{\rm term}$ | all-generator-sum term | status |
|---|---|---|---|---|---|---|---|
| $\sigma$ |  | radial (closed on the vertex) | $\frac{1}{2}$ | 1 | $m_\sigma$ | — | zero $\delta Z$ (exact) |
| $\sigma$ |  | ALP-like (closed on the vertex) | $\frac{1}{2}$ | 1 | $m_a$ | — | zero $\delta Z$ (exact) |
| $\sigma$ |  | meson (closed on the vertex) | $\frac{1}{2}$ | 1 | $m_P$, $m_a$ | exactly zero | zero $\delta Z$ (exact) |
| $\sigma$ | σ1 | radial, radial | $\frac{1}{2}$ | 2 | $m_\sigma$ | — | kept |
| $\sigma$ | σ2 | ALP-like, radial | $1$ | 2 | $\max(m_a, m_\sigma)$ | — | kept |
| $\sigma$ | σ3 | meson, radial | $1$ | 2 | $\max(m_P, m_\sigma)$, $\max(m_a, m_\sigma)$ | — | kept |
| $\sigma$ | σ4 | ALP-like, ALP-like | $\frac{1}{2}$ | 2 | $m_a$ | — | kept |
| $\sigma$ | σ5 | ALP-like, meson | $1$ | 2 | $m_a$, $\max(m_P, m_a)$ | — | kept |
| $\sigma$ | σ6 | meson, meson | $\frac{1}{2}$ | 2 | $m_P$, $m_a$, $\max(m_P, m_a)$ | kept (two $\sigma$ ends) | kept |
| $\pi^+$ |  | radial (closed on the vertex) | $\frac{1}{2}$ | 1 | $m_\sigma$ | — | zero $\delta Z$ (exact) |
| $\pi^+$ | π1 | ALP-like (closed on the vertex) | $\frac{1}{2}$ | 1 | $m_a$ | — | kept |
| $\pi^+$ | π2 | meson (closed on the vertex) | $\frac{1}{2}$ | 1 | $m_a$ | excluded (D4) | kept (D4 term removed) |
| $\pi^+$ |  | radial, radial | $\frac{1}{2}$ | 2 | $m_\sigma$ | — | vanishes identically |
| $\pi^+$ |  | ALP-like, radial | $1$ | 2 | $\max(m_a, m_\sigma)$ | — | vanishes identically |
| $\pi^+$ | π3 | meson, radial | $1$ | 2 | $\max(m_P, m_\sigma)$ | — | kept |
| $\pi^+$ |  | ALP-like, ALP-like | $\frac{1}{2}$ | 2 | $m_a$ | — | vanishes identically |
| $\pi^+$ | π4 | ALP-like, meson | $1$ | 2 | $\max(m_P, m_a)$ | — | kept |
| $\pi^+$ | π5 | meson, meson | $\frac{1}{2}$ | 2 | $\max(m_P, m_a)$ | exactly zero | kept |

## Remarks

- **Meson-only graphs.** For the $(M,M)$ bubbles with a single pion at one vertex and for the meson-only triangle, the all-generator-sum term is exactly zero: a vertex with three meson legs needs $\widehat P_{0,0}$, which vanishes. Their terms have the scale $\max(m_P,m_a)$ (one bare line and one line in a piece at $m_a$). In the $(M,M)$ bubbles (and in π5) only the bare generators and the ALP-field entry $1-1/N(m_a^2)$ contribute, with at least one line in its ALP-field entry; in the triangle the $\pi^-$–$\sigma$ and $\sigma$–$\pi^+$ lines also contribute through the mixed entries ([meson-line piece record](murayama-sigma-pipi-loop-meson-line-pieces.md)). The meson-only $\sigma$ self-energy keeps all its terms, because it has two $\sigma$ ends; its all-generator-sum term has the scale $m_P$.
- **Scales.** The scale of every term comes from its lines only (D3); no scale is tied to $m_\sigma$ in a graph with only meson lines.
- **Routing.** Ten kept graphs depend on the routing (13, 14, 15, 20, 21, 25, 26, σ4, σ5, σ6); they enter through the termwise routing maximum. In 80 products the largest routing depends on the masses (Piecewise coefficients; [1PI trace record](murayama-sigma-pipi-loop-traces-1pi.md)).

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_graphs.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_graphs
```

[Script](../../../scripts/derive_sigma_pipi_loop_graphs.py), [tests](../../../scripts/test_derive_sigma_pipi_loop_graphs.py). The tests check the counts, the absence of boxes and $(2,5)$ graphs, the symmetry factors, that every graph is kept or vanishes exactly for every routing, the excluded terms and the restated D4 criterion, the exactness of the split (all terms = all-generator-sum term + the rest, routing by routing), the per-term scales, and the routing tests of task 4.4 (graphs independent of the routing enter with their single-routing value; $|C_t^{\rho_t^*}|\ge|C_t^\rho|$; the three branch patterns; $\delta Z_{\pi1}+\delta Z_{\pi2}=0$; the field-basis identity).
