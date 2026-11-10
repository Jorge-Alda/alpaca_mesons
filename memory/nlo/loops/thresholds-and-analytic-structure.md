---
title: Thresholds and analytic structure of the one-loop two-point functions
status: derived
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, Gasser:1983yg, Gasser:1984gg, Kaiser:1998ds]
provenance: analytical-derivation
---

# Thresholds and analytic structure

## Inputs and assumptions

The vertices and allowed graphs are those of [vertices and topologies](vertices-and-topologies.md); the tadpole and subtraction are given in [loop functions and subtraction](loop-functions-and-subtraction.md). The scope is one loop with LO meson vertices and first order in $1/f_a$, with no CP-violating sources, constant $\theta$ term or resonance loops. Internal LO masses are taken real and positive. Masses and $\mu$ are in GeV and $s$ and squared thresholds in GeV$^2$ [@Alda:2026nrd; @Gasser:1983yg; @Kaiser:1998ds].

## Retained functions

| Function | Topology and analytic structure |
|---|---|
| Tadpole contraction | Its finite part is $\frac{m_k^2}{16\pi^2}\log(m_k^2/\mu^2)$. It is independent of $s$, so it has no branch cut or threshold in $s$ and no absorptive part. The logarithm is non-analytic as a function of the mass squared, not an $s$-channel threshold. |
| $\Sigma_{ij}(s)$ | The $\Phi^4$ tadpole is the only retained graph. Its dependence on $s$ is a polynomial of degree at most one, with no $s$-channel cut; $\operatorname{Im}\Sigma_{ij}(s)=0$ for real $s$ at this order. |
| $\Sigma_{i\hat a}(s)$ | The $c_q^A$ and $c_q^P$ vertices give tadpoles; the $c_q^V$, $c_q^S$ and $c_G$ tadpole contributions vanish. The result is a polynomial of degree at most one in $s$, with no $s$-channel cut and $\operatorname{Im}\Sigma_{i\hat a}(s)=0$ for real $s$ at this order. |

The generic $c_q^V$ and $c_q^S$ vertices $a\Phi^2$ cannot form a mixed bubble without a $\Phi^3$ vertex, which is absent under these assumptions. Two such vertices can form a bubble in the ALP self-energy, at order $1/f_a^2$, outside the present calculation. Local subtraction terms do not introduce a cut. The statements above therefore hold for the retained contributions at all $s$; they do not determine thresholds or absorptive parts from omitted resonance loops or higher orders.

## Conditional CP-odd extension

If CP-violating sources or a constant $\theta$ term are included, a $\Phi^3$ vertex can be present. Where the relevant couplings and quantum numbers allow it, a bubble with internal mesons of masses $m_k$ and $m_l$ can then have a two-meson threshold at $(m_k+m_l)^2$ GeV$^2$ when the masses are expressed in GeV. This is a conditional extension, not a contribution established in the present calculation; its channels and imaginary parts require a separate calculation.

## Verification

The topology list in [vertices and topologies](vertices-and-topologies.md) identifies the required quartic vertices and excludes both bubbles at the retained order. The explicit [meson](meson-self-energy.md) and [mixed](mixed-alp-meson-self-energy.md) expressions use the momentum-independent tadpole contraction and contain at most one power of $s$.
