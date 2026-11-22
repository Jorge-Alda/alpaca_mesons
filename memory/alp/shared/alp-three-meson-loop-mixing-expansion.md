---
title: Mixing expansion of meson lines in the three-meson loop estimate
status: derived
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Mixing expansion of meson lines

The [dressed inverse propagator](alp-lo-exact-pole-equations.md#step-3-exact-alp-pole-equation) and the [tree pole projections](alp-finite-mixing-three-meson-tree.md#result) are the inputs. Work in the unmixed meson eigenbasis with generators $\lambda_n$, with the dressed kinetic matrix equal to the identity. Keep the ALP-like pole mass $m_a$ fixed. In the following equations, the neutral mesons with nonzero mixing are separated from $m_a$ and from one another. The unmixed charged sectors can be projected at fixed charge, since the ALP has zero charge.

## Pole positions and projectors

The mixed inverse-propagator entry is $\epsilon\langle\lambda_n\mathcal C\rangle$. Solving the ALP pole equation for the quadratic ALP entry gives
$$
D_{aa}(s)=s-m_a^2+
\epsilon^2\sum_j\frac{\langle\lambda_j\mathcal C\rangle^2}{m_a^2-m_{j,0}^2}.
$$
The meson eigenvector equation then gives the mass shift and the ALP component,
$$
m_n^2=m_{n,0}^2+\epsilon^2
\frac{\langle\lambda_n\mathcal C\rangle^2}{m_{n,0}^2-m_a^2}
+O(\epsilon^4),\qquad
v_a^{(n)}=-\epsilon\frac{\langle\lambda_n\mathcal C\rangle}{m_{n,0}^2-m_a^2}
+O(\epsilon^3).
$$
Projecting back onto the meson eigenbasis and normalizing gives
$$
v_j^{(n)}
=\epsilon^2\frac{\langle\lambda_j\mathcal C\rangle\langle\lambda_n\mathcal C\rangle}
{(m_{n,0}^2-m_{j,0}^2)(m_{n,0}^2-m_a^2)}
+O(\epsilon^4),\qquad j\ne n,
$$
$$
v_n^{(n)}
=1-\frac{\epsilon^2}{2}
\frac{\langle\lambda_n\mathcal C\rangle^2}{(m_{n,0}^2-m_a^2)^2}
+O(\epsilon^4).
$$
The off-diagonal meson component comes from inserting the first-order ALP component into the meson row of the eigenvector equation. The parallel meson component follows from unit norm. For an exactly degenerate neutral meson subspace use the spectral projector after diagonalizing the projected mixing perturbation within that subspace.

## Expansion of internal lines

For the full tree propagator $G=D^{-1}$, expanding the inverse gives, in particular,
$$
G_{jk}(s)=\frac{\delta_{jk}}{s-m_{j,0}^2}
+\epsilon^2
\frac{\langle\lambda_j\mathcal C\rangle\langle\lambda_k\mathcal C\rangle}
{(s-m_{j,0}^2)(s-m_{k,0}^2)(s-m_a^2)}
+O(\epsilon^4).
$$
This expression contains the ALP-like pole, including in the meson block. For the contribution with meson-like propagating poles only, expand the meson part of the spectral sum,
$$
\sum_n\frac{v^{(n)}v^{(n)T}}{s-m_n^2},
$$
using the preceding components. The mass shift differentiates the propagator with respect to its squared mass and produces a double denominator. The overlap correction varies both ends of the propagator. The ALP-like spectral term is counted separately in the one-ALP-line cycles. This avoids counting the ALP pole twice.

With one external ALP, the meson-only amplitude starts at order $\epsilon$. To obtain its mixing-dependent part at order $\epsilon^3$, expand the complete on-shell assembly through that order, including the tree vertices with higher ALP derivatives, the changes in the internal and external meson projectors, and the changes in on-shell kinematics. The same procedure applies to exchange, residue and vacuum terms. It uses the [projected cycles and on-shell assembly](alp-three-meson-one-loop-estimate.md#assembly-at-one-loop) before assigning loop sizes.

The [SymPy script](../../../scripts/derive_alp_three_meson_loop_power.py) expands the inverse block and solves the eigenvector equations symbolically. It verifies the propagator expansion and the meson pole and normalization conditions.

```sh
.venv/bin/python scripts/derive_alp_three_meson_loop_power.py
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_three_meson_loop_power
```
