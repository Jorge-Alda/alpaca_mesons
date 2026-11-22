---
title: Crossing condition with the exact finite-mixing tree amplitude
status: derived
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd, Manohar:1983md]
provenance: analytical-derivation
---

# Crossing condition with the exact tree amplitude

Use a nondegenerate ALP-like pole of mass $m_a$ and a kinematic point where the tree amplitude is nonzero. Vary $\epsilon=f_\pi/f_a$ at fixed $m_a$, with the bare mass adjusted by the pole equation. The crossing criterion is $R_{\max}>0$, the order-one value of the estimated loop-to-tree magnitude ratio chosen in [design decision D6](../../../openspec/changes/estimate-loop-tree-crossing-finite-mixing/design.md#decisions).

## One internal ALP-like line

The external ALP-like pole supplies $\epsilon/\sqrt{N(m_a^2)}$ to the tree and to the estimated loop contribution. After extracting this common external factor, a single internal ALP-like line supplies a further $\epsilon^2/N(m_a^2)$ relative to the tree. Here
$$
N(m_a^2)=1+2\epsilon^2\left\langle[\mathbb{M}_{m_a^2}^{-1}(\mathcal C)]^2\right\rangle .
$$

- **Reduced exact tree amplitude** is the [pole-projected contact-plus-exchange amplitude](alp-finite-mixing-three-meson-tree.md#result) divided by the external factor $\epsilon/\sqrt{N(m_a^2)}$. Its limit uses the [leading tree coefficient](alp-finite-mixing-three-meson-tree.md#difference-from-the-leading-expansion).
- **Reduced invariant loop contraction** is the class-dependent contraction after removing the external factor, the internal factor $\epsilon^2/N(m_a^2)$, and the loop scale. For the contraction of a closed line, the [projected quartic and cubic expressions](alp-dressed-vertex-power-counting-invariants.md#result-3-contraction-of-a-closed-alp-like-line-verified) give the vertex factors. The [class table](alp-dressed-vertex-power-counting-invariants.md#result-4-loop-classes-of-the-amplitude) specifies how the other contractions are built from the projected vertices.

Write $\mathcal M^{(0)}$ for the linked exact tree amplitude and $\mathcal M^{(1)}$ for the power-counting estimate of the selected one-loop contribution. For $\epsilon>0$, the two reductions defined above are, respectively,
$$
\frac{\sqrt{N(m_a^2)}}{\epsilon}\mathcal M^{(0)}(\epsilon),
\qquad
\frac{16\pi^2 f_\pi^2}{m^2}
\frac{N(m_a^2)^{3/2}}{\epsilon^3}\mathcal M^{(1)}(\epsilon).
$$
The second reduction uses the [class-dependent loop factor](alp-loop-naturalness-assumption.md). The reduced coefficient is the magnitude ratio of these two reductions (loop over tree), multiplied by that loop factor. In terms of the loop-to-tree magnitude ratio $R(\epsilon)$, the crossing equation and leading power are
$$
R(\epsilon)=
\frac{\left|\mathcal M^{(1)}(\epsilon)\right|}
     {\left|\mathcal M^{(0)}(\epsilon)\right|},
\qquad
R(\epsilon)=R_{\max},
\qquad
R(\epsilon)=O(\epsilon^2).
$$

The class-5 meson-line mixing insertions have the same leading relative power but their counting uses meson pole overlaps. The criterion here applies to a contribution with one internal ALP-like pole. With several contributions, the crossing is determined by their combined size under the chosen estimate, subject to the [cancellation conditions](alp-loop-estimate-cancellations.md).

## Dependence on the crossing value of $\epsilon$

The reduced coefficient still depends on $\epsilon$. Even at fixed $m_a$, the exact daughter masses and overlaps depend on $\epsilon$. Also, eliminating the bare ALP mass in the exchange denominator using $D_a(m_a^2)=0$ gives
$$
D_a(q^2)=q^2-m_a^2+
2\epsilon^2\left\langle
\mathcal C\left[\mathbb M_{q^2}^{-1}(\mathcal C)-
\mathbb{M}_{m_a^2}^{-1}(\mathcal C)\right]\right\rangle.
$$
The resolvent difference is generically nonzero. The extra $\epsilon^2$ term remains after the external ALP pole factor has been divided out of the exchange amplitude. The exact contact terms also contain higher powers, for example through the two-ALP-leg vertex and an ALP overlap on a daughter meson. No symmetry in the stated source assumptions forces these effects to cancel against the loop contraction.

For an explicit estimate, take the limit of the reduced coefficient as $\epsilon\to0$. Since $N(m_a^2)\to1$ in that limit, the leading reduced coefficient is $\lim_{\epsilon\to0^+}R(\epsilon)/\epsilon^2$. Holding it fixed while retaining the internal ALP pole normalization gives the positive branch
$$
\frac{f_\pi}{f_a}
=\left[
\frac{R_{\max}}
{\displaystyle
\lim_{\epsilon\to0^+}\frac{R(\epsilon)}{\epsilon^2}
-2R_{\max}\left\langle[\mathbb{M}_{m_a^2}^{-1}(\mathcal C)]^2\right\rangle}
\right]^{1/2}.
$$
This branch requires a strictly positive denominator. With the exact tree denominator and the full dependence of the reduced coefficient, the crossing condition is implicit. Selecting a channel, kinematics and the leading loop classes would specify the equation to solve.

## Symbolic check

The [SymPy derivation](../../../scripts/derive_alp_loop_crossing_criterion.py) cancels the external factor, expands the remaining ratio in $\epsilon$, solves the condition with the reduced coefficient fixed, and substitutes the root into that condition. On a smooth branch away from zeros and kinematic singularities, the reduced coefficient is a function of $\epsilon^2$, giving a leading $\epsilon^2$ term and corrections starting at $\epsilon^4$. The script then eliminates the bare mass from the inverse propagator using the pole equation and checks the resulting resolvent difference. The projected loop coefficient stays unevaluated. The symbolic tests verify the algebra of the criterion and the pole subtraction.

```sh
.venv/bin/python scripts/derive_alp_loop_crossing_criterion.py
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_loop_crossing_criterion
```
