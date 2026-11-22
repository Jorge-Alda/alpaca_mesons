---
title: Crossing condition with the exact finite-mixing tree amplitude
status: derived
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd, Manohar:1983md]
provenance: analytical-derivation
---

# Crossing condition with the exact tree amplitude

## Accepted crossing equation

For $a\to P_{n_1}P_{n_2}P_{n_3}$, compare the leading mixing-dependent loop-size estimate $\mathcal M^{(1)}$ with the exact finite-mixing tree amplitude $\mathcal M^{(0)}$. At fixed ALP-like pole mass $m_a$ and two independent on-shell kinematic invariants, vary $\epsilon=f_\pi/f_a$ and adjust the bare mass by the pole equation. Use a nondegenerate branch with a nonzero leading tree coefficient and a nonzero exact tree amplitude at the crossing. The crossing satisfies

$$
\epsilon^3
\lim_{\epsilon'\to0^+}
\frac{\left|\mathcal M^{(1)}(\epsilon')\right|}{(\epsilon')^3}
=R_{\max}\left|\mathcal M^{(0)}(\epsilon)\right|.
$$

This equation truncates the loop numerator at its leading mixing-dependent order and keeps the tree denominator exact, including the daughter poles and exchange propagators. The power counting leaves the number of solutions on a given branch open. The criterion $R_{\max}>0$ is a general threshold of order one for the estimated loop-to-tree magnitude ratio. Setting $R_{\max}=1$ compares equal sizes. [Design decision D6](../../../openspec/changes/estimate-loop-tree-crossing-finite-mixing/design.md#decisions).

## Naturalness assumption and contributing classes

The estimate treats the assembled one-loop magnitude as the natural size of the next-order mixing-dependent correction. Combine coefficients of equal scalar moments before taking magnitudes, then add the estimated magnitudes of the distinct moments. This assumes their signed sum has the same order of magnitude. The estimate applies away from large logarithms, threshold enhancements and substantial cancellations with local terms. The reference mass $m$ in $m^2/(16\pi^2 f_\pi^2)$ is class dependent. Tadpoles with an ALP-like pole use $m_a$. For bubbles and triangles use the scale set by the vertices or kinematics after momentum reduction. Meson-loop insertions use the underlying meson scale. [Naturalness prescription](alp-loop-naturalness-assumption.md), [scalar-moment weights](alp-loop-momentum-reduction.md#scalar-moments-and-derivatives), [cancellation conditions](alp-loop-estimate-cancellations.md).

| Loop class | Leading power relative to tree | Role in the crossing |
|---|---:|---|
| 1. Meson-only loops at leading order | $\epsilon^0$ | Ordinary chiral correction, omitted from the mixing-dependent numerator. |
| 2. ALP-like tadpoles | $\epsilon^2$ | Included, with the direct six-leg term and corrections to exchange and pole data. |
| 3. Bubbles and triangles with one ALP-like propagator | $\epsilon^2$ | Included, using the projected kinetic and potential vertices. |
| 4. Loops with two ALP-like propagators | $\epsilon^4$ | Beyond the leading numerator. |
| 5. Meson loops with extra ALP–meson mixing | $\epsilon^2$ | Included through the expansion of the meson-like poles and overlaps. The explicit ALP-like spectral term is counted in its own class. |
| 6. One-point source from an ALP-like tadpole | $\epsilon^2$ | Included through the linearized vacuum displacement when the source is nonzero. |

The powers are relative to a tree term linear in $\epsilon$. The numerator includes the corresponding exchange, external-pole and residue corrections. [Class definitions](alp-dressed-vertex-power-counting-invariants.md#result-4-loop-classes-of-the-amplitude), [on-shell assembly and cubic term](alp-three-meson-one-loop-estimate.md#extracting-and-estimating-the-leading-coefficient).

## One internal ALP-like line

The external ALP-like pole supplies $\epsilon/\sqrt{N(m_a^2)}$ to the tree and to the estimated loop contribution. After extracting this common external factor, a single internal ALP-like line supplies a further $\epsilon^2/N(m_a^2)$ relative to the tree. Here
$$
N(m_a^2)=1+2\epsilon^2\left\langle[\mathbb{M}_{m_a^2}^{-1}(\mathcal C)]^2\right\rangle .
$$

- **Reduced exact tree amplitude** is the [pole-projected contact-plus-exchange amplitude](alp-finite-mixing-three-meson-tree.md#result) divided by the external factor $\epsilon/\sqrt{N(m_a^2)}$. Its limit uses the [leading tree coefficient](alp-finite-mixing-three-meson-tree.md#difference-from-the-leading-expansion).
- **Reduced invariant loop contraction** is the class-dependent contraction after removing the external factor, the internal factor $\epsilon^2/N(m_a^2)$, and the loop scale. The [three-meson loop estimate](alp-three-meson-one-loop-estimate.md#extracting-and-estimating-the-leading-coefficient) specifies the projected one-loop amplitude and its leading coefficient. The [class table](alp-dressed-vertex-power-counting-invariants.md#result-4-loop-classes-of-the-amplitude) gives the underlying vertex counting.

In this subsection take only the contribution with one internal ALP-like pole in $\mathcal M^{(1)}$. For $\epsilon>0$, the two reductions defined above are, respectively,
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

The class-5 meson-line mixing insertions have the same leading relative power but their counting uses meson pole overlaps. The internal normalization factored out in the second reduction applies to a contribution with one internal ALP-like pole. For the combined leading estimate, extract the cubic coefficient of the contributions before this normalization is retained beyond leading order. The [cancellation conditions](alp-loop-estimate-cancellations.md) delimit when the separate sizes can be combined.

## Dependence on the crossing value of $\epsilon$

The reduced coefficient still depends on $\epsilon$. Even at fixed $m_a$, the exact daughter masses and overlaps depend on $\epsilon$. Also, eliminating the bare ALP mass in the exchange denominator using $D_a(m_a^2)=0$ gives
$$
D_a(q^2)=q^2-m_a^2+
2\epsilon^2\left\langle
\mathcal C\left[\mathbb M_{q^2}^{-1}(\mathcal C)-
\mathbb{M}_{m_a^2}^{-1}(\mathcal C)\right]\right\rangle.
$$
The resolvent difference is generically nonzero. The extra $\epsilon^2$ term remains after the external ALP pole factor has been divided out of the exchange amplitude. The exact contact terms also contain higher powers, for example through the two-ALP-leg vertex and an ALP overlap on a daughter meson. No symmetry in the stated source assumptions forces these effects to cancel against the loop contraction.

### Leading-order inversions

For an explicit estimate with one internal ALP-like pole, take the limit of its reduced coefficient as $\epsilon\to0$. Since $N(m_a^2)\to1$ in that limit, the leading reduced coefficient is $\lim_{\epsilon\to0^+}R(\epsilon)/\epsilon^2$ for that contribution. Holding it fixed while retaining the internal ALP pole normalization gives the positive branch
$$
\frac{f_\pi}{f_a}
=\left[
\frac{R_{\max}}
{\displaystyle
\lim_{\epsilon\to0^+}\frac{R(\epsilon)}{\epsilon^2}
-2R_{\max}\left\langle[\mathbb{M}_{m_a^2}^{-1}(\mathcal C)]^2\right\rangle}
\right]^{1/2}.
$$
This branch applies to the estimate with a single explicit internal ALP-like line and requires a strictly positive denominator. For the combined estimate truncated at strict leading order on both sides, the positive branch is instead
$$
\frac{f_\pi}{f_a}
=\left[
\frac{R_{\max}}
{\displaystyle\lim_{\epsilon\to0^+}R(\epsilon)/\epsilon^2}
\right]^{1/2}.
$$
With the exact tree denominator, the crossing condition is implicit. Selecting the daughter channel and fixing its two independent invariants specifies the remaining root equation.

## Notation and inputs

- $a$ and $P_{n_i}$ denote the ALP-like pole and the three meson poles. The tree record defines their masses, $\epsilon=f_\pi/f_a$, the bare mass condition, the fixed kinematic invariants, $q^2=(p_{n_j}+p_{n_k})^2$, and the exchange denominator $D_a(q^2)$. The prime on $\epsilon'$ denotes only the variable in the limit. [Exact pole projections and kinematics](alp-finite-mixing-three-meson-tree.md#result), [fixed-mass expansion](alp-finite-mixing-three-meson-tree.md#difference-from-the-leading-expansion).
- The same tree record defines $\mathbb M_s$, its inverse at $s=m_a^2$ or $q^2$, the pole norm $N(s)$, and $\mathcal C=B_0\widehat P_1-(m_0^2/N_f)\widehat\Theta_1\mathbf 1$. The dressed vertex record defines the invariant source derivatives, $B_0$, $m_0^2$, the flavour trace $\langle\cdot\rangle$, $N_f=3$ and the identity $\mathbf 1$. The unmixed meson masses enter through $\mathbb M_s$. [Pole operator](alp-lo-exact-pole-equations.md), [dressed conventions](alp-dressed-vertex-power-counting-invariants.md#conventions-and-units).
- $\mathcal M^{(0)}$, $\mathcal M^{(1)}$ and $R(\epsilon)$ retain the meanings in the crossing derivation. The scales $m$, $m_a$ and $f_\pi$ have mass dimension one, and $R_{\max}$, $\epsilon$, and $N(s)$ are dimensionless. [Projected loop estimate](alp-three-meson-one-loop-estimate.md), [loop factor](alp-loop-naturalness-assumption.md).

## Symbolic check

The [SymPy derivation](../../../scripts/derive_alp_loop_crossing_criterion.py) cancels the external factor, expands the remaining ratio in $\epsilon$, solves the condition with the reduced coefficient fixed, and substitutes the root into that condition. On a smooth branch away from zeros and kinematic singularities, the reduced coefficient is a function of $\epsilon^2$, giving a leading $\epsilon^2$ term and corrections starting at $\epsilon^4$. The script then eliminates the bare mass from the inverse propagator using the pole equation and checks the resulting resolvent difference. The projected loop coefficient stays unevaluated. The symbolic tests verify the algebra of the criterion and the pole subtraction.

```sh
.venv/bin/python scripts/derive_alp_loop_crossing_criterion.py
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_loop_crossing_criterion
```
