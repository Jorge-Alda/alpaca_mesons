---
title: Naturalness assumption for estimating the loop correction at finite ALP–meson mixing
status: conjectured
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd, Manohar:1983md]
provenance: assumption
---

# Naturalness assumption for the loop-size estimate

For an on-shell amplitude whose leading tree term is first order in $1/f_a$, take the magnitude of the neglected next-order correction to be of the order of the one-loop contribution estimated from the dressed-frame vertices. This is the naturalness assumption of design D3. It assumes that large logarithms, threshold enhancements and accidental cancellations do not control the net correction. A cancellation between contributions of the same order, including a cancellation due to relations among the invariants or between loop and local contributions, can make the correction smaller than this estimate. Use the [exact finite-mixing tree amplitude](alp-finite-mixing-three-meson-tree.md) as the reference for the comparison.

Use the dimensionless loop factor of [@Manohar:1983md] in the normalization of [@Alda:2026nrd],

$$
\frac{m^2}{16\pi^2 f_\pi^2}.
$$

Here $m$ is the mass scale for the loop class, as in the [loop power-counting conventions](alp-loop-power-counting.md). The factor multiplies the invariant contraction and the power from the ALP ends in [Result 4 of the dressed-frame counting](alp-dressed-vertex-power-counting-invariants.md#result-4-loop-classes-of-the-amplitude).

- For ALP-like tadpoles and the closed ALP-like line in a meson one-point function, take $m^2=m_A^2$, the exact tree-level ALP-like pole mass squared. For a meson tadpole use $m^2=m_P^2$, the squared mass of the meson in the loop, evaluated at the unmixed mass for leading meson-line counting.
- For bubbles and triangles, take $m^2$ to be the squared-mass scale carried by the vertices and the kinematics after the loop contribution is put in the same units as the tree amplitude. Use the relevant kinematic invariant $s$ where it sets that scale, or the mass scale from the vertices where it sets the scale. A dimensionless bubble counts as order unity after extracting the one-loop factor $1/(16\pi^2)$. For a triangle, the inverse squared-mass scaling from its three propagators must also be included when counting the net power. The presence of an ALP-like line alone does not set the vertex scale to $m_A^2$.
- For the meson loops with the additional mixing factors of class 5, use the meson scale of the underlying loop and keep those mixing factors in the contraction. The leading meson-only correction has the usual meson chiral scale.

The estimate keeps the class-dependent masses and invariant contractions symbolic. The [original power-counting record](alp-loop-power-counting.md) defines $m$ and the loop normalization. The [dressed-frame counting](alp-dressed-vertex-power-counting-invariants.md) defines $f_a$ and $f_\pi$ through its conventions and identifies $m_P$ and $s$ in Result 4. The [tree record](alp-finite-mixing-three-meson-tree.md) defines the pole $A$ and $m_A^2$. No new notation is introduced.
