---
title: Cancellations in the finite-mixing loop-size estimate
status: derived
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd, Manohar:1983md]
provenance: analytical-derivation
---

# Cancellations in the finite-mixing loop-size estimate

This record uses the class numbering, dressed invariants, and pole vectors of the [dressed-frame counting](alp-dressed-vertex-power-counting-invariants.md) and the class-dependent mass scales of the [naturalness estimate](alp-loop-naturalness-assumption.md). The comparison presupposes a nonzero [tree amplitude](alp-finite-mixing-three-meson-tree.md). The conditions below describe where estimating the correction from individual loop sizes can overestimate the coherent correction.

## Cancellation between classes

Classes 2, 3 and 5 have the same leading relative power, $\epsilon^2$. Include the effect sourced by class 6 at that order if it contributes to the chosen on-shell projection, with its treatment subject to the vacuum prescription tracked in `alp-woa`. At fixed kinematics the leading term cancels when the sum of the signed, mass-weighted contributions at that power vanishes. Each contribution includes its projected invariant contraction and the weight of the corresponding loop function. When absorptive parts are present, both the real and imaginary parts of the sum must vanish. Equal powers of $\epsilon$ permit the cancellation. Its occurrence for a physical amplitude depends on the weights and kinematics.

In the factors displayed in Result 4 of the counting record, classes 2, 3 and 6 carry $\epsilon^2/N(m_a^2)$, while the leading class-5 mixing insertion carries $\epsilon^2$. Here $m_a$ is the tree-level mass of the ALP-like pole and
$$
N(m_a^2)=1+2\epsilon^2\left\langle[\mathbb{M}_{m_a^2}^{-1}(\mathcal{C})]^2\right\rangle.
$$
A cancellation at a given finite $\epsilon$ requires the weighted sum at that point, including the pole factors and the dependence on the exact kinematics, to vanish. The meson-only term of class 1 and the two-ALP-line term of class 4 start at different powers and survive in general. Once the leading ALP-line coefficient cancels, the crossing estimate based on its uncancelled size does not apply.

## Zeros and cancellations in the invariant contractions

- **Quartic contraction in class 2.** For identical meson generators $\lambda$ commuting with $\mathbb{M}_{m_a^2}^{-1}(\mathcal{C})$, Result 3 gives a zero when
  $$
  \left\langle\lambda^2\left(\widehat S_2+4\mathbb{M}_{m_a^2}^{-1}(\mathcal{C})\widehat P_1-4[\mathbb{M}_{m_a^2}^{-1}(\mathcal{C})]^2m_q\right)\right\rangle=0.
  $$
  The condition is on the flavour projection. The full noncommuting contraction also contains commutator terms, and the zero of this quartic projection does not constrain the six-leg tadpole vertex.
- **Cubic contraction sourcing class 6.** The source of the meson one-point function vanishes in the projection $\lambda$ when
  $$
  \left\langle\lambda\left[B_0\left(\widehat P_2-2\{\mathbb{M}_{m_a^2}^{-1}(\mathcal{C}),\widehat S_1\}\right)-\frac{m_0^2}{N_f}\widehat\Theta_2\mathbf 1\right]\right\rangle=0.
  $$
  This relates the pseudoscalar, scalar and topological terms. The same zero removes any class-4 contribution containing that cubic vertex in the same projection.
- **Mixing factors in class 5.** Each displayed product of mixing factors vanishes if either meson projection of $\mathcal C$ vanishes. In particular, $B_0\widehat P_1=(m_0^2/N_f)\widehat\Theta_1\mathbf 1$ cancels the complete mixing matrix $\mathcal C$, eliminating the mixing-insertion terms. A direct ALP vertex can still contribute.

## Symbolic checks

The [SymPy script](../../../scripts/derive_alp_loop_cancellations.py) keeps the quark masses, pole mass, and chiral parameters symbolic in the three-flavour theory, $N_f=3$. The source relations below exhibit the cancellation conditions.

1. **Quartic.** The dressed mass $\widehat M(x)=m_q e^{ix}$ with $x=a/f_a$ gives $\widehat P_1=m_q$ and $\widehat S_2=-m_q$. Choose the first topological coefficient on the symbolic cancellation locus
   $$
   \widehat\Theta_1=\frac{N_f(m_a^2-m_0^2)}{2m_0^2}.
   $$
   The script applies $\mathbb{M}_{m_a^2}$ to $\mathbf 1/2$ and recovers $\mathcal{C}$. Where the operator is invertible, $\mathbb{M}_{m_a^2}^{-1}(\mathcal{C})=\mathbf 1/2$. Inserting the derivatives and this pole relation gives $-m_q+2m_q-m_q=0$ before projection.
2. **Cubic.** With the [general second-order mass sources](alp-dressed-mass-second-order-invariants.md) and the [second-order topological source](alp-dressed-topological-second-order.md), impose
   $$
   \widehat P_2=2\{\mathbb{M}_{m_a^2}^{-1}(\mathcal{C}),\widehat S_1\}+\frac{m_0^2}{N_f B_0}\widehat\Theta_2\mathbf 1.
   $$
   The cubic contraction then vanishes before projection, with $\widehat S_1$ and $\widehat\Theta_2$ still symbolic.
3. **Mixing.** Substitute $\widehat P_1=(m_0^2/(N_fB_0))\widehat\Theta_1\mathbf 1$ into $\mathcal{C}$. The two terms cancel for any $\widehat\Theta_1$, so every meson projection vanishes.
4. **Equal power.** Impose that the signed class-5 coefficient is the negative of the sum of the signed coefficients of classes 2, 3 and 6, after weighting them by their respective loop scales. The leading sum vanishes. With the displayed pole factors, the reduced sum is the class-5 coefficient times $1-1/N(m_a^2)$. It vanishes at leading order, and the residual in the relative loop term starts at $\epsilon^4$ at fixed pole mass and fixed weights. Establishing a physical cancellation at specified kinematics requires the corresponding loop coefficients.

```sh
.venv/bin/python scripts/derive_alp_loop_cancellations.py
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_loop_cancellations
```
