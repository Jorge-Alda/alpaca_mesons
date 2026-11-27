---
title: Validity and limiting cases of the finite-mixing crossing estimate
status: mixed
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd, Manohar:1983md]
provenance: analytical-derivation
---

# Validity and limiting cases of the crossing estimate

Use the [implicit crossing equation](alp-loop-crossing-exact-tree-condition.md#accepted-crossing-equation) with the exact LO tree amplitude and the leading mixing-dependent loop-size numerator. The pole and source conventions are those of the [dressed pole equations](../../../../meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-exact-pole-equations.md) and the [loop estimate](../../three-meson-amplitudes/alp-three-meson-one-loop-estimate.md).

## Domain

- Use stable, isolated LO poles, a smooth on-shell kinematic branch and a nonzero leading tree coefficient. The exact tree amplitude must also be nonzero at the comparison point. Stay away from thresholds and singular internal exchange propagators, with masses and momenta below the chiral scale. [Naturalness assumption](../loop-counting/alp-loop-naturalness-assumption.md), [LO truncation limits](../../../../meson-mixing/finite-mixing/pole-limits/alp-lo-truncation-limits.md).
- The exact tree denominator uses the full LO pole solution. The numerator uses a nondegenerate expansion. In particular, for each coupled meson its leading ALP overlap is small when
  $$
  \frac{\epsilon|\langle\lambda_n\mathcal C\rangle|}
       {|m_a^2-m_{n,0}^2|}\ll1 .
  $$
  The meson pole shifts and induced rotations must also be small compared with the relevant meson mass-squared gaps in the [meson mixing expansion](../../three-meson-amplitudes/alp-three-meson-loop-mixing-expansion.md). Within an exactly degenerate unmixed subspace one must first resolve the mixing with the ALP. The condition applies only to directions with nonzero projected mixing.
- The loop estimate retains the leading mixing-dependent one-loop term. Perturbative external-pole and residue corrections require the appropriate matrix elements of the [loop self-energy](../../../../../nlo/loop-corrections-to-alp-meson-lsz-mixing.md) to be small relative to the pole separations, and its projected derivative to be small relative to the tree pole norm. For the amplitude, $R(\epsilon)\ll1$ is a smallness test under the naturalness assumption, alongside ordinary meson chiral counting. A crossing with $R_{\max}$ of order one locates the loss of that estimated hierarchy. The LO pole solver checks the spectral conditions and does not test the natural size of unknown higher-order terms.

## Decoupling at fixed masses

As $f_a\to\infty$ at fixed $f_\pi$, pole mass $m_a$, dressed invariants and independent kinematic invariants, $\epsilon\to0$. For fixed nonzero mass gaps, the off-diagonal mass entries are linear in $\epsilon$, the shifts of the meson poles start at $\epsilon^2$, and $N(m_a^2)\to1$. The ALP-like pole approaches a pure ALP state, while the mesons approach their unmixed states. The [exact tree amplitude](../../../../decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-finite-mixing-three-meson-tree.md#difference-from-the-leading-expansion) starts at order $\epsilon$. The mixing-dependent loop numerator starts at order $\epsilon^3$, giving
$$
\lim_{\epsilon\to0^+}R(\epsilon)=0,\qquad
R(\epsilon)=O(\epsilon^2).
$$
The ordinary meson-loop contribution to this decay also vanishes in absolute size, starting at order $\epsilon$. Its ratio to the tree can approach a finite chiral correction because both scale linearly in $\epsilon$. This contribution is excluded from the crossing numerator.

## Near an ALP–meson degeneracy

If the separation of the ALP diagonal level and a coupled meson level is comparable to $\epsilon|\langle\lambda_n\mathcal C\rangle|$, their overlap can be of order one. Powers inferred by expanding the two projectors at fixed gaps are then nonuniform. A pole normalization factor by itself cannot restore the expansion used for the loop numerator.

For an isolated pair the exact replacement for the two tree pole masses is the [canonical two-level result](../../../../meson-mixing/finite-mixing/pole-limits/alp-lo-near-degeneracy-two-level.md),
$$
m_\pm^2=
\frac{m_{n,0}^2+M^2_{\hat a\hat a,\rm can}}{2}
\pm\frac12\sqrt{
(M^2_{\hat a\hat a,\rm can}-m_{n,0}^2)^2
+4\epsilon^2\langle\lambda_n\mathcal C\rangle^2}.
$$
Here $M^2_{\hat a\hat a,\rm can}$ is the canonical diagonal squared-mass entry, as defined in the [pole equations](../../../../meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-exact-pole-equations.md). Diagonal equality refers to this entry and the unmixed meson mass, prior to diagonalization. For the canonical $2\times2$ squared-mass matrix $M^2$, the exact projectors and their ALP weights are
$$
v^{(\pm)}v^{(\pm)T}
=\frac{M^2-m_\mp^2\mathbf 1}{m_\pm^2-m_\mp^2},
$$
$$
|v_a^{(\pm)}|^2=
\frac12\left(
1\pm\frac{M^2_{\hat a\hat a,\rm can}-m_{n,0}^2}
{\sqrt{(M^2_{\hat a\hat a,\rm can}-m_{n,0}^2)^2
+4\epsilon^2\langle\lambda_n\mathcal C\rangle^2}}
\right).
$$
At equality of the diagonal entries and nonzero mixing, each pole has ALP weight $1/2$ and the mass-squared splitting is $2\epsilon|\langle\lambda_n\mathcal C\rangle|$. Use both pole projectors there. If the projected mixing vanishes, the isolated levels instead cross. At exact coincidence the two-state spectral projector is well defined, while an individual eigenbasis can be chosen within that subspace.

With other coupled mesons present the exact replacement in the pair is the energy-dependent Schur equation
$$
(s-m_{n,0}^2)[s-\widetilde M_n(s)]
=\epsilon^2\langle\lambda_n\mathcal C\rangle^2,
$$
where $s$ is the squared momentum and
$$
\widetilde M_n(s)
=M^2_{\hat a\hat a,\rm can}
+\epsilon^2\sum_{n'\ne n}
\frac{\langle\lambda_{n'}\mathcal C\rangle^2}{s-m_{n',0}^2}.
$$
The [isolated-pole record](../../../../meson-mixing/finite-mixing/pole-limits/alp-lo-isolated-meson-pole.md#step-3-error-bound-and-validity) bounds the approximation obtained by freezing this function. If several nearby mesons couple, use the full pole problem. These replacements treat the LO poles. A crossing estimate in the near-degenerate regime would require loop counting with the mixed projectors and, for comparable loop shifts and level spacings, a joint pole treatment.

## Symbolic verification

The [limits script](../../../../../../scripts/derive_alp_loop_crossing_limits.py) uses the exact solver's inverse-propagator construction to derive the two eigenvalues and projectors. It verifies their null-space, idempotence and completeness identities, the second-order fixed-gap mass shifts, and the decoupling of the ALP weights. The check uses the exact matrix identities because the spectrum solver's final sorting step requires numerical ordering.

The script checks the joint limit by substituting $M^2_{\hat a\hat a,\rm can}-m_{n,0}^2=\epsilon\langle\lambda_n\mathcal C\rangle$ with $\langle\lambda_n\mathcal C\rangle>0$ into the [two-level splitting and projector formulas](#near-an-alpmeson-degeneracy). Along this path the diagonal gap decreases in proportion to the mixing entry. The splitting vanishes linearly in $\epsilon$, while the ALP weights remain finite and depend on the ratio of the diagonal gap to the mixing entry. With a spectator meson, the Schur equation agrees exactly with the characteristic determinant away from the spectator's unmixed pole. All masses and mixing sources remain symbolic. The checks establish the spectral and power-counting limits. The naturalness hypothesis and the smallness of uncomputed higher-order amplitudes are separate inputs.

```sh
.venv/bin/python scripts/derive_alp_loop_crossing_limits.py
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_loop_crossing_limits
```
