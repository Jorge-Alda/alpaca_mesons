---
title: Projection of neutral-kaon flavor amplitudes onto K_S and K_L
status: mixed
created: 2026-09-23
updated: 2026-09-23
bibtex_keys: [Alda:2026bar, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Projection of Neutral-Kaon Amplitudes

## Scope and inputs

Let $X$ denote all other particles and kinematic labels in an ALP decay. Define

$$
\mathcal A_0=\langle K^0X|T|a\rangle,
\qquad
\overline{\mathcal A}_0=\langle \bar K^0X|T|a\rangle.
$$

These flavor amplitudes may be obtained from the direct leading-order interactions recorded in [the three-meson amplitude derivation](../alp/three-meson/alp-three-meson-amplitudes-widths.md), whose chiral input is documented in [the matrix-level derivation](../alp/three-meson/alp-phi3-lo-chipt-derivation.md) [@Alda:2026bar]. The mixing input is recorded in [the PDG neutral-kaon extract](../extracts/ParticleDataGroup:2024cfk/neutral-kaon-mixing.md) [@ParticleDataGroup:2024cfk]. The executable derivation is [`scripts/derive_neutral_kaon_projection.py`](../../scripts/derive_neutral_kaon_projection.py).

Assumptions and conventions:

- The neutral-kaon mixing convention is
  $$
  |K_S\rangle=p|K^0\rangle+q|\bar K^0\rangle,
  \qquad
  |K_L\rangle=p|K^0\rangle-q|\bar K^0\rangle.
  $$
- The phase convention for $K^0$ and $\bar K^0$ must be used consistently in the states, fields, and flavor couplings.
- Neglecting CPT violation in mixing,
  $$
  p=\frac{1+\epsilon_K}{\sqrt{2(1+|\epsilon_K|^2)}},
  \qquad
  q=\frac{1-\epsilon_K}{\sqrt{2(1+|\epsilon_K|^2)}}.
  $$
- The formulas describe production amplitudes. A measured channel involving the decay products of an unstable neutral kaon can additionally require coherent $K^0$--$\bar K^0$ propagation, decay amplitudes, detector acceptance, and decay-time cuts.

## Physical-state projection

Because amplitudes contain the final-state bra, the mixing coefficients are complex conjugated:

$$
\boxed{
\mathcal A_S=p^*\mathcal A_0+q^*\overline{\mathcal A}_0,
\qquad
\mathcal A_L=p^*\mathcal A_0-q^*\overline{\mathcal A}_0.
}
$$

Thus an exclusive prediction for a tagged $K_S$ or $K_L$ should use $\mathcal A_S$ or $\mathcal A_L$, including interference before taking the absolute square. Computing in the flavor basis first is nevertheless fully equivalent and is often more transparent. The states are not exact CP eigenstates and are not orthogonal:

$$
\langle K_S|K_L\rangle
=|p|^2-|q|^2
=\frac{2\operatorname{Re}\epsilon_K}{1+|\epsilon_K|^2}.
$$

Only in the explicit CP-conserving limit $\epsilon_K=0$, where $p=q=1/\sqrt2$, do the projections reduce to

$$
\boxed{
\mathcal A_S=\frac{\mathcal A_0+\overline{\mathcal A}_0}{\sqrt2},
\qquad
\mathcal A_L=\frac{\mathcal A_0-\overline{\mathcal A}_0}{\sqrt2}.
}
$$

Directly replacing the flavor fields by their $K_S,K_L$ combinations in the Lagrangian gives the same result; it is a change of basis rather than a way to discard one flavor structure. If more than one neutral kaon occurs, apply the mixing transformation to every neutral-kaon leg, retaining all flavor combinations that can contribute.

## What hermiticity does and does not imply

Hermiticity of a flavor-coupling matrix gives, for example,

$$
c_{ds}=c_{sd}^*.
$$

This relates coefficients of conjugate operators. **Hermiticity alone does not imply**

$$
\overline{\mathcal A}_0=\mathcal A_0^*.
$$

Indeed, complex conjugation of $\langle K^0X|T|a\rangle$ produces the matrix element of the reverse transition, not automatically $\langle\bar K^0X|T|a\rangle$. Moreover, CP-conjugate decay amplitudes have conjugated weak couplings but the same final-state strong phases, whereas full complex conjugation would also reverse those phases.

A useful restricted case occurs when all of the following hold:

1. $a$ and the remainder $X$ map back to the channel under the relevant CP/CPT transformation, with phase conventions fixed;
2. the hadronic dynamics used to relate the two flavor channels is CP conserving;
3. each conjugate operator pair has the same hadronic/kinematic factor $F_r$.

Then one may write

$$
\mathcal A_0=\sum_r F_r c_r,
\qquad
\overline{\mathcal A}_0=\eta\sum_r F_r c_r^*,
$$

where the convention-dependent phase $\eta$ contains the CP phases of the external states. This relation lets one calculate the common factors $F_r$ once and obtain both flavor amplitudes from the conjugate couplings. It still is not generally $\overline{\mathcal A}_0=\mathcal A_0^*$ when the $F_r$ contain absorptive or final-state-interaction phases.

For the especially restrictive single-structure tree-level case with $\eta=1$ and real $F$,

$$
\mathcal A_0=Fc,
\qquad
\overline{\mathcal A}_0=Fc^*,
$$

so the CP-limit projections reduce to

$$
\boxed{
\mathcal A_S=\sqrt2F\operatorname{Re}c,
\qquad
\mathcal A_L=i\sqrt2F\operatorname{Im}c.
}
$$

This last shortcut is an additional tree-level/real-matrix-element assumption, not a consequence of hermiticity by itself. If $X$ is not self-conjugate, the conjugation relation instead connects $a\to K^0X$ to the channel with $\bar X$, so it cannot generally supply the missing amplitude for the same final state.

### Specialization to the implemented LO contact amplitudes

For the direct tree-level contact terms of the real LO chiral Lagrangian used here, the kinematic polynomials, masses, neutral-meson generators, and $B_0m_q$ inputs are real in the adopted phase convention. There are no loop absorptive parts or inserted final-state-interaction phases. Hermiticity therefore supplies conjugate $d$--$s$ coefficients with common real kinematic factors. The three conditions above hold for channels that are self-conjugate and hold **pairwise** for channels exchanged by charge conjugation. For example, $\pi^+K^-K^0$ is related to $\pi^-K^+\bar K^0$ rather than to itself.

This conclusion ceases to follow from the contact-term argument after complex form factors or final-state-interaction functions are inserted. The physical-state projection remains valid, but the dressed conjugate amplitudes must then be related using the stated strong-interaction assumptions rather than coupling hermiticity alone.

Electric charge conservation is independent of this argument. The intended flavour-violating channel is $\pi^0\eta'K^0$; the previously stated $\pi^-\eta'K^0$ was a typo and would have total charge $-1$.

## Practical conclusion

For a physical $K_S$ or $K_L$ rate, project amplitudes and then square. One can avoid two independent chiral calculations when the explicit operator structure shows that the two flavor amplitudes differ only by hermitian-conjugate couplings and known external-state phases. One cannot eliminate one amplitude by invoking coupling hermiticity alone.

## Verification

Run

```bash
.venv/bin/python scripts/derive_neutral_kaon_projection.py
```

The script constructs the mixing matrix, derives the bra-level projection and inverse transformation, takes the CP-conserving limit, and derives the real common-factor special case. Assertions verify each relation.
