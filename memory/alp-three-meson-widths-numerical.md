---
title: Numerical implementation of ALP three-meson decay widths
status: numerical
created: 2026-09-22
updated: 2026-09-24
bibtex_keys: [Alda:2026bar, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# Numerical ALP Three-Meson Widths

## Inputs and conventions

The implementation is [`scripts/implementation/three_meson_width.py`](../scripts/implementation/three_meson_width.py). It evaluates the amplitudes and Dalitz integrals derived in [ALP to three-meson amplitudes and decay widths](alp-three-meson-amplitudes-widths.md). The pseudoscalar matrix, generator normalization, $f_\pi=0.093\,\mathrm{GeV}$, fitted products $B_0m_q$, neutral physical generators, and meson masses follow Eqs. (III.1), (C.4), and (C.6) and Table I of `Alda:2026bar` [@Alda:2026bar].

The inputs $\widehat c_q^A$ and $\widehat c_q^P$ are Hermitian $3\times3$ arrays already expressed in the no-mixing basis. Since the fit supplies $B_0m_q$ directly, the symmetry-breaking matrix is evaluated as

$$
B_0\{\widehat c_q^P,m_q\}
=\{\widehat c_q^P,B_0m_q\}.
$$

Natural units are used. All masses and $f_a$ are in GeV, and `decay_width` returns a width in GeV. This leading-order calculation includes only the direct interactions recorded in the analytical derivation; it does not add higher-order chiral corrections or resonance form factors.

Amplitude construction is separated from phase-space integration through `build_amplitude`, which returns a callable `ThreeMesonAmplitude`. For a channel ordered as $(1,2,3)$, its arguments retain the convention $s_{12}=(p_1+p_2)^2$ and $s_{23}=(p_2+p_3)^2$, with $s_{13}=m_a^2+m_1^2+m_2^2+m_3^2-s_{12}-s_{23}$. This ordering is retained when repeated particles occur in positions $(1,3)$ or $(2,3)$.

## Numerical procedure

Definite-flavor channels that violate electric charge, have total $|S|>1$, or lie at or below threshold return zero before coupling validation or integration. A definite-flavor final state with total $|S|=1$ is allowed when at least one of $\widehat c_q^A$ and $\widehat c_q^P$ has a nonzero off-diagonal $d$-$s$ entry; otherwise it returns zero after coupling validation and before integration. Physical $K_S$ and $K_L$ states do not have definite strangeness, so their allowed flavor components are selected by the flavor traces instead of this early filter. Three-identical and two-identical channels use one-dimensional formulas in the repeated-particle pair invariant, with symmetry factors $3!$ and $2!$. Pairwise-distinct channels use nested adaptive quadrature over $s_{12}$ and $s_{23}$ with symmetry factor one. Small negative Kallen products caused by endpoint roundoff are clipped to zero.

The accepted canonical meson identifiers are `pi+`, `pi-`, `pi0`, `K+`, `K-`, `K0`, `K0bar`, `KS`, `KL`, `eta`, and `eta_prime`. The implementation uses

$$
K_S=pK^0+q\bar K^0,
\qquad
K_L=pK^0-q\bar K^0,
$$

where $p=(1+\epsilon_K)/\sqrt{2(1+|\epsilon_K|^2)}$ and $q=(1-\epsilon_K)/\sqrt{2(1+|\epsilon_K|^2)}$. The default is $\epsilon_K=2.228\times10^{-3}e^{i43.52^\circ}$, with provenance and qualifications in [the PDG extract](extracts/ParticleDataGroup:2024cfk/neutral-kaon-mixing.md) [@ParticleDataGroup:2024cfk]. The `epsilon_k` keyword permits an override, including `epsilon_k=0` for the CP-conserving limit. Coherent interference is included at amplitude level through the bra coefficients $p^*$ and $q^*$. A common neutral-kaon mass is used. Time-dependent neutral-kaon propagation is not included.

## Verification

Run

```bash
.venv/bin/python -m unittest scripts.implementation.test_three_meson_width
```

The deterministic tests cover forbidden-channel short circuits, the $d$-$s$ selection rule for unit-strangeness final states, matrix validation, vanishing couplings, quadratic coupling and inverse-$f_a$ scaling, input-order invariance for repeated and distinct states, an allowed net-zero-strangeness channel, one-neutral-kaon projection, and two-neutral-kaon interference with the required identical-particle factors.
