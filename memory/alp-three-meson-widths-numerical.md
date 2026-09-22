---
title: Numerical implementation of ALP three-meson decay widths
status: numerical
created: 2026-09-22
updated: 2026-09-22
bibtex_keys: [Alda:2026bar]
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

## Numerical procedure

Channels that violate electric charge, have total $|S|>1$, or lie at or below threshold return zero before coupling validation or integration. A final state with total $|S|=1$ is allowed when at least one of $\widehat c_q^A$ and $\widehat c_q^P$ has a nonzero off-diagonal $d$-$s$ entry; otherwise it returns zero after coupling validation and before integration. Three-identical and two-identical channels use the one-dimensional $s_{12}$ formulas with symmetry factors $3!$ and $2!$. Pairwise-distinct channels use nested adaptive quadrature over $s_{12}$ and $s_{23}$ with symmetry factor one. Small negative Kallen products caused by endpoint roundoff are clipped to zero.

The accepted canonical meson identifiers are `pi+`, `pi-`, `pi0`, `K+`, `K-`, `K0`, `K0bar`, `eta`, and `eta_prime`.

## Verification

Run

```bash
.venv/bin/python -m unittest scripts.implementation.test_three_meson_width
```

The deterministic tests cover forbidden-channel short circuits, the $d$-$s$ selection rule for unit-strangeness final states, matrix validation, vanishing couplings, quadratic coupling and inverse-$f_a$ scaling, input-order invariance for repeated and distinct states, and an allowed net-zero-strangeness channel.
