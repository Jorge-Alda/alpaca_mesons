---
title: Limits of the LO truncation of the ALP-meson spectrum
status: mixed
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: assumption
---

# Limits of the LO truncation of the ALP-meson spectrum

## What is exact and what is not

The [pole solver](alp-lo-pole-solver.md) diagonalizes the LO quadratic Lagrangian exactly, without expanding in $f_\pi/f_a$. This is not a resummation of the chiral expansion: the result describes the LO model, and it is a diagnostic of that model, not a prediction (assumption, see the finite-mixing record [finite-alp-meson-mixing.md](finite-alp-meson-mixing.md) and [@Alda:2026nrd]).

## Not included (deferred)

- NLO local terms. They enter $K$ and $M^2$ with relative size $m^2/\Lambda_\chi^2$, not suppressed by $f_\pi/f_a$; see the [NLO quadratic Lagrangian](../../nlo/quadratic-lagrangian.md) and [NLO neutral matrices](../../nlo/neutral-meson-matrices.md). The solver accepts general real symmetric $K>0$, $M^2$, so these blocks can be supplied without changing the solver.
- Loop self-energies and renormalized low-energy constants: see [loop corrections to ALP-meson mixing](../../nlo/loop-corrections-to-alp-meson-lsz-mixing.md).
- Resonance mixing: see [resonance saturation](../../nlo/resonance-saturation-and-dynamical-masses.md).
- Momentum-dependent $D(s)$ that is not linear in $s$, complex poles and widths. The solver finds only real poles of a polynomial $D(s)$ and does not describe an unstable state with a width.
- Charge-violating sources (the neutral reduction is rejected for them; see the [decoupling record](alp-lo-charged-meson-decoupling.md)).
- Decay amplitudes and numerical phenomenological inputs: the numbers in the records are illustrative.

## Heuristic validity statements (assumptions, not derived bounds)

The LO result should be flagged as outside the controlled domain in the following cases.

1. $|K_{ia}|\gtrsim1$, so that the kinetic mixing, which is proportional to $f_\pi/f_a$, is not small.
2. An ALP component $O(1)$ in a meson-like pole, or a meson component $O(1)$ in the ALP-like pole. The weak-mixing expansion of the [weak-mixing record](alp-lo-weak-mixing-expansion.md) then does not apply; the residual scales as $\epsilon^2$ to $\epsilon^4$ only when the mixing is small. At the illustrative large-mixing point of the pole solver, the ALP-like pole has $v_a\simeq0.64$ and meson components of order $1$.
3. A pole near an ALP-meson degeneracy ([two-level record](alp-lo-near-degeneracy-two-level.md)): the mixing angle is large however small the mixing entry is, and the label "ALP-like" is ambiguous ([tracking record](alp-lo-pole-tracking.md)).
4. A pole with $m_n^2$ approaching the chiral-expansion scale. The size of the scale is not fixed here; treating it as a number would be an assumption.
5. A near-vanishing of the vacuum-stability margin ([stability record](alp-lo-stability-conditions.md)).

Consequences: the sum $\sum_nv^{(n)}v^{(n)\mathsf T}=K^{-1}$ and the linear algebra remain valid in these cases, but whether the LO Lagrangian describes the physical system there is not decided by this analysis.

## Verification status

No script checks the validity statements; they are heuristic and labelled as assumptions. The statements about what the solver does and does not cover are verified by the interface of [`scripts/solve_alp_meson_spectrum.py`](../../../scripts/solve_alp_meson_spectrum.py) (general $K>0$, $M^2$; rejects non-positive $K$) and its tests.
