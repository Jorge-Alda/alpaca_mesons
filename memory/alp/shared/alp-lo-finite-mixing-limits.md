---
title: Validity and limits of isolated LO finite-mixing pole approximations
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Validity and limits of isolated LO finite-mixing pole approximations

## Inputs and assumptions

Use the exact canonical arrow spectrum of the [arrow-matrix record](alp-lo-arrow-matrix-spectrum.md), the isolated meson approximation of the [isolated-meson record](alp-lo-isolated-meson-pole.md), and the isolated ALP approximation of the [isolated-ALP record](alp-lo-isolated-alp-pole.md), in Alda:2026nrd conventions [@Alda:2026nrd]. The approximations assume the pole is separated from the other unmixed meson masses. They do not require small $f_\pi/f_a$ by itself.

## Spacing condition and level repulsion

For an isolated meson $n$, the relevant small ratios are
$$
\frac{\epsilon^2\operatorname{tr}(\lambda_{n'}\mathcal C)^2}{(m_{n,0}^2-m_{n',0}^2)^2}\ll1\quad(n'\ne n),\qquad |s-m_{n,0}^2|\ll|m_{n,0}^2-m_{n',0}^2|,
$$
where $\epsilon=f_\pi/f_a$. The exact bound and its conditions are given in the isolated-meson record. For an isolated ALP pole, the frozen closed form assumes
$$
q=\epsilon^2\sum_n\frac{\operatorname{tr}(\lambda_n\mathcal C)^2}{(M^2_{\hat a\hat a,\rm can}-m_{n,0}^2)^2}\ll1,
$$
as derived in the isolated-ALP record. These are approximation conditions; the arrow eigenvalues remain the exact result when either condition fails.

For two levels of unmixed squared masses $m^2$ and $M_a^2$ coupled by $g=\epsilon\operatorname{tr}(\lambda_n\mathcal C)$, the exact eigenvalues are
$$
s_\pm=\frac{m^2+M_a^2}{2}\pm\sqrt{\frac{(M_a^2-m^2)^2}{4}+g^2}.
$$
The meson-like level moves away from the ALP level: it shifts downward when $m^2<M_a^2$ and upward when $m^2>M_a^2$. At degeneracy the splitting is $2|g|$. This formula and the criteria for treating close poles together are recorded in [the near-degeneracy record](alp-lo-near-degeneracy-two-level.md). For several levels, the exact interlacing and secular equation in the arrow-matrix record apply.

## Decoupling and weak-mixing regime

As $f_a\to\infty$, $\epsilon\to0$, the canonical mixed entries vanish, $M^2_{\hat a\hat a,\rm can}\to m_{a,0}^2$, and the exact eigenvalues return to the unmixed meson masses $m_{n,0}^2$ and the unmixed ALP mass $m_{a,0}^2$. The exact finite-mixing equations in the [pole-equation record](alp-lo-exact-pole-equations.md) make this limit explicit.

At small mixing the leading displacement of a non-degenerate pole is quadratic in the off-diagonal entry, hence of order $(f_\pi/f_a)^2$; there is no linear mass shift in $f_\pi/f_a$. This follows from the exact secular equation or ordinary non-degenerate eigenvalue perturbation theory and is illustrated in the [second-order shift record](alp-lo-meson-mass-shift-second-order.md). At larger mixing, truncated shifts need not be accurate: use the exact poles as eigenvalues of the canonical mass matrix. All equations here are tree-level LO statements; decay constants, residues, mixing observables, other multiplets, loops and uncertainties are outside this record.
