---
title: Closed form and fixed-point iteration for the LO ALP pole far from the unmixed meson masses, for any $f_\pi/f_a$
status: mixed
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: assumption
---

# Closed form and fixed-point iteration for the LO ALP pole far from the unmixed meson masses, for any $f_\pi/f_a$

## Assumptions, inputs and notation

- Exact ALP pole equation, Step 3 of the [exact pole-equation record](alp-lo-exact-pole-equations.md), with $\mathbb M_s^{-1}$ defined there, and the interlacing of the [arrow-matrix record](alp-lo-arrow-matrix-spectrum.md). Conventions of Alda:2026nrd [@Alda:2026nrd], $\epsilon=f_\pi/f_a$, no expansion in $\epsilon$. Sums over $n$ run over mesons with $\operatorname{tr}(\lambda_n\mathcal C)\neq0$.
- $s_a$ denotes the ALP-like pole defined in each case below.
- **Assumption (approximation, stated explicitly).** For the closed form: $M^2_{\hat a\hat a,\rm can}$ is far from all unmixed meson masses, in the sense $q\ll1$ with
$$
q=\epsilon^2\sum_n\frac{\operatorname{tr}(\lambda_n\mathcal C)^2}{\big(M^2_{\hat a\hat a,\rm can}-m_{n,0}^2\big)^2}.
$$
  The iteration and the bounds below do not need $q\ll1$.

## Step 1. Fixed-point form of the exact equation

The exact equation $D^{\rm eff}_{aa}(s)=0$ is $s=g(s)$ with
$$
g(s)=M^2_{\hat a\hat a,\rm can}-2\epsilon^2\operatorname{tr}\Big(\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\Big)=M^2_{\hat a\hat a,\rm can}+\epsilon^2\sum_n\frac{\operatorname{tr}(\lambda_n\mathcal C)^2}{s-m_{n,0}^2},
$$
and $|g'(s)|=\epsilon^2\sum_n\operatorname{tr}(\lambda_n\mathcal C)^2/(s-m_{n,0}^2)^2$.

## Step 2. Frozen closed form

Evaluating $g$ at $s=M^2_{\hat a\hat a,\rm can}$ (approximation under the assumption above) gives
$$
\boxed{\,s_a\simeq s_a^{(1)}=M^2_{\hat a\hat a,\rm can}-2\epsilon^2\operatorname{tr}\Big(\mathcal C\,\mathbb M^{-1}_{M^2_{\hat a\hat a,\rm can}}(\mathcal C)\Big)=M^2_{\hat a\hat a,\rm can}+\epsilon^2\sum_n\frac{\operatorname{tr}(\lambda_n\mathcal C)^2}{M^2_{\hat a\hat a,\rm can}-m_{n,0}^2}\,}
$$
$\mathbb M^{-1}_{M^2_{\hat a\hat a,\rm can}}(\mathcal C)$ is the solution $Y$ of $\mathbb M_s(Y)=\mathcal C$ at $s=M^2_{\hat a\hat a,\rm can}$: the first-order equation of Alda:2026nrd ([weak-mixing record](alp-lo-weak-mixing-expansion.md), Part B) with $m_{a,0}^2$ replaced by $M^2_{\hat a\hat a,\rm can}$. No expansion in $\epsilon$ is involved in this statement.

## Step 3. Iteration and bounds

The fixed-point iteration $s^{(k+1)}=g(s^{(k)})$, $s^{(0)}=M^2_{\hat a\hat a,\rm can}$, has $s^{(1)}=s_a^{(1)}$ of Step 2.

ALP entry above all unmixed meson masses ($M^2_{\hat a\hat a,\rm can}>m_{n,0}^2$ for all $n$). By interlacing, exactly one pole $s_a$ lies above the largest $m_{n,0}^2$. On that region $g$ is decreasing and convex, and $g(s)>M^2_{\hat a\hat a,\rm can}$. Therefore $g$ maps $[M^2_{\hat a\hat a,\rm can},s^{(1)}]$ into itself, and on this interval $|g'|\le|g'(M^2_{\hat a\hat a,\rm can})|=q$. Exact consequences:
- if $q<1$, $g$ shrinks distances on this interval, the iteration converges to $s_a$, $s_a$ lies between $M^2_{\hat a\hat a,\rm can}$ and $s^{(1)}$, and
$$
\boxed{\,|s^{(1)}-s_a|\le\frac{q}{1+q}\,|s^{(1)}-M^2_{\hat a\hat a,\rm can}|,\qquad |s^{(k)}-s_a|\le q^k\,|s^{(1)}-M^2_{\hat a\hat a,\rm can}|\,}
$$
- for any $q$, $|g'(s_a)|\le(s_a-M^2_{\hat a\hat a,\rm can})/(s_a-m_{\rm max}^2)<1$, with $m_{\rm max}^2$ the largest unmixed mass. So $s_a$ always attracts the iteration from nearby points. Convergence from $s^{(0)}=M^2_{\hat a\hat a,\rm can}$ when $q\ge1$ is verified numerically below, not proved.

The same holds with all inequalities reversed for $M^2_{\hat a\hat a,\rm can}$ below all unmixed meson masses ($s_a$ below the smallest $m_{n,0}^2$, $g$ decreasing and concave there).

ALP entry between meson masses. The iteration can converge to a pole only if $|g'|<1$ there. A pole with $|g'(s_a)|>1$ repels the iteration. Then use the two-level form of the [isolated meson-pole record](alp-lo-isolated-meson-pole.md) with the nearest meson as the other level, or the eigenvalues of the arrow matrix.

The frozen form is the first iterate; its error is of order $q$ times the shift $|s^{(1)}-M^2_{\hat a\hat a,\rm can}|$. This is small when the shift is small compared with the distance to the meson masses, for any $\epsilon$.

## Verification

[`scripts/derive_alp_lo_arrow_matrix_poles.py`](../../../scripts/derive_alp_lo_arrow_matrix_poles.py), exact poles from the 50-digit eigenvalues of the arrow matrix. The sign of the $\mathbb M_s^{-1}$ form is checked exactly in [`derive_alp_meson_mass_shifts_exact.py`](../../../scripts/derive_alp_meson_mass_shifts_exact.py).

Arrow model (test values, GeV$^2$): mesons at $40,\dots,44$, $\operatorname{tr}(\lambda_n\mathcal C)=(1,-\tfrac12,\tfrac7{10},\tfrac3{10},-\tfrac65)$; ALP entry $44+d$ (above), $40-d$ (below) or $42.5$ (between). Iterations counted to $10^{-45}$.

| side | $\epsilon$ | $d$ | $\lvert s^{(1)}-s_a\rvert$ | $q$ | $\lvert g'(s_a)\rvert$ | iterations |
|---|---|---|---|---|---|---|
| above | 0.1 | 1 | $2.9\times10^{-4}$ | 0.016 | 0.015 | 25 |
| above | 0.1 | 32 | $2.8\times10^{-8}$ | $2.9\times10^{-5}$ | $2.9\times10^{-5}$ | 11 |
| above | 1 | 1 | 0.85 | 1.57 | 0.41 | 120 |
| above | 1 | 32 | $2.8\times10^{-4}$ | $2.9\times10^{-3}$ | $2.9\times10^{-3}$ | 19 |
| above | 2 | 1 | 5.0 | 6.29 | 0.60 | 206 |
| above | 2 | 32 | $4.4\times10^{-3}$ | 0.012 | 0.011 | 24 |
| below | 2 | 32 | $4.2\times10^{-3}$ | 0.011 | 0.011 | 24 |
| between | 0.1 | 0.5 | $1.3\times10^{-4}$ | 0.032 | 0.032 | 30 |
| between | 0.5 | 0.5 | $4.4\times10^{-2}$ | 0.81 | 0.75 | 348 |
| between | 1 | 0.5 | 0.31 | 3.2 | 2.9 | no convergence |

Checks: the bounds of Step 3 at all points with $q<1$ above and below; convergence to $s_a$ from $M^2_{\hat a\hat a,\rm can}$ at all points above and below, including $q>1$; $|g'(s_a)|\le(s_a-M^2_{\hat a\hat a,\rm can})/(s_a-m_{\rm max}^2)<1$; and, in the between case, convergence only to poles with $|g'|\le1$.

LO blocks (neutral test sources of the exact record): $M^2_{\hat a\hat a,\rm can}$ lies above all unmixed masses for $\epsilon=1/2,1,2$. For $\epsilon=1/2$, $M^2_{\hat a\hat a,\rm can}=1.08$ GeV$^2$ is $0.05$ GeV$^2$ above the heaviest unmixed meson, so $q\gg1$ and $s^{(1)}=2.79$ GeV$^2$ is far from $s_a=1.42$ GeV$^2$. The iteration still converges ($|g'(s_a)|=0.63$). For $\epsilon=1$: $s^{(1)}=3.01$, $s_a=2.93$ GeV$^2$. For $\epsilon=2$: $s^{(1)}=8.80$, $s_a=8.77$ GeV$^2$.

```bash
.venv/bin/python scripts/derive_alp_lo_arrow_matrix_poles.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_arrow_matrix_poles.py
```

All checks pass.
