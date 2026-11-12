---
title: LO neutral poles as eigenvalues of an arrow matrix, and interlacing with the unmixed meson masses
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# LO neutral poles as eigenvalues of an arrow matrix, and interlacing with the unmixed meson masses

## Assumptions, inputs and notation

- Canonical frame of the [exact pole-equation record](alp-lo-exact-pole-equations.md) (Step 2): $T^{\mathsf T}D(s)T=s\mathbb 1-M^2_{\rm can}$ with $\det T=1$, $M^2_{{\rm can},PP}=M^2_{PP}$, mixed entry $-\epsilon\operatorname{tr}(\lambda_i\mathcal C)$, ALP entry $M^2_{\hat a\hat a,\rm can}$. Conventions of Alda:2026nrd [@Alda:2026nrd], $\epsilon=f_\pi/f_a$, no expansion in $\epsilon$.
- $\lambda_n$, $m_{n,0}^2$, $\mathcal C$, $\mathbb M_s$ as in the [second-order shift record](alp-lo-meson-mass-shift-second-order.md) (Steps 1-3 there, which do not use an expansion in $\epsilon$).
- The poles are $\mu_1\le\dots\le\mu_6$. $\operatorname{tr}(\lambda_n\mathcal C)$ (GeV$^2$) is the coupling of the meson $n$ to the ALP in the canonical frame; its sign follows the free sign of $\lambda_n$, and only $\operatorname{tr}(\lambda_n\mathcal C)^2$ enters the poles.
- Arrow matrix: a square matrix whose only nonzero entries are on the diagonal, in the last row and in the last column.

## Step 1. Poles as eigenvalues

Because $\det T=1$ and the canonical kinetic matrix is $\mathbb 1$, the poles of the LO neutral system are exactly the eigenvalues of the real symmetric $6\times6$ matrix $M^2_{\rm can}$. Rotating the meson block to the mass basis, $\lambda_i\to\lambda_n$, is orthogonal and does not change the eigenvalues; the mixed entry becomes $-\epsilon\operatorname{tr}(\lambda_n\mathcal C)$ (Step 2 of the second-order record, applied to $\mathcal C$). The result is the arrow matrix
$$
\boxed{\,M^2_{\rm can}\;\simeq\;\begin{pmatrix}\operatorname{diag}\big(m_{n,0}^2\big)&\big(-\epsilon\operatorname{tr}(\lambda_n\mathcal C)\big)_n\\\big(-\epsilon\operatorname{tr}(\lambda_n\mathcal C)\big)_n^{\mathsf T}&M^2_{\hat a\hat a,\rm can}\end{pmatrix}\,}
$$
Its characteristic polynomial is
$$
\det(s\mathbb 1-M^2_{\rm can})=\det(s\mathbb 1-M^2_{PP})\,D^{\rm eff}_{aa}(s)=(s-M^2_{\hat a\hat a,\rm can})\prod_n(s-m_{n,0}^2)-\epsilon^2\sum_n\operatorname{tr}(\lambda_n\mathcal C)^2\prod_{n'\neq n}(s-m_{n',0}^2),
$$
with $D^{\rm eff}_{aa}(s)=s-M^2_{\hat a\hat a,\rm can}-\epsilon^2\sum_n\operatorname{tr}(\lambda_n\mathcal C)^2/(s-m_{n,0}^2)$, the arrow form of the secular equation (Step 3 of the exact record). The last expression is a polynomial identity and holds also at $s=m_{n,0}^2$.

## Step 2. Interlacing

In this step the unmixed masses are labelled by their order, $m_{1,0}^2\le\dots\le m_{5,0}^2$. $M^2_{PP}$ is the upper-left $5\times5$ block of $M^2_{\rm can}$; by the Cauchy interlacing theorem for a real symmetric matrix and its block with one row and column removed,
$$
\boxed{\,\mu_1\le m_{1,0}^2\le\mu_2\le m_{2,0}^2\le\mu_3\le\dots\le\mu_5\le m_{5,0}^2\le\mu_6\,}
$$
i.e. $\mu_k\le m_{k,0}^2\le\mu_{k+1}$ for $k=1,\dots,5$, for any $\epsilon$.

If the $m_{k,0}^2$ are distinct and $\operatorname{tr}(\lambda_n\mathcal C)\neq0$ for every meson $n$, the inequalities are strict and each interval contains exactly one pole: on every interval between neighbouring $m_{k,0}^2$ the function $D^{\rm eff}_{aa}$ is continuous, has $dD^{\rm eff}_{aa}/ds=1+\epsilon^2\sum_n\operatorname{tr}(\lambda_n\mathcal C)^2/(s-m_{n,0}^2)^2>0$, and goes from $-\infty$ to $+\infty$. So one pole lies below $m_{1,0}^2$, one in each gap $(m_{k,0}^2,m_{k+1,0}^2)$, and one above $m_{5,0}^2$. If $\operatorname{tr}(\lambda_n\mathcal C)=0$ for a meson $n$, $\mu=m_{n,0}^2$ is an exact pole (Step 7 of the exact record). For the $K^0$ doublet of diagonal $\chi_0$, one combination has $\operatorname{tr}(\lambda\,\mathcal C)=0$ and gives an exact pole at $m_{K^0,0}^2$, where two neighbouring inequalities become equalities.

Interlacing places every pole in a known interval without approximation. It does not tell which pole is "ALP-like"; that label is discussed in the [pole-tracking record](alp-lo-pole-tracking.md). The sum of the poles is $\operatorname{tr}M^2_{\rm can}=\operatorname{tr}M^2_{PP}+M^2_{\hat a\hat a,\rm can}$.

## Verification

[`scripts/derive_alp_lo_arrow_matrix_poles.py`](../../../scripts/derive_alp_lo_arrow_matrix_poles.py), LO blocks of `canonical_mass_matrix_from_sources` ([`derive_alp_meson_mass_shifts_exact.py`](../../../scripts/derive_alp_meson_mass_shifts_exact.py)) with neutral and flavour-diagonal test sources, $m_{a,0}^2=3/5$ GeV$^2$, $B_0=27/10$ GeV, $\operatorname{tr}S_2=1/20$ GeV, $\epsilon=1/2,1,2$, 50 digits:

- the rotated $M^2_{\rm can}$ equals the arrow matrix with border $-\epsilon\operatorname{tr}(\lambda_n\mathcal C)$; its eigenvalues equal those of $M^2_{\rm can}$ and the generalized eigenvalues of the bare $(K,M^2)$; the determinant formula holds at $s=1/3$ GeV$^2$ and at every $m_{n,0}^2$ (Step 1); the identity $\det(s\mathbb 1-M^2_{\rm can})=\det(sK-M^2)$ is also checked exactly in the exact-record script;
- $\mu_k\le m_{k,0}^2\le\mu_{k+1}$ holds in all LO cases and in arrow models with $M^2_{\hat a\hat a,\rm can}$ below, inside and above the meson masses and $\epsilon$ up to 2, while the opposite order $m_{k,0}^2\le\mu_k$ fails in every case (Step 2).

Units: GeV$^2$ for masses squared, $\operatorname{tr}(\lambda_n\mathcal C)$ and $M^2_{\hat a\hat a,\rm can}$; $\epsilon$ dimensionless. Inputs are test values, not fits.

```bash
.venv/bin/python scripts/derive_alp_lo_arrow_matrix_poles.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_arrow_matrix_poles.py
```

All checks pass.
