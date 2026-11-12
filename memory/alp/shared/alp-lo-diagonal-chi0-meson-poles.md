---
title: Explicit exact LO pole equations of $\pi^0$, $\eta$, $\eta'$ and the ALP for diagonal $\chi_0$
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Explicit exact LO pole equations of $\pi^0$, $\eta$, $\eta'$ and the ALP for diagonal $\chi_0$

## Assumptions, inputs and notation

- Exact pole equations of the [exact pole-equation record](alp-lo-exact-pole-equations.md): the ALP equation of Step 3, $D^{\rm eff}_{aa}(s)=s-M^2_{\hat a\hat a,\rm can}+2\epsilon^2\operatorname{tr}(\mathcal C\,\mathbb M_s^{-1}(\mathcal C))$ with $-2\operatorname{tr}(\mathcal C\,\mathbb M_s^{-1}(\mathcal C))=\sum_n\operatorname{tr}(\lambda_n\mathcal C)^2/(s-m_{n,0}^2)$, and the meson equation of Step 5. Notation $\mathbb M_s$, $\mathcal C$, $\lambda_n$, $m_{n,0}^2$, $M^2_{\hat a\hat a,\rm can}$, $\epsilon=f_\pi/f_a$ as there; conventions of Alda:2026nrd [@Alda:2026nrd]. No expansion in $\epsilon$.
- $\chi_0=\operatorname{diag}(x_u,x_d,x_s)$ with general $x_u$, $x_d$, $x_s$ ($m_u\neq m_d$ allowed); sources Hermitian and commuting with $Q$. No isospin-limit rotation is used: $\pi^0$, $\eta$, $\eta'$ are the labels of the three flavour-diagonal mass eigenstates in increasing order of $m_{n,0}^2$, whatever their flavour content.
- For diagonal $\chi_0$ the flavour-diagonal mesons and the $K^0$ doublet do not mix directly ([neutral-block record](alp-lo-neutral-block-reduction.md), (a)). The $K^0$ doublet is treated in the [$K^0$ doublet record](alp-lo-k0-doublet-exact-pole.md); its contribution to the secular sum, $4m_{K^0,0}^4|(c_q^{AP})_{ds}|^2/(s-m_{K^0,0}^2)$, is taken from there.
- Entries of $\mathcal C$: $(\mathcal C)_{qq}=-x_q(c_q^{AP})_{qq}-\tfrac{m_0^2}{N_f}c_{GA}$.
- $\mathcal N(s)$: name for the quadratic numerator defined in Step 2, inside this record and the exact record only; it has no separate physical meaning. Script function `mixing_numerator_polynomial`.

## Step 1. Unmixed masses

$\mathbb M_0$ maps diagonal matrices to diagonal matrices, as $X\mapsto\operatorname{diag}(x_qX_{qq})+\tfrac{m_0^2}{N_f}\operatorname{tr}(X)\mathbb 1$. Its eigenvalues on this three-dimensional space, $m_{\pi^0,0}^2<m_{\eta,0}^2<m_{\eta',0}^2$, are the roots of
$$
(s-m_{\pi^0,0}^2)(s-m_{\eta,0}^2)(s-m_{\eta',0}^2)=\prod_q(s-x_q)-\frac{m_0^2}{N_f}\sum_q\prod_{r\neq q}(s-x_r),
$$
and $\det(s\mathbb 1-M^2_{PP})=(s-m_{K^0,0}^2)^2(s-m_{\pi^0,0}^2)(s-m_{\eta,0}^2)(s-m_{\eta',0}^2)$.

## Step 2. Secular sum

Solving $\mathbb M_s(Y)=\mathcal C$ for $Y=\mathbb M_s^{-1}(\mathcal C)$ in flavour space (diagonal part: Sherman-Morrison for the rank-one $m_0^2$ term; $d$-$s$ entries: division by $m_{K^0,0}^2-s$) gives
$$
-2\operatorname{tr}\Big(\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\Big)=\frac{2\,\mathcal N(s)}{(s-m_{\pi^0,0}^2)(s-m_{\eta,0}^2)(s-m_{\eta',0}^2)}+\frac{4m_{K^0,0}^4|(c_q^{AP})_{ds}|^2}{s-m_{K^0,0}^2},
$$
$$
\mathcal N(s)=\sum_q(\mathcal C)_{qq}^2\prod_{r\neq q}(s-x_r)-\frac{m_0^2}{N_f}\sum_{q<r}\big[x_q(c_q^{AP})_{qq}-x_r(c_q^{AP})_{rr}\big]^2(s-x_{p}),
$$
where $p$ is the flavour different from $q$ and $r$. $c_{GA}$ enters only through $(\mathcal C)_{qq}$.

## Step 3. Residues

$\mathcal N(s)$ has degree two in $s$, one less than the cubic of Step 1, so the first term of Step 2 is a sum of three simple fractions. Comparing with $\sum_n\operatorname{tr}(\lambda_n\mathcal C)^2/(s-m_{n,0}^2)$:
$$
\boxed{\,\operatorname{tr}(\lambda_n\mathcal C)^2=\frac{2\,\mathcal N(m_{n,0}^2)}{\prod_{n'\neq n}(m_{n,0}^2-m_{n',0}^2)},\qquad n,n'\in\{\pi^0,\eta,\eta'\}\,}
$$
For example $\operatorname{tr}(\lambda_{\pi^0}\mathcal C)^2=2\mathcal N(m_{\pi^0,0}^2)/[(m_{\pi^0,0}^2-m_{\eta,0}^2)(m_{\pi^0,0}^2-m_{\eta',0}^2)]$.

## Step 4. Secular polynomial

With the canonical kinetic matrix equal to $\mathbb 1$, $\det D(s)=\det(s\mathbb 1-M^2_{\rm can})$ is the polynomial of degree 6 with leading coefficient 1
$$
\boxed{\,\begin{aligned}\det D(s)=(s-m_{K^0,0}^2)\Big\{&(s-m_{K^0,0}^2)\Big[(s-M^2_{\hat a\hat a,\rm can})(s-m_{\pi^0,0}^2)(s-m_{\eta,0}^2)(s-m_{\eta',0}^2)-2\epsilon^2\mathcal N(s)\Big]\\&-4\epsilon^2m_{K^0,0}^4|(c_q^{AP})_{ds}|^2(s-m_{\pi^0,0}^2)(s-m_{\eta,0}^2)(s-m_{\eta',0}^2)\Big\}\end{aligned}\,}
$$
where the product of the three factors is the cubic of Step 1. It contains the exact $K^0$ pole as an overall factor, depends on the $d$-$s$ sources only through $|(c_q^{AP})_{ds}|^2$, and at $\epsilon=0$ reduces to $(s-m_{K^0,0}^2)^2(s-m_{a,0}^2)(s-m_{\pi^0,0}^2)(s-m_{\eta,0}^2)(s-m_{\eta',0}^2)$.

## Step 5. Explicit per-meson equations

Step 5 of the exact record with the residues of Step 3 and the $K^0$ term reads, for $n\in\{\pi^0,\eta,\eta'\}$,
$$
\boxed{\,(s-m_{n,0}^2)\Big[s-M^2_{\hat a\hat a,\rm can}-\epsilon^2\sum_{n'\neq n}\frac{2\,\mathcal N(m_{n',0}^2)}{\prod_{n''\neq n'}(m_{n',0}^2-m_{n'',0}^2)\,(s-m_{n',0}^2)}-\frac{4\epsilon^2m_{K^0,0}^4|(c_q^{AP})_{ds}|^2}{s-m_{K^0,0}^2}\Big]=\epsilon^2\,\frac{2\,\mathcal N(m_{n,0}^2)}{\prod_{n'\neq n}(m_{n,0}^2-m_{n',0}^2)}\,}
$$
where $n'$ and $n''$ run over $\pi^0,\eta,\eta'$. Multiplied by $(s-m_{K^0,0}^2)\prod_{n'\neq n}(s-m_{n',0}^2)$ it is the polynomial of Step 4 divided by $(s-m_{K^0,0}^2)$. All three equations have the same roots; each is written so that its root near $m_{n,0}^2$ is the one of interest. The $K^0$ equation is in the $K^0$ doublet record.

## Verification

[`scripts/derive_alp_diagonal_chi0_meson_poles.py`](../../../scripts/derive_alp_diagonal_chi0_meson_poles.py), which reuses `flavour_diagonal_meson_polynomial`, `mixing_numerator_polynomial` and `explicit_secular_polynomial` of [`derive_alp_meson_mass_shifts_exact.py`](../../../scripts/derive_alp_meson_mass_shifts_exact.py) (the latter also checks Step 2 symbolically and Step 4 against $\det(s\mathbb 1-M^2_{\rm can})$):

- Symbolic: $\mathcal N(s)$ has degree at most 2 (symbolic $x_q$, $c_q^{AP}$, $c_{GA}$, $m_0^2$); the partial fractions of Step 3 for symbolic roots and a generic quadratic numerator; the equation of Step 5 times $(s-m_{K^0,0}^2)\prod_{n'\neq n}(s-m_{n',0}^2)$ equals the polynomial of Step 4 divided by $s-m_{K^0,0}^2$, for each $n$.
- Step 3 against $\operatorname{tr}(\lambda_n\mathcal C)^2$ with 50-digit $\lambda_n$ from the eigenvectors of $M^2_{PP}$: largest difference $2\times10^{-60}$ GeV$^4$ (neutral, CP-conserving and flavour-diagonal test sources).
- Exact rationals, the same three source sets with $m_{a,0}^2=3/5$ GeV$^2$, $B_0=27/10$ GeV, $\operatorname{tr}S_2=1/20$ GeV, $\epsilon=1/2$ and $2$: the polynomial of Step 4 equals $\det(sK-M^2)$ exactly. Exact real-root isolation (to $10^{-50}$) gives six real roots, which agree with [`solve_alp_meson_spectrum.py`](../../../scripts/solve_alp_meson_spectrum.py) within $1.8\times10^{-14}$ GeV$^2$. The three equations of Step 5 vanish at every root other than $m_{K^0,0}^2$, with largest residual $10^{-43}$ GeV$^4$ at 50 digits. At $\epsilon=2$ one root is negative ($M^2$ not positive); that case tests the algebra only.

Units: $x_q,m_0^2,M^2_{\hat a\hat a,\rm can},s$ and masses squared in GeV$^2$; $\mathcal N$ and $\operatorname{tr}(\lambda_n\mathcal C)^2$ in GeV$^6$ and GeV$^4$; couplings and $\epsilon$ dimensionless. Inputs are test values, not fits.

```bash
.venv/bin/python scripts/derive_alp_diagonal_chi0_meson_poles.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_diagonal_chi0_meson_poles.py
```

All checks pass.
