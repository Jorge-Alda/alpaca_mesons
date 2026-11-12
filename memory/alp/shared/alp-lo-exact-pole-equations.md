---
title: Exact LO pole equations of the neutral mesons and the ALP for finite $f_\pi/f_a$
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Exact LO pole equations of the neutral mesons and the ALP for finite $f_\pi/f_a$

## Assumptions, inputs and notation

- LO neutral system on $(P_3,P_8,P_0,P_6,P_7,a)$: blocks of the [finite-mixing record](finite-alp-meson-mixing.md), charged mesons decoupled ([charged-decoupling record](alp-lo-charged-meson-decoupling.md)). Conventions of Alda:2026nrd [@Alda:2026nrd], $N_f=3$, $\epsilon=f_\pi/f_a$, $D(s)=sK-M^2$. Kinetic blocks: $K_{PP}=\mathbb 1$, $K_{ia}=-\epsilon\operatorname{tr}(\lambda_ic_q^A)$ for the neutral Gell-Mann generators $\lambda_i$, $K_{aa}=1+2\epsilon^2\operatorname{tr}(c_q^Ac_q^A)$.
- Quark spurion linear in $a$ plus the $a^2$ term of the [joint-vertex record](alp-lo-joint-vertex-expansion.md): $M^2_{aa}=m_{a,0}^2+\epsilon^2\big(\tfrac{2m_0^2}{N_f}c_G^2-B_0\operatorname{tr}S_2\big)$, with $\operatorname{tr}S_2$ one additive input. Its dressed form $B_0\operatorname{tr}\widehat S_2=B_0\operatorname{tr}S_2+2\operatorname{tr}(c_q^A\{c_q^P,\chi_0\})-2\operatorname{tr}((c_q^A)^2\chi_0)$ is (4.1) of the [sum-rule record](alp-lo-mass-sum-rules.md).
- From the [second-order shift record](alp-lo-meson-mass-shift-second-order.md): the neutral space, the operator $\mathbb M_s(X)=\tfrac12\{\chi_0,X\}+\tfrac{m_0^2}{N_f}\operatorname{tr}(X)\mathbb 1-sX$, the tadpole source $\mathcal C=-\tfrac12\{\chi_0,c_q^{AP}\}-\tfrac{m_0^2}{N_f}c_{GA}\mathbb 1$, the matrix $\Sigma(s)$ with $D_{ia}(s)=\epsilon\operatorname{tr}(\lambda_i\Sigma(s))$ and $\Sigma(s)=\mathbb M_s(c_q^A)+\mathcal C$, and the generators $\lambda_n$ of the meson mass eigenstates, $\mathbb M_0(\lambda_n)=m_{n,0}^2\lambda_n$, $\operatorname{tr}(\lambda_n\lambda_{n'})=2\delta_{nn'}$ (Steps 1-3 there). $M^2_{PP}$ is the matrix of $\mathbb M_0$ in the Gell-Mann basis. Only the statements of that record that do not involve an expansion in $\epsilon$ are used here.
- Inverse of $\mathbb M_s$: for $s$ different from every $m_{n,0}^2$, $\mathbb M_s^{-1}(X)$ denotes the neutral matrix $Y$ with $\mathbb M_s(Y)=X$. It is a linear map of the neutral space into itself, and in the eigenbasis
$$
\mathbb M_s^{-1}(X)=\sum_n\lambda_n\,\frac{\operatorname{tr}(\lambda_nX)}{2\,(m_{n,0}^2-s)} .
$$
  Since $\mathbb M_s(Y)=\mathbb M_0(Y)-sY$, the inverse of $s\mathbb 1-M^2_{PP}$ corresponds to $-\mathbb M_s^{-1}$. This is a rewriting, not an additional input.
- $D^{\rm eff}_{aa}(s)$ (approved notation of the [weak-mixing record](alp-lo-weak-mixing-expansion.md)): ALP entry of $D$ after eliminating the mesons. $m_{K^0,0}^2=\tfrac12(x_d+x_s)$ ([neutral-block record](alp-lo-neutral-block-reduction.md)).
- No expansion in $\epsilon$ and no assumption on its size. Sources Hermitian and commuting with $Q$; Steps 1-7 allow a $d$-$s$ entry in $\chi_0$, Step 8 assumes $\chi_0=\operatorname{diag}(x_u,x_d,x_s)$.

## Step 1. Schur factorization

Eliminating the meson rows of $D(s)\Phi=0$ gives, for $\det(s\mathbb 1-M^2_{PP})\neq0$,
$$
\det D(s)=\det\big(s\mathbb 1-M^2_{PP}\big)\,D^{\rm eff}_{aa}(s),\qquad
D^{\rm eff}_{aa}(s)=D_{aa}(s)-D_{aP}(s)\big(s\mathbb 1-M^2_{PP}\big)^{-1}D_{Pa}(s),
$$
with $D_{aa}(s)=sK_{aa}-M^2_{aa}$. In flavour space, using the mixed entry and $\mathbb M_s^{-1}$,
$$
D^{\rm eff}_{aa}(s)=sK_{aa}-M^2_{aa}+2\epsilon^2\operatorname{tr}\Big(\Sigma(s)\,\mathbb M_s^{-1}\big(\Sigma(s)\big)\Big).
$$

## Step 2. Canonical frame: all blocks in invariants

Take the field change $P_i=P_i'+\epsilon\operatorname{tr}(\lambda_ic_q^A)\,a$ for every neutral Gell-Mann generator $\lambda_i$, $a=a'$ (the canonical coordinates of the finite-mixing record). Its matrix $T$ has $\det T=1$. It is a congruence $D\to T^{\mathsf T}DT$ and leaves poles and $\det D$ unchanged ([congruence record](alp-inverse-propagator-congruence.md)). Using $\Sigma(s)=\mathbb M_s(c_q^A)+\mathcal C$ and (4.1),
$$
T^{\mathsf T}D(s)T=s\mathbb 1-M^2_{\rm can},\qquad
M^2_{{\rm can},PP}=M^2_{PP},\quad
(M^2_{\rm can})_{ia}=-\epsilon\operatorname{tr}(\lambda_i\mathcal C),
$$
$$
M^2_{\hat a\hat a,\rm can}=m_{a,0}^2+\epsilon^2\Big[\frac{2m_0^2}{N_f}c_{GA}^2-B_0\operatorname{tr}\widehat S_2\Big].
$$
The kinetic matrix becomes $\mathbb 1$ and the mixed entry no longer depends on $s$; in the mass basis it is $-\epsilon\operatorname{tr}(\lambda_n\mathcal C)$. In $M^2_{\hat a\hat a,\rm can}$ the $c_G^2$, $(\operatorname{tr}c_q^A)^2$ and $c_G\operatorname{tr}c_q^A$ terms combine into $c_{GA}^2$, and the $\operatorname{tr}((c_q^A)^2\chi_0)$ and $\operatorname{tr}(c_q^A\{c_q^P,\chi_0\})$ terms are absorbed by $\operatorname{tr}\widehat S_2$. This is why the $\operatorname{tr}\widehat S_2$ form is used below: with $\operatorname{tr}S_2$ the ALP entry is not invariant by itself.

## Step 3. Exact ALP pole equation

The Schur complement is unchanged by $T$ (adding multiples of meson rows and columns to the ALP row and column). Hence
$$
\boxed{\,D^{\rm eff}_{aa}(s)=s-M^2_{\hat a\hat a,\rm can}+2\epsilon^2\operatorname{tr}\Big(\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\Big)\,}
$$
identical as a function of $s$ to the form of Step 1. By the eigenbasis form of $\mathbb M_s^{-1}$,
$$
-2\operatorname{tr}\Big(\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\Big)=\sum_n\frac{\big[\operatorname{tr}(\lambda_n\mathcal C)\big]^2}{s-m_{n,0}^2}.
$$
The factor 2 comes from the normalization $\operatorname{tr}(\lambda_n\lambda_{n'})=2\delta_{nn'}$, and the sign from $\mathbb M_s=\mathbb M_0-s$. The ALP-like pole, and every pole that is not an unmixed meson mass, solves $D^{\rm eff}_{aa}(m^2)=0$.

## Step 4. Invariance of the full spectrum

$M^2_{\rm can}$ depends only on $\chi_0$, $m_0^2$, $m_{a,0}^2$, $\epsilon$, $c_q^{AP}$, $c_{GA}$ and $\operatorname{tr}\widehat S_2$. Under the constant axial redefinition ($c_q^{A,P}\to c_q^{A,P}+\kappa$, $c_G\to c_G-\operatorname{tr}\kappa$, and $S_2$ transformed by $M'=g_RMg_L^\dagger$, which gives $\operatorname{tr}S_2'=\operatorname{tr}S_2-4\operatorname{tr}(\kappa P_1)-4\operatorname{tr}(\kappa^2m_q)$ with $P_1=\{c_q^P,m_q\}$) all of them are unchanged. Therefore every pole is invariant for any value of $\epsilon$. With the bare $\operatorname{tr}S_2$ held fixed instead, the spectrum would change.

## Step 5. Exact meson pole equation

For a non-degenerate $m_{n,0}^2$, multiply $D^{\rm eff}_{aa}(s)=0$ by $s-m_{n,0}^2$:
$$
\boxed{\,(s-m_{n,0}^2)\Big[s-M^2_{\hat a\hat a,\rm can}-\epsilon^2\sum_{n'\neq n}\frac{\big[\operatorname{tr}(\lambda_{n'}\mathcal C)\big]^2}{s-m_{n',0}^2}\Big]=\epsilon^2\big[\operatorname{tr}(\lambda_n\mathcal C)\big]^2\,}
$$
In flavour space the sum is $-2\operatorname{tr}\big(\mathcal C\,\mathbb M_s^{-1}(\mathcal C-\tfrac12\operatorname{tr}(\lambda_n\mathcal C)\,\lambda_n)\big)$: the term $\tfrac12\operatorname{tr}(\lambda_n\mathcal C)\,\lambda_n$ is subtracted from $\mathcal C$ before applying $\mathbb M_s^{-1}$. The bare form
$$
(s-m_{n,0}^2)\Big[sK_{aa}-M^2_{aa}-\sum_{n'\neq n}\frac{D_{n'a}(s)^2}{s-m_{n',0}^2}\Big]=D_{na}(s)^2,
$$
with the mass-basis entries $D_{na}(s)=\epsilon\operatorname{tr}(\lambda_n\Sigma(s))$, is the same equation: both sides differ by terms that cancel identically. The root near $m_{n,0}^2$ is the meson pole.

## Step 6. Poles as eigenvalues of the canonical mass matrix

Since $\det T=1$ and the canonical kinetic matrix is $\mathbb 1$, $\det D(s)=\det(s\mathbb 1-M^2_{\rm can})=\det(s\mathbb 1-M^2_{PP})\,D^{\rm eff}_{aa}(s)$: the six poles are exactly the eigenvalues of the real symmetric matrix $M^2_{\rm can}$. In the meson mass basis this matrix has the arrow form, from which follow the interlacing of poles and unmixed meson masses ([arrow-matrix record](alp-lo-arrow-matrix-spectrum.md)) and closed-form approximations for a meson pole far from the other meson masses ([isolated meson-pole record](alp-lo-isolated-meson-pole.md)) and for an ALP pole far from the meson masses ([isolated ALP-pole record](alp-lo-isolated-alp-pole.md)).

## Step 7. Domain of the Schur form

$D^{\rm eff}_{aa}$ has poles at $s=m_{n,0}^2$, so the Schur form says nothing at these points; the polynomial $\det D(s)$ is valid everywhere. For a non-degenerate $m_{n,0}^2$,
$$
\det D(m_{n,0}^2)=-\epsilon^2\big[\operatorname{tr}(\lambda_n\mathcal C)\big]^2\prod_{n'\neq n}(m_{n,0}^2-m_{n',0}^2).
$$
Consequences: if $\operatorname{tr}(\lambda_n\mathcal C)\neq0$, $m_{n,0}^2$ is not a pole for any $\epsilon\neq0$ and every pole is a zero of $D^{\rm eff}_{aa}$. If $\operatorname{tr}(\lambda_n\mathcal C)=0$, the meson $n$ does not mix in the canonical frame and stays at $m_{n,0}^2$ for any $\epsilon$; $D^{\rm eff}_{aa}$ is then regular there and does not show this pole. For an eigenvalue of multiplicity $g$ (the $K^0$ doublet, $g=2$), $g-1$ poles stay at the unmixed value, and all $g$ if $\operatorname{tr}(\lambda\,\mathcal C)=0$ for every neutral $\lambda$ in the eigenspace. When $m^2$ approaches $m_{n,0}^2$ the Schur form is numerically ill-conditioned; use the polynomial or the [two-level record](alp-lo-near-degeneracy-two-level.md).

## Step 8. Explicit form for diagonal $\chi_0$

Let $\chi_0=\operatorname{diag}(x_u,x_d,x_s)$, so the entries of $\mathcal C$ are $(\mathcal C)_{qq}=-x_q(c_q^{AP})_{qq}-\tfrac{m_0^2}{N_f}c_{GA}$ and $(\mathcal C)_{ds}=-m_{K^0,0}^2(c_q^{AP})_{ds}$. $\mathbb M_0$ maps diagonal matrices to diagonal matrices and acts on the $d$-$s$ entries as multiplication by $m_{K^0,0}^2$, so the $\pi^0,\eta,\eta'$ and the $K^0$ doublet do not mix. The unmixed $\pi^0$, $\eta$, $\eta'$ masses are the roots of
$$
\prod_{n=\pi^0,\eta,\eta'}(s-m_{n,0}^2)=\prod_q(s-x_q)-\frac{m_0^2}{N_f}\sum_q\prod_{r\neq q}(s-x_r),
$$
and $\det(s\mathbb 1-M^2_{PP})=(s-m_{K^0,0}^2)^2\prod_n(s-m_{n,0}^2)$. Solving $\mathbb M_s(Y)=\mathcal C$ for $Y=\mathbb M_s^{-1}(\mathcal C)$ in flavour space (diagonal part: Sherman-Morrison for the rank-one $m_0^2$ term; $d$-$s$ entries: division by $m_{K^0,0}^2-s$) gives
$$
-2\operatorname{tr}\Big(\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\Big)=\frac{2\,\mathcal N(s)}{\prod_n(s-m_{n,0}^2)}+\frac{4m_{K^0,0}^4|(c_q^{AP})_{ds}|^2}{s-m_{K^0,0}^2},
$$
$$
\mathcal N(s)=\sum_q(\mathcal C)_{qq}^2\prod_{r\neq q}(s-x_r)-\frac{m_0^2}{N_f}\sum_{q<r}\big[x_q(c_q^{AP})_{qq}-x_r(c_q^{AP})_{rr}\big]^2(s-x_{p}),
$$
where $p$ is the flavour different from $q$ and $r$, and the products and sums over $n$ run over $\pi^0,\eta,\eta'$. Here $\mathcal N(s)$ is only a name for this quadratic numerator inside this record (it has no separate physical meaning); the script function is `mixing_numerator_polynomial`. $c_{GA}$ enters only through $(\mathcal C)_{qq}$.

The secular polynomial, of degree 6 with leading coefficient $\det K=1$, is
$$
\boxed{\,\det D(s)=(s-m_{K^0,0}^2)\Big\{(s-m_{K^0,0}^2)\Big[(s-M^2_{\hat a\hat a,\rm can})\prod_n(s-m_{n,0}^2)-2\epsilon^2\mathcal N(s)\Big]-4\epsilon^2m_{K^0,0}^4|(c_q^{AP})_{ds}|^2\prod_n(s-m_{n,0}^2)\Big\}\,}
$$
It contains the exact $K^0$ pole as an overall factor, depends on the $d$-$s$ sources only through $|(c_q^{AP})_{ds}|^2$, and at $\epsilon=0$ reduces to $(s-m_{K^0,0}^2)^2(s-m_{a,0}^2)\prod_n(s-m_{n,0}^2)$.

Per meson: comparing residues with Step 3, for $n=\pi^0,\eta,\eta'$
$$
\big[\operatorname{tr}(\lambda_n\mathcal C)\big]^2=\frac{2\mathcal N(m_{n,0}^2)}{\prod_{n'\neq n}(m_{n,0}^2-m_{n',0}^2)},
$$
and for the $K^0$ doublet the sum of $[\operatorname{tr}(\lambda_n\mathcal C)]^2$ over the two generators is $4m_{K^0,0}^4|(c_q^{AP})_{ds}|^2$. Inserted in Step 5 these give the explicit meson pole equations.

## Verification

[`scripts/derive_alp_meson_mass_shifts_exact.py`](../../../scripts/derive_alp_meson_mass_shifts_exact.py):

- Kinetic blocks in trace form, and $\sum_i[\operatorname{tr}(\lambda_iX)]^2=2\operatorname{tr}(XX)$ over the neutral Gell-Mann generators for neutral $X$ (symbolic).
- Step 1 for a generic symmetric four-meson block (symbolic), and for the exact LO blocks with symbolic $s$. The flavour-space form of Step 1 against the Gell-Mann-basis form, with $\mathbb M_s^{-1}$ obtained by solving the matrix equation on the neutral space (exact rationals, symbolic $s$).
- Step 3 in trace form: the sign of the rewrite, $(s-\mathbb M_0)^{-1}(\mathcal C)=-\mathbb M_s^{-1}(\mathcal C)$, and $-2\operatorname{tr}(\mathcal C\,\mathbb M_s^{-1}(\mathcal C))$ against the Gell-Mann-basis form (exact, symbolic $s$); the boxed $D^{\rm eff}_{aa}$ against the bare Schur complement (exact); the sum over $\lambda_n$ and the eigenbasis form of $\mathbb M_s^{-1}$ (50 digits at $s=1/3$ and $3/5$ GeV$^2$); neutral and flavour-diagonal sources.
- Steps 2 and 4 with symbolic $Q$-commuting sources, $\chi_0$ diagonal and with a complex $d$-$s$ entry, symbolic $B_0$, $\operatorname{tr}S_2$, $\kappa$; $S_2'$ from `constant_rotation` of [`derive_alp_lo_mass_sum_rules.py`](../../../scripts/derive_alp_lo_mass_sum_rules.py). Also $\det(sK-M^2)$ unchanged for a rational $\kappa$.
- Steps 3, 5, 7 in a generic symbolic mass-basis model (three mesons); the flavour-space form of the sum in Step 5 at 50 digits.
- Step 6: $\det(s\mathbb 1-M^2_{\rm can})=\det(sK-M^2)=\det(s\mathbb 1-M^2_{PP})\,D^{\rm eff}_{aa}(s)$ (exact rationals, symbolic $s$). The arrow form and the approximations are checked in [`derive_alp_lo_arrow_matrix_poles.py`](../../../scripts/derive_alp_lo_arrow_matrix_poles.py).
- Step 8 symbolically against $\det(s\mathbb 1-M^2_{\rm can})$ ($x_q$, $m_0^2$, $c_q^{AP}$, $c_{GA}$, $\epsilon$, $M^2_{\hat a\hat a,\rm can}$ symbolic), and exactly against $\det(sK-M^2)$ of `closed_form_blocks` with rational sources ($m_{a,0}^2=3/5$ GeV$^2$, $B_0=27/10$ GeV, $\operatorname{tr}S_2=1/20$ GeV) at $\epsilon=1/2$ and $\epsilon=2$. The second value gives one negative root ($M^2$ not positive); it tests the algebra only.
- Roots: exact real-root isolation (six real roots, $m_{K^0,0}^2$ of multiplicity one for complex $d$-$s$ sources and two for flavour-diagonal sources) against [`solve_alp_meson_spectrum.py`](../../../scripts/solve_alp_meson_spectrum.py): largest difference $2.2\times10^{-16}$ GeV$^2$ ($\epsilon=1/2$) and $1.8\times10^{-14}$ GeV$^2$ ($\epsilon=2$); $|D^{\rm eff}_{aa}|<2.5\times10^{-34}$ at all roots other than $m_{K^0,0}^2$.
- Per-meson closed forms of Step 8 against $[\operatorname{tr}(\lambda_n\mathcal C)]^2$ with 50-digit $\lambda_n$ (difference below $10^{-40}$).

Units: $f_\pi,f_a,B_0,\operatorname{tr}S_2$ in GeV; $\chi_0,m_0^2,m_{a,0}^2,M^2,s$ in GeV$^2$; couplings and $\epsilon$ dimensionless. Rational inputs are test values, not fits.

```bash
.venv/bin/python scripts/derive_alp_meson_mass_shifts_exact.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_meson_mass_shifts_exact.py
```

All checks pass (about 35 s).
