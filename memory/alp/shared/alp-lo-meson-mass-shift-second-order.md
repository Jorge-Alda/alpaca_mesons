---
title: $O(f_\pi^2/f_a^2)$ shift of the LO neutral-meson poles in invariant form
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# $O(f_\pi^2/f_a^2)$ shift of the LO neutral-meson poles in invariant form

## Assumptions, inputs and notation

- LO neutral system on $(P_3,P_8,P_0,P_6,P_7,a)$ with the blocks $K$, $M^2$ of the [finite-mixing record](finite-alp-meson-mixing.md), reduced as in the [charged-decoupling record](alp-lo-charged-meson-decoupling.md) and the [neutral-block record](alp-lo-neutral-block-reduction.md). Conventions of Alda:2026nrd [@Alda:2026nrd]: $\operatorname{tr}(\lambda_i\lambda_k)=2\delta_{ik}$, $N_f=3$, $D(s)=sK-M^2$.
- Sources $c_q^A$, $c_q^P$, $\chi_0=2B_0m_q$ Hermitian and commuting with $Q$ ($\chi_0$ may have a $d$-$s$ entry). $\epsilon=f_\pi/f_a$ (approved notation of the [weak-mixing record](alp-lo-weak-mixing-expansion.md)), varied at fixed sources.
- Neutral space: the Hermitian $3\times3$ matrices that commute with $Q$. It is spanned by the generators $\lambda_3,\lambda_8,\lambda_0,\lambda_6,\lambda_7$, which are orthonormal under the inner product $\tfrac12\operatorname{tr}(XY)$.
- Invariants $c_q^{AP}=c_q^A-c_q^P$, $c_{GA}=c_G+\operatorname{tr}c_q^A$, and Alda:2026nrd's tadpole source ([dressed topological record](alp-dressed-topological-invariant.md))
$$
\mathcal C=-\tfrac12\{\chi_0,c_q^{AP}\}-\tfrac{m_0^2}{N_f}c_{GA}\mathbb 1 .
$$
- Flavour-space operator (the operator of the first-order equation of Alda:2026nrd, see Part B of the weak-mixing record):
$$
\mathbb M_s(X)=\tfrac12\{\chi_0,X\}+\tfrac{m_0^2}{N_f}\operatorname{tr}(X)\mathbb 1-sX ,
$$
  so that $\mathbb M_0$ is the unmixed meson operator.
- Mass eigenstates: $n$ labels the neutral mesons ($\pi^0,\eta,\eta',K^0,\bar K^0$ for diagonal $\chi_0$; for a $d$-$s$ entry in $\chi_0$ the five neutral states mix and $n$ is a generic label). $\lambda_n$ is the neutral Hermitian generator of the meson $n$, defined by
$$
\mathbb M_0(\lambda_n)=m_{n,0}^2\,\lambda_n,\qquad \operatorname{tr}(\lambda_n\lambda_{n'})=2\delta_{nn'} .
$$
  The $\lambda_n$ form an orthonormal eigenbasis of $\mathbb M_0$ on the neutral space under $\tfrac12\operatorname{tr}(XY)$. They are distinct from the Gell-Mann generators $\lambda_3,\lambda_8,\lambda_0,\lambda_6,\lambda_7$. The overall sign of each $\lambda_n$ is free; only squares of traces with $\lambda_n$ enter the poles. Assumption: $m_{n,0}^2\neq m_{a,0}^2$, and $m_{n,0}^2$ is non-degenerate; the degenerate $K^0$ doublet is treated at the end.
- $\Sigma(s)=-s\,c_q^A+\tfrac12\{c_q^P,\chi_0\}-\tfrac{m_0^2c_G}{N_f}\mathbb 1$: the flavour matrix that gives the mixed entry of $D(s)$ (Step 2). It is a matrix of the neutral space.
- Constant axial redefinition: $c_q^{A}\to c_q^A+\kappa$, $c_q^P\to c_q^P+\kappa$, $c_G\to c_G-\operatorname{tr}\kappa$, with $\kappa$ Hermitian and commuting with $Q$ ([general redefinition record](alp-general-chiral-redefinition-invariants.md), [source-transformation record](alp-chiral-source-transformations.md); the $c_G$ law is the Jacobian law of the latter with $g_R=g_L^\dagger=e^{ix\kappa}$).

## Step 1. Meson mass operator in flavour space

$\mathbb M_0$ maps the neutral space into itself, and the meson block of $M^2$ is its matrix in the Gell-Mann basis:
$$
(M^2_{PP})_{ik}=\tfrac12\operatorname{tr}\big(\lambda_i\,\mathbb M_0(\lambda_k)\big),\qquad i,k\in\{3,8,0,6,7\}.
$$
This is the only normalization used below: with the inner product $\tfrac12\operatorname{tr}(XY)$, under which the Gell-Mann generators are orthonormal, $M^2_{PP}$ and $\mathbb M_0$ are the same linear map. Hence the eigenvalues of $M^2_{PP}$ are the $m_{n,0}^2$, and the eigenvectors of $M^2_{PP}$ correspond to the $\lambda_n$ defined above.

$\mathbb M_s$ is self-adjoint with respect to the trace: for any $3\times3$ matrices $X$, $Y$,
$$
\operatorname{tr}\big(Y\,\mathbb M_s(X)\big)=\operatorname{tr}\big(\mathbb M_s(Y)\,X\big).
$$
With $Y=\lambda_n$ this gives, for every $s$ and every $X$,
$$
\operatorname{tr}\big(\lambda_n\,\mathbb M_s(X)\big)=(m_{n,0}^2-s)\operatorname{tr}(\lambda_nX).
$$

## Step 2. Mixed entry

From the blocks, the entry of $D(s)$ between the meson along the neutral generator $\lambda_i$ and the ALP is
$$
D_{ia}(s)=sK_{ia}-M^2_{ia}=\epsilon\operatorname{tr}\big(\lambda_i\Sigma(s)\big).
$$
The change from the Gell-Mann basis to the mass eigenstates is orthogonal, and $D_{ia}$ is linear in $\lambda_i$; therefore the entry between the meson $n$ and the ALP is
$$
D_{na}(s)=\epsilon\operatorname{tr}\big(\lambda_n\Sigma(s)\big).
$$

## Step 3. Splitting off the redundant part

Substitute $c_q^P=c_q^A-c_q^{AP}$ and $c_G=c_{GA}-\operatorname{tr}c_q^A$ in $\Sigma(s)$:
$$
\Sigma(s)=\mathbb M_s(c_q^A)+\mathcal C .
$$
All dependence on $c_q^A$ beyond $c_q^{AP}$ and $c_{GA}$ sits in the first term. With Step 1,
$$
D_{na}(s)=\epsilon\Big[(m_{n,0}^2-s)\operatorname{tr}(\lambda_nc_q^A)+\operatorname{tr}(\lambda_n\mathcal C)\Big].
$$

## Step 4. Value at the unmixed meson mass

At $s=m_{n,0}^2$ the first term of Step 3 vanishes:
$$
D_{na}(m_{n,0}^2)=\epsilon\operatorname{tr}\big(\lambda_n\mathcal C\big).
$$
This holds to all orders in $\epsilon$, because $D_{na}$ is linear in $\epsilon$.

## Step 5. Shift

The exact meson-pole equation of the [exact pole-equation record](alp-lo-exact-pole-equations.md) (also Step 2 of the weak-mixing record), expanded in $\epsilon$, gives
$$
\boxed{\,m_n^2-m_{n,0}^2=\frac{D_{na}(m_{n,0}^2)^2}{m_{n,0}^2-m_{a,0}^2}+O(\epsilon^4)=\epsilon^2\,\frac{\big[\operatorname{tr}(\lambda_n\mathcal C)\big]^2}{m_{n,0}^2-m_{a,0}^2}+O(\epsilon^4)\,}
$$
with $\mathcal C=-\tfrac12\{\chi_0,c_q^{AP}\}-\tfrac{m_0^2}{N_f}c_{GA}\mathbb 1$. Neither $M^2_{aa}-m_{a,0}^2$ nor $K_{aa}-1$ nor the $\operatorname{tr}S_2$ term enters at this order. The sign is that of $m_{n,0}^2-m_{a,0}^2$: the meson is pushed away from the bare ALP mass.

For diagonal $\chi_0=\operatorname{diag}(x_u,x_d,x_s)$ the numerator for $\pi^0,\eta,\eta'$ has a closed form in terms of the unmixed masses; see Step 8 of the [exact pole-equation record](alp-lo-exact-pole-equations.md).

## Step 6. Invariance

Under the constant axial redefinition, $\Sigma(s)\to\Sigma(s)+\mathbb M_s(\kappa)$, so $D_{na}(s)\to D_{na}(s)+\epsilon\,(m_{n,0}^2-s)\operatorname{tr}(\lambda_n\kappa)$ by Step 1, while $c_q^{AP}$, $c_{GA}$ and therefore $\mathcal C$ are unchanged. The added term vanishes at $s=m_{n,0}^2$. Hence $D_{na}(m_{n,0}^2)$ is invariant exactly, and the $O(\epsilon^2)$ shift is invariant. The invariance of the full spectrum to all orders is shown in the [exact pole-equation record](alp-lo-exact-pole-equations.md).

## Step 7. Relation with the first-order equation of Alda:2026nrd

Let $Y$ solve Alda:2026nrd's first-order equation $\mathbb M_{m_{a,0}^2}(Y)=\mathcal C$. Taking the trace with $\lambda_n$ and using Step 1,
$$
\operatorname{tr}(\lambda_n\mathcal C)=(m_{n,0}^2-m_{a,0}^2)\operatorname{tr}(\lambda_nY),
$$
so that
$$
m_n^2-m_{n,0}^2=\epsilon^2\,(m_{n,0}^2-m_{a,0}^2)\,\big[\operatorname{tr}(\lambda_nY)\big]^2+O(\epsilon^4).
$$
In the weak-mixing record the meson part of the first-order ALP pole vector is the matrix $2\epsilon(c_q^A+Y)$; its projection on the meson $n$, obtained with $\tfrac12\operatorname{tr}(\lambda_n\,\cdot)$, is $\epsilon\operatorname{tr}\big(\lambda_n(c_q^A+Y)\big)=D_{na}(m_{a,0}^2)/(m_{n,0}^2-m_{a,0}^2)$, by Step 3. The meson shift is fixed by the $Y$ part only.

## $K^0$ doublet

For diagonal $\chi_0$ the eigenvalue $m_{K^0,0}^2=\tfrac12(x_d+x_s)$ of $\mathbb M_0$ is doubly degenerate, with eigenspace the neutral matrices with only $d$-$s$ entries. Choose the two generators $\lambda_n$ of this eigenspace such that one of them has $\operatorname{tr}(\lambda_n\mathcal C)=0$. That meson stays at $m_{K^0,0}^2$ exactly ([neutral-block record](alp-lo-neutral-block-reduction.md)). The other one is shifted by $\epsilon^2[\operatorname{tr}(\lambda_n\mathcal C)]^2/(m_{K^0,0}^2-m_{a,0}^2)$, where, for any orthonormal choice of the two generators,
$$
\sum_{n\in K^0\ \rm doublet}\big[\operatorname{tr}(\lambda_n\mathcal C)\big]^2=4|(\mathcal C)_{ds}|^2=4m_{K^0,0}^4|(c_q^{AP})_{ds}|^2 .
$$
The fuller treatment is task 2.3 of the openspec change.

## Verification

[`scripts/derive_alp_meson_mass_shifts_second_order.py`](../../../scripts/derive_alp_meson_mass_shifts_second_order.py):

- Symbolic (Hermitian $Q$-commuting $c_q^A$, $c_q^P$, $\chi_0$ with complex $d$-$s$ entries, symbolic $c_G$, $m_0^2$, $f_\pi$, $f_a$, $\kappa$): orthonormality of the Gell-Mann generators, $(M^2_{PP})_{ik}=\tfrac12\operatorname{tr}(\lambda_i\mathbb M_0(\lambda_k))$, self-adjointness for generic complex $X$, $Y$, $\mathbb M_0$ preserving the neutral space (Step 1); Step 2 in the Gell-Mann basis; $\Sigma(s)=\mathbb M_s(c_q^A)+\mathcal C$ as a matrix identity, and no $c_q^A$ in $\mathcal C$ after the substitution (Step 3); Step 6.
- Exact rationals: $Y$ obtained by solving $\mathbb M_{m_{a,0}^2}(Y)=\mathcal C$ as a matrix equation on the neutral space, equal to the solution obtained with $M^2_{PP}-m_{a,0}^2$ in the Gell-Mann basis (Step 7).
- 50 digits, $\lambda_n$ built from the eigenvectors of $M^2_{PP}$: eigen-equation $\mathbb M_0(\lambda_n)=m_{n,0}^2\lambda_n$, $\operatorname{tr}(\lambda_n\lambda_{n'})=2\delta_{nn'}$ (Step 1); $D_{na}(s)$ from the blocks equal to $\epsilon\operatorname{tr}(\lambda_n\Sigma(s))$ and to the split of Step 3 at $s=m_{n,0}^2$, $m_{a,0}^2$, $1/3$ GeV$^2$; Step 4; both equations of Step 7. Largest difference $10^{-59}$.
- Comparison with [`solve_alp_meson_spectrum.py`](../../../scripts/solve_alp_meson_spectrum.py) (float64), inputs `illustrative_lo_inputs` of [`derive_alp_lo_neutral_blocks.py`](../../../scripts/derive_alp_lo_neutral_blocks.py) with $m_{a,0}^2=3/5$ GeV$^2$ (test value away from the $\eta$ and $K^0$ masses), $\operatorname{tr}S_2=0$, $\epsilon=1/20,\dots,1/320$. Largest residual over the five meson poles, fitted exponents: neutral sources (complex $d$-$s$ entries) $4.08$ for the boxed formula and $1.96$ for $m_n^2-m_{n,0}^2$; flavour-diagonal sources $3.99$ and $1.99$.

Units: $f_\pi,f_a$ in GeV; $\chi_0,m_0^2,m_{a,0}^2,s$ in GeV$^2$; couplings and $\epsilon$ dimensionless. Rational inputs are test values, not fits.

```bash
.venv/bin/python scripts/derive_alp_meson_mass_shifts_second_order.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_meson_mass_shifts_second_order.py
```

All checks pass.
