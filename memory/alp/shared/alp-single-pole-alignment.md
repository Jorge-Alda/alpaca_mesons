---
title: ALP field aligned with one selected pole
status: derived
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# ALP field aligned with one selected pole

## Inputs and assumptions

- $D(s)=sK-M^2$ on $\Phi=(P_i,a)$, its congruence $D'=T^{\mathsf T}DT$ and the pole vector $v^{(k)}$: [congruence record](alp-inverse-propagator-congruence.md).
- The pole equation and $w_i^{(a)}(s)=[(M^2-sK)^{-1}]_{ij}(sK_{j\hat a}-M^2_{j\hat a})$: [finite-mixing record](finite-alp-meson-mixing.md), notation of `Alda:2026nrd` [@Alda:2026nrd]. Here $w_i^{(k)}$ is the same object for a selected pole $s=m_k^2$.
- $K$ positive definite, $K,M^2$ real symmetric. No expansion in the mixing or in $f_\pi/f_a$.

Write $D_{PP}$, $D_{P\hat a}$, $D_{\hat a\hat a}$ for the meson, mixed and ALP blocks of $D$.

## Regular pole

1. The shift $P_i=P_i'+w_ia'$, $a=a'$ is $\Phi=T\Phi'$ with $T=\begin{pmatrix}\mathbf 1&w\\0&1\end{pmatrix}$. Then $D'_{PP}=D_{PP}$ and
$$
D'_{i\hat a}(s)=D_{i\hat a}(s)+D_{ij}(s)w_j,\qquad
D'_{\hat a\hat a}(s)=D_{\hat a\hat a}+2w_iD_{i\hat a}+w_iD_{ij}w_j .
$$

2. If $D_{PP}(m_k^2)=m_k^2K_{PP}-M^2_{PP}$ is invertible, choose $w=w^{(k)}(m_k^2)=-D_{PP}(m_k^2)^{-1}D_{P\hat a}(m_k^2)$. Then
$$
m_k^2K'_{i\hat a}-M^{2\prime}_{i\hat a}=0 .
$$
Because $\det D=\det D_{PP}\,(D_{\hat a\hat a}-D_{\hat aP}D_{PP}^{-1}D_{P\hat a})$, also $D'_{\hat a\hat a}(m_k^2)=0$: the primed pole vector is $(0,\dots,0,1)$, and the pole normalization becomes $K'_{\hat a\hat a}$. The amplitude is $\mathcal V'_{\hat a}/\sqrt{K'_{\hat a\hat a}}$, as in the finite-mixing record.

3. What the shift does not do: $K'_{i\hat a}$ and $M^{2\prime}_{i\hat a}$ are each nonzero, and $D'_{i\hat a}(s)\neq0$ for $s\neq m_k^2$. The ALP coordinate is decoupled from the mesons only at the selected pole.

## Regularity condition

The shift is unique if and only if $D_{PP}(m_k^2)$ is invertible. Given that $m_k^2$ is a pole, this is equivalent to: no pole vector at $m_k^2$ has a vanishing ALP component.

- If $D_{PP}(m_k^2)$ is invertible, a null vector $(u,0)$ would need $D_{PP}u=0$, so $u=0$.
- If $D_{PP}(m_k^2)$ has a one-dimensional null space spanned by $u$, rotate the mesons orthogonally so that $D_{PP}=\operatorname{diag}(0,d_2,\dots)$ with $d_j\neq0$. Then $\det D(m_k^2)=-b_1^2\,d_2\cdots$, with $b_1=u^{\mathsf T}D_{P\hat a}$. A pole therefore requires $b_1=0$, so $(u,0)$ is a pole vector.
- If the null space of $D_{PP}$ has dimension two or more, some combination $u$ obeys $D_{\hat aP}u=0$, and $(u,0)$ is again a pole vector.

Near such points $w^{(k)}$ is still exact but can be large; no small-mixing expansion is involved.

## Exceptional poles: when the shift fails

Use the generalized eigenproblem directly: compute the null space of $m_k^2K-M^2$.

1. **Pole without ALP component.** The null space is spanned by $(u,0)$. The meson rows of the pole equation can be solved for $w$ (with free components), but the ALP row then fails, so $(w,1)$ is not a pole vector. No shift aligns this pole with $a'$; the state does not contain the ALP coordinate at $s=m_k^2$. A general $T$ with a column proportional to $(u,0)$ aligns it with a meson coordinate.
2. **Degenerate pole.** The null space has dimension at least two and contains a vector with nonzero ALP component. Any such vector, rescaled to $(w,1)$, aligns the pole; $w$ is fixed only up to vectors $u$ with $(u,0)$ in the null space. The residue of $D^{-1}$ is $\sum_rv^{(k,r)}v^{(k,r)\mathsf T}$ over a $K$-orthonormal basis of the null space. It does not depend on the choice of basis, so the sum of squared amplitudes over the degenerate states does not either.

If $K$ is not positive definite, a pole vector with $v^{(k)\mathsf T}Kv^{(k)}=0$ can occur and the pole of $D^{-1}$ is not simple; this lies outside the physical case and is discussed in the [several-pole record](alp-multi-pole-alignment.md).

## Verification

[`derive_alp_single_pole_alignment.py`](../../../scripts/derive_alp_single_pole_alignment.py) uses the exact rational $N_f=3$ data of the congruence record. For each of the four poles it checks: invertible meson block; $w^{(k)}$ equals the null vector of the generalized eigenproblem; $m_k^2K'_{i\hat a}-M^{2\prime}_{i\hat a}=0$; $K'_{i\hat a}\neq0$, $M^{2\prime}_{i\hat a}\neq0$; $D'_{i\hat a}(s)\neq0$ as a function of $s$; $D'_{\hat a\hat a}(m_k^2)=0$; $K'_{\hat a\hat a}$ equals the pole normalization; and agreement with the general expressions of [`derive_alp_finite_mixing.py`](../../../scripts/derive_alp_finite_mixing.py) (meson overlap and ALP block after the shift, with its $\rho=-w$). It checks the determinant $-b_1^2d_2d_3$ symbolically, and two exceptional examples with a decoupled $P_3$: a pole without ALP component (shift impossible), and a degenerate pole (any $w_3$ aligns; rank-two residue equal to the orthonormal sum). All checks pass.

```sh
.venv/bin/python scripts/derive_alp_single_pole_alignment.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_single_pole_alignment.py
```

The shift is the linear part of a chiral redefinition, see the [chiral-generator record](alp-pole-alignment-chiral-generator.md).
