---
title: Finite-mixing ALP couplings to the baryon octet at low energy
status: derived
created: 2026-11-03
updated: 2026-11-03
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Finite-mixing ALP couplings to the baryon octet

## Assumptions and inputs

Use the LO meson theory and canonical bare ALP quadratic term in the [finite-mixing setup](../shared/finite-alp-meson-mixing.md), the CCWZ baryon Lagrangian and notation of [Alda:2026nrd, Sec. V](../../extracts/Alda:2026nrd/baryon-couplings.md), and a stable physical ALP pole below the range where baryon form factors are needed. Keep the one-ALP baryon vertex at first order in the external ALP field, but do not expand the meson component $w_i^{(a)}$ of the physical pole. The baryon parameters are the LO constants $D_B,F_B,S_B,\widetilde S_B$.

Set $s=m_a^2$ in

$$
(sK_{ij}-M^2_{ij})w_j^{(a)}(s)+(sK_{i\hat a}-M^2_{i\hat a})=0.
$$

Use

$$
c_q^{AP}=c_q^A-c_q^P,\qquad c_q^{VS}=c_q^V-c_q^S,\qquad
c_{GA}=c_G+\operatorname{tr}c_q^A,\qquad \chi_0=2B_0m_q,
$$

$$
\mathcal C=-\frac12\{\chi_0,c_q^{AP}\}
-\frac{m_0^2}{N_f}c_{GA}\mathbf1,
\qquad
\mathbb M_a(X)=B_0\{m_q,X\}
+\frac{m_0^2}{N_f}\operatorname{tr}(X)\mathbf1-m_a^2X.
$$

The inverse below requires that $m_a^2$ is not an unmixed meson pole; otherwise use the full degenerate pole subspace.

## Exact meson-source relation

For the complete normalized $U(N_f)$ basis, completeness gives

$$
\sum_iK_{i\hat a}\lambda_i=-\frac{2f_\pi}{f_a}c_q^A,
\qquad
\sum_iM^2_{i\hat a}\lambda_i
=-\frac{2f_\pi}{f_a}\left[\mathcal C+
B_0\{m_q,c_q^A\}+\frac{m_0^2}{N_f}\operatorname{tr}(c_q^A)\mathbf1\right].
$$

Multiplying the pole equation by $\lambda_i$ and summing therefore yields

$$
\mathbb M_a\!\left(\sum_iw_i^{(a)}(m_a^2)\lambda_i\right)
=\frac{2f_\pi}{f_a}
\left[\mathcal C+\mathbb M_a(c_q^A)\right].
$$

Define the ALP-only axial source in the pole-aligned CCWZ coordinates by

$$
\mathbb M_a^{-1}(\mathcal C)=\frac{f_a}{2f_\pi}\sum_iw_i^{(a)}(m_a^2)\lambda_i-c_q^A.
$$

The previous equation proves, without expanding in $w_i^{(a)}$ or $f_\pi/f_a$,

$$
\boxed{\mathbb M_a\!\left(\frac{f_a}{2f_\pi}\sum_iw_i^{(a)}(m_a^2)\lambda_i-c_q^A\right)=\mathcal C,
\qquad \frac{f_a}{2f_\pi}\sum_iw_i^{(a)}(m_a^2)\lambda_i-c_q^A=\mathbb M_a^{-1}(\mathcal C),
\qquad \sum_iw_i^{(a)}(m_a^2)\lambda_i
=\frac{2f_\pi}{f_a}\left[c_q^A+\mathbb M_a^{-1}(\mathcal C)\right].}
$$

The inverse is evaluated at the exact physical $m_a^2$. Thus the invariant matrix appearing in Section V retains the same form, but it is not a small-mixing approximation here: it determines the exact pole eigenvector through the last equality.

## CCWZ baryon building blocks and pole normalization

The linear CCWZ expansions are $u_\mu=-\partial_\mu\Phi/f_\pi+2(\partial_\mu a)c_q^A/f_a$ and $X_{TC}=-\operatorname{tr}\Phi/f_\pi-2c_Ga/f_a$. On the ALP pole, with $P_i=w_i^{(a)}a$, they become

$$
\left.u_\mu\right|_{a,\mathrm{pole}}
=-\frac{2\partial_\mu a}{f_a}\mathbb M_a^{-1}(\mathcal C),
\qquad
\left.\partial_\mu X_{TC}\right|_{a,\mathrm{pole}}
=-\frac{2\partial_\mu a}{f_a}
\left[c_{GA}+\operatorname{tr}(\mathbb M_a^{-1}(\mathcal C))\right].
$$

The vector connection gives

$$
\left.\nabla_\mu\mathcal B\right|_{a,\mathrm{pole}}
=-i\frac{\partial_\mu a}{f_a}[c_q^{VS},\mathcal B].
$$

The pole vector has been normalized here with its original $a$ component equal to one. Its kinetic norm is

$$
K_{\hat a\hat a}+2K_{i\hat a}w_i^{(a)}+K_{ij}w_i^{(a)}w_j^{(a)}
=1+\sum_i\left(w_i^{(a)}-v_i\right)^2
=1+\frac{2f_\pi^2}{f_a^2}\operatorname{tr}\!\left([\mathbb M_a^{-1}(\mathcal C)]^2\right),
\qquad
v_i=\frac{f_\pi}{f_a}\operatorname{tr}(c_q^A\lambda_i).
$$

The last equality follows from $w_i^{(a)}-v_i=(f_\pi/f_a)\operatorname{tr}(\mathbb M_a^{-1}(\mathcal C)\lambda_i)$ and completeness of the normalized Hermitian generators. Thus $a_{\mathrm{phys}}=\sqrt{1+2f_\pi^2\operatorname{tr}([\mathbb M_a^{-1}(\mathcal C)]^2)/f_a^2}\,a$ in the pole-aligned coordinates. Every coupling to a canonically normalized external ALP is divided by this square root.

## Effective baryon couplings

Substitution into the Section V baryon Lagrangian gives

$$
\mathcal L_{aB}=\frac{\partial_\mu a_{\mathrm{phys}}}{f_a}
\sum_{B_1,B_2}\bar B_1\gamma^\mu
\left(c_{B_1B_2}^{V,\mathrm{pole}}+
c_{B_1B_2}^{A,\mathrm{pole}}\gamma_5\right)B_2,
$$

with

$$
c_{B_1B_2}^{V,\mathrm{pole}}
=-\frac{\operatorname{tr}\!\left(
c_q^{VS}[\lambda_{B_1}^\dagger,\lambda_{B_2}]\right)}
{\sqrt{1+2f_\pi^2\operatorname{tr}([\mathbb M_a^{-1}(\mathcal C)]^2)/f_a^2}}.
$$

$$
\begin{aligned}
c_{B_1B_2}^{A,\mathrm{pole}}=-\frac{1}{\sqrt{1+2f_\pi^2\operatorname{tr}([\mathbb M_a^{-1}(\mathcal C)]^2)/f_a^2}}\Bigl\{&
(D_B+F_B)\operatorname{tr}\!\left(
\lambda_{B_1}^\dagger\mathbb M_a^{-1}(\mathcal C)\lambda_{B_2}\right)
+(D_B-F_B)\operatorname{tr}\!\left(
\lambda_{B_1}^\dagger\lambda_{B_2}\mathbb M_a^{-1}(\mathcal C)\right)\\
&+\left[(S_B+\widetilde S_B)\operatorname{tr}(\mathbb M_a^{-1}(\mathcal C))
+\widetilde S_Bc_{GA}\right]
\operatorname{tr}(\lambda_{B_1}^\dagger\lambda_{B_2})\Bigr\}.
\end{aligned}
$$

The algebraic numerator is the Section V expression, but the finite-mixing answer is not generally just that expression times an external normalization. The operator $\mathbb M_a$ contains $-m_a^2X$, and $m_a^2$ must now be the root of the full mixed determinant $\det(sK-M^2)=0$. Thus the axial numerator changes when the physical pole moves. For comparison with a perturbative evaluation at the bare input mass $m_{a,0}$, define $\mathbb M_s(X)=B_0\{m_q,X\}+m_0^2\operatorname{tr}(X)\mathbf1/N_f-sX$. Then

$$
\mathbb M_{m_a^2}^{-1}(\mathcal C)-\mathbb M_{m_{a,0}^2}^{-1}(\mathcal C)
=(m_a^2-m_{a,0}^2)
\mathbb M_{m_a^2}^{-1}\!\left(\mathbb M_{m_{a,0}^2}^{-1}(\mathcal C)\right).
$$

The axial baryon numerator therefore changes through the exact pole mass as well as through the normalization denominator. The vector numerator $c_q^{VS}$ has no $m_a^2$ dependence at this order, so for the pointlike low-energy vertex it changes only through the external-pole normalization. If the earlier perturbative formula was already evaluated at the physical $m_a^2$, its algebraic numerator agrees with the finite-pole result; the additional finite-mixing correction is then the pole normalization. At first order, the mass shift and normalization corrections are higher order in the original small-mixing counting, which is why they did not appear in Section V. No baryon form factors or higher chiral orders are included here.

## Reproducibility

The executable derivation is [`scripts/derive_alp_finite_baryon_couplings.py`](../../../scripts/derive_alp_finite_baryon_couplings.py). It derives the flavor-matrix pole equation, the CCWZ source projection, and the Section V baryon trace coefficients.
