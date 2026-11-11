---
title: Pole-alignment shift as a constant-generator chiral redefinition
status: derived
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Pole-alignment shift as a constant-generator chiral redefinition

## Inputs and conventions

- $x=a/f_a$, $g_\eta=e^{ixG_\eta}$ for constant generators, $C_A=(C_R-C_L)/2$: [partonic record](alp-general-chiral-redefinition-invariants.md).
- $U'=g_RUg_L^\dagger$, $\chi'=g_R\chi g_L^\dagger$, $r'_\mu=g_Rr_\mu g_R^\dagger+ig_R\partial_\mu g_R^\dagger$ (same for $\ell_\mu$), $\Theta'=\Theta-\frac1{2i}\log(\det g_R/\det g_L)$: [source-transformation record](alp-chiral-source-transformations.md).
- $U=\exp(i\Phi/f_\pi)$, $\Phi=P_i\lambda_i$ with the nine $U(3)$ generators, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, $\lambda_0=\sqrt{2/3}\,\mathbf 1$: [CCWZ expansion](alp-ccwz-arbitrary-order-expansion.md) and `Alda:2026nrd` [@Alda:2026nrd].
- Group element of the pole alignment, $L_a=\exp[-ia\,w_i^{(a)}\lambda_i/(2f_\pi)]$, $R_a=L_a^\dagger$, and the LO blocks of $K,M^2$: [finite-mixing record](finite-alp-meson-mixing.md).
- LO Lagrangian as in `scripts/derive_alp_lo_finite_mixing.py`: $u_\mu=-\partial_\mu\Phi/f_\pi+2(\partial_\mu x)C_A$ at linear order, $\chi(x)=\chi_0+ix\{c_q^P,\chi_0\}$, $X_{TC}=-\langle\Phi\rangle/f_\pi-2\Theta$ with $\Theta=c_Gx$. This sign of the $c_q^P$ term is the one that reproduces $M^2_{i\hat a}$ of the finite-mixing record; it is a consistency check within the project conventions, not an independent check of `Alda:2026nrd`.
- Tree level, LO; vacuum $\Phi=0$, $g_{R,L}(0)=\mathbf 1$.

## Result

1. Choose
$$
G_R=-G_L=-\frac{f_a}{2f_\pi}\,w_i\lambda_i .
$$
Then $g_R=e^{ixG_R}=L_a$ and $g_L=L_a^\dagger=R_a$. Note the labels: the element called $L_a$ in the finite-mixing record acts on right-handed quarks in the convention $U'=g_RUg_L^\dagger$. With the opposite assignment the shift below would change sign.

2. The linear part of $U'=g_RUg_L^\dagger$ is
$$
P_i'=P_i-w_ia,\qquad a'=a,
$$
i.e. $\Phi=T\Phi'$ with the shift $T$ of the [single-pole record](alp-single-pole-alignment.md). With $w_i=w_i^{(k)}(m_k^2)$ this is the pole-aligned shift. The generator-to-$w$ relation is
$$
w_i=-\frac{f_\pi}{f_a}\operatorname{tr}(\lambda_i\,G_R)=-\frac{f_\pi}{2f_a}\operatorname{tr}\bigl[\lambda_i(G_R-G_L)\bigr].
$$

3. For any $g_\eta(x)$ with $g_\eta(0)=\mathbf 1$ and $G_{\eta,1}=-i\,\partial_xg_\eta(0)$, the linear part is again a shift with $w_i=-\frac{f_\pi}{2f_a}\operatorname{tr}[\lambda_i(G_{R,1}-G_{L,1})]$. The vector part and the higher Taylor terms first enter at order $aP$ or $a^2$. Therefore chiral redefinitions of this class produce only shift-type linear maps; the general $T$ of the [several-pole record](alp-multi-pole-alignment.md) is outside this class.

4. Transformed sources at the order needed for the quadratic terms:
$$
C_A'=c_q^A+G_R,\qquad c_q^{P\prime}=c_q^P+G_R,\qquad
c_G'=c_G+\frac{f_a}{2f_\pi}w_i\operatorname{tr}\lambda_i ,
$$
so $c_q^A-c_q^P$ is unchanged, in agreement with the invariant $c_q^{AP}$ of the [dressed-mass record](alp-dressed-mass-first-order-invariants.md). The anomaly shift uses $\operatorname{tr}(G_R-G_L)=-(f_a/f_\pi)\sqrt6\,w_0$. The transformed $\chi'(x)$ also contains an $x^2$ term, which contributes to $M^{2\prime}_{\hat a\hat a}$.

5. The Hessian of the LO Lagrangian built from $C_A'$, $\chi'(x)$ (to order $x^2$) and $\Theta'$ equals $T^{\mathsf T}KT$ and $T^{\mathsf T}M^2T$ exactly. The congruence of the [congruence record](alp-inverse-propagator-congruence.md) is therefore what the chiral redefinition does to the quadratic terms.

The nonlinear terms of the same map (first term of order $aP_iP_j$ for the axial rotation) are treated in the [finite CCWZ redefinition record](alp-ccwz-finite-redefinitions.md); they do not affect the quadratic terms.

## Verification

[`derive_alp_pole_alignment_chiral_generator.py`](../../../scripts/derive_alp_pole_alignment_chiral_generator.py) works with all nine $U(3)$ coordinates plus $a$, keeping the nine $w_i$ symbolic and using exact rational source data (Hermitian $c_q^A,c_q^P$ with complex flavour-changing entries, diagonal $\chi_0$, values not physical). It checks: the original Hessian reproduces the $K$ and $M^2$ blocks of the finite-mixing record; $ixG_R$ equals the exponent of $L_a$; the linear map $P_i'=P_i-w_ia$; the transformed Hessian equals the congruence with the shift; the transformed $\chi'$ has an $x^2$ term; $c_q^{P\prime}=c_q^P+G_R$, $c_q^A-c_q^P$ unchanged and the $c_G$ shift; and, for arbitrary noncommuting $G_{R,1},G_{L,1}$, the shift form of the linear map. All checks pass (about 15 s).

```sh
.venv/bin/python scripts/derive_alp_pole_alignment_chiral_generator.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_pole_alignment_chiral_generator.py
```
