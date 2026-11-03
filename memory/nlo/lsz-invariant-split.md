---
title: Invariant split of the ALP–meson LSZ overlap
status: derived
created: 2026-11-02
updated: 2026-11-03
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds]
provenance: analytical-derivation
---

# Invariant split of the ALP–meson overlap

## Assumptions and inputs

Use the full-field notation of Alda:2026nrd, with $\Phi_I=(P_i,a)$ and $P=(\pi_3,\eta_8,\eta_0)^T$. Write the local quadratic Lagrangian as

$$
\mathcal L_{\rm quad}=
\frac12(\partial_\mu\Phi_I)K_{IJ}(\partial^\mu\Phi_J)
-\frac12\Phi_I M^2_{IJ}\Phi_J.
$$

The real symmetric matrices $K_{IJ}$ and $M^2_{IJ}$ are divided into meson blocks $K_{ij}$, $M^2_{ij}$ and ALP–meson blocks $K_{i\hat a}$, $M^2_{i\hat a}$. The meson kinetic matrix and the pole operator $M^2_{ij}-m_a^2K_{ij}$ are assumed invertible at the pole. Retain only first order in ALP–meson mixing. Use Alda:2026nrd's trace convention $\langle Y\rangle=\operatorname{tr}(Y)$. The local NLO meson blocks are in [the neutral meson matrix record](neutral-meson-matrices.md); the mixed blocks follow from [the quadratic NLO Lagrangian](quadratic-lagrangian.md). The operator basis and source conventions are those of [@Kaiser:1998ds; @Alda:2026nrd]. Loops and the displayed higher-order remainder in the NLO operator basis are omitted.

## LSZ components and invariant source

Let $s=p^2$. The LSZ eigenvector equation of Alda:2026nrd, with $w_{\hat a}^{(a)}=1$ at this order, gives

$$
(sK_{ij}-M^2_{ij})w_j^{(a)}(s)
+(sK_{i\hat a}-M^2_{i\hat a})=0.
$$

Define the field-redefinition-invariant meson source in this quadratic basis by

$$
\mathcal C_i^{\rm pole}
=M^2_{ij}(K_{\rm meson}^{-1})_{jk}K_{k\hat a}-M^2_{i\hat a},
\qquad (K_{\rm meson})_{ij}=K_{ij}.
$$

This source has no dependence on $s$. The identity

$$
sK_{i\hat a}-M^2_{i\hat a}
=\mathcal C_i^{\rm pole}+(sK_{ij}-M^2_{ij})(K_{\rm meson}^{-1})_{jk}K_{k\hat a}
$$

then splits the LSZ meson components into a pole contribution and the term associated with the ALP–meson kinetic block:

$$
\mathbb M_a^{\rm NLO}\equiv M^2_{ij}-m_a^2K_{ij},
\qquad
w_i^{(a)}(s)
=[(M^2-sK)^{-1}]_{ij}\mathcal C_j^{\rm pole}
-(K_{\rm meson}^{-1})_{ij}K_{j\hat a}.
$$

For CCWZ derivative terms, write $K_{i\hat a}=-K_{ij}v_j+\Delta_i$, where $v_i=(f_\pi/f_a)\operatorname{tr}(c_q^A\lambda_i)$ and $\Delta_i$ contains contributions involving the pure-ALP building block $D_\mu\theta$. Then

$$
w^{(a)}(s)=(M^2-sK)^{-1}\mathcal C^{\rm pole}+v-K^{-1}\Delta.
$$

The retained $u_\mu$-quadratic basis has $\Delta=0$. If source-dependent operators such as $O_{47}$ are included, alignment need not hold; the $O_{47}$ contribution to $\Delta$ is singlet before multiplication by $K^{-1}$. See the [CCWZ kinetic-alignment record](kinetic-block-alignment.md) for its explicit form.

At the ALP pole the first term is $(\mathbb M_a^{\rm NLO})^{-1}_{ij}\mathcal C_j^{\rm pole}$. At LO, this operator reduces in flavor-matrix notation to Alda:2026nrd's Eq. (III.27),

$$
\mathbb M_a(\bm{\mathcal C})
=B_0\{m_q,\bm{\mathcal C}\}
+\frac{m_0^2}{3}\langle\bm{\mathcal C}\rangle\mathbf 1
-m_a^2\bm{\mathcal C}.
$$

For the ALP-dependent field shift $P_i=P_i'-\rho_i a$, the mixed blocks become

$$
K'_{i\hat a}=K_{i\hat a}-K_{ij}\rho_j,
\qquad
M_{i\hat a}^{2\prime}=M^2_{i\hat a}-M^2_{ij}\rho_j.
$$

It follows that $\mathcal C_i^{\rm pole}$ and $\Delta_i=K_{i\hat a}+K_{ij}v_j$ are unchanged, while the second term in the split shifts by $+\rho_i$. Thus $w_i^{(a)\prime}=w_i^{(a)}+\rho_i$, as required by $P_i=P_i'-\rho_i a$. For the ALP-dependent axial transformation $c_q^A\to c_q^A+\kappa$, $c_G\to c_G-\operatorname{tr}\kappa$, one has $\rho_i=\delta v_i=(f_\pi/f_a)\operatorname{tr}(\kappa\lambda_i)$. Hence $w_i^{(a)}$ is an affine field-coordinate component, whereas $w_i^{(a)}-v_i$ is invariant, even when $\Delta\ne0$. Under a regular change of meson coordinates, $\mathcal C_i^{\rm pole}$ transforms as a source covector and $w_i^{(a)}$ as a field vector.

For a linear source $\mathcal J=V_{\hat a}a+V_iP_i$, the field shift gives $V'_{\hat a}=V_{\hat a}-V_i\rho_i$. Its direct-plus-pole coefficient is invariant:

$$
\begin{aligned}
V_{\hat a}+V_iw_i^{(a)}(s)
&=\underbrace{V_{\hat a}-V_i(K_{\rm meson}^{-1})_{ij}K_{j\hat a}}_{\text{shift invariant}}\\
&\quad+V_i[(M^2-sK)^{-1}]_{ij}\mathcal C_j^{\rm pole}.
\end{aligned}
$$

The invariant source and pole operator are the NLO generalization of the leading-order LSZ organization in Alda:2026nrd, Secs. III.C–D.

## U(3) invariant source through the retained local NLO terms

Use a complete normalized $U(3)$ generator basis and arbitrary Hermitian flavor matrices $c_q^A$ and $c_q^P$, without assuming they commute with the diagonal quark-mass spurion. Set $\chi_0=2B_0m_q$, $c_q^{AP}=c_q^A-c_q^P$, $c_{GA}=c_G+\operatorname{tr}c_q^A$, and $N_f=3$. Alda:2026nrd's invariant flavor matrix is

$$
\bm{\mathcal C}=-\frac12\{\chi_0,c_q^{AP}\}-\frac{m_0^2}{3}c_{GA}\mathbf 1.
$$

Including the retained local NLO operators gives

$$
\bm{\mathcal C}_{\rm NLO}=\bm{\mathcal C}+\Delta\bm{\mathcal C}_{\rm NLO},
$$

where, to first order in each retained local LEC,

$$
\begin{aligned}
\Delta\bm{\mathcal C}_{\rm NLO}={}&-
\frac{4L_8}{f_\pi^2}\{\chi_0,\{\chi_0,c_q^{AP}\}\}
-\frac{8L_6}{f_\pi^2}\operatorname{tr}(\chi_0)\{\chi_0,c_q^{AP}\}
-\frac{16L_7}{f_\pi^2}\operatorname{tr}(\chi_0c_q^{AP})\chi_0
\\
&-\frac{\Lambda_2}{3}\left[
c_{GA}\chi_0+\operatorname{tr}(\chi_0c_q^{AP})\mathbf 1
\right]
+\frac{8L_{25}}{f_\pi^2}\left[
c_{GA}\chi_0^2+\operatorname{tr}(\chi_0^2c_q^{AP})\mathbf 1
\right].
\end{aligned}
$$

Projection onto the meson basis gives the source in the LSZ equation:

$$
\mathcal C_i^{\rm pole}
=\frac{f_\pi}{f_a}\left\langle\bm{\mathcal C}_{\rm NLO}\lambda_i^\dagger\right\rangle.
$$

For diagonal spurions, the anticommutators reduce to the corresponding products. For general $c_q^{AP}$, the leading term is not $-\chi_0c_q^{AP}$; the $L_8$ term is the nested anticommutator shown, and the $L_6$ term also requires an anticommutator. The $L_7$, $\Lambda_2$, and $L_{25}$ trace structures retain the displayed orderings. The axial matrix cancels from $\bm{\mathcal C}_{\rm NLO}$ when the mixed kinetic and mass blocks are combined. The derivative LECs $L_4,L_5,L_{18},\Lambda_1$ change these blocks together and do not enter this invariant source.

For these local terms, the derivative bilinears give the kinetic-block relation derived in [the kinetic-alignment record](kinetic-block-alignment.md):

$$
-(K_{\rm meson}^{-1})_{ij}K_{j\hat a}
=\frac{f_\pi}{f_a}\left\langle c_q^A\lambda_i^\dagger\right\rangle.
$$

The ALP-pole overlap is therefore

$$
\boxed{
\langle0|P_i|a\rangle=w_i^{(a)}(m_a^2)
=\frac{f_\pi}{f_a}\left\langle c_q^A\lambda_i^\dagger\right\rangle
+[(\mathbb M_a^{\rm NLO})^{-1}]_{ij}\mathcal C_j^{\rm pole}.
}
$$

At LO this is Eq. (III.25) of Alda:2026nrd, with $K_{ij}=\delta_{ij}$ and diagonal $M^2_{ij}=m_i^2\delta_{ij}$. The displayed NLO expression retains the local operators above for general, noncommuting flavor matrices. It does not include loops or omitted higher-order operators.

## Reproducibility

The general field-shift derivation is [`scripts/derive_nlo_lsz_invariant_split.py`](../../scripts/derive_nlo_lsz_invariant_split.py); the noncommuting flavor-matrix derivation is [`scripts/derive_nlo_lsz_source_projection.py`](../../scripts/derive_nlo_lsz_source_projection.py). Run

```text
.venv/bin/python scripts/derive_nlo_lsz_invariant_split.py
.venv/bin/python scripts/derive_nlo_lsz_source_projection.py
```

The first script constructs the inverse-propagator block matrix and derives the invariant source split and direct-plus-pole cancellation. The second differentiates the LO and NLO trace operators using noncommuting $3\times3$ matrix symbols and tests the ordered formulas on a non-diagonal Hermitian example. Its tests establish the axial-block cancellation, recovery of Alda:2026nrd's leading invariant matrix $\bm{\mathcal C}$, and the NLO spurion structures.
