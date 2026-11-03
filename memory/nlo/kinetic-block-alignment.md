---
title: Condition for ALP–meson kinetic-block alignment
status: derived
created: 2026-11-02
updated: 2026-11-03
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds]
provenance: analytical-derivation
---

# Condition for ALP–meson kinetic-block alignment

## Assumptions and inputs

Use $\Phi=\sum_i\lambda_iP_i$, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, and Alda:2026nrd's linear vielbein

$$
u_\mu^{(1)}=-\frac{\partial_\mu\Phi}{f_\pi}+2\frac{\partial_\mu a}{f_a}c_q^A
$$

from [the Appendix A expansion](../extracts/Alda:2026nrd/appendix-a-expansion.md) [@Alda:2026nrd]. The retained local derivative operators are the LO kinetic term and $L_4,L_5,L_{18},\Lambda_1$ in [the NLO Lagrangian](lagrangian.md) and [its quadratic expansion](quadratic-lagrangian.md) [@Kaiser:1998ds; @Alda:2026nrd]. Work at quadratic field order and first order in ALP–meson mixing. The meson kinetic matrix $K_{ij}$ is symmetric and invertible.

## Derivation

Expand the axial matrix in the complete generator basis. With the stated normalization, its coefficient along $\lambda_i$ is $\tfrac12\operatorname{tr}(c_q^A\lambda_i)$, so the linear vielbein can be written

$$
u_\mu^{(1)}=-\frac{1}{f_\pi}\sum_i\lambda_i\left(\partial_\mu P_i-v_i\partial_\mu a\right),\qquad
v_i=\frac{f_\pi}{f_a}\operatorname{tr}(c_q^A\lambda_i).
$$

The dagger in the projection may be included for a general Hermitian basis; for the real Hermitian Gell-Mann basis used here, $\lambda_i^\dagger=\lambda_i$. Every retained derivative operator is a quadratic form in $u_\mu^{(1)}$ or its trace. Consequently its entire quadratic derivative Lagrangian has the form

$$
\mathcal L_{\rm kin}^{(2)}=\frac12(\partial_\mu P-v\partial_\mu a)^T K(\partial^\mu P-v\partial^\mu a)
$$

for the same meson kinetic matrix $K$ that appears in the ALP-free meson term. Differentiating with respect to the meson and ALP derivatives gives

$$
K_{i\hat a}=-K_{ij}v_j.
$$

Since $K$ is invertible, this proves

$$
-(K^{-1})_{ij}K_{j\hat a}=v_i=\frac{f_\pi}{f_a}\operatorname{tr}(c_q^A\lambda_i^\dagger).
$$

The argument applies term by term: the LO, $L_5$, $\Lambda_1$, $L_4$, and $L_{18}$ quadratic derivative terms use the same replacement and each mixed coefficient is the polarization of its meson kinetic quadratic form. Thus the relation holds for the displayed LO and local NLO Lagrangians, through the retained order. It is not a generic identity for arbitrary kinetic matrices.

For example, add an independent allowed quadratic mixing $\Delta_i\partial_\mu P_i\partial^\mu a$ to a theory whose meson kinetic matrix is still $K$. Then $K_{i\hat a}=-K_{ij}v_j+\Delta_i$, and

$$
-(K^{-1})_{ij}K_{j\hat a}=v_i-(K^{-1})_{ij}\Delta_j.
$$

In one dimension this is $v-\Delta/K$, displaying the failure without the common-derivative assumption. More generally, the necessary and sufficient condition for the stated overlap is precisely $K_{i\hat a}=-K_{ij}v_j$. A symmetry or operator construction that enforces the common derivative replacement can guarantee it; quadratic block algebra alone cannot.

## CCWZ restriction on the independent mixing

At quadratic field order, a meson derivative enters through $u_\mu$. The other relevant one-derivative building block is

$$
D_\mu\theta=\langle u_\mu\rangle-\partial_\mu X_{TC}
=2c_{GA}\frac{\partial_\mu a}{f_a},
$$

which is purely ALP and already $O(1/f_a)$. Consequently, terms quadratic in $u_\mu$ obey the common derivative replacement used above. An independent mixed column $\Delta_i$ can arise from terms containing $D_\mu\theta$ paired with the meson-derivative part of $u_\mu$ (and any allowed non-derivative spurions).

For $O_{47}=D_\mu\hat\theta\langle C^\mu\rangle\langle M\rangle$, expansion about the real-mass vacuum gives

$$
\mathcal L_{47}^{(a\Phi)}=-\frac{4g_{47}c_{GA}\langle\chi_0\rangle}{f_af_\pi}\,\partial_\mu a\,\langle\partial^\mu\Phi\rangle,
\qquad
\Delta_i^{(47)}=-\frac{4g_{47}c_{GA}\langle\chi_0\rangle}{f_af_\pi}\,\operatorname{tr}(\lambda_i),
$$

where $g_{47}$ is the operator coefficient and $\Phi=\lambda_iP_i$. Thus $O_{47}$ contributes only to the $U(1)$ singlet component in a trace-orthogonal generator basis. A non-diagonal meson kinetic matrix can spread $K^{-1}\Delta^{(47)}$ across components, while the operator's own mixed column remains singlet. The single-trace $O_{46}$ gives a mass-spurion-weighted direction and is removed by the higher-order field redefinition recorded in [the singlet-source record](singlet-source-derivative-basis.md). For the retained $L_4,L_5,L_{18},\Lambda_1$ basis, no independent $D_\mu\theta$ term is present, so $\Delta_i=0$ and the alignment result applies.

More generally, at this order the source-dependent derivative structures listed in [the singlet-source record](singlet-source-derivative-basis.md) give the following candidate directions:

$$
\Delta_i^{(V_4)}=\frac{2V_4c_{GA}}{f_af_\pi}\operatorname{tr}(\lambda_i),\qquad
\Delta_i^{(46)}=-\frac{4g_{46}c_{GA}}{f_af_\pi}\operatorname{tr}(\lambda_i\chi_0),\qquad
\Delta_i^{(47)}=-\frac{4g_{47}c_{GA}\langle\chi_0\rangle}{f_af_\pi}\operatorname{tr}(\lambda_i).
$$

The $V_4$ term is removed in the chosen reduced convention, and $O_{46}$ is removable by the cited higher-order field redefinition. $O_{50}$ contains two factors of $D_\mu\theta$ and therefore has no $a\Phi$ kinetic term at first order in $1/f_a$. These directions are the restriction on the possible $\Delta_i$ at this derivative and field order; they do not permit an arbitrary independent flavor vector.

## Transformation of the LSZ components

Under an ALP-dependent axial chiral field transformation, $c_q^A\to c_q^A+\kappa$ and $c_G\to c_G-\operatorname{tr}\kappa$ [@Alda:2026nrd]. Hence $c_{GA}=c_G+\operatorname{tr}c_q^A$ and $D_\mu\theta$ are invariant. The induced meson-coordinate shift is $P_i=P_i'-\rho_i a$, with

$$
\rho_i=\delta v_i=\frac{f_\pi}{f_a}\operatorname{tr}(\kappa\lambda_i).
$$

Since $K'_{i\hat a}=K_{i\hat a}-K_{ij}\rho_j$ and $v_i'=v_i+\rho_i$, the residual column $\Delta_i=K_{i\hat a}+K_{ij}v_j$ is invariant. The LSZ component therefore transforms affinely,

$$
w_i^{(a)\prime}(s)=w_i^{(a)}(s)+\rho_i,
\qquad
w_i^{(a)\prime}(s)-v_i'=w_i^{(a)}(s)-v_i.
$$

In general,

$$
w^{(a)}(s)=\bigl[(M^2-sK)^{-1}\bigr]\mathcal C^{\rm pole}+v-K^{-1}\Delta.
$$

Thus $w_i^{(a)}$ is a field-coordinate component rather than an invariant. Its affine shift is fixed by $\kappa$, while the departure from the common-derivative alignment is controlled by the allowed $D_\mu\theta$ operators and is invariant under this transformation. For $O_{47}$, that departure starts in the singlet direction before multiplication by $K^{-1}$.

## Scope

The claim is limited to the displayed local derivative operators, with loops and omitted higher-order operators excluded as in [the invariant-split record](lsz-invariant-split.md). The general source-dependent derivative expansion includes $V_4D_\mu\psi D^\mu\theta$; Kaiser:1998ds removes it by a meson-field redefinition, as recorded in [the source extract](../extracts/Kaiser:1998ds/singlet-source-derivative-terms.md). If retained explicitly, it contributes to the independent $D_\mu\theta$ column and alignment must include it. The algebraic factorization itself is exact for any symmetric $K$ obeying it; when $K$ is truncated at NLO, it is understood to that truncation order.

## Reproducibility

The executable SymPy derivation constructs the common quadratic form, differentiates the LO and each retained NLO derivative term, and adds an independent mixing coefficient to obtain a counterexample. The CCWZ restriction, including the $O_{47}$ singlet column and affine LSZ covariance, is in [the CCWZ derivation script](../../scripts/derive_ccwz_kinetic_covariance.py). Run .venv/bin/python scripts/derive_alp_kinetic_alignment.py and .venv/bin/python scripts/derive_ccwz_kinetic_covariance.py.
