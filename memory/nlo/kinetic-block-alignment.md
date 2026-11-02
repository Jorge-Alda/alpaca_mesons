---
title: Condition for ALP–meson kinetic-block alignment
status: derived
created: 2026-11-02
updated: 2026-11-02
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

## Scope

The claim is limited to the displayed local derivative operators, with loops and omitted higher-order operators excluded as in [the invariant-split record](lsz-invariant-split.md). Extra derivative interactions containing independent ALP sources can change the mixed kinetic block. The algebraic factorization itself is exact for any symmetric $K$ obeying it; when $K$ is truncated at NLO, it is understood to that truncation order.

## Reproducibility

The executable SymPy derivation constructs the common quadratic form, differentiates the LO and each retained NLO derivative term, and adds an independent mixing coefficient to obtain a counterexample. Run .venv/bin/python scripts/derive_alp_kinetic_alignment.py.
