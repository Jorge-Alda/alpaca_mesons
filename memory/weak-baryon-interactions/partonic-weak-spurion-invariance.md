---
title: Partonic weak spurion and ALP field-redefinition invariance
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, Borasoy:2000ki]
provenance: analytical-derivation
---

# Partonic weak spurion and ALP field-redefinition invariance

## Inputs and assumptions

Use the quark field redefinition and coupling shifts of Alda:2026nrd, Eqs. (II.2)–(II.3), and the left-handed weak spurion of the [lowest-order weak Lagrangian](lowest-order-weak-lagrangian.md). Consider the illustrative partonic operator $\bar q_L\gamma^\mu\bm w(a)q_L J_\mu+\mathrm{h.c.}$, with $J_\mu$ unaffected by the quark rotation. The weak matrix is allowed a general explicit ALP dependence at first order. For a current containing rotated quarks, its transformation must also be included. The derivation uses natural units and retains first order in $a/f_a$.

## Derivation

At $a=0$ the weak matrix is the fixed $\bm w$ of the weak Lagrangian. The coefficient of the one-ALP weak interaction is $f_a\left.\partial\bm w(a)/\partial a\right|_{a=0}$. Under the left-handed field redefinition in Alda:2026nrd,

$$
\bm g_L=\exp\!\left[\frac{ia}{f_a}(\bm\delta-\bm\kappa)\right],\qquad
\bm w(a)\to\bm g_L\bm w(a)\bm g_L^\dagger,
$$

and therefore

$$
f_a\left.\frac{\partial\bm w(a)}{\partial a}\right|_{a=0}
\to f_a\left.\frac{\partial\bm w(a)}{\partial a}\right|_{a=0}
+i[\bm\delta-\bm\kappa,\bm w].
$$

Alda:2026nrd gives $\bm{c_q^V}-\bm{c_q^A}\to\bm{c_q^V}-\bm{c_q^A}+\bm\delta-\bm\kappa$. Hence the following matrix is invariant at the retained order:

$$
f_a\left.\frac{\partial\bm w(a)}{\partial a}\right|_{a=0}
-i[\bm{c_q^V}-\bm{c_q^A},\bm w].
$$

Its adjoint gives the coefficient for the Hermitian-conjugate current. The explicit one-ALP weak coupling is independent input; the value $\partial\bm w(a)/\partial a=0$ in one parametrization generally does not remain zero after a field redefinition.

## Relation to the baryon result

In the no-mixing and vector parametrization used in the [one-ALP weak baryon calculation](one-alp-baryon-couplings.md), $\widehat{\bm c_q^V}-\widehat{\bm c_q^A}=\bm{c_q^{VS}}+\mathbb M_a^{-1}(\bm{\mathcal C})$. The transformed one-ALP coefficient is therefore the invariant above plus

$$
i[\bm{c_q^{VS}}+\mathbb M_a^{-1}(\bm{\mathcal C}),\bm w].
$$

When the explicit one-ALP weak coupling vanishes in the starting parametrization, this reduces to $i[\bm{c_q^A}-\bm{c_q^S}+\mathbb M_a^{-1}(\bm{\mathcal C}),\bm w]$, which is the term used there. This reduction does not make that starting assumption invariant; it expresses an invariant result in terms of the parameters of the selected starting parametrization.

## Reproduction

The symbolic derivation is [weak_partonic_invariant.py](../../scripts/weak_partonic_invariant.py). Run `.venv/bin/python scripts/weak_partonic_invariant.py` and `.venv/bin/python scripts/weak_partonic_invariant.py --test`.

## Expansion in weak matrices and diagonal ALP couplings

Using the requested expansion

$$
\bm w(a)=\bm w_0+\bm w_1\frac{a}{f_a}+\mathcal O\!\left(\frac{a^2}{f_a^2}\right),
\qquad (\bm w_0)^a{}_b=\delta^a_2\delta^3_b,
$$

the invariant obtained above is

$$
\bm w_1-i[\bm c_q^V-\bm c_q^A,\bm w_0].
$$

The notation $\bm w_{1,\mathrm{inv}}$ is proposed for this one-ALP weak coefficient; it requires user approval before adoption. For diagonal Hermitian $\bm c_q^V$ and $\bm c_q^A$, their diagonal entries are real. Substituting the specified $\bm w_0$ gives

$$
\bm w_1-i\bigl((c_q^V)_{dd}-(c_q^A)_{dd}-(c_q^V)_{ss}+(c_q^A)_{ss}\bigr)\bm w_0.
$$

In particular, its $(2,3)$ element is $(\bm w_1)_{23}-i\bigl((c_q^V)_{dd}-(c_q^A)_{dd}-(c_q^V)_{ss}+(c_q^A)_{ss}\bigr)$; all other elements equal those of $\bm w_1$. The Hermitian-conjugate matrix has the corresponding conjugate $(3,2)$ element. If $\bm w_1=0$ in the chosen starting parametrization, the invariant vanishes precisely when the left-handed diagonal couplings of $d$ and $s$ are equal. These statements retain the assumption that the other current does not transform.

[weak_partonic_diagonal_couplings.py](../../scripts/weak_partonic_diagonal_couplings.py) calculates the matrix result and its adjoint. Run `.venv/bin/python scripts/weak_partonic_diagonal_couplings.py` or pass `--test` to run its checks.
