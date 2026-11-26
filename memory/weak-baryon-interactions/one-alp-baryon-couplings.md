---
title: One-ALP terms from the lowest-order weak meson-baryon Lagrangian
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, Borasoy:2000ki]
provenance: analytical-derivation
---

# One-ALP terms from the lowest-order weak meson-baryon Lagrangian

## Inputs and assumptions

Use the [lowest-order weak meson-baryon Lagrangian and spurion convention](lowest-order-weak-lagrangian.md), from Borasoy:2000ki, and the ALP–meson mixing and field redefinitions of Alda:2026nrd, Eqs. (II.2), (III.26)–(III.30), and (IV.1). The weak spurion is assumed independent of the ALP in the starting quark-field parametrization. Its induced ALP dependence is retained when changing parametrization. Work in natural units, at first order in $1/f_a$, with one ALP, two baryon fields, and no external meson fields. The mixing expansion applies away from pseudoscalar poles.

The baryon expansion is $\bm{\mathcal B}=\sum_i\bm{\lambda_{B_i}} B_i$, with $\langle\bm{\lambda_{B_i}}^\dagger\bm{\lambda_{B_j}}\rangle=\delta_{ij}$; the matrices are left generic. The definitions of $\mathbb M_a^{-1}(\bm{\mathcal C})$ and $\bm{c_q^{VS}}$ are recorded in the [Alda:2026nrd baryon-coupling extract](../extracts/Alda:2026nrd/baryon-couplings.md); the operator $\mathbb M_a$ is given below.

## Field redefinition and expansion

Choose the no-mixing parametrization and vector convention of Alda:2026nrd,

$$
\bm\kappa=-\bm{c_q^A}-\mathbb M_a^{-1}(\bm{\mathcal C}),\qquad
\bm\delta=-\bm{c_q^S},\qquad
\bm g_L=\exp\!\left[\frac{ia}{f_a}(\bm\delta-\bm\kappa)\right].
$$

The weak spurion transforms as $\bm w\to\bm g_L\bm w\bm g_L^\dagger$. The same applies to $\bm w^\dagger$. Setting the transformed meson field to zero, expand the transformed spurion to first order:

$$
\left.\bm w_+\right|_{a,\,\Phi=0}=
\frac{ia}{f_a}\left[\bm{c_q^A}-\bm{c_q^S}+\mathbb M_a^{-1}(\bm{\mathcal C}),\bm w+\bm w^\dagger\right].
$$

Here $\left.\cdots\right|_{a,\,\Phi=0}$ means the part linear in $a$ with no meson fields. The result follows by differentiating the left-handed rotation and inserting the two field-redefinition matrices. The trace and commutator order are retained when inserting the baryon expansion.

## Result

$$
\begin{aligned}
\left.\mathcal L_W^{(0)}\right|_{a\bar BB}={}&\frac{ia}{f_a}\sum_{B_1,B_2}\bar B_1 B_2\Bigg[
 d_w\left\langle\bm{\lambda_{B_1}}^\dagger\left\{\left[\bm{c_q^A}-\bm{c_q^S}+\mathbb M_a^{-1}(\bm{\mathcal C}),\bm w+\bm w^\dagger\right],\bm{\lambda_{B_2}}\right\}\right\rangle\\
&\hspace{30mm}+f_w\left\langle\bm{\lambda_{B_1}}^\dagger\left[\left[\bm{c_q^A}-\bm{c_q^S}+\mathbb M_a^{-1}(\bm{\mathcal C}),\bm w+\bm w^\dagger\right],\bm{\lambda_{B_2}}\right]\right\rangle\Bigg].
\end{aligned}
$$

The mass operator and its input are those of Alda:2026nrd:

$$
\mathbb M_a(\bm A)=B_0\{\bm m_q,\bm A\}+\frac{m_0^2}{N_f}\langle\bm A\rangle\bm 1-m_a^2\bm A,\qquad
\bm{\mathcal C}=-B_0\{\bm m_q,\bm{c_q^{AP}}\}-\frac{m_0^2}{N_f}c_{GA}\bm 1.
$$

This direct weak vertex is generally nonzero. It vanishes when the displayed commutator vanishes. The combination $\bm{c_q^A}-\bm{c_q^S}$ is not invariant under an arbitrary initial quark-field redefinition by itself: the assumption of an ALP-independent weak spurion selects the starting parametrization. For a different starting parametrization, the induced ALP dependence of the weak spurion must also be retained. Other diagrams, including baryon-pole diagrams, are not included in this contact term.

## Reproduction

The SymPy calculation [weak_baryon_one_alp.py](../../scripts/weak_baryon_one_alp.py) expands the rotation, differentiates to extract the ALP coefficient, substitutes the no-mixing redefinition, and projects on generic baryon matrices. Run `.venv/bin/python scripts/weak_baryon_one_alp.py`; its checks run with `.venv/bin/python scripts/weak_baryon_one_alp.py --test`.

The [partonic weak spurion invariance record](partonic-weak-spurion-invariance.md) derives the general field-redefinition invariant and explains how this result follows when the explicit one-ALP weak coupling vanishes in the starting parametrization.
