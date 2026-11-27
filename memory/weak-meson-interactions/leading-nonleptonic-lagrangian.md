---
title: Leading non-leptonic weak meson Lagrangian
status: mixed
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, Gerard:2005yk, Borasoy:2000ki, Tandean:1999nv]
provenance: analytical-derivation
---

# Leading non-leptonic weak meson Lagrangian

## Inputs and scope

The $U(3)$ operator content is taken from Gerard:2005yk, Eqs. (6)–(9). The convention for $\bm u$, $\bm u_\mu$, $\bm\chi_+$ and $f_\pi$ follows Alda:2026nrd. The weak octet spurion and 27-plet tensor are specified in [the weak meson-baryon Lagrangian](../weak-baryon-interactions/lowest-order-weak-lagrangian.md) and [the 27-plet note](../weak-baryon-interactions/nonleptonic-27-plet.md). The coefficients $G_8$, $G_8^s$, $G_8^m$ and $G_{27}$ use the names of Gerard:2005yk; their numerical normalization and phase here depend on the displayed spurion and tensor conventions and are not identified with tabulated values without conversion.

Consider even-intrinsic-parity, $|\Delta S|=1$ operators at chiral order $p^2$, first order in the weak interaction. The displayed terms use the $SU(3)_L$ octet and 27-plet sources, with the singlet pseudoscalar retained in $U(3)$. They are separate from the order-$p^0$ weak baryon operators. Electroweak-penguin structures and a complete basis with arbitrary ALP-dependent weak sources are outside this record.

## Octet

A convenient form for the leading terms is

$$
\begin{aligned}
\mathcal L_{W,8}^{(2)}=f_\pi^4\bigg[&
G_8\left\langle\bm u\bm w\bm u^\dagger\,\bm u_\mu\bm u^\mu\right\rangle
+G_8^s\left\langle\bm u\bm w\bm u^\dagger\,\bm u_\mu\right\rangle\left\langle\bm u^\mu\right\rangle\\
&+G_8^m\left\langle\bm u\bm w\bm u^\dagger\,\bm\chi_+\right\rangle\bigg]+\mathrm{h.c.}
\end{aligned}
$$

The first term is the usual octet derivative interaction. The second couples the weak octet to the singlet current and vanishes in the $SU(3)$ theory without singlet external axial sources; it need not vanish in the $U(3)$ theory. The third is the weak mass operator. For real $G_8$, the first term and its Hermitian conjugate can equivalently be written $G_8 f_\pi^4\langle\bm w_+\bm u_\mu\bm u^\mu\rangle$, using the definition of $\bm w_+$ in the weak baryon note. The Hermitian-conjugate form above also permits complex coefficients.

Gerard:2005yk discusses the removal of the weak mass term for its specified hadronic processes and assumptions. This does not establish that it may be discarded in an ALP calculation with ALP-dependent sources or after a different field redefinition; the full transformed Lagrangian must be used consistently.

## 27-plet

Using the $9\times9$ tensor-product matrix $\bm T$ from the [27-plet note](../weak-baryon-interactions/nonleptonic-27-plet.md), the leading derivative term is

$$
\mathcal L_{W,27}^{(2)}=G_{27}f_\pi^4\,
\operatorname{Tr}_{\!\mathrm{flavor}\otimes\mathrm{flavor}}
\left[\bm T\left((\bm u^\dagger\bm u_\mu\bm u)\otimes(\bm u^\dagger\bm u^\mu\bm u)\right)\right]+\mathrm{h.c.}
$$

Each $\bm u^\dagger\bm u_\mu\bm u$ transforms in the left-handed flavor space; the contraction is invariant when $\bm T$ transforms as specified in the 27-plet note. A tensor with both isospin components represents both; a pure $|\Delta I|=3/2$ tensor gives only that projection. The source choice, normalization and weak coefficient must be matched consistently. $G_{27}$ is independent of $G_8$ and of the baryon coefficient $\beta_{27}$.

## First terms in the meson expansion

With external sources set to zero, Alda:2026nrd gives $\bm u=\exp[i\bm\Phi/(2f_\pi)]$ and

$$
\bm u_\mu=-\frac{\partial_\mu\bm\Phi}{f_\pi}+\mathcal O(\bm\Phi^3),\qquad
\bm u^\dagger\bm u_\mu\bm u=-\frac{\partial_\mu\bm\Phi}{f_\pi}+\mathcal O(\bm\Phi^2).
$$

Thus the two displayed derivative structures start at quadratic order as

$$
\begin{aligned}
\left.\mathcal L_{W,8}^{(2)}\right|_{\Phi^2}&=G_8 f_\pi^2\left\langle\bm w\,\partial_\mu\bm\Phi\,\partial^\mu\bm\Phi\right\rangle+\mathrm{h.c.},\\
\left.\mathcal L_{W,27}^{(2)}\right|_{\Phi^2}&=G_{27} f_\pi^2\operatorname{Tr}_{\!\mathrm{flavor}\otimes\mathrm{flavor}}\left[\bm T\bigl(\partial_\mu\bm\Phi\otimes\partial^\mu\bm\Phi\bigr)\right]+\mathrm{h.c.}
\end{aligned}
$$

These are only the derivative structures proportional to $G_8$ and $G_{27}$; the singlet and weak mass terms are not expanded in this calculation. For nonzero external sources and ALP-dependent spurions, use the full covariant building blocks before expanding.

## Reproduction

[weak_pseudoscalar_leading_terms.py](../../scripts/weak_pseudoscalar_leading_terms.py) expands the two derivative operators with vanishing external sources. Run `.venv/bin/python scripts/weak_pseudoscalar_leading_terms.py`; pass `--test` to run its checks. The script does not derive the weak mass term, the singlet term, or an exhaustive basis in the presence of ALP sources.

## One ALP and the lowest meson multiplicities

Use $\bm w(a)=\bm w_0+(a/f_a)\bm w_1+\mathcal O(a^2/f_a^2)$, as in the [partonic spurion record](../weak-baryon-interactions/partonic-weak-spurion-invariance.md), and allow the 27-plet tensor to have its own ALP dependence. All couplings and spurions in the following expressions are in the same arbitrary starting quark-field parametrization. Set external sources other than the ALP to zero and work to first order in $1/f_a$, in natural units. The fields $\bm\Phi$ are the meson fields, not physical mass eigenstate amplitudes. The expressions are the terms for the specific weak operators above, not a complete basis of independent ALP-dependent weak operators.

The required inputs from [Alda:2026nrd Appendix A](../extracts/Alda:2026nrd/appendix-a-expansion.md) are $\bm\chi_0=2B_0\bm m_q$,

$$
\begin{aligned}
\bm u_\mu&=-\frac{\partial_\mu\bm\Phi}{f_\pi}+2\frac{\partial_\mu a}{f_a}\bm c_q^A+\cdots,\\
\bm\chi_+&=2\bm\chi_0+2i\frac{a}{f_a}[\bm c_q^S,\bm\chi_0]
+\frac{a}{f_af_\pi}\{\bm\Phi,\{\bm c_q^P,\bm\chi_0\}\}+\cdots,\\
\bm u\bm w(a)\bm u^\dagger&=\bm w_0+\frac{a}{f_a}\bm w_1
+\frac{i}{2f_\pi}[\bm\Phi,\bm w_0]
+\frac{ia}{2f_af_\pi}[\bm\Phi,\bm w_1]+\cdots.
\end{aligned}
$$

Here the omitted terms include higher field orders and terms irrelevant to the coefficients below; the displayed part of $\bm\chi_+$ is restricted to the orders required here. The derivative operators contain no term with one ALP and no mesons. Their lowest terms with mesons are

$$
\begin{aligned}
\left.\mathcal L_{W,8}^{(2)}\right|_{a\Phi,\,G_8}
={}&-2G_8 f_\pi^3\frac{\partial_\mu a}{f_a}
\left\langle\bm w_0\{\bm c_q^A,\partial^\mu\bm\Phi\}\right\rangle+\mathrm{h.c.},\\
\left.\mathcal L_{W,8}^{(2)}\right|_{a\Phi,\,G_8^s}
={}&-2G_8^s f_\pi^3\frac{\partial_\mu a}{f_a}
\left[\langle\bm w_0\bm c_q^A\rangle\langle\partial^\mu\bm\Phi\rangle
+\langle\bm w_0\partial^\mu\bm\Phi\rangle\langle\bm c_q^A\rangle\right]+\mathrm{h.c.},\\
\left.\mathcal L_{W,27}^{(2)}\right|_{a\Phi}
={}&-2G_{27}f_\pi^3\frac{\partial_\mu a}{f_a}
\operatorname{Tr}_{\!\mathrm{flavor}\otimes\mathrm{flavor}}\!\left[\bm T\left(
\bm c_q^A\otimes\partial^\mu\bm\Phi+\partial^\mu\bm\Phi\otimes\bm c_q^A\right)\right]+\mathrm{h.c.}
\end{aligned}
$$

The tensor $\bm T$ in this expression is its value at $a=0$. An explicit one-ALP coefficient of $\bm T(a)$, like $\bm w_1$ in the derivative octet terms, first enters these derivative operators at one ALP and two meson fields. The terms displayed above are weak ALP–meson bilinears; their presence in a particular field parametrization is not itself a physical amplitude.

The weak mass operator may contain a term with one ALP and no meson fields. Keeping also its first term with a meson gives

$$
\begin{aligned}
\left.\mathcal L_{W,8}^{(2)}\right|_{a,\,G_8^m}
={}&2G_8^m f_\pi^4\frac{a}{f_a}\left\langle
\bm w_1\bm\chi_0+i\bm w_0[\bm c_q^S,\bm\chi_0]\right\rangle+\mathrm{h.c.},\\
\left.\mathcal L_{W,8}^{(2)}\right|_{a\Phi,\,G_8^m}
={}&G_8^m f_\pi^3\frac{a}{f_a}\big\langle
 i[\bm\Phi,\bm w_1]\bm\chi_0
 +\bm w_0\{\bm\Phi,\{\bm c_q^P,\bm\chi_0\}\}
 -[\bm\Phi,\bm w_0][\bm c_q^S,\bm\chi_0]
\big\rangle+\mathrm{h.c.}
\end{aligned}
$$

The ALP-only term or particular bilinears can vanish for specific couplings, spurion choices, or phases. The Hermitian-conjugate contribution must be included before drawing that conclusion.

### Alda:2026nrd no-mixing parametrization

After the strong-sector ALP–meson no-mixing field redefinition of Alda:2026nrd, substitute

$$
\widehat{\bm c_q^A}=-\mathbb M_a^{-1}(\bm{\mathcal C}),\qquad
\widehat{\bm c_q^P}=-\bm c_q^{AP}-\mathbb M_a^{-1}(\bm{\mathcal C}),\qquad
\widehat{\bm c_q^S}=0,
$$

and use the weak spurion in that same parametrization. At first order its coefficient is

$$
\widehat{\bm w_1}=\bm w_1+i\left[\bm c_q^A-\bm c_q^S+\mathbb M_a^{-1}(\bm{\mathcal C}),\bm w_0\right].
$$

In particular, the two weak-mass expressions become

$$
\begin{aligned}
\left.\mathcal L_{W,8}^{(2)}\right|_{a,\,G_8^m}
={}&2G_8^m f_\pi^4\frac{a}{f_a}\langle\widehat{\bm w_1}\bm\chi_0\rangle+\mathrm{h.c.},\\
\left.\mathcal L_{W,8}^{(2)}\right|_{a\Phi,\,G_8^m}
={}&G_8^m f_\pi^3\frac{a}{f_a}\left\langle
 i[\bm\Phi,\widehat{\bm w_1}]\bm\chi_0
 -\bm w_0\{\bm\Phi,\{\bm c_q^{AP}+\mathbb M_a^{-1}(\bm{\mathcal C}),\bm\chi_0\}\}
\right\rangle+\mathrm{h.c.}
\end{aligned}
$$

For the derivative terms replace $\bm c_q^A$ by its hatted value above. Removing strong-sector ALP–meson mixing does not remove the weak bilinears calculated here. The inverse operator requires the non-degenerate perturbative treatment of Alda:2026nrd, away from meson poles. An ALP-dependent 27-plet tensor must also be transformed for calculations at one ALP and two or more mesons.

The ordered-matrix expansion of all four terms is calculated by [weak_meson_one_alp_minimal.py](../../scripts/weak_meson_one_alp_minimal.py). Run `.venv/bin/python scripts/weak_meson_one_alp_minimal.py`; pass `--test` for its checks. The script produces the ordered products; the commutators and anticommutators above group those products without assuming flavor matrices commute.
