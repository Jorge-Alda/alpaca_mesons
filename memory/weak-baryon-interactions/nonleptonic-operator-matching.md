---
title: Non-leptonic weak operators and the baryon octet spurion
status: mixed
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Salone:2026byq, Borasoy:2000ki, Alda:2026nrd]
provenance: analytical-derivation
---

# Non-leptonic weak operators and the baryon octet spurion

## Inputs and scope

The partonic operator is from Salone:2026byq, Eq. (12); the octet baryon Lagrangian is from Borasoy:2000ki and is recorded in [lowest-order weak Lagrangian](lowest-order-weak-lagrangian.md). The field redefinition follows Alda:2026nrd, Eq. (II.2). The calculation uses natural units, three light flavors ordered $(u,d,s)$, and retains terms linear in $a/f_a$. Color and Lorentz indices are spectators in the flavor calculation. The example is one current-current operator; it is not the complete weak Hamiltonian.

## Charged-current operator

After integrating out the $W$ boson, a representative $s\to d$ term is

$$
\mathcal L^{\Delta S=1}=\frac{4G_F}{\sqrt2}V_{ud}V_{us}^*
(\bar u\gamma^\mu P_Ls)(\bar d\gamma_\mu P_Lu)+\mathrm{h.c.}
$$

The two flavor matrices in $\bar q_L\gamma^\mu(\cdots)q_L$ are respectively

$$
\begin{pmatrix}0&0&1\\0&0&0\\0&0&0\end{pmatrix},\qquad
\begin{pmatrix}0&0&0\\1&0&0\\0&0&0\end{pmatrix}.
$$

Their tensor product carries four flavor indices. The prefactor is displayed explicitly; QCD evolution, further Wilson coefficients, other color contractions and penguin operators are not calculated here.

## Field redefinition of both currents

In the convention used for spurion transformations in the [partonic invariance record](partonic-weak-spurion-invariance.md), each matrix transforms as $\bm g_L(\cdots)\bm g_L^\dagger$, with

$$
\bm g_L=\exp\!\left[\frac{ia}{f_a}(\bm\delta-\bm\kappa)\right].
$$

To first order, the change in the tensor product is

$$
\begin{aligned}
\frac{ia}{f_a}\bigg(&
\left[\bm\delta-\bm\kappa,
\begin{pmatrix}0&0&1\\0&0&0\\0&0&0\end{pmatrix}\right]
\otimes\begin{pmatrix}0&0&0\\1&0&0\\0&0&0\end{pmatrix}\\
&+\begin{pmatrix}0&0&1\\0&0&0\\0&0&0\end{pmatrix}\otimes
\left[\bm\delta-\bm\kappa,
\begin{pmatrix}0&0&0\\1&0&0\\0&0&0\end{pmatrix}\right]\bigg).
\end{aligned}
$$

Both terms are required. For diagonal $\bm\delta-\bm\kappa$, the first is proportional to $(\delta-\kappa)_{uu}-(\delta-\kappa)_{ss}$ and the second to $(\delta-\kappa)_{dd}-(\delta-\kappa)_{uu}$. Their sum is $(\delta-\kappa)_{dd}-(\delta-\kappa)_{ss}$: the two up-quark contributions cancel.

## Octet contraction and baryon matching

One flavor-equivariant contraction of the two matrices gives

$$
\begin{pmatrix}0&0&0\\1&0&0\\0&0&0\end{pmatrix}
\begin{pmatrix}0&0&1\\0&0&0\\0&0&0\end{pmatrix}
=\begin{pmatrix}0&0&0\\0&0&1\\0&0&0\end{pmatrix}=\bm w_0.
$$

It is traceless and transforms as an octet matrix. Contracting the first-order result gives

$$
\bm w_0\to\bm w_0+\frac{ia}{f_a}[\bm\delta-\bm\kappa,\bm w_0]
+\mathcal O\!\left(\frac{a^2}{f_a^2}\right).
$$

For diagonal $\bm\delta-\bm\kappa$ this becomes

$$
\bm w_0\to\left[1+\frac{ia}{f_a}\bigl((\delta-\kappa)_{dd}-(\delta-\kappa)_{ss}\bigr)\right]\bm w_0.
$$

This calculation confirms the octet transformation law. The contraction is not a complete irreducible decomposition and does not determine the octet Wilson coefficient, its normalization, or the hadronic low-energy constants. Different operators, color structures, and octet components may have independent matching coefficients.

Salone:2026byq, Sec. II and App. A, also identifies a 27-plet; it includes $\Delta I=3/2$ contributions and requires separate chiral operators. The two-coupling Lagrangian of Borasoy:2000ki retains the octet component. At the baryon level it reads

$$
\mathcal L_W^{(0)}=d_w\langle\bar{\bm{\mathcal B}}\{\bm w_+,\bm{\mathcal B}\}\rangle
+f_w\langle\bar{\bm{\mathcal B}}[\bm w_+,\bm{\mathcal B}]\rangle,\qquad
\bm w_+=\bm u\bm w\bm u^\dagger+\bm u\bm w^\dagger\bm u^\dagger.
$$

Matching a general explicit ALP-dependent weak interaction requires its own operator coefficients and potentially independent baryon low-energy constants. The earlier [one-ALP baryon calculation](one-alp-baryon-couplings.md) assumes that the ALP dependence can be included in the same spurion with the same $d_w,f_w$.

## Other currents and limitations

Products involving leptonic currents give semileptonic operators when leptons are kept explicitly. Heavy-quark operators contribute through threshold matching and evolution when heavy fields are integrated out; operators with explicit heavy fields require a theory containing those fields. These contributions are not evaluated above. The illustrative construction in [partonic weak spurion invariance](partonic-weak-spurion-invariance.md) assumes that the second current is unaffected by the light-quark rotation; it is not by itself a matching calculation for the two light currents above.

## Reproduction

[nonleptonic_weak_octet_rotation.py](../../scripts/nonleptonic_weak_octet_rotation.py) expands both rotated current matrices, forms their tensor product, contracts the flavor indices, and evaluates the diagonal limit. Run `.venv/bin/python scripts/nonleptonic_weak_octet_rotation.py`; pass `--test` for the three checks.
