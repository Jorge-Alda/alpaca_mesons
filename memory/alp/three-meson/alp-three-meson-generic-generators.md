---
title: Generic-generator expansion of one-ALP three-meson interactions
status: derived
created: 2026-09-22
updated: 2026-09-22
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Generic-Generator Expansion of One-ALP Three-Meson Interactions

## Input and conventions

The starting result is the derived leading-order one-ALP, three-pseudoscalar Lagrangian in [the matrix-level derivation](alp-phi3-lo-chipt-derivation.md), whose CCWZ inputs originate from [the Appendix A extract](../../extracts/Alda:2026nrd/appendix-a-expansion.md) [@Alda:2026nrd]:

$$
\left.\mathcal L\right|_{a\Phi^3}
=K\,\partial_\mu a\,
\left\langle c_q^A[\Phi,[\Phi,\partial^\mu\Phi]]\right\rangle
+M\,a\,\left\langle\Phi^3P\right\rangle,
$$

where

$$
K=\frac{1}{6f_af_\pi},\qquad
M=-\frac{B_0}{6f_af_\pi},\qquad
P=\{c_q^P,m_q\}.
$$

Expand the pseudoscalar matrix as

$$
\Phi=\sum_r\phi_r\lambda_r,
\qquad
\partial^\mu\Phi=\sum_r(\partial^\mu\phi_r)\lambda_r.
$$

Assumptions and conventions:

- The component fields $\phi_r$ are commuting Lorentz scalars.
- The flavour generators $\lambda_r$ are generic noncommuting matrices.
- No normalization, Hermiticity, tracelessness, or particular flavour group is assumed.
- The expressions are Lagrangian coefficients. Factorials from functional differentiation into vertex rules are not included.
- Flavour brackets $\langle\cdots\rangle$ denote a trace.

The executable derivation is [`scripts/derive_alp_three_meson_generators.py`](../../../scripts/derive_alp_three_meson_generators.py). It substitutes three abstract matrix generators into both matrix traces and uses SymPy polynomial coefficient extraction to isolate each scalar-field monomial.

## Ordered generator structures

For later use, define

$$
T_{iij}^{A}
=\left\langle
c_q^A[\lambda_i,[\lambda_i,\lambda_j]]
\right\rangle,
$$

$$
T_{ij;k}^{A}
=\left\langle c_q^A\left(
[\lambda_i,[\lambda_j,\lambda_k]]
+[\lambda_j,[\lambda_i,\lambda_k]]
\right)\right\rangle,
$$

where $T_{ij;k}^{A}=T_{ji;k}^{A}$, and

$$
T_{iij}^{P}
=\left\langle
(\lambda_i^2\lambda_j+\lambda_i\lambda_j\lambda_i
+\lambda_j\lambda_i^2)P
\right\rangle,
$$

$$
T_{ijk}^{P}
=\left\langle
\sum_{\sigma\in S_3}
\lambda_{\sigma(i)}\lambda_{\sigma(j)}\lambda_{\sigma(k)}P
\right\rangle.
$$

The sums over orderings in $T^P$ must not be replaced by ordinary numerical multiplicities unless additional generator commutation relations justify it.

## Three identical mesons

For $\Phi=\phi_i\lambda_i$, the nested commutator vanishes because $[\lambda_i,\lambda_i]=0$. Therefore

$$
\boxed{
\mathcal L_{a\phi_i\phi_i\phi_i}
=M\,a\,\phi_i^3
\left\langle\lambda_i^3P\right\rangle.
}
$$

There is no derivative contribution from the kinetic term.

## Two identical and one different meson

Keep the monomials containing two fields $\phi_i$ and one field $\phi_j$, with $i\ne j$. Direct expansion gives

$$
[\Phi,[\Phi,\partial^\mu\Phi]]\big|_{iij}
=[\lambda_i,[\lambda_i,\lambda_j]]
\left(
\phi_i^2\partial^\mu\phi_j
-\phi_i\phi_j\partial^\mu\phi_i
\right),
$$

and

$$
\Phi^3\big|_{iij}
=\phi_i^2\phi_j
(\lambda_i^2\lambda_j+\lambda_i\lambda_j\lambda_i
+\lambda_j\lambda_i^2).
$$

Thus

$$
\boxed{
\mathcal L_{a\phi_i\phi_i\phi_j}
=K\,\partial_\mu a\,T_{iij}^{A}
\left(
\phi_i^2\partial^\mu\phi_j
-\phi_i\phi_j\partial^\mu\phi_i
\right)
+M\,a\,\phi_i^2\phi_j\,T_{iij}^{P},
\qquad i\ne j.
}
$$

If there is no flavour violation, the only processes allowed by charge and strangeness conservation with two identical mesons are those where all mesons have a diagonal generator ($\pi^0$, $\eta$ or $\eta^\prime$), but in that case $T_{iij}^A=0$. The axial part only contributes in $|\Delta S|=1$ decays, like $a \to \pi^0 \pi^0 K^0$.

## Three different mesons

For pairwise distinct $i,j,k$, group the kinetic expansion by the differentiated field:

$$
\left.[\Phi,[\Phi,\partial^\mu\Phi]]\right|_{ijk}
=\phi_i\phi_j(\partial^\mu\phi_k)\,
\mathcal K_{ij;k}
+\phi_i\phi_k(\partial^\mu\phi_j)\,
\mathcal K_{ik;j}
+\phi_j\phi_k(\partial^\mu\phi_i)\,
\mathcal K_{jk;i},
$$

where

$$
\mathcal K_{ij;k}
=[\lambda_i,[\lambda_j,\lambda_k]]
+[\lambda_j,[\lambda_i,\lambda_k]].
$$

The symmetry-breaking expansion contains all six generator orderings:

$$
\left.\Phi^3\right|_{ijk}
=\phi_i\phi_j\phi_k
\sum_{\sigma\in S_3}
\lambda_{\sigma(i)}\lambda_{\sigma(j)}\lambda_{\sigma(k)}.
$$

The Lagrangian is therefore

$$
\boxed{
\begin{aligned}
\mathcal L_{a\phi_i\phi_j\phi_k}
={}&K\,\partial_\mu a\left[
\phi_i\phi_j(\partial^\mu\phi_k)T_{ij;k}^{A}
+\phi_i\phi_k(\partial^\mu\phi_j)T_{ik;j}^{A}
+\phi_j\phi_k(\partial^\mu\phi_i)T_{jk;i}^{A}
\right]\\
&+M\,a\,\phi_i\phi_j\phi_k\,T_{ijk}^{P},
\qquad i,j,k\ \text{pairwise distinct}.
\end{aligned}
}
$$

## No-mixing parametrization

The same formulas apply after the replacements from the matrix-level derivation,

$$
c_q^A\longrightarrow-\mathbb M_a^{-1}(C),
\qquad
P\longrightarrow
-\{c_q^{AP}+\mathbb M_a^{-1}(C),m_q\}.
$$

Equivalently, these substitutions reproduce the signs in the no-mixing matrix Lagrangian before the generator expansion.
