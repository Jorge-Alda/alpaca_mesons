---
title: Explicit finite-mixing neutral two-pseudoscalar amplitude
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Explicit finite-mixing neutral two-pseudoscalar amplitude

## Inputs and assumptions

Use natural units and the CCWZ Lagrangian, stationary aligned vacuum, generators and source derivatives of the [joint expansion](../shared/alp-lo-joint-vertex-expansion.md), following [@Alda:2026nrd]. The calculation is tree-level LO with no expansion in mixing. Put $x=a/f_a$, $\chi_n=2B_0\partial_x^n M(0)$ and $\chi_0=2B_0\operatorname{diag}(m_u,m_d,m_s)$. These quark masses are Lagrangian inputs, not physical mixed-pole masses. Current and topological derivatives are $C_{R,L,n}$ and $\Theta_n$. Higher derivatives are independent inputs. Include the $\chi_2$ correction to the $aa$ mass block from the joint expansion.

The labels $a,1,2$ denote chosen simple physical poles with squared masses $s_a,s_1,s_2$. A flavour $K^0$ is not generally a mass eigenstate. Initially all three poles have nonzero ALP components and invertible meson blocks. The parent alignment and exceptional poles are described in the [pole-aligned record](../shared/alp-lo-pole-aligned-tree-vertices.md). Momenta are incoming and conserve momentum; the result excludes the Feynman-rule factor $i$.

## Pole overlaps

The [operator inverse](../shared/alp-lo-finite-mixing-operator-inverse.md) extends to general source functions by using the matrix
$$
\frac{f_\pi}{f_a}\left[-2sC_{A,0}+2B_0P_1-\frac{2m_0^2}{3}\Theta_1\mathbf1\right]
$$
on the right-hand side of its meson operator equation. This matrix is also $\sum_iD_{ia}(s)\lambda_i$, where $D(s)=sK-M^2$. The component solution is explicitly
$$
\left[\sum_iw_i(s)\lambda_i\right]_{jk}
=\frac{f_\pi}{f_a}\frac{-2s(C_{A,0})_{jk}+2B_0(P_1)_{jk}}{B_0(m_j+m_k)-s}\qquad(j\ne k),
$$
$$
\left[\sum_iw_i(s)\lambda_i\right]_{jj}
=\frac{f_\pi/f_a}{2B_0m_j-s}\left[
-2s(C_{A,0})_{jj}+2B_0(P_1)_{jj}-\frac{2m_0^2}{3}\Theta_1
-\frac{m_0^2}{3}
\frac{\displaystyle\sum_k\frac{-2s(C_{A,0})_{kk}+2B_0(P_1)_{kk}-2m_0^2\Theta_1/3}{2B_0m_k-s}}
{\displaystyle 1+\frac{m_0^2}{3}\sum_k\frac{1}{2B_0m_k-s}}
\right].
$$
Finally $w_i(s)=\frac12\operatorname{tr}(\lambda_i\sum_jw_j(s)\lambda_j)$. These expressions require the displayed denominators and the trace denominator to be nonzero. A singularity in this representation can be removable; take a limit or solve the full operator equation in that case. Exceptional poles with zero ALP component require a separate eigenvector. The roots include the general-source $aa$ curvature. No denominator is identified with a physical meson mass.

The $K$-normalized ALP component at a pole is
$$
\left[1+\frac12\operatorname{tr}\left[\left(\sum_i w_i(s)\lambda_i-\frac{f_\pi}{f_a}(C_{R,0}-C_{L,0})\right)^2\right]\right]^{-1/2}.
$$
The positive square root fixes the state sign. The source-transformed, pole-aligned daughter directions have meson components $w_i(s_n)-w_i(s_a)$ and ALP component one, before this normalization. The exact symmetric chiral map has no quadratic coordinate term, so its linear term suffices for cubic vertices at a stationary vacuum.

For the neutral $uu,dd,ss,ds,sd$ block, the same normalization is explicitly
$$
\left[1+\frac12\sum_{j=u,d,s}\left(\left[\sum_iw_i(s)\lambda_i\right]_{jj}-\frac{2f_\pi}{f_a}(C_{A,0})_{jj}\right)^2
+\left|\left[\sum_iw_i(s)\lambda_i\right]_{ds}-\frac{2f_\pi}{f_a}(C_{A,0})_{ds}\right|^2\right]^{-1/2}.
$$
Hermiticity is assumed; use the two off-diagonal components separately if
complexified directions are used to represent crossed flavour fields.

## Explicit aligned cubic terms

The following trace formulas are intermediate steps in the derivation of the component expression above. The three terms are evaluated at $p_r^2=s_r$ and $p_a+p_1+p_2=0$. Their sum, multiplied by the three normalization factors above, is the physical three-pole amplitude. Primed source derivatives are evaluated at the origin after the exact chiral map in the [source-transformation record](../shared/alp-chiral-source-transformations.md). They are fixed explicitly by the original sources and parent overlap:
$$
\chi'_n=\sum_{k+l+m=n}\frac{n!}{k!l!m!}
\left(-\frac{i f_a}{2f_\pi}\right)^{k+m}
\left(\sum_iw_i(s_a)\lambda_i\right)^k\chi_l
\left(\sum_iw_i(s_a)\lambda_i\right)^m,
$$
$$
\begin{aligned}
C'_{R,n}={}&\sum_{k+l+m=n}\frac{n!}{k!l!m!}
\left(-\frac{i f_a}{2f_\pi}\right)^k
\left(\frac{i f_a}{2f_\pi}\right)^m
\left(\sum_iw_i(s_a)\lambda_i\right)^k C_{R,l}
\left(\sum_iw_i(s_a)\lambda_i\right)^m
-\delta_{n0}\frac{f_a}{2f_\pi}\sum_iw_i(s_a)\lambda_i,\\
C'_{L,n}={}&\sum_{k+l+m=n}\frac{n!}{k!l!m!}
\left(\frac{i f_a}{2f_\pi}\right)^k
\left(-\frac{i f_a}{2f_\pi}\right)^m
\left(\sum_iw_i(s_a)\lambda_i\right)^k C_{L,l}
\left(\sum_iw_i(s_a)\lambda_i\right)^m
+\delta_{n0}\frac{f_a}{2f_\pi}\sum_iw_i(s_a)\lambda_i.
\end{aligned}
$$
Here $k,l,m$ are nonnegative integers, $\Theta'_1=\Theta_1+\frac{f_a}{2f_\pi}\operatorname{tr}\sum_iw_i(s_a)\lambda_i$, and $\Theta'_2=\Theta_2$. Only $n\le3$ for $\chi'$ and $n\le1$ for the currents occur below. These are Taylor derivatives, including the factorials shown.
$$
\begin{aligned}
\mathcal M_{aPP}={}&\frac{i(p_a\cdot p_1-p_a\cdot p_2)}{4f_a}
\operatorname{tr}\!\left[\sum_i(w_i(s_1)-w_i(s_a))\lambda_i
 [C'_{R,0}+C'_{L,0},\sum_j(w_j(s_2)-w_j(s_a))\lambda_j]\right]\\
&-\frac{1}{8f_a}\operatorname{tr}\!\left[(\chi'_1+\chi_1'^\dagger)
 \left\{\sum_i(w_i(s_1)-w_i(s_a))\lambda_i,
 \sum_j(w_j(s_2)-w_j(s_a))\lambda_j\right\}\right].
\end{aligned}
$$
$$
\begin{aligned}
\mathcal M_{aaP}={}&-\frac{i f_\pi}{2f_a^2}\operatorname{tr}\!\left[(C'_{R,0}-C'_{L,0})
 \left(p_a\cdot p_2[C'_{R,0}+C'_{L,0},\sum_i(w_i(s_1)-w_i(s_a))\lambda_i]
 +p_a\cdot p_1[C'_{R,0}+C'_{L,0},\sum_i(w_i(s_2)-w_i(s_a))\lambda_i]\right)\right]\\
&+\frac{f_\pi}{2f_a^2}\operatorname{tr}\!\left[\left((p_a\cdot p_1+p_1\cdot p_2)\sum_i(w_i(s_1)-w_i(s_a))\lambda_i
 +(p_a\cdot p_2+p_1\cdot p_2)\sum_i(w_i(s_2)-w_i(s_a))\lambda_i\right)(C'_{R,1}-C'_{L,1})\right]\\
&+\frac{i f_\pi}{4f_a^2}\operatorname{tr}\!\left[\sum_i(w_i(s_1)+w_i(s_2)-2w_i(s_a))\lambda_i(\chi_2'^\dagger-\chi'_2)\right]\\
&-\frac{f_\pi m_0^2\Theta'_2}{3f_a^2}\operatorname{tr}\!\left[\sum_i(w_i(s_1)+w_i(s_2)-2w_i(s_a))\lambda_i\right].
\end{aligned}
$$
$$
\mathcal M_{aaa}=
-\frac{f_\pi^2(p_a\cdot p_1+p_a\cdot p_2+p_1\cdot p_2)}{f_a^3}
\operatorname{tr}[(C'_{R,0}-C'_{L,0})(C'_{R,1}-C'_{L,1})]
+\frac{f_\pi^2}{4f_a^3}\operatorname{tr}(\chi'_3+\chi_3'^\dagger)
-\frac{2f_\pi^2m_0^2}{f_a^3}\Theta'_1\Theta'_2.
$$
The three scalar products are $(s_2-s_a-s_1)/2$, $(s_1-s_a-s_2)/2$ and $(s_a-s_1-s_2)/2$, respectively. These are pole squared masses used in kinematics. The terms include higher source derivatives and both daughter ALP components.

## Explicit decoupled-kaon component

For the exact zero-ALP kaon pole of the [neutral-block record](../shared/alp-lo-neutral-block-reduction.md), choose its normalized components as $(-D_{7a},D_{6a})/\sqrt{D_{6a}^2+D_{7a}^2}$ in the $(P_6,P_7)$ basis. The entries are evaluated at its own physical pole. The normalization follows from the unit meson kinetic block. These entries are explicitly determined by the source matrix above: $D_{6a}-iD_{7a}$ is its $ds$ element. This real pole is not, in general, the flavour state $K^0$.

This direction requires $D_{6a}^2+D_{7a}^2>0$. If both entries vanish, choose a basis in the decoupled degenerate kaon subspace; the displayed normalized direction is undefined. Let the kaon pole be daughter 1 and let daughter 2 be a pion-like pole with nonzero ALP component. After the two nonzero-ALP normalizations at $s_a$ and $s_2$, the explicit reduced amplitude is
$$
\begin{aligned}
&\frac{-D_{7a}-iD_{6a}}{\sqrt{D_{6a}^2+D_{7a}^2}}\Bigg\{
\frac{i(s_2-s_1)}{4f_a}
\sum_{k=u,d,s}\sum_j(w_j(s_2)-w_j(s_a))
 \left[(C'_{R,0}+C'_{L,0})_{sk}(\lambda_j)_{kd}
 -(\lambda_j)_{sk}(C'_{R,0}+C'_{L,0})_{kd}\right]\\
&-\frac{1}{8f_a}\sum_{k=u,d,s}\sum_j(w_j(s_2)-w_j(s_a))
 \left[(\lambda_j)_{sk}(\chi'_1+\chi_1'^\dagger)_{kd}
 +(\chi'_1+\chi_1'^\dagger)_{sk}(\lambda_j)_{kd}\right]\\
&-\frac{i f_\pi(s_1-s_a-s_2)}{4f_a^2}
\sum_{k=u,d,s}\left[(C'_{R,0}-C'_{L,0})_{sk}(C'_{R,0}+C'_{L,0})_{kd}
 -(C'_{R,0}+C'_{L,0})_{sk}(C'_{R,0}-C'_{L,0})_{kd}\right]\\
&-\frac{f_\pi s_1}{2f_a^2}
 (C'_{R,1}-C'_{L,1})_{sd}
 +\frac{i f_\pi}{4f_a^2}(\chi_2'^\dagger-\chi'_2)_{sd}
\Bigg\}
\quad+\quad\text{the same braces with }d\leftrightarrow s
\text{ and complex-conjugate prefactor}.
\end{aligned}
$$
The indicated second contribution interchanges $d$ and $s$ in the matrix entries and uses the complex-conjugate prefactor; it does not conjugate the coefficients outside that prefactor. $D_{6a}-iD_{7a}$ is the $ds$ entry of the explicitly displayed general-source matrix, since $(\lambda_6)_{ds}=1$ and $(\lambda_7)_{ds}=-i$. The prefactor includes the kaon's normalization. Multiply the displayed reduced amplitude by the parent and daughter-2 factors from the previous section. The $aaa$ term vanishes because the kaon has no ALP component. The singlet $aaP$ term vanishes because this kaon direction is traceless. The scalar products use the three physical pole masses in the kinematic identities above; no unmixed-mass substitution is used.

## Verification and limitations

The executable [derivation](../../../scripts/derive_alp_explicit_neutral_amplitude.py) evaluates the three cubic terms, applies the general-source operator inverse and performs the vector/axial and scalar/pseudoscalar substitution. The [unit tests](../../../scripts/test_derive_alp_explicit_neutral_amplitude.py) check them against the independent CCWZ generator, including noncommuting complex sources, the general-source quadratic inverse, Bose exchange and a complex Hermitian vector source supported only on the $d$–$s$ entries. A further test obtains the mixed parent and daughters from the same source-derived $U(3)$ quadratic matrices, and compares the original and parent-aligned amplitudes at the corresponding on-shell scalar products. This is an algebraic on-shell comparison; no claim is made that the chosen masses admit real decay phase space. The zero-ALP branch is tested algebraically. Widths and weak kaon dynamics are outside this calculation.

Run from the repository root:

```sh
.venv/bin/python scripts/derive_alp_explicit_neutral_amplitude.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_explicit_neutral_amplitude.py
```
