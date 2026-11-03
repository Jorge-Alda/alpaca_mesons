---
title: Finite-mixing LO amplitudes for ALP to three pions
status: derived
created: 2026-11-03
updated: 2026-11-03
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Finite-mixing LO amplitudes for $a\to3\pi$

## Inputs and conventions

This derivation uses Alda:2026nrd's CCWZ building blocks and the notation of
the [small-mixing three-meson amplitude](alp-three-meson-amplitudes-widths.md)
[@Alda:2026nrd]. Flavor matrices are general and need not commute or be
diagonal. Hermiticity is required for a real Lagrangian. Use natural units,
the mostly-minus metric, $S=1+iT$,
$\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$,
$\phi=\Phi/(2f_\pi)$, and $\chi_0=2B_0m_q$.

The starting LO terms are

$$
\mathcal L_{\rm LO}=\frac{f_\pi^2}{4}
\operatorname{tr}(u_\mu u^\mu+\chi_+)+\mathcal L_{TC},
\qquad
\mathcal L_{TC}=-\frac{f_\pi^2m_0^2}{4N_f}X_{TC}^2.
$$

The ALP mass is included in the quadratic pole matrix. The amplitudes below
are the tree contributions from the LO chiral terms. No numerical source
matrices or pole spectrum are chosen.

## 1. Invariant couplings and transformed CCWZ building blocks

Alda:2026nrd's invariant coupling combinations are

$$
c_q^{AP}=c_q^A-c_q^P,\qquad
c_q^{VS}=c_q^V-c_q^S,\qquad
c_{GA}=c_G+\operatorname{tr}c_q^A.
$$

The invariant matrix entering the ALP--meson mass block is

$$
\mathcal C=-\frac12\{\chi_0,c_q^{AP}\}
-\frac{m_0^2}{N_f}c_{GA}\mathbf1
=-B_0\{m_q,c_q^{AP}\}
-\frac{m_0^2}{N_f}c_{GA}\mathbf1.
$$

Make the finite chiral field redefinition that removes the ALP derivative
sources. Transform the external sources and anomalous singlet source with
the meson fields. Then

$$
u_\mu'=h u_\mu h^\dagger,\qquad
\chi_\pm'=h\chi_\pm h^\dagger,\qquad X_{TC}'=X_{TC},
$$

so all traces in $\mathcal L_{\rm LO}$ stay unchanged. The result below uses
only these CCWZ building blocks.

For Alda's linear source convention, at zero meson field,

$$
\left.\chi_+\right|_{\Phi=0}
=2\chi_0+\frac{2ia}{f_a}[c_q^S,\chi_0],
\qquad
\left.\chi_-\right|_{\Phi=0}
=\frac{2ia}{f_a}\{c_q^P,\chi_0\}.
$$

The finite transformed building blocks at the same point are

$$
\begin{aligned}
\chi_+'(a,0)=\frac12\Big[&
e^{-ia(c_q^V+c_q^A)/f_a}(\chi_++\chi_-)
e^{ia(c_q^V-c_q^A)/f_a}\\
&+e^{-ia(c_q^V-c_q^A)/f_a}(\chi_+-\chi_-)
e^{ia(c_q^V+c_q^A)/f_a}\Big]_{\Phi=0},
\end{aligned}
$$

$$
\begin{aligned}
\chi_-'(a,0)=\frac12\Big[&
e^{-ia(c_q^V+c_q^A)/f_a}(\chi_++\chi_-)
e^{ia(c_q^V-c_q^A)/f_a}\\
&-e^{-ia(c_q^V-c_q^A)/f_a}(\chi_+-\chi_-)
e^{ia(c_q^V+c_q^A)/f_a}\Big]_{\Phi=0}.
\end{aligned}
$$

Their first derivatives give the invariant one-ALP sources,

$$
\left.\partial_a\chi_+'\right|_0
=-\frac{4iB_0}{f_a}[c_q^{VS},m_q],\qquad
\left.\partial_a\chi_-'\right|_0
=-\frac{4iB_0}{f_a}\{c_q^{AP},m_q\}.
$$

For the linear source convention, every higher derivative needed below is
fixed by the exact ordered products. Explicitly, for $n=1,2,3,4$,

$$
\begin{aligned}
\left.\partial_a^n\chi_\pm'\right|_0={}&
\frac12\sum_{r+s=n}\frac{n!}{r!s!}
\left(-\frac{i(c_q^V+c_q^A)}{f_a}\right)^r
(2\chi_0)
\left(\frac{i(c_q^V-c_q^A)}{f_a}\right)^s\\
&+\frac12\sum_{r+s=n-1}\frac{n!}{r!s!}
\left(-\frac{i(c_q^V+c_q^A)}{f_a}\right)^r
\frac{2i}{f_a}\bigl([c_q^S,\chi_0]+\{c_q^P,\chi_0\}\bigr)
\left(\frac{i(c_q^V-c_q^A)}{f_a}\right)^s\\
&\mathbin{\pm}\frac12\sum_{r+s=n}\frac{n!}{r!s!}
\left(-\frac{i(c_q^V-c_q^A)}{f_a}\right)^r
(2\chi_0)
\left(\frac{i(c_q^V+c_q^A)}{f_a}\right)^s\\
&\mathbin{\pm}\frac12\sum_{r+s=n-1}\frac{n!}{r!s!}
\left(-\frac{i(c_q^V-c_q^A)}{f_a}\right)^r
\frac{2i}{f_a}\bigl([c_q^S,\chi_0]-\{c_q^P,\chi_0\}\bigr)
\left(\frac{i(c_q^V+c_q^A)}{f_a}\right)^s.
\end{aligned}
$$

The upper sign is for $\chi_+'$, the lower sign for $\chi_-'$; an empty sum
is zero. Matrix order is as displayed. The higher derivatives are therefore
fixed by the original source convention, not independent couplings.

In these coordinates,

$$
X_{TC}'=-\frac{\operatorname{tr}\Phi}{f_\pi}
-2c_{GA}\frac{a}{f_a},
\qquad K_{IJ}=\delta_{IJ}.
$$

## 2. Pole components from Alda's quadratic matrices

Use Alda's field vector $\Phi_I=(P_i,a)$ and

$$
\mathcal L_2=\frac12\partial_\mu\Phi_I K_{IJ}\partial^\mu\Phi_J
-\frac12\Phi_I M^2_{IJ}\Phi_J.
$$

In the derivative-source-free CCWZ basis,

$$
M^2_{ij}=\frac14\operatorname{tr}\!\left(
\{\lambda_i,\lambda_j\}\chi_0\right)
+\frac{m_0^2}{2N_f}\operatorname{tr}(\lambda_i)\operatorname{tr}(\lambda_j),
\qquad
M^2_{i\hat a}=-\frac{f_\pi}{f_a}\operatorname{tr}(\mathcal C\lambda_i),
$$

$$
M^2_{\hat a\hat a}=m_{a,0}^2
-\frac{f_\pi^2}{4}\operatorname{tr}\!\left(
\left.\partial_a^2\chi_+'\right|_0\right)
+\frac{2f_\pi^2m_0^2c_{GA}^2}{N_ff_a^2}.
$$

For an ALP-overlap pole at $s=m_n^2$, use Alda's convention
$w_{\hat a}^{(n)}=1$. Then

$$
w_i^{(n)}=\frac{f_\pi}{f_a}
\operatorname{tr}\!\left[\lambda_i\mathbb M_a^{-1}(\mathcal C)\right]
\bigg|_{m_a^2=m_n^2},
$$

$$
K_{\hat a\hat a}+2K_{i\hat a}w_i^{(n)}+K_{ij}w_i^{(n)}w_j^{(n)}
=1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}\!\left[\bigl(\mathbb M_a^{-1}(\mathcal C)\bigr)^2\right]
\bigg|_{m_a^2=m_n^2}.
$$

Here Alda's pole operator is

$$
\mathbb M_a(X)=B_0\{m_q,X\}
+\frac{m_0^2}{N_f}\operatorname{tr}(X)\mathbf1-m_a^2X.
$$

The pole equation is

$$
M^2_{\hat a\hat a}-m_n^2
-\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}\!\left[\mathcal C\mathbb M_a^{-1}(\mathcal C)\right]
\bigg|_{m_a^2=m_n^2}=0.
$$

In the quark-mass basis, $m_q$ is diagonal; no coupling matrix is assumed
diagonal. The matrix inverse needed in these expressions is

$$
\begin{aligned}
\mathbb M_a^{-1}(\mathcal C)={}&
\int_0^\infty d\tau\,
e^{-\tau(B_0m_q-\frac12m_a^2\mathbf1)}\mathcal C
e^{-\tau(B_0m_q-\frac12m_a^2\mathbf1)}\\
&-\frac{m_0^2}{N_f}\,
\frac{\operatorname{tr}\!\left[\int_0^\infty d\tau\,
e^{-\tau(B_0m_q-\frac12m_a^2\mathbf1)}\mathcal C
e^{-\tau(B_0m_q-\frac12m_a^2\mathbf1)}\right]}
{1+\frac{m_0^2}{N_f}\operatorname{tr}
[(2B_0m_q-m_a^2\mathbf1)^{-1}]}
(2B_0m_q-m_a^2\mathbf1)^{-1}.
\end{aligned}
$$

The integral is defined where it converges and continued to the pole. A
neutral pure-meson pole with zero ALP overlap instead obeys
$w_{\hat a}^{(n)}=0$,
$M^2_{ij}w_j^{(n)}=m_n^2K_{ij}w_j^{(n)}$,
$\operatorname{tr}(\mathcal C\lambda_i)w_i^{(n)}=0$, and
$K_{ij}w_i^{(n)}w_j^{(n)}=1$. At a singular pole operator or a degeneracy,
solve the complete generalized eigenproblem. The physical charged-pion poles
are $(\lambda_{\pi^\pm},0)$:
electric charge conservation forbids bilinear mixing of a neutral ALP with
either charged pion.

For compactness in the amplitudes only, define the normalized neutral pole
column

$$
\mathcal P(s)=\left(
\frac{\frac{2f_\pi}{f_a}\mathbb M_a^{-1}(\mathcal C)|_{m_a^2=s}}
{\sqrt{1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}[(\mathbb M_a^{-1}(\mathcal C)|_{m_a^2=s})^2]}},
\frac{1}{\sqrt{1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}[(\mathbb M_a^{-1}(\mathcal C)|_{m_a^2=s})^2]}}
\right).
$$

It is the pair consisting of the meson matrix and ALP component, after LSZ
normalization. It introduces no new coupling or overlap input.

## 3. Explicit cubic and quartic vertices

The cubic and quartic interactions through the required field orders are

$$
\begin{aligned}
\mathcal L_3={}&
-\frac a8\operatorname{tr}\!\left[
\left.\partial_a\chi_+'\right|_0\Phi^2\right]
+\frac{f_\pi a^2}{8i}\operatorname{tr}\!\left[
\left.\partial_a^2\chi_-'\right|_0\Phi\right]
+\frac{f_\pi^2a^3}{24}\operatorname{tr}\!\left[
\left.\partial_a^3\chi_+'\right|_0\right],\\
\mathcal L_4={}&
\frac{\operatorname{tr}([\Phi,\partial_\mu\Phi]
[\Phi,\partial^\mu\Phi])}{48f_\pi^2}
+\frac{\operatorname{tr}\!\left[\left.\chi_+'\right|_0\Phi^4\right]}{96f_\pi^2}
-\frac{a}{24if_\pi}\operatorname{tr}\!\left[
\left.\partial_a\chi_-'\right|_0\Phi^3\right]
-\frac{a^2}{16}\operatorname{tr}\!\left[
\left.\partial_a^2\chi_+'\right|_0\Phi^2\right]\\
&+\frac{f_\pi a^3}{24i}\operatorname{tr}\!\left[
\left.\partial_a^3\chi_-'\right|_0\Phi\right]
+\frac{f_\pi^2a^4}{96}\operatorname{tr}\!\left[
\left.\partial_a^4\chi_+'\right|_0\right].
\end{aligned}
$$

For pole directions $(X,x),(Y,y),(Z,z)$, the cubic vertex is the equation

$$
\boxed{\begin{aligned}
\mathcal V^{(3)}[(X,x),(Y,y),(Z,z)]={}&
-\frac18\operatorname{tr}\!\left[
\left.\partial_a\chi_+'\right|_0
\bigl(x\{Y,Z\}+y\{X,Z\}+z\{X,Y\}\bigr)\right]\\
&+\frac{f_\pi}{4i}\operatorname{tr}\!\left[
\left.\partial_a^2\chi_-'\right|_0(xyZ+xzY+yzX)\right]
+\frac{f_\pi^2xyz}{4}\operatorname{tr}\!\left[
\left.\partial_a^3\chi_+'\right|_0\right].
\end{aligned}}
$$

These cubic vertices are not identically zero. For the charge-allowed
Hermitian source matrix
$c_q^{VS}=|d\rangle\langle s|+|s\rangle\langle d|$ and $m_s\ne m_d$,
take $\lambda_{K^-}=\sqrt2|s\rangle\langle u|$ and
$\lambda_{K^+}=\lambda_{K^-}^\dagger$. The two vertices are

$$
\begin{aligned}
\mathcal V^{(3)}[(0,1),(\lambda_{\pi^+},0),(\lambda_{K^-},0)]
&=\frac{iB_0}{2f_a}\operatorname{tr}\!\left(
[c_q^{VS},m_q]\{\lambda_{\pi^+},\lambda_{K^-}\}\right)
=\frac{iB_0(m_s-m_d)}{f_a},\\
\mathcal V^{(3)}\!\left[\left(0,
\frac{1}{\sqrt{1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}[(\mathbb M_a^{-1}(\mathcal C))^2]_{m_{\pi^0}^2}}}\right),
(\lambda_{\pi^-},0),(\lambda_{K^+},0)\right]
&=\frac{iB_0}{2f_a\sqrt{1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}[(\mathbb M_a^{-1}(\mathcal C))^2]_{m_{\pi^0}^2}}}
\operatorname{tr}\!\left(
[c_q^{VS},m_q]\{\lambda_{\pi^-},\lambda_{K^+}\}\right)\\
&=-\frac{iB_0(m_s-m_d)}{f_a\sqrt{1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}[(\mathbb M_a^{-1}(\mathcal C))^2]_{m_{\pi^0}^2}}}.
\end{aligned}
$$

They give a $K^\pm$ exchange contribution in the $t$ channel,

$$
\frac{B_0^2(m_s-m_d)^2}
{f_a^2(m_K^2-t)\sqrt{1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}[(\mathbb M_a^{-1}(\mathcal C))^2]_{m_{\pi^0}^2}}},
$$

before the incoming-ALP LSZ factor. It vanishes if the pion has no
ALP component. This explicit allowed example shows why the exchange is
generically present at finite mixing, while other source choices can make it
zero.

For four directions $(X_r,x_r)$ and incoming momenta $k_r$, the tensor
vertex is

$$
\boxed{\mathcal V^{(4)}=\mathcal V^{(4)}_{\rm kin}
+\mathcal V^{(4)}_{\chi_+}.}
$$

Its kinetic part is

$$
\boxed{\mathcal V^{(4)}_{\rm kin}=
-\frac{1}{48f_\pi^2}\sum_{\sigma\in S_4}
(k_{\sigma_2}\!\cdot k_{\sigma_4})
\operatorname{tr}\!\left([X_{\sigma_1},X_{\sigma_2}]
[X_{\sigma_3},X_{\sigma_4}]\right).}
$$

Its mass-source part is

$$
\boxed{\begin{aligned}
\mathcal V^{(4)}_{\chi_+}={}&
\frac{1}{96f_\pi^2}\sum_{\sigma\in S_4}
\operatorname{tr}\!\left[\left.\chi_+'\right|_0
X_{\sigma_1}X_{\sigma_2}X_{\sigma_3}X_{\sigma_4}\right]\\
&-\frac{1}{24if_\pi}\sum_{r=1}^4x_r
\sum_{\sigma\in S_3^{(r)}}\operatorname{tr}\!\left[
\left.\partial_a\chi_-'\right|_0
X_{\sigma_1}X_{\sigma_2}X_{\sigma_3}\right]\\
&-\frac18\sum_{r<s}x_rx_s
\operatorname{tr}\!\left[\left.\partial_a^2\chi_+'\right|_0
\{X_k,X_l\}\right]\\
&+\frac{f_\pi}{4i}\sum_{r=1}^4
\left(\prod_{j\ne r}x_j\right)
\operatorname{tr}\!\left[\left.\partial_a^3\chi_-'\right|_0X_r\right]
+\frac{f_\pi^2}{4}\left(\prod_{r=1}^4x_r\right)
\operatorname{tr}\!\left[\left.\partial_a^4\chi_+'\right|_0\right].
\end{aligned}}
$$

Here $S_4$ is the set of all four-leg orderings, $S_3^{(r)}$ the six
orderings excluding leg $r$, and $k,l$ the legs left after choosing $r,s$.
The components $\mathcal V^{(3)}_{IJK}$ and $\mathcal V^{(4)}_{IJKL}$ are
these equations with meson directions $(\lambda_i,0)$ and ALP direction
$(0,1)$. The Feynman rules are $i\mathcal V^{(3)}$ and $i\mathcal V^{(4)}$.

## 4. The two amplitudes and the three exchanges

For $a(p)\to\pi^0(p_0)\pi^+(p_+)\pi^-(p_-)$ let

$$
s=(p_++p_-)^2,\qquad t=(p_0+p_-)^2,\qquad
u=(p_0+p_+)^2.
$$

For $a(p)\to\pi^0(p_1)\pi^0(p_2)\pi^0(p_3)$ let
$s=(p_2+p_3)^2$, $t=(p_1+p_3)^2$, $u=(p_1+p_2)^2$.
Use $\lambda_{\pi^0}=\operatorname{diag}(1,-1,0)$,
$\lambda_{\pi^+}=\sqrt2\,|u\rangle\langle d|$, and
$\lambda_{\pi^-}=\lambda_{\pi^+}^\dagger$.

The two amplitudes are

$$
\boxed{\begin{aligned}
\mathcal M_{000}={}&
\mathcal V^{(4)}[\mathcal P(m_a^2),\mathcal P(m_{\pi^0}^2),
\mathcal P(m_{\pi^0}^2),\mathcal P(m_{\pi^0}^2)]\\
&+\mathcal V^{(3)}[\mathcal P(m_a^2),\mathcal P(m_{\pi^0}^2),I]
(M^2-sK)^{-1}_{IJ}\mathcal V^{(3)}[J,\mathcal P(m_{\pi^0}^2),\mathcal P(m_{\pi^0}^2)]\\
&+\mathcal V^{(3)}[\mathcal P(m_a^2),\mathcal P(m_{\pi^0}^2),I]
(M^2-tK)^{-1}_{IJ}\mathcal V^{(3)}[J,\mathcal P(m_{\pi^0}^2),\mathcal P(m_{\pi^0}^2)]\\
&+\mathcal V^{(3)}[\mathcal P(m_a^2),\mathcal P(m_{\pi^0}^2),I]
(M^2-uK)^{-1}_{IJ}\mathcal V^{(3)}[J,\mathcal P(m_{\pi^0}^2),\mathcal P(m_{\pi^0}^2)].
\end{aligned}}
$$

$$
\boxed{\begin{aligned}
\mathcal M_{0+-}={}&
\mathcal V^{(4)}[\mathcal P(m_a^2),\mathcal P(m_{\pi^0}^2),
(\lambda_{\pi^+},0),(\lambda_{\pi^-},0)]\\
&+\mathcal V^{(3)}[\mathcal P(m_a^2),\mathcal P(m_{\pi^0}^2),I]
(M^2-sK)^{-1}_{IJ}\mathcal V^{(3)}[J,(\lambda_{\pi^+},0),(\lambda_{\pi^-},0)]\\
&+\mathcal V^{(3)}[\mathcal P(m_a^2),(\lambda_{\pi^+},0),I]
(M^2-tK)^{-1}_{IJ}\mathcal V^{(3)}[J,\mathcal P(m_{\pi^0}^2),(\lambda_{\pi^-},0)]\\
&+\mathcal V^{(3)}[\mathcal P(m_a^2),(\lambda_{\pi^-},0),I]
(M^2-uK)^{-1}_{IJ}\mathcal V^{(3)}[J,\mathcal P(m_{\pi^0}^2),(\lambda_{\pi^+},0)].
\end{aligned}}
$$

The internal indices run over all meson and ALP fields. The three pairings
are $(a,0)|(0,0)$ at $s,t,u$ and $(a,0)|(+,-)$,
$(a,+)|(0,-)$, $(a,-)|(0,+)$, respectively. The charged poles are exactly
$(\lambda_{\pi^\pm},0)$: the neutral ALP cannot mix with a state of nonzero
electric charge. If the neutral pion has zero ALP overlap, replace
$\mathcal P(m_{\pi^0}^2)$ by its pure-meson pole column from the meson block.

The contact is $\mathcal V^{(4)}_{\rm kin}+\mathcal V^{(4)}_{\chi_+}$, not
just the kinetic term. For three identical neutral pole matrices,
$\mathcal V^{(4)}_{\rm kin}=0$; for $\pi^0\pi^+\pi^-$ it is the explicit
commutator trace in the displayed $\mathcal V^{(4)}_{\rm kin}$ equation.

## 5. Explicit comparison with $f_a\gg f_\pi$

At first order in the ALP--meson mixing, Alda's invariant no-mixing
coefficients are

$$
\widehat c_q^A=-\mathbb M_a^{-1}(\mathcal C),\qquad
\widehat c_q^P=-c_q^{AP}-\mathbb M_a^{-1}(\mathcal C),
$$

and the small-mixing amplitudes from the linked record are

$$
\mathcal M_{000}^{\rm linear}
=6M\operatorname{tr}\!\left[
\lambda_{\pi^0}^3\{\widehat c_q^P,m_q\}\right],
\qquad M=-\frac{B_0}{6f_af_\pi},
$$

$$
\begin{aligned}
\mathcal M_{0+-}^{\rm linear}={}&K\left[
T_{+-,0}^A\,p\!\cdot p_0+
T_{0-,+}^A\,p\!\cdot p_++
T_{0+,-}^A\,p\!\cdot p_-\right]+M T_{0+-}^P,\\
K={}&\frac{1}{6f_af_\pi},\\
T_{ij;k}^A={}&\operatorname{tr}\!\left(\widehat c_q^A
([\lambda_i,[\lambda_j,\lambda_k]]+
[\lambda_j,[\lambda_i,\lambda_k]])\right),\\
T_{0+-}^P={}&\operatorname{tr}\!\left[
\left(\sum_{\sigma\in S_3}
\lambda_{\sigma(0)}\lambda_{\sigma(+)}\lambda_{\sigma(-)}\right)
\{\widehat c_q^P,m_q\}\right].
\end{aligned}
$$

The finite pole differs from its linear incoming column by its LSZ
normalization. The finite neutral-pion column differs from the reference
$(\lambda_{\pi^0},0)$ by

$$
\mathcal P(m_{\pi^0}^2)-(\lambda_{\pi^0},0).
$$

Thus, for the charged channel, the exact difference is

$$
\boxed{\begin{aligned}
\mathcal M_{0+-}-\mathcal M_{0+-}^{\rm linear}={}&
\left[\frac{1}{\sqrt{1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}[(\mathbb M_a^{-1}(\mathcal C))^2]_{m_a^2}}}-1\right]
\mathcal M_{0+-}^{\rm linear}\\
&+\frac{\mathcal V^{(4)}\!\left[
(\frac{2f_\pi}{f_a}\mathbb M_a^{-1}(\mathcal C)|_{m_a^2},1),
\mathcal P(m_{\pi^0}^2)-(\lambda_{\pi^0},0),
(\lambda_{\pi^+},0),(\lambda_{\pi^-},0)\right]}
{\sqrt{1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}[(\mathbb M_a^{-1}(\mathcal C))^2]_{m_a^2}}}\\
&+\mathcal V^{(3)}[\mathcal P(m_a^2),\mathcal P(m_{\pi^0}^2),I]
(M^2-sK)^{-1}_{IJ}\mathcal V^{(3)}[J,(\lambda_{\pi^+},0),(\lambda_{\pi^-},0)]\\
&+\mathcal V^{(3)}[\mathcal P(m_a^2),(\lambda_{\pi^+},0),I]
(M^2-tK)^{-1}_{IJ}\mathcal V^{(3)}[J,\mathcal P(m_{\pi^0}^2),(\lambda_{\pi^-},0)]\\
&+\mathcal V^{(3)}[\mathcal P(m_a^2),(\lambda_{\pi^-},0),I]
(M^2-uK)^{-1}_{IJ}\mathcal V^{(3)}[J,\mathcal P(m_{\pi^0}^2),(\lambda_{\pi^+},0)].
\end{aligned}}
$$

For $3\pi^0$, the comparison is explicitly

$$
\boxed{\begin{aligned}
\mathcal M_{000}-\mathcal M_{000}^{\rm linear}={}&
\left[\frac{1}{\sqrt{1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}[(\mathbb M_a^{-1}(\mathcal C))^2]_{m_a^2}}}-1\right]
\mathcal M_{000}^{\rm linear}\\
&+\frac{3\mathcal V^{(4)}[(\frac{2f_\pi}{f_a}\mathbb M_a^{-1}(\mathcal C)|_{m_a^2},1),
\mathcal P(m_{\pi^0}^2)-(\lambda_{\pi^0},0),
(\lambda_{\pi^0},0),(\lambda_{\pi^0},0)]}{
\sqrt{1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}[(\mathbb M_a^{-1}(\mathcal C))^2]_{m_a^2}}}\\
&+\frac{3\mathcal V^{(4)}[(\frac{2f_\pi}{f_a}\mathbb M_a^{-1}(\mathcal C)|_{m_a^2},1),
\mathcal P(m_{\pi^0}^2)-(\lambda_{\pi^0},0),
\mathcal P(m_{\pi^0}^2)-(\lambda_{\pi^0},0),
(\lambda_{\pi^0},0)]}{
\sqrt{1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}[(\mathbb M_a^{-1}(\mathcal C))^2]_{m_a^2}}}\\
&+\frac{\mathcal V^{(4)}[(\frac{2f_\pi}{f_a}\mathbb M_a^{-1}(\mathcal C)|_{m_a^2},1),
\mathcal P(m_{\pi^0}^2)-(\lambda_{\pi^0},0),
\mathcal P(m_{\pi^0}^2)-(\lambda_{\pi^0},0),
\mathcal P(m_{\pi^0}^2)-(\lambda_{\pi^0},0)]}{
\sqrt{1+\frac{2f_\pi^2}{f_a^2}
\operatorname{tr}[(\mathbb M_a^{-1}(\mathcal C))^2]_{m_a^2}}}\\
&+\mathcal V^{(3)}[\mathcal P(m_a^2),\mathcal P(m_{\pi^0}^2),I]
(M^2-sK)^{-1}_{IJ}\mathcal V^{(3)}[J,\mathcal P(m_{\pi^0}^2),\mathcal P(m_{\pi^0}^2)]\\
&+\mathcal V^{(3)}[\mathcal P(m_a^2),\mathcal P(m_{\pi^0}^2),I]
(M^2-tK)^{-1}_{IJ}\mathcal V^{(3)}[J,\mathcal P(m_{\pi^0}^2),\mathcal P(m_{\pi^0}^2)]\\
&+\mathcal V^{(3)}[\mathcal P(m_a^2),\mathcal P(m_{\pi^0}^2),I]
(M^2-uK)^{-1}_{IJ}\mathcal V^{(3)}[J,\mathcal P(m_{\pi^0}^2),\mathcal P(m_{\pi^0}^2)].
\end{aligned}}
$$

The three contact corrections retain the full $\mathcal V^{(4)}$, including
its kinetic part. Their sum preserves the neutral kinetic zero when the three
exact pion columns are identical. Away from a mixed-pole degeneracy, these
finite corrections start at relative order $(f_\pi/f_a)^2$.

In the strict no-mixing limit the outgoing-pion ALP components vanish, and
the pure-meson cubic vertex is zero by the displayed $\mathcal V^{(3)}$
formula. Thus the exchange terms vanish in that limit. At large but finite
$f_a$, they can be nonzero but start beyond the linear $1/f_a$ amplitude;
the example above is suppressed by the small $\pi^0$ ALP component.

The script checks the CCWZ source derivatives, the pole inverse, the explicit
vertices and contacts, the charge-protected charged poles, and the linear
limit. Run:

```text
.venv/bin/python scripts/derive_alp_finite_three_pion.py
.venv/bin/python scripts/derive_alp_finite_three_pion.py --test
```
