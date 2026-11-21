---
title: Pole projection for three pseudoscalars
status: derived
created: 2026-11-12
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Pole projection for three pseudoscalars

## Inputs and assumptions

Use the conventions of `Alda:2026nrd` [@Alda:2026nrd]: natural units, mostly-minus metric, $S=1+iT$, $\Phi=\lambda_IP_I$ with $\langle\lambda_I\lambda_J\rangle=2\delta_{IJ}$, $U=\exp(i\Phi/f_\pi)$, $\chi_0=2B_0m_q$. The parent is labelled $a$, the daughters $i,j,k$, and internal fields $\ell$ or $m$; the ALP is one value of an internal label. Squared masses, Mandelstam invariants and entries of $M^2$ are in one common unit, for example GeV$^2$. The calculation is tree level in LO chiral perturbation theory.

- **Invariant couplings.** $c_q^{AP}=c_q^A-c_q^P$, $c_q^{VS}=c_q^V-c_q^S$, $c_{GA}=c_G+\langle c_q^A\rangle$, and $\mathcal C=-B_0\{m_q,c_q^{AP}\}-\frac{m_0^2}{N_f}c_{GA}\mathbf1$. The pole operator at squared mass $s$ is $\mathbb M_{s}(X)=B_0\{X,m_q\}+\frac{m_0^2}{N_f}\langle X\rangle\mathbf 1-sX$, with inverse $\mathbb M_{s}^{-1}(\mathcal C)$. It is used at $s=m_a^2$ for the parent and at $s=m_r^2$ for daughter $r$. The normalization is
$$
N(s)=1+2\,\frac{f_\pi^2}{f_a^2}\,\big\langle\mathbb M_{s}^{-1}(\mathcal C)^2\big\rangle ,
$$
the same object as in the [operator-inverse record](../shared/alp-lo-finite-mixing-operator-inverse.md), which writes the prefactor as $2(f_\pi/f_a)^2$. Alda's first-order coefficients at the parent pole are $\widehat c_q^A=-\mathbb M_{m_a^2}^{-1}(\mathcal C)$ and $\widehat c_q^P=-c_q^{AP}-\mathbb M_{m_a^2}^{-1}(\mathcal C)$, as in the [generic-generator expansion](alp-three-meson-generic-generators.md).

- **Starting coordinates.** Begin with generic unprimed source functions $C_V(x)$, $C_A(x)$, $M(x)=S(x)+iP(x)$, and $\Theta(x)$, with Taylor derivatives denoted by subscripts. Transport them to derivative-source-free coordinates using $\partial_xW_R=i(C_V+C_A)W_R$, $\partial_xW_L=i(C_V-C_A)W_L$, $W_{R,L}(0)=\mathbf1$. The invariant sources are $\widehat M=W_R^\dagger M W_L$ and $\widehat\Theta=\Theta+\int_0^x\operatorname{tr}C_A(x')\,dx'$, with the singlet normalization of the [dressed topological record](../shared/alp-dressed-topological-invariant.md); in particular $\widehat\Theta_n=\Theta_n+\operatorname{tr}C_{A,n-1}$. Thus derivative-source-free describes the transported coordinates and does not require the original unprimed $C_{V,n}$ or $C_{A,n}$ to vanish. Use the quadratic matrices and the vertices of sections 2 and 3 of the [finite-mixing three-pion record](alp-finite-mixing-three-pion-amplitudes.md), with $\chi=2B_0\widehat M$, and include the topological-potential vertices from $\widehat\Theta$ given in step 3. With $\widehat M=\widehat S+i\widehat P$, at zero meson field $\chi_+=4B_0\widehat S$ and $\chi_-=4iB_0\widehat P$:
$$
\partial_a^n\chi_+|_0=\frac{4B_0}{f_a^n}\widehat S_n,\qquad \partial_a^n\chi_-|_0=\frac{4iB_0}{f_a^n}\widehat P_n,
$$
with the Taylor coefficients $S_n$, $P_n$ (derivatives with respect to $x=a/f_a$) of the [general-source record](../shared/alp-general-chiral-redefinition-invariants.md). In the derivative-source-free coordinates the currents vanish and the mass source is the dressed one, $M=\widehat M$, so there $S_n=\widehat S_n$ and $P_n=\widehat P_n$. In the Alda:2026nrd parametrization ($S_0=m_q$, $P_0=0$, $S_1=i[c_q^S,m_q]$, $P_1=\{c_q^P,m_q\}$, $C_{V,0}=c_q^V$, $C_{A,0}=c_q^A$) the [first-order record](../shared/alp-dressed-mass-first-order-invariants.md) gives
$$
\widehat S_1=-i[c_q^{VS},m_q],\qquad \widehat P_1=-\{c_q^{AP},m_q\},
$$
and the [second-order record](../shared/alp-dressed-mass-second-order-invariants.md) gives $\widehat S_2$, $\widehat P_2$ in closed form; for Alda's linear sources ($S_2=P_2=C_{V,1}=C_{A,1}=0$)
$$
\widehat S_2=\big[c_q^V,[2c_q^S-c_q^V,m_q]\big]+\big\{c_q^A,\{2c_q^P-c_q^A,m_q\}\big\},\qquad
\widehat P_2=-i\big[c_q^V,\{2c_q^P-c_q^A,m_q\}\big]-i\big\{c_q^A,[2c_q^S-c_q^V,m_q]\big\}.
$$
Higher orders follow from the transport recursion $\widehat M_n=\sum_{r+s+t=n}\frac{n!}{r!s!t!}W_{R,r}^\dagger M_tW_{L,s}$ of the general-source record. These are the same derivatives as the ordered-product formula of section 1 of the three-pion record (checked by the script up to second order).

- **Ordered tensors.** For matrices $Y_i,Y_j,Y_k$ define
$$
T^A_{ij;k}(Y_i,Y_j,Y_k)=\big\langle\widehat c_q^A\big([Y_i,[Y_j,Y_k]]+[Y_j,[Y_i,Y_k]]\big)\big\rangle,\qquad T^A_{ij;k}=T^A_{ji;k},
$$
$$
T^P_{ijk}(Y_i,Y_j,Y_k)=\Big\langle\sum_{\sigma\in S_3}Y_{\sigma(i)}Y_{\sigma(j)}Y_{\sigma(k)}\,\widehat P\Big\rangle,\qquad \widehat P=\{\widehat c_q^P,m_q\},
$$
where the sum over $\sigma\in S_3$ runs over the six orderings of the three arguments. These are the tensors $T^A$ and $T^P$ of the [generic-generator record](alp-three-meson-generic-generators.md), with $c_q^A\to\widehat c_q^A$ and $\{c_q^P,m_q\}\to\widehat P$, and with generators $Y=\lambda$ as arguments. The arguments may be any matrices, in particular the meson columns of step 2. The tensors are multilinear, so a factor multiplying a column can be pulled out.

- Nothing is expanded in $f_\pi/f_a$. Daughters are not assumed distinct, and the neutral daughter pole columns retain their $\xi_r$ entries in the pole-aligned coordinates. These entries are coordinates of the pole vectors, not basis-independent properties of the states. The selected parent pole is simple, has nonzero overlap with the ALP interpolating field in the starting coordinates, and has an invertible meson block at its pole, as required by the [single-pole record](../shared/alp-single-pole-alignment.md).

## Step 1: pole-aligned coordinates

Apply the pole alignment of the [pole-aligned vertex record](../shared/alp-lo-pole-aligned-tree-vertices.md) to the derivative-source-free coordinates, with the parent overlap of section 2 of the three-pion record, $w_I\lambda_I=\frac{2f_\pi}{f_a}\mathbb M_{m_a^2}^{-1}(\mathcal C)$:
$$
e^{i\Phi/f_\pi}=e^{ia'\mathbb M_{m_a^2}^{-1}(\mathcal C)/f_a}\,e^{i\Phi'/f_\pi}\,e^{ia'\mathbb M_{m_a^2}^{-1}(\mathcal C)/f_a},\qquad a=a'.
$$
This is the constant-generator axial rotation $G_R=-G_L=-\mathbb M_{m_a^2}^{-1}(\mathcal C)$ of the [chiral-generator record](../shared/alp-pole-alignment-chiral-generator.md); its linear part is the shift $P_I=P_I'+w_Ia'$ of the [single-pole record](../shared/alp-single-pole-alignment.md). From that record: the meson–ALP entry of $m_a^2K'-M^{2\prime}$ vanishes, so the parent pole is purely along $a'$, with normalization $K'_{\hat a\hat a}=N(m_a^2)$; the mixing is removed only at $s=m_a^2$, and $K'_{I\hat a}=w_I$ and $M^{2\prime}_{I\hat a}$ remain nonzero.

The generic unprimed currents have already been absorbed into $\widehat M$ and $\widehat\Theta$ by transport. The [source-transformation record](../shared/alp-chiral-source-transformations.md) and the exact Taylor laws of the [general-source record](../shared/alp-general-chiral-redefinition-invariants.md), applied with $g_R=e^{ix\widehat c_q^A}$, $g_L=g_R^\dagger$ (that is, $\delta=0$ and $\kappa=\widehat c_q^A$ in the general record) to these derivative-source-free sources, give the sources in the pole-aligned coordinates (primed, as in the general record):
$$
C_{V,n}'=0,\qquad C_{A,0}'=\widehat c_q^A=-\mathbb M_{m_a^2}^{-1}(\mathcal C),\qquad C_{A,n}'=0\ (n\ge1),\qquad
\Theta_1'=\widehat\Theta_1-\langle\widehat c_q^A\rangle=c_{GA}+\langle\mathbb M_{m_a^2}^{-1}(\mathcal C)\rangle,
$$
$$
M_n'=\sum_{r+s+t=n}\frac{n!}{r!s!t!}\,(i\widehat c_q^A)^r\,\widehat M_t\,(i\widehat c_q^A)^s ,
$$
with $S_n'$ and $P_n'$ the Hermitian and anti-Hermitian parts of $M_n'$. The vertices below use $\partial_a^n\chi_+|_0=4B_0S_n'/f_a^n$ and $\partial_a^n\chi_-|_0=4iB_0P_n'/f_a^n$. Equivalently, the closed forms of the [second-order record](../shared/alp-dressed-mass-second-order-invariants.md) with $C_{V,n}'=0$, $C_{A,0}'=\widehat c_q^A$, $C_{A,1}'=0$ express the invariants through these coefficients, $\widehat S_1=S_1'$, $\widehat P_1=P_1'-\{\widehat c_q^A,m_q\}$, $\widehat S_2=S_2'+2\{\widehat c_q^A,P_1'\}-\{\widehat c_q^A,\{\widehat c_q^A,m_q\}\}$ and $\widehat P_2=P_2'-2\{\widehat c_q^A,S_1'\}$. Solved for the primed coefficients, the first two orders are
$$
S_1'=\widehat S_1=-i[c_q^{VS},m_q],\qquad
P_1'=\widehat P_1+\{\widehat c_q^A,m_q\}=-\{c_q^{AP}+\mathbb M_{m_a^2}^{-1}(\mathcal C),m_q\}=\{\widehat c_q^P,m_q\},
$$
$$
S_2'=\widehat S_2-2\{\widehat c_q^A,\widehat P_1\}-\{\widehat c_q^A,\{\widehat c_q^A,m_q\}\},\qquad
P_2'=\widehat P_2+2\{\widehat c_q^A,\widehat S_1\}.
$$
With the first-order values inserted, $S_2'=\widehat S_2-2\{\mathbb M_{m_a^2}^{-1}(\mathcal C),\{c_q^{AP},m_q\}\}-\{\mathbb M_{m_a^2}^{-1}(\mathcal C),\{\mathbb M_{m_a^2}^{-1}(\mathcal C),m_q\}\}$ and $P_2'=\widehat P_2+2i\{\mathbb M_{m_a^2}^{-1}(\mathcal C),[c_q^{VS},m_q]\}$. Orders three and four follow from the same Taylor law with $\widehat M_3$, $\widehat M_4$ of the transport recursion. Thus in the pole-aligned coordinates the one-ALP couplings are exactly Alda's no-mixing coefficients $\widehat c_q^A$ and $\widehat c_q^P$ of the [generic-generator expansion](alp-three-meson-generic-generators.md), with $\mathbb M_{m_a^2}^{-1}(\mathcal C)$ at the exact parent mass, and $c_q^{VS}$ is unchanged.

## Step 2: external columns

Each leg $x$ carries a meson column $\lambda_x$, a matrix in flavour space, and a scalar coordinate $\xi_x$ along the ALP field. The table gives these pole-vector coordinates in the pole-aligned basis of step 1. In particular, whether a column has $\xi_x=0$ is basis dependent; only the complete pole-projected amplitude is invariant under nonsingular field redefinitions.

| leg | meson column $\lambda_x$ | ALP-field coordinate $\xi_x$ |
|---|---|---|
| parent $a$ | $0$ | $\dfrac{1}{\sqrt{N(m_a^2)}}$ |
| neutral daughter pole with nonzero $\xi_r$ in this basis | $\dfrac{2\,(f_\pi/f_a)\,(m_a^2-m_r^2)}{\sqrt{N(m_r^2)}}\,\mathbb M_{m_r^2}^{-1}(\widehat c_q^A)$ | $\dfrac{1}{\sqrt{N(m_r^2)}}$ |
| daughter pole with $\xi_r=0$ in this basis (including charged poles) | its LO meson matrix $\lambda_r$ | $0$ |
| internal meson field $\ell$ | its generator $\lambda_\ell$ | $0$ |
| internal ALP field | $0$ | $1$ |

The daughter entry is the aligned meson component $2(f_\pi/f_a)\big[\mathbb M_{m_r^2}^{-1}(\mathcal C)-\mathbb M_{m_a^2}^{-1}(\mathcal C)\big]/\sqrt{N(m_r^2)}$ of the [operator-inverse record](../shared/alp-lo-finite-mixing-operator-inverse.md). The resolvent identity gives $\mathbb M_{m_r^2}^{-1}(\mathcal C)-\mathbb M_{m_a^2}^{-1}(\mathcal C)=(m_a^2-m_r^2)\,\mathbb M_{m_r^2}^{-1}(\widehat c_q^A)$, because $\mathbb M_s$ depends linearly on $s$ and $\widehat c_q^A=-\mathbb M_{m_a^2}^{-1}(\mathcal C)$. The columns depend on the ALP couplings only through $\mathcal C$, that is through $c_q^{AP}$ and $c_{GA}$. The coupling $c_q^{VS}$ enters only through the vertices.

The invariant inputs are $\mathcal C$, $c_q^{VS}$ and the dressed coefficients $\widehat S_n$, $\widehat P_n$ with $n\ge2$. The LO meson-spectrum inputs are $m_q$, $B_0$, $m_0^2$, $f_\pi$, the bare ALP mass, the flavour matrices $\lambda_r$ of daughters without ALP component, and the pole masses $m_a^2$, $m_r^2$. The pole masses are roots of the pole equation of the operator-inverse record, with the second-derivative correction to $M^2_{\hat a\hat a}$ of the [joint-vertex record](../shared/alp-lo-joint-vertex-expansion.md) for nonlinear sources.

## Step 3: vertices in the pole-aligned coordinates

The pole-aligned Lagrangian is the LO CCWZ Lagrangian with the sources of step 1. Each leg enters through its column $(\lambda_x,\xi_x)$ of step 2. The vertices below are the LO vertices evaluated on these columns. The one-ALP kinetic vertex is the one of the [generic-generator expansion](alp-three-meson-generic-generators.md) with $c_q^A\to\widehat c_q^A$. The two-ALP kinetic term is not written in the existing notes; it is checked against the joint CCWZ expansion. The mass-source terms are those of section 3 of the [three-pion record](alp-finite-mixing-three-pion-amplitudes.md) with $\partial_a^n\chi_+|_0=4B_0S_n'/f_a^n$ and $\partial_a^n\chi_-|_0=4iB_0P_n'/f_a^n$. The generic unprimed currents enter through $\widehat M_n$ and $\widehat\Theta_n$, not as residual derivative currents: the derivative-source-free transport precedes the pole alignment. In particular $\widehat\Theta_n=\Theta_n+\operatorname{tr}C_{A,n-1}$ and the pole alignment changes only its linear coefficient, $\Theta_1' = \widehat\Theta_1-\langle\widehat c_q^A\rangle$, while $\Theta_n'=\widehat\Theta_n$ for $n\ge2$.

### Quartic vertex

For the parent $a$ and the three daughters $i,j,k$:
$$
\boxed{\begin{aligned}
\mathcal V^{(4)}_{aijk}={}&\xi_a\Big[\frac{1}{6f_\pi f_a}\Big(T^A_{jk;i}(\lambda_j,\lambda_k,\lambda_i)\,p\cdot p_i+T^A_{ik;j}(\lambda_i,\lambda_k,\lambda_j)\,p\cdot p_j+T^A_{ij;k}(\lambda_i,\lambda_j,\lambda_k)\,p\cdot p_k\Big)-\frac{B_0}{6f_\pi f_a}\,T^P_{ijk}(\lambda_i,\lambda_j,\lambda_k)\Big]\\
&+\xi_a\Big[\frac{1}{f_a^2}\sum_{r}\xi_r\,p\cdot p_r\,\big\langle[\widehat c_q^A,\lambda_s][\widehat c_q^A,\lambda_t]\big\rangle
-\frac{B_0}{2f_a^2}\sum_{r}\xi_r\big\langle S_2'\{\lambda_s,\lambda_t\}\big\rangle
+\frac{B_0f_\pi}{f_a^3}\sum_{r}\xi_s\xi_t\big\langle P_3'\,\lambda_r\big\rangle\Big]\\
&+\frac{B_0f_\pi^2}{f_a^4}\,\xi_a\xi_i\xi_j\xi_k\big\langle S_4'\big\rangle
 -\frac{f_\pi m_0^2\Theta_3'}{N_ff_a^3}\xi_a\sum_{r\in\{i,j,k\}}\operatorname{tr}(\lambda_r)\prod_{s\in\{i,j,k\},\,s\ne r}\xi_s
 -\frac{f_\pi^2m_0^2}{N_ff_a^4}\bigl(8\Theta_1'\Theta_3'+6(\Theta_2')^2\bigr)\xi_a\xi_i\xi_j\xi_k .
\end{aligned}}
$$
In the sums, $r$ runs over $\{i,j,k\}$ and $\{s,t\}$ are the two daughters other than $r$. A daughter pole with $\xi=0$ in this basis makes the terms requiring that coordinate vanish, so no separate case is needed. The $P_3'$ term is nonzero only when both other daughter columns have nonzero $\xi$ in this basis.

### Cubic vertex

For any three legs $x,y,z$, each taken from $\{a,i,j,k\}$ or the internal fields:
$$
\boxed{\begin{aligned}
\mathcal V^{(3)}_{xyz}={}&\frac{B_0}{f_a}\Big[\frac{i}{2}\Big\langle[c_q^{VS},m_q]\big(\xi_x\{\lambda_y,\lambda_z\}+\xi_y\{\lambda_z,\lambda_x\}+\xi_z\{\lambda_x,\lambda_y\}\big)\Big\rangle\\
&\qquad+\frac{f_\pi}{f_a}\big\langle P_2'\big(\xi_x\xi_y\lambda_z+\xi_x\xi_z\lambda_y+\xi_y\xi_z\lambda_x\big)\big\rangle+\frac{f_\pi^2}{f_a^2}\,\xi_x\xi_y\xi_z\,\langle S_3'\rangle\Big]
 -\frac{f_\pi m_0^2\Theta_2'}{N_ff_a^2}\sum_{r\in\{x,y,z\}}\operatorname{tr}(\lambda_r)\prod_{s\in\{x,y,z\},\,s\ne r}\xi_s
 -\frac{6f_\pi^2m_0^2\Theta_1'\Theta_2'}{N_ff_a^3}\xi_x\xi_y\xi_z.
\end{aligned}}
$$
The internal fields follow the table of step 2: a meson $\ell$ has $\lambda_\ell$ and $\xi_\ell=0$; the ALP has $\lambda=0$ and $\xi=1$. The vertex is symmetric in its three legs.

## Step 4: amplitude

The exact LO amplitude is the pole-projected sum of the quartic contact and the three exchanges:
$$
\boxed{\begin{aligned}
\mathcal M(a_{\rm phys}\to P_iP_jP_k)={}&\mathcal V^{(4)}_{aijk}\\
&+\sum_{\ell,m}\mathcal V^{(3)}_{ai\ell}\,G_{\ell m}(s_{jk})\,\mathcal V^{(3)}_{mjk}
+\sum_{\ell,m}\mathcal V^{(3)}_{aj\ell}\,G_{\ell m}(s_{ik})\,\mathcal V^{(3)}_{mik}
+\sum_{\ell,m}\mathcal V^{(3)}_{ak\ell}\,G_{\ell m}(s_{ij})\,\mathcal V^{(3)}_{mij},
\end{aligned}}
$$
with $G(s)=(M^{2\prime}-sK')^{-1}$ the full inverse of the pole-aligned quadratic form of step 1, not a sum over a few poles. The internal labels $\ell,m$ run over the nine meson fields and the ALP. The pairings are $(a,i)|(j,k)$ at $s_{jk}$, $(a,j)|(i,k)$ at $s_{ik}$ and $(a,k)|(i,j)$ at $s_{ij}$. The plus sign is the minus sign of the [pole-projected construction](../shared/alp-lo-pole-projected-tree-amplitudes.md) with $(M^2-sK)^{-1}=-D^{-1}$.

### Identical daughters

The amplitude uses the fully symmetrized vertex and has no factor $1/n!$. Identical daughters give repeated arguments in the tensors. With $\lambda_\alpha$ and $\lambda_\beta$ the generators of two species,
$$
T^A_{\alpha\alpha;\alpha}=0,\qquad T^A_{\alpha\alpha;\beta}=2\big\langle\widehat c_q^A[\lambda_\alpha,[\lambda_\alpha,\lambda_\beta]]\big\rangle,
$$
$$
T^P_{\alpha\alpha\alpha}=6\big\langle\lambda_\alpha^3\widehat P\big\rangle,\qquad
T^P_{\alpha\alpha\beta}=2\big\langle(\lambda_\alpha^2\lambda_\beta+\lambda_\alpha\lambda_\beta\lambda_\alpha+\lambda_\beta\lambda_\alpha^2)\widehat P\big\rangle .
$$
The factors $6=3!$ and $2=2!$ are the multiplicities of the small-mixing amplitudes. The identical-particle factor $S$ enters only the width, through the Dalitz integral of the [amplitude record](alp-three-meson-amplitudes-widths.md).

## Cross-check: contraction in the derivative-source-free coordinates

Route 1 of the [finite-mixing construction](../shared/alp-lo-finite-mixing-amplitude-construction.md) contracts the three-pion record's $\mathcal V^{(4)}$ and $\mathcal V^{(3)}$ with the derivative-source-free columns. In these coordinates the parent has the meson column $2(f_\pi/f_a)\mathbb M_{m_a^2}^{-1}(\mathcal C)/\sqrt{N(m_a^2)}$, so the pure-meson kinetic quartic term carries the ALP dependence. The same three exchanges appear at $(M^2-s_{xy}K)^{-1}$. Reducing the kinetic term with the decay kinematics, for three daughters with ALP component,

$$
\begin{aligned}
\mathcal V^{(4)}_{\rm kin}={}&-\frac{2(f_\pi/f_a)^4}{3f_\pi^2\sqrt{N(m_a^2)N(m_i^2)N(m_j^2)N(m_k^2)}}\Big\{
\bigl\langle[\mathbb M_{m_a^2}^{-1}(\mathcal C),\mathbb M_{m_i^2}^{-1}(\mathcal C)][\mathbb M_{m_j^2}^{-1}(\mathcal C),\mathbb M_{m_k^2}^{-1}(\mathcal C)]\bigr\rangle(s_{ik}-s_{ij})\\
&+\bigl\langle[\mathbb M_{m_a^2}^{-1}(\mathcal C),\mathbb M_{m_j^2}^{-1}(\mathcal C)][\mathbb M_{m_i^2}^{-1}(\mathcal C),\mathbb M_{m_k^2}^{-1}(\mathcal C)]\bigr\rangle(s_{jk}-s_{ij})\\
&+\bigl\langle[\mathbb M_{m_a^2}^{-1}(\mathcal C),\mathbb M_{m_k^2}^{-1}(\mathcal C)][\mathbb M_{m_i^2}^{-1}(\mathcal C),\mathbb M_{m_j^2}^{-1}(\mathcal C)]\bigr\rangle(s_{jk}-s_{ik})\Big\},
\end{aligned}
$$

and for daughter $i$ with ALP component and daughters $j$, $k$ without,

$$
\begin{aligned}
\mathcal V^{(4)}_{\rm kin}={}&-\frac{(f_\pi/f_a)^2}{6f_\pi^2\sqrt{N(m_a^2)N(m_i^2)}}\Big\{
\bigl\langle[\mathbb M_{m_a^2}^{-1}(\mathcal C),\mathbb M_{m_i^2}^{-1}(\mathcal C)][\lambda_j,\lambda_k]\bigr\rangle(s_{ik}-s_{ij})
+\bigl\langle[\mathbb M_{m_a^2}^{-1}(\mathcal C),\lambda_j][\mathbb M_{m_i^2}^{-1}(\mathcal C),\lambda_k]\bigr\rangle(s_{jk}-s_{ij})\\
&+\bigl\langle[\mathbb M_{m_a^2}^{-1}(\mathcal C),\lambda_k][\mathbb M_{m_i^2}^{-1}(\mathcal C),\lambda_j]\bigr\rangle(s_{jk}-s_{ik})\Big\}.
\end{aligned}
$$

The mass-source part of route 1 is the $\mathcal V^{(4)}_{\chi_+}$ of the three-pion record on the same columns. The kinetic and mass-source contacts of the two routes differ term by term. Their sums differ by the nonlinear contact of the [pole-aligned vertex record](../shared/alp-lo-pole-aligned-tree-vertices.md), which vanishes when all external legs are on their poles.

**Three-pion amplitudes.** Substitute the choices of section 4 of the [three-pion record](alp-finite-mixing-three-pion-amplitudes.md) in route 1: the parent column $\mathcal P(m_a^2)$; for $3\pi^0$ all daughters are $\mathcal P(m_{\pi^0}^2)$; for $\pi^0\pi^+\pi^-$ the daughter columns are $\mathcal P(m_{\pi^0}^2)$ for $i$, and the charged generators $\lambda_{\pi^+}=\sqrt2|u\rangle\langle d|$ and $\lambda_{\pi^-}=\lambda_{\pi^+}^\dagger$ for $j$ and $k$, both with no ALP component. The Mandelstam variables are $s=s_{jk}$, $t=s_{ik}$, $u=s_{ij}$. The pairings then become those of the record. The generic result equals the boxed $\mathcal M_{000}$ and $\mathcal M_{0+-}$ of the three-pion record as rational functions of $s,t,u$.

## Weak-mixing limit

### Power counting

Hold the Taylor derivatives of the source functions with respect to $x=a/f_a$ fixed as $f_a$ varies. A vertex with $n$ ALP legs is $O\big((f_\pi/f_a)^n\big)$. The pure-meson cubic vertex is zero. The meson–ALP entries of $M^2$, and hence of $(M^2-sK)^{-1}$, are $O(f_\pi/f_a)$. For nondegenerate poles, the meson components of the parent and the ALP components of meson-like daughters are $O(f_\pi/f_a)$, and the other components are $1+O\big((f_\pi/f_a)^2\big)$. Expanding the amplitude under these rules gives no $O(f_a^0)$ and no $O\big((f_\pi/f_a)^2\big)$ term. Every exchange starts at $O\big((f_\pi/f_a)^3\big)$. The counting is not uniform near a degeneracy between poles, or when an exchanged $s_{xy}$ approaches a pole of $M^2-sK$.

### Leading term

With unmixed daughters $\lambda_i$, $\lambda_j$, $\lambda_k$ (no ALP component) and $N(m_a^2)^{-1/2}$ set to one,
$$
\mathcal M=\frac{1}{6f_af_\pi}\Big(T^A_{jk;i}\,p\cdot p_i+T^A_{ik;j}\,p\cdot p_j+T^A_{ij;k}\,p\cdot p_k\Big)-\frac{B_0}{6f_af_\pi}\,T^P_{ijk}+O\big(f_a^{-3}\big),
$$
with $\widehat c_q^A=-\mathbb M_{m_a^2}^{-1}(\mathcal C)$ and $\widehat P=\{\widehat c_q^P,m_q\}$. This is the small-mixing $\mathcal M_{ijk}$ of the [amplitude record](alp-three-meson-amplitudes-widths.md) for any generators, without any on-shell condition.

For identical daughters, the reductions of the previous section give
$$
\mathcal M_{\alpha\alpha\alpha}=-\frac{B_0}{6f_af_\pi}\,T^P_{\alpha\alpha\alpha},\qquad
\mathcal M_{\alpha\alpha\beta}(s_{jk})=\frac{T^A_{\alpha\alpha;\beta}}{24f_af_\pi}\,(m_a^2+3m_\beta^2-3s_{jk})-\frac{B_0}{6f_af_\pi}\,T^P_{\alpha\alpha\beta},
$$
where in the second case the two identical daughters are $j$ and $k$ and the third daughter $i$ has species $\beta$. The Mandelstam sum rule $s_{ij}+s_{ik}+s_{jk}=m_a^2+m_i^2+m_j^2+m_k^2$ is used. With $\lambda_{\pi^0}$ and $\lambda_{\pi^\pm}$ these are the $\mathcal M^{\rm linear}_{000}$ and $\mathcal M^{\rm linear}_{0+-}$ of section 5 of the [three-pion record](alp-finite-mixing-three-pion-amplitudes.md). The remainder, of relative order $(f_\pi/f_a)^2$, contains the parent normalization $N(m_a^2)^{-1/2}-1$, the daughter corrections of step 2, the contact terms with daughter ALP components and all exchanges.

In route 1 the same order is $\mathcal V^{(4)}$ on the derivative-source-free columns. It contains the squared-mass differences of its kinetic term and orderings of $\mathbb M_{m_a^2}^{-1}(\mathcal C)$ between daughter generators. It equals the forms above only when each daughter obeys its LO on-shell condition $B_0\{m_q,\lambda_r\}+\frac{m_0^2}{N_f}\langle\lambda_r\rangle\mathbf1=m_r^2\lambda_r$, which is the vanishing of the nonlinear contact at this order. For $\alpha\alpha\alpha$ this follows in two lines: the condition implies $[\lambda_\alpha^2,m_q]=0$, and the extra orderings combine into $\langle\mathbb M_{m_a^2}^{-1}(\mathcal C)(\lambda_\alpha^2m_q\lambda_\alpha+\lambda_\alpha m_q\lambda_\alpha^2-\lambda_\alpha^3m_q-m_q\lambda_\alpha^3)\rangle=0$.

## Reproducibility and limits

[`scripts/derive_alp_finite_mixing_three_pseudoscalars.py`](../../../scripts/derive_alp_finite_mixing_three_pseudoscalars.py) imports $\mathcal V^{(3)}$, $\mathcal V^{(4)}$, the source derivatives and the matrix inverse from [`scripts/derive_alp_finite_three_pion.py`](../../../scripts/derive_alp_finite_three_pion.py), the joint CCWZ expansion from [`scripts/derive_alp_lo_joint_vertices.py`](../../../scripts/derive_alp_lo_joint_vertices.py) and the pole solver from [`scripts/solve_alp_meson_spectrum.py`](../../../scripts/solve_alp_meson_spectrum.py). It checks the quantities used here: the congruence and the vanishing meson–ALP entry at the parent pole; the aligned source derivatives; the aligned contact and cubic vertices against the joint CCWZ expansion; the exact three-pion substitutions; the route comparison at LO poles of one source choice; the multiplicities; and the weak limits.

The nonlinear singlet-potential additions in step 3 are derived by [`scripts/derive_alp_three_pseudoscalar_topological_vertices.py`](../../../scripts/derive_alp_three_pseudoscalar_topological_vertices.py) and compared with the joint CCWZ vertices by [`scripts/test_derive_alp_three_pseudoscalar_topological_vertices.py`](../../../scripts/test_derive_alp_three_pseudoscalar_topological_vertices.py). For the singlet potential $-f_\pi^2m_0^2[\langle\Phi\rangle/f_\pi+2\Theta'(a/f_a)]^2/(4N_f)$, the script also differentiates the general source through fourth order. If $\Theta'_0\ne0$, add $2\Theta'_0\Theta'_4$ to the last parenthesized factor in the quartic expression. The test verifies its coefficient. Generic unprimed $C_{V,n}$ and $C_{A,n}$ enter the amplitude through the transport recursion for $\widehat M_n$ and $\widehat\Theta_n$; the zero primed higher derivative currents are specific to the intermediate derivative-source-free representation followed by pole alignment.

Limits:

- Tree level, LO. No loops, widths, resonance exchange or NLO terms. External poles are stable simple poles of the LO spectrum.
- If an exchanged $s_{xy}$ in the physical region reaches a pole of $M^2-sK$ (for example a kaon or $\eta$ pole for a heavy enough parent), the tree propagator diverges and a width prescription beyond this calculation is needed.
- The alignment needs a parent with nonzero ALP component and an invertible meson block at $m_a^2$. Near a degeneracy, or with $O(1)$ admixtures, the pole vectors depend on the basis and the weak-mixing counting fails, as in the [pole-projected construction](../shared/alp-lo-pole-projected-tree-amplitudes.md).
- For Alda's linear sources the dressed coefficients $\widehat M_n$ with $n\ge2$ depend on $c_q^V$ and $c_q^A$ separately, not only on $c_q^{VS}$ and $c_q^{AP}$, as stated in the [general-source record](../shared/alp-general-chiral-redefinition-invariants.md). They enter only at relative order $(f_\pi/f_a)^2$ and beyond.
- The route comparison at LO poles is numerical for one source choice. The three-pion reproduction is exact for rational inputs, not for general symbolic sources.

```sh
.venv/bin/python scripts/derive_alp_finite_mixing_three_pseudoscalars.py
.venv/bin/python scripts/derive_alp_finite_mixing_three_pseudoscalars.py --test
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_finite_mixing_three_pseudoscalars.py
.venv/bin/python scripts/derive_alp_three_pseudoscalar_topological_vertices.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_three_pseudoscalar_topological_vertices.py
```
