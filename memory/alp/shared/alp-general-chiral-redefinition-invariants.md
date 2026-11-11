---
title: General chiral field redefinitions and invariant partonic ALP couplings
status: derived
created: 2026-11-03
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd, Bonnefoy:2022rik]
provenance: analytical-derivation
---

# General chiral field redefinitions and invariant partonic ALP couplings

## Scope and inputs

For the current working truncation, retain the QCD gauge and quark kinetic
terms, the quark mass source, ALP-dependent mass couplings, and derivative
currents with at most one derivative acting on $a$. The free ALP kinetic
term is included separately. Do not include dipoles, higher-derivative
quark bilinears, electroweak fields, or interactions with two or more ALP
derivatives. Set $x=a/f_a$ and allow general ALP-dependent unitary flavor redefinitions
$$
q_R'=g_R(x)q_R,\qquad q_L'=g_L(x)q_L,\qquad g_R(0)=g_L(0)=\mathbf1.
$$
The current and chirality-changing source conventions follow
`Alda:2026nrd` [@Alda:2026nrd]. The broader partonic operator generators are
recorded separately in the [operator-generation note](alp-shift-symmetric-quark-eft-operators.md).
The distinction between shift-symmetric derivative couplings and
spurion-dependent fermion couplings is also discussed in
`Bonnefoy:2022rik` [@Bonnefoy:2022rik]. This record stays in the partonic
theory.

## Restricted partonic Lagrangian

Write the QCD gauge and massless quark kinetic terms as
$\mathcal L_{\mathrm{QCD}}^{\mathrm{gauge+kin}}$. The retained Lagrangian is
$$
\begin{aligned}
\mathcal L_{\mathrm{partonic}}={}&\mathcal L_{\mathrm{QCD}}^{\mathrm{gauge+kin}}
+\frac12\partial_\mu a\partial^\mu a-V(a)
-\bar q_R X_{RL}(x)q_L-\bar q_LX_{RL}^\dagger(x)q_R\\
&+\frac{\partial_\mu a}{f_a}
\left(\bar q_R\gamma^\mu C_R(x)q_R+
\bar q_L\gamma^\mu C_L(x)q_L\right)
+\Phi(x)G^A_{\mu\nu}G^{A\mu\nu}
+\Theta(x)G^A_{\mu\nu}\widetilde G^{A\mu\nu}.
\end{aligned}
$$
Here $X_{RL}(x)=M(x)=S(x)+iP(x)$ is the full chirality-changing mass
function, with $M(0)=m_q$. The $C_{R,L}(x)$ are Hermitian current
coefficient functions. $\Phi(x)$ and $\Theta(x)$ are the ALP-dependent
coefficients of the scalar and topological gluon terms, respectively. Every
ALP interaction in this truncation has at most one factor of $\partial_\mu a$.
Absorb $\Phi(0)$ into the QCD gauge coupling; $\Theta(0)$ contributes to the
physical QCD vacuum angle together with the phase of $M(0)$.
Since $g_{L,R}$ depend only on $x$, transforming the QCD kinetic term
generates only one such derivative; this restricted operator set closes
under the field redefinitions considered here. The anomaly shift and CP
selection rules for the gluonic coefficient functions are in the separate
[QCD topological-sector note](alp-qcd-topological-cp.md).

## Chiral transformations of the retained functions

The only derivative-current coefficient functions are $C_R(x)$ and
$C_L(x)$; they transform as connections in $x$ because of the quark kinetic
terms:
$$
C_\eta'(x)=g_\eta C_\eta g_\eta^\dagger
+i g_\eta\partial_xg_\eta^\dagger,\qquad \eta=R,L.
$$

The retained mass function $M(x)=X_{RL}(x)$ transforms homogeneously,
$$
M'(x)=g_R(x)M(x)g_L^\dagger(x).
$$
If $M_n=d^nM/dx^n|_0$ and
$g_{\eta,r}=d^rg_\eta/dx^r|_0$, its exact Taylor law is
$$
M_n'=\sum_{r+s+t=n}\frac{n!}{r!s!t!}
g_{R,r}M_tg_{L,s}^\dagger.
$$
Here $M=S+iP$ contains the scalar/pseudoscalar mass and ALP source
coefficients.

For their Taylor derivatives,
$$
C_{\eta,n}'=
\sum_{r+s+t=n}\frac{n!}{r!s!t!}
\left.\frac{d^rg_\eta}{dx^r}\right|_0
C_{\eta,t}
\left.\frac{d^sg_\eta^\dagger}{dx^s}\right|_0
+i\sum_{r+s=n}\binom nr
\left.\frac{d^rg_\eta}{dx^r}\right|_0
\left.\frac{d^{s+1}g_\eta^\dagger}{dx^{s+1}}\right|_0.
$$
The inhomogeneous sum produces the higher powers of $x$ multiplying the
single-derivative ALP current retained here.

For the usual constant vector and axial generators,
$$
g_R=e^{ixG_R},\qquad g_L=e^{ixG_L},\qquad
G_R=\delta+\kappa,\quad G_L=\delta-\kappa,
$$
the current coefficients are $C_R=C_V+C_A$ and $C_L=C_V-C_A$. Their
all-order derivative laws are
$$
C_{V,n}'=\frac12\sum_{m=0}^n\binom nm
\left[(i\operatorname{ad}_{G_R})^m C_{R,n-m}
+(i\operatorname{ad}_{G_L})^m C_{L,n-m}\right]
+\delta_{n0}\delta,
$$
$$
C_{A,n}'=\frac12\sum_{m=0}^n\binom nm
\left[(i\operatorname{ad}_{G_R})^m C_{R,n-m}
-(i\operatorname{ad}_{G_L})^m C_{L,n-m}\right]
+\delta_{n0}\kappa,
$$
where $\operatorname{ad}_G(X)=[G,X]$. In particular,
$$
\begin{aligned}
C_{V,1}'&=C_{V,1}+i[\delta,C_{V,0}]+i[\kappa,C_{A,0}],\\
C_{A,1}'&=C_{A,1}+i[\delta,C_{A,0}]+i[\kappa,C_{V,0}].
\end{aligned}
$$
For $C_V(x)=C_{V,0}+ixc_q^{V,2}+\cdots$ and an axial leading coefficient
$C_{A,0}=c_q^A$, this gives
$$
c_q^{V,2\prime}=c_q^{V,2}+[\delta,C_{V,0}]+[\kappa,c_q^A].
$$
With a pure axial redefinition, $\delta=0$, this is the example
$c_q^{V,2\prime}=c_q^{V,2}+[\kappa,c_q^A]$. The corresponding operator is
$i(a\partial_\mu a/f_a^2)\bar q\gamma^\mu c_q^{V,2}q$; Hermiticity requires
$c_q^{V,2}$ to be anti-Hermitian in this explicit-$i$ convention.

For scalar and pseudoscalar sources, decompose the transformed mass function
as $M'(x)=S'(x)+iP'(x)$. At first order for the constant generators,
$$
\begin{aligned}
S_1'&=S_1+i[\delta,S_0]-\{\kappa,P_0\},\\
P_1'&=P_1+i[\delta,P_0]+\{\kappa,S_0\}.
\end{aligned}
$$
At higher orders, take the Hermitian and anti-Hermitian parts of the Taylor
coefficients $M_n'$ above. These transformations apply to every scalar and
pseudoscalar mass coupling retained in this Lagrangian. For a flavor-singlet
axial rotation, include its anomalous shift of $\Theta(x)$ separately; it is
a QCD topological coefficient, not another quark-bilinear coefficient
[@Alda:2026nrd].

## Invariants at every order

The connection coefficients $C_R,C_L$ by themselves are basis coordinates.
Remove their field-redefinition dependence with transport matrices
$$
\partial_xW_\eta=iC_\eta W_\eta,\qquad W_\eta(0)=\mathbf1.
$$
Hermitian $C_\eta$ gives unitary $W_\eta$. Under the inhomogeneous connection
law, $W_\eta'=g_\eta W_\eta$. The retained chirality-changing mass function
has the invariant dressing
$$
\widehat M(x)=W_R^\dagger(x)M(x)W_L(x),
$$
which is invariant for $g_\eta(0)=\mathbf1$. Its Hermitian and
anti-Hermitian parts define $\widehat S=(\widehat M+\widehat M^\dagger)/2$
and $\widehat P=(\widehat M-\widehat M^\dagger)/(2i)$. If constant flavor
rotations at $x=0$ are also allowed, these dressed matrices transform by that
residual constant rotation; their invariant traces and physical contractions
are unchanged.

This gives a recursion for invariant Taylor coefficients at arbitrary order.
Writing $W_{\eta,n}=d^nW_\eta/dx^n|_0$,
$$
W_{\eta,0}=\mathbf1,\qquad
W_{\eta,n+1}=i\sum_{r=0}^n\binom nrC_{\eta,r}W_{\eta,n-r},
$$
and
$$
\widehat M_n=\sum_{r+s+t=n}\frac{n!}{r!s!t!}
W_{R,r}^\dagger M_tW_{L,s}.
$$
Thus, within this truncation, the invariant mass/source coefficients at
order $1/f_a^n$ are the $n$th Taylor coefficients of $\widehat M(x)$. A subtraction
such as $c_q^{V,2}-[\kappa,c_q^A]$ reconstructs the coefficient before one
specified rotation, but it is not invariant under arbitrary higher-order
$g_L,g_R$. The dressed mass-source function supplies the general invariant
combinations in this truncation. At first order, applying
Alda:2026nrd's source
normalization reduces these dressed combinations to the familiar
$c_q^{VS}=c_q^V-c_q^S$ and $c_q^{AP}=c_q^A-c_q^P$
[@Alda:2026nrd]; the derivation, including the kernel of the commutator map, is in the [first-order record](alp-dressed-mass-first-order-invariants.md). At higher orders they are not obtained by taking the same
pairwise differences of coefficients at that order: each dressed Taylor
coefficient contains products and commutators of lower-order current and
spurion coefficients. For example, $\widehat M_2$ contains $C_{R,1}$ and
$C_{L,1}$; their vector part $C_{V,1}=ic_q^{V,2}$ contributes
$[c_q^{V,2},M_0]$ together with the other lower-order terms.


## Matching to the $U(3)$ CCWZ theory

The invariant partonic functions specify ALP-dependent chiral sources. They
do not specify the complete meson Lagrangian or its higher-order low-energy
constants. For the leading meson theory, use the $U(3)$ conventions of
[Alda:2026nrd](../../extracts/Alda:2026nrd/appendix-a-expansion.md):
$$
U=u^2=\exp(i\Phi/f_\pi),\qquad \langle\lambda_a\lambda_b\rangle=2\delta_{ab}.
$$
Take the chirality-changing source and its CCWZ spurion to be
$$
\widehat M(x)=\widehat S(x)+i\widehat P(x),\qquad
\chi(x)=2B_0\widehat M(x).
$$
Here $M(0)=m_q$ has the partonic mass normalization in the preceding
sections; $B_0$ converts it to the conventional dimension-two chiral source.
If the input mass function already includes $2B_0$, use it directly as
$\chi$.
$$
\chi_+=u^\dagger\chi u^\dagger+u\chi^\dagger u,\qquad
\chi_-=u^\dagger\chi u^\dagger-u\chi^\dagger u.
$$
The leading even-parity meson Lagrangian is
$$
\mathcal L_2=\frac{f_\pi^2}{4}\langle u_\mu u^\mu+\chi_+\rangle+
\mathcal L_{\rm singlet},
$$
where the singlet potential and topological-source dependence must be fixed
in the selected $U(3)$ power counting.

For the derivative source, decompose $C_R=C_V+C_A$ and $C_L=C_V-C_A$ and
set
$$
r_\mu=(\partial_\mu x)C_R(x),\qquad
\ell_\mu=(\partial_\mu x)C_L(x),\qquad x=a/f_a.
$$
The CCWZ blocks, in the convention of [Alda:2026nrd], are
$$
u_\mu=i\left[u^\dagger(\partial_\mu-ir_\mu)u-u(\partial_\mu-i\ell_\mu)u^\dagger\right],\qquad
\Gamma_\mu=\frac12\left[u^\dagger(\partial_\mu-ir_\mu)u+u(\partial_\mu-i\ell_\mu)u^\dagger\right].
$$
The covariant derivative uses $\Gamma_\mu$; $u_\mu$ contains the Goldstone
vielbein and axial source. At $\Phi=0$,
$u_\mu=(\partial_\mu x)(C_R-C_L)=2(\partial_\mu x)C_A$, fixing the
normalization. Under the CCWZ compensator, $u_\mu\to h u_\mu h^\dagger$
and $\Gamma_\mu\to h\Gamma_\mu h^\dagger+h\partial_\mu h^\dagger$.

There is also a useful exact covariant form of the field-redefinition
invariant construction. Define
$$
\widehat U=W_R^\dagger U W_L,\qquad
D_\mu U=\partial_\mu U-ir_\mu U+iU\ell_\mu,
$$
with $r_\mu=(\partial_\mu x)C_R$ and $\ell_\mu=(\partial_\mu x)C_L$.
The transport equations imply
$$
\partial_\mu\widehat U=W_R^\dagger(D_\mu U)W_L,
$$
and the mass spurion becomes $\widehat M=W_R^\dagger M W_L$. Thus the
representative with vanishing ALP-dependent quark connections has
$\widehat U$ as its Goldstone coordinate, zero $r_\mu,\ell_\mu$, and mass
spurion $2B_0\widehat M$. The pullback connections have zero field strengths
$F_R=F_L=0$, even for noncommuting $C_{R,L}(x)$, since they depend on the
single coordinate $x$ only. This is a source-level change of variables;
the anomalous Wess-Zumino and QCD topological terms must still be transformed
with their anomaly.

The coefficient of $x^n/n!$ in $\chi(x)$ is $2B_0\widehat M_n$. At first
order, specialize to the Alda:2026nrd source parametrization
$$
S(x)=m_q+ix[c_q^S,m_q]+\cdots,\quad
P(x)=x\{c_q^P,m_q\}+\cdots,\quad
C_V(x)=c_q^V+\cdots,\quad C_A(x)=c_q^A+\cdots.
$$
Then the first-order invariant combinations are
$c_q^{VS}=c_q^V-c_q^S$ and $c_q^{AP}=c_q^A-c_q^P$. The singlet axial
anomaly is carried separately by the invariant topological source. The
selected partonic record writes the raw coefficient $\Theta_{\rm raw}$ of
$G\widetilde G$; before using the normalized convention of
[Alda:2026nrd](../../extracts/Alda:2026nrd/appendix-a-expansion.md),
convert it as $\Theta_{\rm Alda}=(4\pi/\alpha_s)\Theta_{\rm raw}$.
At higher order the CCWZ spurion uses the dressed $\widehat M_n$ and is not
given by same-order pairwise differences. Physical amplitudes also depend on
the Goldstone expansion, loops, and independent higher-order CCWZ operators;
the leading $p^2$ matching does not determine those corrections.

## Inputs and reproducible calculation

The endpoint transformations and source conventions follow
`Alda:2026nrd` [@Alda:2026nrd]. The all-order coefficient transformations,
transport invariance, and Taylor recursion are derived in
[`derive_alp_invariant_coefficient_functions.py`](../../../scripts/derive_alp_invariant_coefficient_functions.py).

Run through third order with:

```sh
.venv/bin/python scripts/derive_alp_invariant_coefficient_functions.py --order 3
```

The SymPy calculation retains noncommuting flavor matrices and verifies
transport covariance by a zero product-rule residual.
