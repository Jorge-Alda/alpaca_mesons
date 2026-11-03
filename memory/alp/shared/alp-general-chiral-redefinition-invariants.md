---
title: General chiral field redefinitions and invariant partonic ALP couplings
status: derived
created: 2026-11-03
updated: 2026-11-03
bibtex_keys: [Alda:2026nrd, Bonnefoy:2022rik]
provenance: analytical-derivation
---

# General chiral field redefinitions and invariant partonic ALP couplings

## Scope and inputs

Work with partonic quark bilinears in a fixed unbroken gauge theory. Set
$x=a/f_a$ and allow general ALP-dependent unitary flavor redefinitions
$$
q_R'=g_R(x)q_R,\qquad q_L'=g_L(x)q_L,\qquad g_R(0)=g_L(0)=\mathbf1.
$$
The current and chirality-changing source conventions follow
`Alda:2026nrd` [@Alda:2026nrd]. The operator generators at fixed ALP grade
and canonical dimension are in the [partonic operator-generation record](alp-shift-symmetric-quark-eft-operators.md).
The distinction between shift-symmetric derivative couplings and
spurion-dependent fermion couplings is also discussed in
`Bonnefoy:2022rik` [@Bonnefoy:2022rik]. This record stays in the partonic
theory.

## Full operator families and their transformation

For every allowed partonic bilinear structure $\mathcal O^I_{\eta\zeta}$,
include its full ALP-dependent coefficient function
$X^I_{\eta\zeta}(x)$. The label $I$ covers all allowed Dirac, gauge-field,
and derivative structures; $\eta,\zeta=L,R$ label the quark endpoints.
Expanding these functions in $x$ supplies the operator coefficient at every
power of $1/f_a$. A local field redefinition also differentiates $g_{L,R}$
when it acts on quark derivatives, generating the corresponding ALP-gradient
operators. Thus the operator family must be closed under both endpoint
transformations and these derivative-generated terms.

At total order $1/f_a^n$, the coefficient slots include, for every allowed
ALP-free bilinear, its $x^n$ coefficient; for every one-ALP derivative
current, the term $x^{n-1}(\partial_\mu a/f_a)$ multiplying that current;
and all higher ALP-derivative jets generated when derivatives act on
$g_{L,R}$. For example, $c_q^{V,2}$ labels the $x(\partial_\mu a/f_a)$
vector-current slot at total order $1/f_a^2$. Every other allowed partonic
bilinear has its own Taylor coefficient at each order. In a manifestly
shift-symmetric parametrization, these slots are restricted to derivative
jets and assigned spurion functions. In a general redefined parametrization,
undifferentiated powers of $x$ also occur; coefficients induced from lower
orders are correlated by the transformation laws below. Independent explicit
breaking functions add independent coefficients only as permitted by their
spurion transformations.

A homogeneous coefficient transforms as
$$
X_{\eta\zeta}^{I\prime}(x)=
g_\eta(x)X^I_{\eta\zeta}(x)g_\zeta^\dagger(x).
$$
If $X_n=d^nX/dx^n|_0$, the exact Taylor law is
$$
X_{\eta\zeta,n}^{I\prime}=
\sum_{r+s+t=n}\frac{n!}{r!s!t!}
\left.\frac{d^rg_\eta}{dx^r}\right|_0
X^I_{\eta\zeta,t}
\left.\frac{d^sg_\zeta^\dagger}{dx^s}\right|_0.
$$
This applies to every left-left, right-right, left-right, and right-left
Wilson matrix, including the higher Taylor coefficients of mass, scalar,
pseudoscalar, tensor, and other spurion sources.

The right- and left-handed derivative-current coefficients are connections
because the quark kinetic terms generate an inhomogeneous term:
$$
C_\eta'(x)=g_\eta C_\eta g_\eta^\dagger
+i g_\eta\partial_xg_\eta^\dagger,\qquad \eta=R,L.
$$
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
The inhomogeneous sum is essential: it produces the lower-order shifts and
the higher ALP-gradient terms from the field redefinition.

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

For scalar and pseudoscalar sources, combine them as $M(x)=S(x)+iP(x)$.
The chirality-changing coefficient transforms as
$$
M'(x)=g_R(x)M(x)g_L^\dagger(x),
$$
so its arbitrary-order Taylor law is the homogeneous formula above. At first
order for the constant generators,
$$
\begin{aligned}
S_1'&=S_1+i[\delta,S_0]-\{\kappa,P_0\},\\
P_1'&=P_1+i[\delta,P_0]+\{\kappa,S_0\}.
\end{aligned}
$$
At higher orders, take the Hermitian and anti-Hermitian parts of $M_n'$.
These transformations apply to every higher scalar/pseudoscalar coefficient,
not just the leading mass source.

## Invariants at every order

The connection coefficients $C_R,C_L$ by themselves are basis coordinates.
Remove their field-redefinition dependence with transport matrices
$$
\partial_xW_\eta=iC_\eta W_\eta,\qquad W_\eta(0)=\mathbf1.
$$
Hermitian $C_\eta$ gives unitary $W_\eta$. Under the inhomogeneous connection
law, $W_\eta'=g_\eta W_\eta$. Therefore, for every homogeneous Wilson matrix
with endpoints $\eta,\zeta$,
$$
\widehat X^I_{\eta\zeta}(x)=
W_\eta^\dagger(x)X^I_{\eta\zeta}(x)W_\zeta(x)
$$
is invariant under all such redefinitions with $g_\eta(0)=\mathbf1$. For a
chirality-changing source this is
$$
\widehat M(x)=W_R^\dagger(x)M(x)W_L(x),
$$
and $\widehat S=(\widehat M+\widehat M^\dagger)/2$,
$\widehat P=(\widehat M-\widehat M^\dagger)/(2i)$. The same dressing applies
to every higher-dimensional Wilson matrix. If constant flavor rotations at
$x=0$ are also allowed, the dressed matrices transform by that residual
constant rotation; their invariant traces and physical contractions are
unchanged.

This gives a recursion for invariant Taylor coefficients at arbitrary order.
Writing $W_{\eta,n}=d^nW_\eta/dx^n|_0$,
$$
W_{\eta,0}=\mathbf1,\qquad
W_{\eta,n+1}=i\sum_{r=0}^n\binom nrC_{\eta,r}W_{\eta,n-r},
$$
and
$$
\widehat X^I_{\eta\zeta,n}=
\sum_{r+s+t=n}\frac{n!}{r!s!t!}
W_{\eta,r}^\dagger X^I_{\eta\zeta,t}W_{\zeta,s}.
$$
Thus the new Wilson coefficients required at order $1/f_a^n$ are the $n$th
Taylor coefficients of the complete operator functions; their invariant
combinations are the coefficients of the dressed functions. A subtraction
such as $c_q^{V,2}-[\kappa,c_q^A]$ reconstructs the coefficient before one
specified rotation, but it is not invariant under arbitrary higher-order
$g_L,g_R$. The full dressed coefficient functions supply the general
invariant combinations of $c_q^{V,A}$, $c_q^{S,P}$, and every higher
partonic Wilson coefficient. At first order, applying Alda:2026nrd's source
normalization reduces these dressed combinations to the familiar
$c_q^{VS}=c_q^V-c_q^S$ and $c_q^{AP}=c_q^A-c_q^P$
[@Alda:2026nrd]. At higher orders they are not obtained by taking the same
pairwise differences of coefficients at that order: each dressed Taylor
coefficient contains products and commutators of lower-order current and
spurion coefficients. For example, $\widehat M_2$ contains $C_{R,1}$ and
$C_{L,1}$; their vector part $C_{V,1}=ic_q^{V,2}$ contributes
$[c_q^{V,2},M_0]$ together with the other lower-order terms.

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
