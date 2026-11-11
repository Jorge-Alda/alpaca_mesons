---
title: Arbitrary-order ALP expansion of the CCWZ building blocks
status: derived
created: 2026-11-03
updated: 2026-11-03
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Arbitrary-order ALP expansion of the CCWZ building blocks

## Assumptions and notation

Use natural units, $x=a/f_a$, and the Appendix A convention of
`Alda:2026nrd`,
$$
U=u^2=\exp(i\Phi/f_\pi),\qquad
\phi=\frac{\Phi}{2f_\pi},\qquad u=e^{i\phi}.
$$
Flavor matrices need not commute. Expand the coefficient functions at $x=0$
using derivatives (so the Taylor factorial is explicit):
$$
C_\eta(x)=\sum_{m\geq0}\frac{x^m}{m!}C_{\eta,m},\qquad \eta=R,L,
$$
$$
X_{RL}(x)=\sum_{m\geq0}\frac{x^m}{m!}X_{RL,m},\qquad
X_{LR}(x)=X_{RL}^\dagger(x)=\sum_{m\geq0}\frac{x^m}{m!}X_{LR,m}.
$$
Here $X_{RL}=M$ is the full chirality-changing mass function in the
[general partonic-source record](alp-general-chiral-redefinition-invariants.md);
$X_{LR}=M^\dagger$. This note interprets the requested $X_{\eta\zeta,m}$ as
these two chirality-changing source coefficients. There are no independent
$X_{RR}$ or $X_{LL}$ mass sources in the stated truncation. The source
normalization is $\chi=2B_0X_{RL}$, as in Appendix A of `Alda:2026nrd`
[@Alda:2026nrd].

Let $k$ count pseudoscalar fields, with a differentiated field counted as
one meson. Introduce the linear maps
$$
\mathcal C_\phi[A]=[\phi,A],\qquad
\mathcal A_\phi[A]=\{\phi,A\}.
$$
The coefficient with fixed $(m,k)$ below is the homogeneous term in $x^m$
and degree $k$ in meson fields. To obtain a labeled $k$-meson vertex,
differentiate that matrix polynomial with respect to the component fields;
this generates the Bose permutations.

## Inputs and provenance

The CCWZ definitions, field normalization, source conventions, and
one-ALP benchmark are those in the [Appendix A extract](../../extracts/Alda:2026nrd/appendix-a-expansion.md)
[@Alda:2026nrd]. The generalized Taylor coefficients and field-redefinition
invariants are defined in the [partonic source record](alp-general-chiral-redefinition-invariants.md).
The calculation is implemented in
[`derive_alp_ccwz_arbitrary_order.py`](../../../scripts/derive_alp_ccwz_arbitrary_order.py).
The [finite-redefinition record](alp-ccwz-finite-redefinitions.md) relates these same building blocks between chiral coordinates and links the corresponding vertex comparisons.

## Derivation from the defining products

The derivative-source building blocks are
$$
u_\mu=i\left[u^\dagger(\partial_\mu-ir_\mu)u-u(\partial_\mu-i\ell_\mu)u^\dagger\right],\qquad
\Gamma_\mu=\frac12\left[u^\dagger(\partial_\mu-ir_\mu)u+u(\partial_\mu-i\ell_\mu)u^\dagger\right],
$$
with $r_\mu=(\partial_\mu x)C_R(x)$ and
$\ell_\mu=(\partial_\mu x)C_L(x)$. The mass blocks are
$$
\chi_\pm=u^\dagger\chi u^\dagger\pm u\chi^\dagger u,
\qquad \chi=2B_0X_{RL}(x).
$$
Insert $u=e^{i\phi}$ and expand each exponential. Each meson insertion in
$u^\dagger C_Ru$ produces $(-i)\mathcal C_\phi$, and each insertion in
$u C_Lu^\dagger$ produces $i\mathcal C_\phi$. In the mass blocks, the
same-side exponentials produce $(-i)\mathcal A_\phi$ for $X_{RL}$ and
$i\mathcal A_\phi$ for $X_{LR}$. The exponential Taylor coefficient supplies
$1/k!$, while the source expansion supplies $1/m!$.

For $m\geq1$, the term of order $x^m$ in $u_\mu$ and $\Gamma_\mu$ is
$$
(u_\mu)_{m,k}=
\frac{x^{m-1}\partial_\mu x}{(m-1)!k!}
\left[(-i)^k\mathcal C_\phi^k[C_{R,m-1}]
-i^k\mathcal C_\phi^k[C_{L,m-1}]\right],
$$
$$
(\Gamma_\mu)_{m,k}=
-\frac{i x^{m-1}\partial_\mu x}{2(m-1)!k!}
\left[(-i)^k\mathcal C_\phi^k[C_{R,m-1}]
+i^k\mathcal C_\phi^k[C_{L,m-1}]\right].
$$
Since $\partial_\mu x=\partial_\mu a/f_a$, these terms scale as
$1/f_a^m$. Therefore strict order $1/f_a^n$ uses
$C_{\eta,0},\ldots,C_{\eta,n-1}$. The requested coefficient
$C_{\eta,n}$ is also displayed in the input series; it first contributes at
$1/f_a^{n+1}$.

The mass blocks at order $x^m$ are, for $m\geq0$,
$$
(\chi_\pm)_{m,k}=
\frac{2B_0x^m}{m!k!}
\left[(-i)^k\mathcal A_\phi^k[X_{RL,m}]
\pm i^k\mathcal A_\phi^k[X_{LR,m}]\right].
$$
Thus all $X_{RL,m}$ and $X_{LR,m}$ through $m=n$ enter the mass blocks at
orders through $1/f_a^n$. If the sources are first put in the dressed
representative of the partonic record, substitute $\widehat X_{RL,m}$ and
its conjugate; those coefficients already contain products of lower-order
$C_{\eta,r}$ from the transport recursion.

There are also source-independent terms in the CCWZ connection and vielbein.
They follow from the same exponential products. For $r\geq0$,
$$
(u_\mu)_{0,2r+1}=-\frac{2(-1)^r}{(2r+1)!}
\mathcal C_\phi^{2r}[\partial_\mu\phi],\qquad
(\Gamma_\mu)_{0,2r+2}=\frac{(-1)^r}{(2r+2)!}
\mathcal C_\phi^{2r+1}[\partial_\mu\phi],
$$
and the other meson degrees in these source-independent pieces vanish.
There is no source-independent contribution to $\chi_\pm$ beyond the
$X_{RL,0}$ mass source already included above.

## Checks, truncation, and use

The formulas are derived in the script by forming the noncommuting
exponential products, differentiating the meson powers by the product rule,
and extracting coefficients in independent counting variables for $x$ and
$\Phi$. Its assertions compare the direct products with the nested
commutator and anticommutator expressions for each retained $(m,k)$, and
check the source-free $u_\mu$ and $\Gamma_\mu$ terms. The $m=1$ derivative
formula reduces to the arbitrary-meson, linear-in-$a/f_a$ Appendix A
expansion of `Alda:2026nrd` [@Alda:2026nrd].

Run an example through $1/f_a^2$ and two mesons with

```sh
.venv/bin/python scripts/derive_alp_ccwz_arbitrary_order.py --order 2 --mesons 2
```

Run the symbolic checks with

```sh
.venv/bin/python scripts/derive_alp_ccwz_arbitrary_order.py --test
```

This is an expansion of the leading CCWZ building blocks, not a finite
higher-order chiral Lagrangian. Products of building blocks require a further
convolution over their ALP orders, and chiral loops or independent higher
chiral operators are not fixed by these source expansions.
