---
title: Higher-order ALP dependence of leading quark sources under chiral transformations
status: derived
created: 2026-11-03
updated: 2026-11-03
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Higher-order ALP dependence of leading quark sources

## Assumptions

Work in natural units and define $x=a/f_a$. The flavor matrices are
noncommuting. The local chiral transformation is parameterized as
$$
g_R(x)=e^{ixG_R},\qquad g_L(x)=e^{ixG_L},\qquad
G_R=\delta+\kappa,\quad G_L=\delta-\kappa,
$$
where $\delta$ and $\kappa$ are constant Hermitian flavor matrices. This
parameterization follows the source transformation convention of
`Alda:2026nrd`; in applications, $\delta$ and $\kappa$ are set by the chosen
field redefinition and ALP source coefficients [@Alda:2026nrd]. In the
vector/axial source decomposition, the induced connections shift the vector
and axial source coefficients by $\delta$ and $\kappa$, respectively.

The scope here is the renormalizable QCD quark Lagrangian with scalar,
pseudoscalar, vector, and axial sources. Higher-dimensional quark bilinears,
such as tensor or chromomagnetic terms, require extra sources or independent
Wilson coefficients and are not implied by the starting Lagrangian. If such
terms are included, each must be transformed in its own chiral
representation.

## Inputs and provenance

The source transformation laws and the linear ALP source convention are in
the [Alda:2026nrd Appendix-A record](../../extracts/Alda:2026nrd/appendix-a-expansion.md).
The existing [finite-mixing three-pion record](../three-pseudoscalar/alp-finite-mixing-three-pion-amplitudes.md)
uses the same $\chi_\pm$ convention and exhibits derivatives beyond first
order. The calculation is implemented in
[`derive_alp_higher_quark_bilinears.py`](../../../scripts/derive_alp_higher_quark_bilinears.py).

## Bilinears in the leading quark-source Lagrangian

Before integrating out quarks, the source-coupled bilinears are
$$
\bar q_L i\!\not\!D q_L,\qquad
\bar q_R i\!\not\!D q_R,\qquad
\bar q_R\chi q_L+\bar q_L\chi^\dagger q_R.
$$
The left and right kinetic terms contain the vector and axial currents;
there is not a separate renormalizable current bilinear. The last line
contains the scalar and pseudoscalar quark bilinears. Under bosonization these
are represented by covariant-derivative building blocks such as $u_\mu$ and
spurions such as $\chi_\pm$. At any fixed chiral order, every meson operator
made from these blocks can inherit ALP dependence. Products of blocks
contribute at total order $x^n$ whenever their Taylor degrees sum to $n$.

This is the leading quark-source sector, not a complete higher-dimensional
quark EFT basis. The separate [EFT operator classification](alp-shift-symmetric-quark-eft-operators.md)
gives the larger generating set and explains what additional truncation is
needed for a finite basis.

Thus $1/f_a$ counting by itself does not select a finite set of chiral
operators: it counts powers of the external field in the transformed sources.
The derivative expansion and chiral order are separate truncations. Repeated
powers at a fixed chiral order are fixed by the source transformation when
the original source is specified; they are not new independent couplings.

## Method: transform the source bilinears

Use the standard source laws
$$
\chi'(x)=g_R(x)\chi(x)g_L^\dagger(x),
$$
$$
r_\mu'(x)=g_R r_\mu g_R^\dagger+i g_R\partial_\mu g_R^\dagger,
\qquad
l_\mu'(x)=g_L l_\mu g_L^\dagger+i g_L\partial_\mu g_L^\dagger.
$$
Here $\chi(x)$ is the input scalar/pseudoscalar source, and $r_\mu,l_\mu$
are the right and left sources coupled to the kinetic currents. Write
$$
\chi_n=\left.\frac{d^n\chi}{dx^n}\right|_{x=0},\qquad
c_{R,n}=\left.\frac{d^n c_R}{dx^n}\right|_{x=0},\qquad
c_{L,n}=\left.\frac{d^n c_L}{dx^n}\right|_{x=0},
$$
where $r_\mu=(\partial_\mu a/f_a)c_R(x)$ and
$l_\mu=(\partial_\mu a/f_a)c_L(x)$ for the ALP-dependent parts. The ordered
mass-source derivatives are
$$
\left.\frac{d^n\chi'}{dx^n}\right|_0
=\sum_{r+s+t=n}\frac{n!}{r!s!t!}
(iG_R)^r\chi_t(-iG_L)^s.
$$
The matrix order is fixed as written. Since $x=a/f_a$,
$\partial_a^n\chi'|_0=f_a^{-n}(d^n\chi'/dx^n)|_0$.

For Alda:2026nrd's linear scalar/pseudoscalar source convention, the input
mass derivatives at zero meson field are
$$
\chi_0=2B_0m_q,\qquad
\chi_1=i\left([c_q^S,\chi_0]+\{c_q^P,\chi_0\}\right),\qquad
\chi_{n\geq2}=0.
$$
These are derivatives with respect to $x$ before the finite transformation;
the resulting transformed derivatives are generally nonzero for every $n$.

For the current sources, the exact induced connections are
$$
i g_R\partial_\mu g_R^\dagger
=G_R\frac{\partial_\mu a}{f_a},\qquad
i g_L\partial_\mu g_L^\dagger
=G_L\frac{\partial_\mu a}{f_a}.
$$
They are linear in $\partial_\mu a$ for constant generators. Higher powers
in $x$ in the current coefficients can still arise from the conjugations
$g_{R,L}c_{R,L}(x)g_{R,L}^\dagger$ if the input current source depends on
$a$.

For example, at zero meson field, retaining the input derivatives
$\chi_0,\chi_1,\chi_2$, the first two transformed mass derivatives are
$$
\chi'_1=iG_R\chi_0-i\chi_0G_L+\chi_1,
$$
$$
\chi'_2=-(G_R)^2\chi_0-\chi_0(G_L)^2+2G_R\chi_0G_L
+2iG_R\chi_1-2i\chi_1G_L+\chi_2.
$$
These expressions retain noncommuting flavor order. Higher derivatives
follow from the same finite product formula. With the meson fields restored,
the transformed spurion feeds the corresponding $\chi_\pm$ blocks and all
operators containing them.

More generally, if a chiral operator contains a product of $k$ transformed
building blocks $B_j(x)$, its nth derivative is
$$
\left.\frac{d^n}{dx^n}\operatorname{tr}\!\left(
B_1(x)\cdots B_k(x)\right)\right|_0
=\sum_{n_1+\cdots+n_k=n}\frac{n!}{n_1!\cdots n_k!}
\operatorname{tr}\!\left(B_{1,n_1}\cdots B_{k,n_k}\right),
$$
where $B_{j,n_j}=d^{n_j}B_j/dx^{n_j}|_0$ and the product order is unchanged.
This gives the higher ALP terms in any retained higher-chiral-order
operator, such as products of $u_\mu$ and $\chi_\pm$, without introducing
new spurion coefficients.

## Result and implication for finite mixing

At a fixed chiral order, the quark bilinear sources that can feed higher
powers of $1/f_a$ are the scalar/pseudoscalar mass spurion and the
left/right kinetic-current sources. A constant-generator chiral
transformation exponentiates the mass spurion into arbitrarily high ordered
products. It induces only a first-power derivative connection in the
kinetic sources, although source conjugation and products of chiral building
blocks can produce higher total powers. For generic flavor matrices, these
terms must be kept to the desired power of $a/f_a$ or included through the
exact exponential unless a specified symmetry removes them.

The derivatives $\partial_a^n\chi_\pm$ required by finite ALP-meson pole
projections therefore follow from the transformed spurion rather than from
the $\mathcal O(1/f_a)$ truncation in `Alda:2026nrd`. The latter fixes the
source conventions and the linear terms, but does not itself supply a
higher-order truncation [@Alda:2026nrd].

## Verification

Run:

```sh
.venv/bin/python scripts/derive_alp_higher_quark_bilinears.py
```

The SymPy script expands noncommuting matrix products and differentiates
them at $x=0$. It prints the transformed mass derivatives through $n=4$ and
the right/left current-source derivatives through $n=3$, with all products
in their derived order. It ran successfully with the repository's `.venv`.
