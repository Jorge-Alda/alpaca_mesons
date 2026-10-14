---
title: Off-shell Lorentz and integration-by-parts reduction of odd vertex candidates
status: derived
created: 2026-10-13
updated: 2026-10-15
bibtex_keys: [Ruiz-Femenia:2003jdx, Roig:2013baa, Alda:2026nrd]
provenance: analytical-derivation
---

Canonical tensor and odd-coupling conventions follow [the Levi-Civita audit](levi-civita-convention.md). The linked original scripts and saved artifacts use the explicitly identified reference backend; apply that record's signed coefficient and tensor adapter when interpreting them in the canonical convention.

# Local tensor-vertex reduction

## Inputs and assumptions

Use the [flavor enumeration](odd-breaking-flavor-enumeration.md),
the [Ruiz derivative seeds](../../extracts/Ruiz-Femenia:2003jdx/vector-odd-parity.md)
and [VVP seeds](../../extracts/Ruiz-Femenia:2003jdx/vvp-three-pseudoscalar.md).
Four-dimensional metric $(+---)$, $\epsilon^{0123}=+1, \epsilon_{0123}=-1$, canonical
antisymmetric tensor vectors and all-incoming momentum conservation apply.
No on-shell mass relation, transversality of a generic tensor vector,
Goldstone EOM, OPE prior or fitted coefficient is used in the primary reduction.
The independent axial block $u_\mu$ is not replaced by a scalar momentum
until the separately labelled projection.

The declared sector has two vector objects and either one $u_\mu$ block
with one additional derivative, optionally one $\chi_+$, or one $i\chi_-$
block with zero/two additional derivatives. Also retain $Vf_-\chi_+$
with zero/two additional derivatives. No second spurion is permitted.
This fills the baseline derivative ceiling, including the higher-derivative
mass/source terms; it is not only a dressing of the original derivative seeds.

## Executable derivation

Run [the script](../../../scripts/enumerate_odd_breaking_vertices.py).
The [JSON](../../../data/vector_sector/odd_vertex_enumeration.json) contains
raw labels, independent labels and exact rational coordinates for every
discarded candidate. The [transcript](../../../data/vector_sector/odd_vertex_enumeration.md)
contains the executed Lorentz polynomials and the literature dictionary audit.

For two 2-forms $A^{ab},B^{cd}$ and two vectors $U^e,D^f$, there are six
index slots $(a,b,c,d,e,f)$. Enumerate all 15 choices of the metric-contracted
pair; the four remaining slots enter $\epsilon$ in their original order.
Contractions of the two indices of a single antisymmetric form vanish
automatically. This enumeration includes divergences and scalar dot products.

1. Compute each contraction by explicit finite Lorentz-component sums with
   SymPy `LeviCivita` and the metric. Flavor matrices stay abstract. This
   enforces antisymmetry and all four-dimensional Schouten relations exactly.
2. For derivative-$u$ vertices, initially place the derivative on $A,B$
   and, when present, $X$. Verify the omitted derivative-on-$U$ contraction
   equals their negative sum under total momentum conservation. This is
   the linear-vertex integration-by-parts reduction. In particular, the
   derivative-on-$X$ term is not discarded.
3. For VJP $f_+$ at one-current order, use
   $B^{\mu\nu}=k_B^\mu e_B^\nu-k_B^\nu e_B^\mu$, implementing its
   linear Bianchi identity. Neither internal VVP tensor is constrained this
   way. For $f_-$, retain a general antisymmetric form: a homogeneous Bianchi
   reduction would lose the $u f_+$ terms relevant to its Goldstone expansion.
4. For two-derivative scalar blocks, distribute the remaining two derivatives
   as $AA,AB,BB$ after integration by parts eliminates the scalar momentum.
   Enumerate the same metric-pair contractions. Second-covariant-derivative
   order and its nonlinear curvature terms must be retained in the chosen
   action representatives when stage B expands beyond this vertex sector.
5. Expand with `Poly`, build the exact monomial-coefficient matrix and use
   rational `rref`. The selected columns multiplied by stored coordinates
   reproduce the entire matrix, giving a symbolic certificate for every
   reduction. No random-rank test or floating-point tolerance defines a basis.
6. For VVP, exchange both vector tensors and their momenta. Decompose into
   symmetric and antisymmetric Lorentz subspaces and pair only equal-parity
   flavor tensors so the complete vertex is Bose-even.

## Results

| Derivative sector | Lorentz rank | VVP symmetric | VVP antisymmetric |
|---|---:|---:|---:|
| $Vf_+u\nabla$ | 5 | — | — |
| $Vf_+u\chi_+\nabla$ | 8 | — | — |
| $VVu\nabla$ | 6 | 3 | 3 |
| $VVu\chi_+\nabla$ | 9 | 4 | 5 |
| $Vf_+(i\chi_-)\nabla\nabla$ | 5 | — | — |
| $VV(i\chi_-)\nabla\nabla$ | 7 | 4 | 3 |
| $Vf_-\chi_+\nabla\nabla$ | 7 | — | — |

The zero-additional-derivative mass/$f_-$ sectors have the unique
$\epsilon_{\mu\nu\rho\sigma}A^{\mu\nu}B^{\rho\sigma}$ Lorentz tensor.
The script reconstructs the published $c_1,c_2,c_5,c_6,c_7$ derivative
seeds and finds rank 5, spanning the VJP space. The symmetrized $d_1,d_3,d_4$
derivative seeds have rank 3 and span its Bose-symmetric VVP subspace.
Their exact coordinates are stored separately from flavor normalization.

## Constant-background and scalar projections

Setting the spurion momentum to zero reduces the derivative-$\chi_+$
Lorentz ranks $8\to5$ for VJP and $9\to6$ for VVP. The lost directions
are source-gradient terms, not full-action redundancies.

Then specialize the tensor legs to canonical wedge tensors and replace
$u_\mu$ by a scalar-pseudoscalar momentum, keeping virtualities general:

| Sector | Projected Lorentz rank | Symmetric | Antisymmetric |
|---|---:|---:|---:|
| VJP, with/without a constant $\chi_+$ | 3 | — | — |
| VVP, with/without a constant $\chi_+$ | 3 | 2 | 1 |

These exact polynomial relations explain why off-shell coefficients can
merge in two-boson targets. They do not determine measurement-weighted rank
or source-marginalized uncertainty, and they are not used to delete axial
source vertices from the catalogue. The $f_-$ Goldstone/photon expansion,
mass-source derivative terms and dynamical ALP insertions still require
normalized covariant-block expansion in stage B.

## Scope of equivalence

The reduction is exact for the specified local external-block vertices.
No nonlinear contact remainder generated by a chiral Bianchi identity,
covariant-derivative commutator or Goldstone EOM is declared zero in the
full action. A canonical covariant representative must retain or explicitly
classify those higher-field/source terms in stage B. The output is not a
claim of an EOM-minimal complete $U(3)$ generating functional, nor a basis
for general multi-pseudoscalar contact amplitudes.
