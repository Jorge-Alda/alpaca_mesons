---
title: Finite graded Hilbert series of the bounded covariant odd action
status: derived
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Graf:2020yxt, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and evidence

User-selected scope: the existing bounded VJP/VVP action, with one vector
nonet, at most two flavor traces, the fixed derivative/spurion ceilings,
independent external sources and LO Goldstone EOM; no resonance EOM.
Conventions are Alda:2026nrd's, with $X=\chi_+$, $Y=i\chi_-$ and
$Z_{TC}=m_0^2X_{TC}/3$. Natural units are used.

This calculation independently derives the **raw multiplicities** by finite
character projections, and constructs the quotient's finite Hilbert series
from the existing exact relation presentation. The latter is an inherited
input, not an independent Molien-Weyl proof that no further covariant action
relations exist. The 228-dimensional completeness claim still relies on
the previous covariant/three-flavor certificates.

Inputs:

- [Source-character method](../../extracts/Graf:2020yxt/hilbert-series-with-sources.md).
- [Flavor trace and C convention](odd-breaking-flavor-enumeration.md).
- [Covariant closure](odd-covariant-closure.md) and [EOM ceiling](odd-goldstone-eom-multiplier-reduction.md).
- [Leading relations](odd-leading-covariant-operator-count.md) and [higher relations/completeness](odd-higher-covariant-operator-count.md).
- [Power-counting grading and filtration](odd-hilbert-power-counting.md).

# Character derivation

The executable [SymPy derivation](../../../scripts/derive_odd_bounded_hilbert_series.py)
uses exact Laurent constant terms. For a U(3) covariant matrix, the central
U(1) acts trivially under conjugation, so the flavor representation is
$\mathbf1\oplus\mathbf8$ of SU(3). With $z_1z_2z_3=1$,

$$
\chi_{\rm nonet}=\sum_{a,b=1}^3z_a/z_b,
\qquad
\int d\mu_{SU(3)}F=\operatorname{CT}\left[
F\prod_{a<b}(1-z_b/z_a)\right]
$$

for Weyl-invariant finite Laurent polynomials. The C-twisted character for
plus-transpose matrices is $1-z-z^{-1}$ on the folded branch; pairs of
minus-transpose slots leave the overall product sign unchanged. Integrating
its powers with $1-z^2$ gives the C insertion trace.

| Distinct nonet slots | SU(3) singlets | C insertion trace | C-even singlets |
|---|---:|---:|---:|
| 3 | 6 | 4 | 5 |
| 4 | 23 | 9 | 16 |

These are **unrestricted trace** counts. The model's at-most-two-trace
Reynolds projector gives four C-even three-slot structures, ten four-slot
structures, and one C-odd three-slot structure for $Vf_-X$. Its exchange
multiplicities are $(3,1)$ and $(7,3)$ for the C-even three/four-slot spaces.
The existing finite-three-flavor Gram certificate establishes independence
within that trace ceiling; a full Haar projection alone cannot impose it.

For SO(4), write the vector character as
$Q=(a+a^{-1})(b+b^{-1})$ and the two-form character as
$T=a^2+1+a^{-2}+b^2+1+b^{-2}$. The two SU(2) Haar constant-term factors are
$(1-a^2)(1-b^2)$. Two distinct two-forms and two ordered vectors have three
pseudoscalar invariants; replacing the latter by their symmetric square
gives two. The disconnected O(4) branch vanishes for each distinct two-form,
so the pseudoscalar projection in this specific calculation divides the
SO(4) count by two.

After IBP eliminates the derivative on the odd slot, there are two momentum
placements in the leading sector and three in the one-$X$ sector. The ordered
two-derivative sector has AA, AB, BB placements. Its AA/BB commutators are
retained: symmetric free jets would give $2+3+2=7$ instead of nine and lose
curvature images of the declared covariant span.

Bose symmetry requires an exchange insertion, not dividing by two again.
The script derives exchange traces zero, minus one and plus one in the
leading, one-$X$ and ordered-dd spaces, respectively. Hence their Lorentz
symmetric/antisymmetric multiplicities are $(3,3)$, $(4,5)$ and $(5,4)$.
Multiplying equal exchange parities in flavor and Lorentz gives:

| Block | Raw count |
|---|---:|
| $Vf_+u\nabla$ | $4\times6=24$ |
| $VVu\nabla$ | $3\times3+1\times3=12$ |
| $Vf_+Y$, $VVY$, $Vf_-X$ | $4,3,1$ |
| $Vf_+uX\nabla$ | $10\times9=90$ |
| $VVuX\nabla$ | $7\times4+3\times5=43$ |
| $Vf_+Y\nabla\nabla$ | $4\times9=36$ |
| $VVY\nabla\nabla$ | $3\times5+1\times4=19$ |
| $Vf_-X\nabla\nabla$ | $1\times9=9$ |

# Quotient presentation and basis

Mark the VJP sector with $jv$ and VVP with $v^2$, and let $t$ count chiral
order. These are **sector markers**, not a field-content grading of nonlinear
Bianchi remainders. The raw space $F$ and independent relation space $R$ give
the finite exact sequence $0\to R\to F\to Q\to0$. All thirteen relation
vectors are homogeneous in these grades and have exact rank thirteen:

$$
H_F=jv(29t^4+135t^6)+v^2(15t^2+62t^4),
$$
$$
H_R=jv(3t^4+8t^6)+v^2(t^2+t^4),
$$
$$
H_Q=H_F-H_R=jv(26t^4+127t^6)+v^2(14t^2+61t^4).
$$

Thus $H_Q(1,1,1)=228=40+188$. Preserving the original rank-eleven action
gives 217 additions as established by the inherited embedding certificate.
This is a finite vector-space Hilbert polynomial, not an unrestricted
invariant-ring Hilbert series: multiplication generally leaves the seed
space or its trace ceiling.

The script reconstructs a filtration-adapted basis by exact row reduction,
eliminating expensive parent columns first. It does not use the saved basis
list to select those representatives. The new JSON contains their raw IDs;
the original normalized explicit representatives remain in
[the action manifest](../../../data/vector_sector/odd_covariant_action_basis.json).
The two basis choices are equivalent within the stated relation presentation.

# Reproduction and limits

```text
.venv/bin/python scripts/derive_odd_bounded_hilbert_series.py
.venv/bin/python scripts/verify_joint_odd_quotient.py
```

[Executed transcript](../../../data/vector_sector/odd_bounded_hilbert_series.md)
and [hashed result/basis](../../../data/vector_sector/odd_bounded_hilbert_series.json)
contain intermediate characters, block counts, ranks and the filtered basis.
The existing joint verifier has been rerun successfully.

The derivation asserts agreement between the independently character-derived
raw polynomial and the covariant seed metadata. It checks relation rank,
grade homogeneity, basis cardinality in each grade and nested filtration
images. It does not assert an independently computed IBP cohomology for a
full resonance jet algebra. No physical SM/ALP kernel or fit rank follows
from this operator count.
