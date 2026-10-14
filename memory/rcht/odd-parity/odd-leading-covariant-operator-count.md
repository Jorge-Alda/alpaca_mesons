---
title: Complete leading covariant odd-action quotient and 29 new operators
status: derived
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988te]
provenance: analytical-derivation
---

# Result and precise scope

The complete **leading seed sector** has **40 independent covariant action
operators** under the selected explicit-vector, source-preserving LO
Goldstone-EOM convention. The original seven VJP and four VVP operators embed
with rank eleven. Hence this sector contains **29 new independent operators**.

This sector includes $Vf_+u\nabla$, $VVu\nabla$, $Vf_+Y$, $VVY$, and
$Vf_-X$, where $X=\chi_+$ and $Y=i\chi_-$. All allowed single/double flavor
traces are included. The one-$X$ derivative operators and two-extra-derivative
$Y$/$f_-X$ operators are **not included in this count**. The total count for
the planned larger action remains open; 29 is not that total.

# Indexed inputs and assumptions

Use [the catalogue](odd-breaking-vertex-catalogue.md),
[the Lorentz construction](odd-breaking-lorentz-reduction.md),
[the flavor construction](odd-breaking-flavor-enumeration.md),
[covariant identities](odd-covariant-closure.md),
[Goldstone multiplier images](odd-goldstone-eom-multiplier-reduction.md), and
[Alda's leading topological action](../../extracts/Alda:2026nrd/lo-topological-action.md).
Preserve the explicit tensor vectors, at-most-double traces, independent
external sources and Alda field conventions. Use the LO Goldstone EOM only;
even resonance interactions are not inserted into that EOM. No resonance EOM,
physical-pole condition, fitted mass substitution or OPE prior is used.

# Executable derivation

Run the following from the repository root:

```text
.venv/bin/python scripts/lift_odd_leading_covariant.py
.venv/bin/python scripts/count_odd_leading_covariant.py
.venv/bin/python scripts/verify_odd_leading_action_count.py
```

The first script records the [lift](../../../data/vector_sector/odd_leading_covariant_lift.json).
The second saves [exact quotient results and basis IDs](../../../data/vector_sector/odd_leading_action_count.json)
and a [step-by-step transcript](../../../data/vector_sector/odd_leading_action_count.md).
The third saves [independent concrete-flavor checks](../../../data/vector_sector/odd_leading_action_verification.json).
All calculations use the repository SymPy environment and exact arithmetic.

## 1. Restore covariant directions lost by the vertex projection

Let $O_0,\ldots,O_4$ denote the five saved VJP Lorentz tensors, now with
an unconstrained $f_+$ tensor before applying Bianchi. Add

$$
N=\epsilon_{ijk\ell}V^{i\alpha}u_\alpha\nabla^\ell f_+^{jk}.
$$

These six tensors have exact Lorentz rank six. With four charge-even flavor
functions they supply 24 VJP tensors. The full chiral identity gives

$$
N=-\frac{i}{2}\epsilon_{ijk\ell}V^{i\alpha}u_\alpha
[u^\ell,f_-^{jk}],
$$

with the matrix ordering determined separately by each flavor trace function.
Thus $N$ is not generally zero. Its four flavor images have rank three.
The fourth, with $f_+$ alone in a flavor trace, vanishes exactly because the
trace of the commutator on the right is zero.

The leading raw Lorentz/flavor list therefore contains

| Seed sector | Raw tensors |
| --- | ---: |
| Generic covariant $Vf_+u\nabla$ | 24 |
| $VVu\nabla$, including equal-parity Bose sectors | 12 |
| $Vf_+Y$ | 4 |
| $VVY$ | 3 |
| $Vf_-X$ | 1 |
| Total | 44 |

## 2. Lift the EOM kernel instead of deleting projected vertices

In the formal derivative convention of the script, integration by parts and
Schouten give the generic divergence tensor

$$
D=2O_1+2O_2+2O_4-2N.
$$

There are no Fourier factors of $i$ in these formal derivative symbols.
The independent-source EOM gives

$$
Y=2\nabla\cdot u+2Z_{TC}I,
\qquad Z_{TC}=\frac{m_0^2X_{TC}}{3}.
$$

Each mass column is projected to its divergence and topological-contact
columns. The residual $-2N$ is retained as the explicit nonlinear Bianchi
column. The matrix thus tests the full remainder of every candidate relation
in the linear EOM kernel. This is a kernel-lift construction, not a claim that
each whole nonlinear operator equals its projected polynomial.

The raw-to-quotient matrix has rank 40 and nullity four. Its kernel consists
of the one traced-Bianchi zero and the three contact-free EOM combinations:
two VJP and one VVP. All four have explicit exact certificates, recorded with
the raw basis IDs in the JSON. The induced topological contacts cancel in
these combinations; no independent topological coefficient is introduced or
silently removed.

The surviving nonlinear Bianchi directions are necessary. Restricting the
lift to the five saved VJP representatives instead gives no contact-free
VJP relation and one VVP relation, as the first lift script demonstrates.
That restricted lift is a different covariant span. Its result must not be
applied to the complete leading seed space.

## 3. Embed the original action before counting additions

The five published VJP derivative seeds have generic rank five, but their
union with the five saved representatives has rank six. For example, the
second published derivative seed has a nonzero $N$ component. Recovery after
the homogeneous one-current Bianchi projection therefore did not establish
recovery of the original covariant action.

The complete six-dimensional VJP Lorentz span resolves this gap. The script
embeds the original $c_1,c_2,c_5,c_6,c_7$ with their generic coordinates,
$c_3,c_4$, the original symmetric $d_1,d_3,d_4$ seeds, and $d_2$. Their
quotient rank is eleven. Row reduction with those columns first produces a
40-element basis that preserves all original operators and lists 29 new
representatives. Changes in nonzero coefficient normalization do not affect
this rank statement; normalization for physical predictors remains a separate
kernel calculation.

# Verification and implications

Exact row reduction reproduces all columns from the selected pivot columns.
The script records all four kernel vectors and the original-seed embedding.
Concrete rationally scaled Gell-Mann/identity probe fields independently give
Bianchi-image rank three and verify the traced-Bianchi zero, ruling out an
extra hidden three-flavor trace relation in this remainder subspace.
Input/script hashes are checked by the verifier.

The numerical equality between 40 here and the previous leading vertex count
is coincidental: restoring three nonzero covariant directions and reducing
three EOM relations changes the representatives. This count establishes
action independence, not the number of scalar kernels or measurable fit
parameters. In particular, the nonlinear Bianchi directions can enter
multi-pseudoscalar or source amplitudes despite disappearing from the saved
linear one-odd-block projection.

The remaining global calculation must restore and reduce the corresponding
Bianchi and derivative-order remainders in the one-$X$ derivative and
two-extra-derivative mass/source blocks, then test cross-block relations.
Do not add their old vertex dimensions to 40 and call that the total action
count.
