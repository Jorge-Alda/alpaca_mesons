---
title: Exact three-flavor trace Gram closure of the covariant operator quotient
status: derived
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Ecker:1988te]
provenance: analytical-derivation
---

# Inputs and purpose

Use the [source-functional normal form](odd-covariant-source-gauge-reduction.md)
and [Alda:2026nrd's canonical neutral generators](alda-neutral-amplitude-inputs.md).
All flavor matrices have dimension three. The purpose is to test genuine
three-flavor trace identities rather than assume abstract trace words remain
independent at arbitrary field degree.

# Executable derivation

Run [the closure script](../../../scripts/close_odd_three_flavor_quotient.py):

```text
.venv/bin/python scripts/close_odd_three_flavor_quotient.py
```

The script constructs the nine Hermitian identity/Gell-Mann directions and
normalizes each using its actual SymPy trace norm. This includes
$T_0=I/\sqrt3$, $T_3=\operatorname{diag}(1,-1,0)/\sqrt2$ and
$T_8=\operatorname{diag}(1,1,-2)/\sqrt6$. Summing their matrix entries gives
the exact completeness tensor

$$
\sum_{a=0}^{8}(T_a)_{ij}(T_a)_{kl}=\delta_{il}\delta_{jk},
\qquad \operatorname{tr}(T_aT_b)=\delta_{ab}.
$$

## Polarized traces

Assign distinct color/momentum labels to each primitive-field occurrence.
Sum the permutations of occurrences of the same field, and impose total
momentum conservation. A product of traces is represented by a permutation
$\sigma$ whose cycles are its trace words. Scalar $Z_{TC}$ factors carry
momentum but no color leg.

Let $\mathcal T_\sigma$ denote the corresponding trace tensor. Contracting
each pair of color legs with the computed completeness tensor closes the
index lines into cycles of $\sigma^{-1}\tau$. Thus the exact Gram matrix is

$$
\langle\mathcal T_\sigma,\mathcal T_\tau\rangle
=N_f^{c(\sigma^{-1}\tau)},
$$

where $c$ includes fixed-point cycles. The script generates these permutations
from the actual operator words and evaluates their cycles with SymPy.
For each independent source/momentum monomial, contraction of its coefficient
vector with this Gram matrix gives a nonnegative squared norm at $N_f=3$.
The sum over monomials is zero if and only if the polarized source functional
vanishes. This incorporates all finite-flavor trace identities at the field
degrees encountered, without replacing the analytical matrix derivation by
generic component matrices.

## Closing the remaining rank gap

The [finite polarized-source evaluations](../../../scripts/certify_odd_higher_rank.py)
give an exact lower bound of 24 for the corrected remainder space. The
abstract trace-word calculation gives upper bound 25. These evaluations use
deterministic concrete U(3) directions and exact arithmetic; their failure to
saturate the upper bound is not evidence for a new operator relation.

The closure script takes the kernel of the accepted probe matrix within the
25-dimensional symbolic span. One direction remains. Its exact squared norm,
in the normalization of the saved candidate vector, is

$$
236544N_f^6-795648N_f^5+439296N_f^4
+1311744N_f^3-675840N_f^2-516096N_f.
$$

At three flavors it is **42,467,328**, strictly positive. The direction is
therefore independent. The remainder rank is exactly **25**, and no additional
three-flavor trace relation lowers the count. The integer norm is a
normalization-dependent proof value, not a physical observable.

Results, probe coordinates and the Gram contraction are pinned in
[the closed certificate](../../../data/vector_sector/odd_three_flavor_action_quotient.json).
The earlier [probe artifact](../../../data/vector_sector/odd_higher_rank_certificate.json)
is a lower-bound certificate, not the final closure.

# Independent checks

[The joint verifier](../../../scripts/verify_joint_odd_quotient.py) computes the
three-flavor Gram ranks for all color trace tensors with at most two traces:
two, five and seventeen for two, three and four color legs respectively.
They have full rank, protecting the low-source projection from hidden flavor
identities. It also computes the full integrated source-functional Gram
matrix of the four leading Bianchi flavor families and obtains rank three,
independently confirming the leading-sector remainder count.

The same verifier checks the original rank, complete basis cardinalities and
all recorded provenance hashes. See
[the verification artifact](../../../data/vector_sector/odd_joint_quotient_verification.json).
