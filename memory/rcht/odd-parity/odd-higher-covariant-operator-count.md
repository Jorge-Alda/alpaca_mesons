---
title: Complete bounded covariant odd-action count and 217 new operators
status: derived
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988te]
provenance: analytical-derivation
---

# Result

The declared covariant VJP/VVP seed space has **228 independent action
operators**, with the original seven VJP/four VVP operators retaining rank
eleven. It therefore has **217 new independent operators**:

$$
228=40+188,\qquad 217=(40-11)+188=29+188.
$$

The 40-operator [leading quotient](odd-leading-covariant-operator-count.md)
is unchanged. The one-spurion derivative and two-extra-derivative mass/source
blocks jointly contribute 188 additional operators.

This is the count for the **specified covariant seed classes**, not for an
arbitrary complete odd U(3) action, a set of fit parameters, or a physical
source/kinematic projection. Independent multi-Goldstone seeds and independent
topological coefficient functions remain outside the declared model.

# Inputs and retained choices

Use [the bounded source-resolved catalogue](odd-breaking-vertex-catalogue.md),
[the leading quotient and original-action embedding](odd-leading-covariant-operator-count.md),
[the Goldstone EOM ceiling audit](odd-goldstone-eom-multiplier-reduction.md),
[the covariant source-frame reduction](odd-covariant-source-gauge-reduction.md), and
[the exact three-flavor Gram closure](odd-three-flavor-trace-gram.md).
Preserve Alda:2026nrd's canonical fields, one vector nonet, the fixed derivative
ceiling, at most one quark spurion and two flavor traces. Only chiral
identities, exact IBP/Schouten and source-preserving LO Goldstone EOM are
used. No resonance EOM, fitted physical masses, on-shell vector constraint,
OPE prior, loop correction or zero-source simplification enters this count.

# Step-by-step reduction

1. Polarize ordered derivative slots before doing the Lorentz reduction.
   The higher raw covariant list has 197 tensors, rather than the 165 higher
   tensors in the earlier vertex projection. In particular, commuting second
   covariant derivatives prematurely would lose their curvature terms.
2. The initial linear source projection has rank 163 and 34 candidate kernel
   vectors. Evaluate their actual nonlinear source functionals with all
   scalar/axial/vector/topological derivatives retained.
3. Make that projection common to all blocks. The four $f_-X\nabla\nabla$
   candidates overlap the lowest source grade of $Vf_+uX\nabla$. Subtract
   the matching full covariant derivative operators. Two candidates vanish
   exactly; two retain higher-field remainders. The exact subtraction maps
   are stored in the joint reduction JSON.
4. The corrected higher-field map has upper rank 25. Exact finite source
   evaluations establish lower rank 24; the full color Gram norm of the
   remaining direction is positive at three flavors. This closes the rank
   at 25 without relying on sampled rank saturation.
5. Hence the higher quotient dimension is $163+25=188$. Its nine exact
   relations are saved with raw operator IDs. Selecting representatives by
   eliminating later blocks first produces the basis below.

| Retained higher block | Raw covariant tensors | Chosen basis representatives |
| --- | ---: | ---: |
| $Vf_+uX\nabla$ | 90 | 89 |
| $VVuX\nabla$ | 43 | 43 |
| $Vf_+Y\nabla\nabla$ | 36 | 31 |
| $VVY\nabla\nabla$ | 19 | 18 |
| $Vf_-X\nabla\nabla$ | 9 | 7 |
| Total | 197 | 188 |

Here $X=\chi_+$ and $Y=i\chi_-$. The per-block representative counts refer
to this declared elimination order; cross-block identities mean they are not
five unrelated quotient dimensions. The joint total is independent of basis
ordering.

# Basis and exact certificates

The [basis manifest](../../../data/vector_sector/odd_covariant_action_basis.json)
preserves the 40-element leading basis, including all original operators,
and lists all 188 higher representatives. Each higher entry records flavor
trace words, trace topology, spurion count, derivative placement/order,
Bose sector, mass dimension and its explicit covariant contraction.

VJP higher operators have dimension seven and coefficient normalization
$1/M_V^3$; VVP higher operators have dimension six and normalization
$1/M_V^2$. The generated representative convention carries the required
Hermitian $i$ for the $f_-$ family and its stated Bose normalization. These
are mathematical basis conventions; normalized SM/ALP predictor matching is
still a separate calculation. Existing physical $c_i,d_i$ implementations
have not been redefined by this count.

See the [readable higher basis](../../../data/vector_sector/odd_covariant_higher_basis.md),
[executed derivation transcript](../../../data/vector_sector/odd_covariant_action_derivation.md),
[joint relation matrix results](../../../data/vector_sector/odd_covariant_higher_reduction.json),
[closed three-flavor certificate](../../../data/vector_sector/odd_three_flavor_action_quotient.json), and
[final verification](../../../data/vector_sector/odd_joint_quotient_verification.json).

## Reproduction

Run the calculation chain listed in
[the source-frame method](odd-covariant-source-gauge-reduction.md), then

```text
.venv/bin/python scripts/verify_joint_odd_quotient.py
.venv/bin/python scripts/export_odd_covariant_basis.py
```

The finite source-probe script saves a lower bound even when it does not
saturate the symbolic bound. The subsequent Gram closure supplies the final
certificate. Replays verify accepted exact probes rather than repeat every
unsuccessful evaluation. All scripts use the repository `.venv` and SymPy.

# Interpretation and next step

The earlier **205** remains a vertex-tensor count, not an action count.
Restoring covariant Bianchi and ordered-derivative directions adds 36 raw
tensors, giving 241 raw tensors including the leading sector. Thirteen
joint chiral/EOM relations then leave 228 independent operators. This is why
subtracting a few EOM relations directly from 205 could not give the answer.

The count is source-general. Restricted photons, constant vacuum spurions,
canonical scalar legs or low-energy measurements identify fewer combinations.
Higher-field or multiple-source directions need not appear in the saved
radiative/TFF projections. Neither 217 nor 228 is a proposed data-fit dimension.

The operator-count milestone is complete. Next derive the normalized physical
SM and independent-source ALP kernels for this basis, determine observable
projections and blind directions, then calibrate only identifiable low-energy
combinations. No calibration or physical ALP prediction has been performed
by the reduction itself.

The [bounded Hilbert-series presentation](odd-bounded-hilbert-series.md)
independently reconstructs the raw multiplicities by character projections
and reproduces this quotient from its saved relation certificates. Its
[parent filtration](odd-hilbert-power-counting.md) has nested dimensions
12, 101 and 228. The quotient part inherits this record's completeness proof;
it is not a second independent Molien-Weyl proof of 228.
