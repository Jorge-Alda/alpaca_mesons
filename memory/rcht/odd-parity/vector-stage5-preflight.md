---
title: Stage-5 data-only sensitivities and single-virtual asymptotics
status: mixed
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Executable preflight

## Inputs and assumptions

Use [stage-4 amplitudes](vector-stage4-amplitudes.md) and the checksum-pinned
[stage-4 rank result](vector-stage4-rank.md). Assumptions are fixed central
sectors, tree-level single-nonet action and independent scalar flavor traces
$q_i=q(T_i,Q)$, $h_i=h(T_i,Q)$, $q_{ij}=q(T_i,T_j)$,
$h_{ij}=h(T_i,T_j)$ as defined in that amplitude record. No commutativity of
abstract matrices is assumed; scalar traces encode the established bilinear
flavor calculation. Precision/source treatment is in the
[joint specification](vector-stage5-joint-fit-specification.md).

## Method and reproduction

Run [the executable](../../../scripts/diagnose_vector_stage5.py):

```text
.venv/bin/python scripts/diagnose_vector_stage5.py
```

The script checks all stage-4 input/kernel hashes, removes every inherited
OPE row, and performs unit-row-normalized SVD with tolerance $10^{-11}$.
It preserves all 64 branch/source/input scenarios. These are signed-amplitude
sensitivities, not width Jacobians, whitened measurement information, or
source-marginalized identifiability.

For the asymptotic derivation, insert the three explicit spectral denominators
into both VJP placements and every VVP pair in $F_P(-Q^2,0)$. SymPy computes
successive limits after subtracting each leading term, producing
$F_P(-Q^2,0)=a_PQ^2+b_P+k_P/Q^2+O(Q^{-4})$.
An independent series at $z=1/Q^2=0$ verifies all three coefficients.
The full input and all three coefficients are saved using `sympy.latex` in
[the result JSON](../../../data/vector_sector/stage5_preflight.json), together
with input/script hashes. Terminal output is restricted to the result path.

## Results

All 64 scenarios have data-only diagnostic ranks:

| Amplitude rows | Rank out of 11 |
|---|---:|
| Ten radiative channels | 7 |
| Above plus three real-photon TFFs | 7 |
| Above plus six single-virtual TFFs | 8 |

Thus virtual rows add one direction without inherited OPE sensitivities.
No claim about actual data errors or fitted uncertainty follows.

The script yields

$$
a_P=-\frac{4(c_1-c_2+c_5)}{FM_V}
\sum_i\frac{R_iq_i}{M_i^2},
$$

and, writing the saved constant coefficient compactly,

$$
b_P=\frac{\sqrt2N_Cq(Q,Q)}{8\pi^2F}
-\frac4{FM_V}(-c_1+c_2+c_5-2c_6)\sum_i R_iq_i
+\frac4{FM_V}\sum_i\frac{R_i}{M_i^2}
\left[q_i(c_1+c_2-c_5)m_P^2+4c_3h_i\right]
+\frac{4\sqrt2d_3}{F}\sum_{ij}\frac{R_iR_jq_{ij}}{M_j^2}.
$$

Vanishing $a_P$ and $b_P$ are necessary for this rational model to fall as
$1/Q^2$. They do not determine the physical leading coefficient, establish
QCD logarithms or constitute full VVP short-distance matching. In particular,
setting only $c_1-c_2+c_5=0$ removes the growing term but leaves a separate
channel-dependent constant condition. Neither condition is imposed as a
new precision prior in this execution.

## Verification and limit

The executable completed with all hash and symbolic-series assertions passing.
Both mass branches and all signs remain separate. No likelihood minimization
was run: the [joint specification](vector-stage5-joint-fit-specification.md)
identifies missing high-energy matching, radiative data and dispersive predictors.
