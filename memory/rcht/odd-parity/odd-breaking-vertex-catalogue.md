---
title: Stage-A bounded odd vertex catalogue and handoff to normalized amplitudes
status: mixed
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Roig:2013baa, Kaiser:2000gs]
provenance: analytical-derivation
---

# Stage-A output

## Inputs and explicit truncation

This executes stage A of the [selected plan](odd-breaking-low-energy-calibration-plan.md).
Its two executable components are
[flavor enumeration](odd-breaking-flavor-enumeration.md) and
[Lorentz/IBP reduction](odd-breaking-lorentz-reduction.md).
The scope is one vector nonet, one odd $u_\mu$ or $i\chi_-$ block (or
$f_-$), at most one quark spurion, and at most two flavor traces.
The derivative ceiling is two derivatives supplied by $u_\mu$ plus an
additional derivative, or by two explicit derivatives on a scalar source
block; the $f_+$ field strength carries its usual extra source derivative.

The script enumerates and reduces this local source-resolved vertex space.
It does not assert completeness for arbitrary nonlinear $U(3)$ source
functions, multiple pseudoscalars, triple traces, extra resonance multiplets
or an EOM-minimal generating functional. No Goldstone/source EOM or
short-distance constraint is imposed. The one-spurion mass terms with two
additional derivatives are retained, since omitting them while keeping
derivative-$\chi_+$ terms at the same ceiling would be an incomplete
enumeration of this chosen block space.

## Reproduction and inventory

```text
.venv/bin/python scripts/enumerate_odd_breaking_vertices.py
```

The [script](../../../scripts/enumerate_odd_breaking_vertices.py) produces
[a step-by-step transcript](../../../data/vector_sector/odd_vertex_enumeration.md)
and [the machine-readable inventory](../../../data/vector_sector/odd_vertex_enumeration.json).
Each entry has a stable block/flavor/Lorentz identifier, explicit abstract
trace, metric-pair and derivative placement, Bose parity, Hermitian prefactor,
operator mass dimension, coefficient normalization, spurion count and source
signature. Lorentz polynomials and all raw-candidate reduction coefficients
are referenced from the common block tables rather than duplicated.
The original eleven structures are included; their coefficient convention
is recovered using the stored derivative-seed coordinates and the trace
normalization stated in the flavor record.

## Reduced source-resolved counts

| Block | Single trace | Double trace | Total | Mass dimension | Dimensionless-coefficient normalization |
|---|---:|---:|---:|---:|---|
| Derivative VJP, no spurion | 5 | 15 | 20 | 5 | $1/M_V$ |
| Derivative VVP, no spurion | 3 | 9 | 12 | 4 | $1$ |
| VJP $i\chi_-$ | 1 | 3 | 4 | 5 | $1/M_V$ |
| VVP $i\chi_-$ | 1 | 2 | 3 | 4 | $1$ |
| VJP $f_-\chi_+$ | 1 | 0 | 1 | 5 | $1/M_V$ |
| Derivative VJP with $\chi_+$ | 24 | 56 | 80 | 7 | $1/M_V^3$ |
| Derivative VVP with $\chi_+$ | 13 | 30 | 43 | 6 | $1/M_V^2$ |
| VJP $i\chi_-$, two extra derivatives | 5 | 15 | 20 | 7 | $1/M_V^3$ |
| VVP $i\chi_-$, two extra derivatives | 4 | 11 | 15 | 6 | $1/M_V^2$ |
| VJP $f_-\chi_+$, two extra derivatives | 7 | 0 | 7 | 7 | $1/M_V^3$ |
| Total | 64 | 141 | 205 | — | — |

These are coefficient tensors of the specified source-resolved block space,
not 205 independent physical fit parameters or a count for every nonlinear
$U(3)$ operator. The first five single-trace rows give the original eleven
coefficients. Zero-spurion derivative and leading mass double traces add
29 tensors; the higher-dimension one-spurion blocks contribute 165.
Derivative-$\chi_+$ double traces have leading hard derivative numerators;
the leading $i\chi_-$ mass double trace alone does not cancel a hard TFF
constant. No cancellation prior is imposed on either family.

## Target/channel and source audit

- Derivative $u_\mu$ blocks enter SM radiative/TFF amplitudes and derivative
  ALP insertions, subject to channel flavor traces. Extra traces access
  singlet slots that the original common single trace cannot distinguish.
- $i\chi_-$ blocks carry mass and independent pseudoscalar-density/ALP source
  insertions. Their derivative extensions cannot be replaced by a measured
  on-shell mass before differentiating the external sources.
- $f_-\chi_+$ blocks require their own Goldstone/photon expansion; the old
  charged $c_4$ term is recovered. Direct diagonal-neutral $f_-$ terms vanish
  under the electromagnetic commuting-generator specialization, while
  charged channels need not vanish.
- Spurion-gradient directions are invisible for a constant mass background;
  their source-sensitive role remains open. Photon trace identities and
  scalar/canonical-tensor projections further collapse the catalogue.
- No independent multi-pseudoscalar contact basis has been derived here.
  The $a\to3P$ audit must classify the retained covariant representatives and
  any nonlinear contact remainders instead of extrapolating these counts.

No data-only rank, numerical prior, calibrated coupling or unique ALP
prediction is produced. Nominal extra-trace OZI suppression is bookkeeping,
not a Gaussian penalty [@Kaiser:2000gs].

## Verification

The enumeration's exact matrix-reconstruction, IBP, exchange, published-seed,
singlet and scalar-projection assertions pass. In addition, run

```text
.venv/bin/python scripts/verify_odd_vertex_enumeration.py
```

This [verification script](../../../scripts/verify_odd_vertex_enumeration.py)
tests all 15 metric-pair contractions with an exact rational proper-Lorentz
boost, an orientation-changing reflection, and off-shell photon Ward
replacement. It tests the polarized trace identity with concrete noncommuting
$U(3)$ flavor directions and checks catalogue hashes/counts. The
[verification result](../../../data/vector_sector/odd_vertex_verification.json)
pins its script and enumeration hashes. These are algebra/invariance checks,
not evidence of data adequacy or a normalized ALP prediction.

## Stage-B handoff and completion limits

The first stage-B subsection is now executed in the
[source/EOM audit](odd-trace-source-eom-audit.md) and
[leading normalized scalar kernels](odd-trace-scalar-kernels.md).
The audit verifies three neutral-linear divergence/mass relations but does
not establish a full EOM-minimal count. The leading scalar EM kernel space
adds 15 formal directions; this is distinct from action-operator independence.
The complete bounded one-spurion action count is the next pre-fit milestone.

The local vertex enumeration/reduction requested in stage A is executed.
Before fitting, stage B must choose covariant representatives, expand the
full blocks with independent Alda:2026nrd sources and account for nonlinear contact
remainders not fixed by equality in this vertex sector. EOM conversion of
singlet/double traces must retain the Alda:2026nrd anomalous singlet and source terms;
the literature basis dictionary is an audit, not permission to discard them.

Begin stage B by deriving exact normalized radiative/TFF kernels for the
zero-spurion single/double-trace sector, then add one-spurion derivative,
mass/source and $f_-$ extensions with the same source accounting. Form a
combined sensitivity matrix and identify observable combinations before
assembling a fit. Constant-singlet rescalings and photon-blind traces must
not be fitted independently from the lower-order coefficients to the same
fixed-input dataset. New SM-blind but ALP-sensitive combinations remain
explicitly uncalibrated until additional information or a stated assumption
is supplied. The old blind-direction result is not inherited unchanged.
