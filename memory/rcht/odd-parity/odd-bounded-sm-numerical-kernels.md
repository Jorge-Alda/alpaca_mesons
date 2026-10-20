---
title: Numerical SM gradients of all bounded odd-action parents
status: numerical
created: 2026-10-14
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx]
provenance: numerical-computation
---

Canonical tensor and odd-coupling conventions follow [the Levi-Civita audit](levi-civita-convention.md). The linked original scripts and saved artifacts use the explicitly identified reference backend; apply that record's signed coefficient and tensor adapter when interpreting them in the canonical convention.

# Inputs and conventions

Use the exact-derived [compact cubic catalogue](../../../data/vector_sector/odd_cubic_parent_kernels.json),
the [charged expansion](odd-charged-cubic-kernels.md),
[normalized neutral source maps](odd-power-counted-channel-maps.md),
the [photon/action phase dictionary](../../alp/photon/alp-photon-odd-action-phase-matching.md),
and the [signed spectral response](vector-stage4-amplitudes.md).
The retained action and leading nonet WZW term are assumptions; their
normalization and source signs follow the linked records and Alda:2026nrd conventions.
Natural units are GeV, metric $(+---)$ and source $\epsilon^{0123}=+1, \epsilon_{0123}=-1$.
All raw parent coefficients are dimensionless after their declared $M_V$
powers; 241 columns represent the bounded 228-dimensional covariant quotient.
No OPE or coefficient prior is inserted.

# Implementation

[`odd_bounded_amplitudes.py`](../../../scripts/implementation/odd_bounded_amplitudes.py)
evaluates the saved rational flavor words and Lorentz polynomials using
concrete physical matrices. It validates the exact-input artifact hashes.
`radiative_gradients` returns $\partial G/\partial g_j$ for each raw parent,
including composite $f_-$ momentum and canonical external $1/M_{\rm ext}$.
`tff_gradients` returns $\partial F/\partial g_j$ with both photon placements,
the VVP assignment factors and the signed $W(s)$ response. `wzw_offset`
supplies the inherited leading nonet anomaly in the same scalar phase.
The explicit photon $Q$ is never replaced by $W$.

`original_embedding` reads the saved exact original $c_1,\ldots,c_7,d_1,
\ldots,d_4$ coordinate map without eliminating any coupling. The two parent
filtrations are applied before amplitude reduction; `retained_columns`
implements soft counting and conservative hard-vector derivative promotion.

The reusable input is the established `VectorAmplitudeInputs`, with explicit
$T_P,\chi,F,m_P,M_V$ and sector. Local kernels must use the LO source-EOM
mass to represent the same covariant quotient. A physical-mass replacement
inside the polynomial requires a separate EOM-remainder audit. Physical
phase-space factors may be replaced separately only as a stated truncation
assumption. No physical-mass replacement changes the quark spurion.

# Independent verification

Run

```text
.venv/bin/python scripts/verify_odd_bounded_amplitudes.py
```

The [verification artifact](../../../data/vector_sector/odd_bounded_amplitude_verification.json)
records 42 original-action channel checks against the independent legacy
amplitudes, including two nonideal angles, split poles, opposite source signs,
charged legs and asymmetric two-virtual-photon kinematics. It also checks
723 charged raw-column values against the prior SymPy artifact, 4338
photon-parent Bose identities, original embedding rank eleven, and nested
parent-filter inclusion. The executable also reconstructs the direct PrimEx
pion-width variance from its stated statistical and systematic errors.
These are implementation checks, not data ranks,
model adequacy or uncertainty estimates.
