---
title: Input contract for the VPPP finite-width dispersive profile
status: mixed
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Hoferichter:2019mqg, Hoferichter:2023bjm, Hoferichter:2025lcz, BABAR:2021cde]
provenance: numerical-computation
---

# Implemented profile layer

The [two-body width builder](../../../scripts/build_vppp_partial_widths.py)
calculates conditional $V\to P\gamma$ and $V\to P_1P_2$ channel shapes from
caller-supplied on-pole widths. The [three-pion width builder](../../../scripts/build_vppp_three_pion_width.py)
converts a caller-supplied full KT scalar grid to a normalized $3\pi$ width.
The [spectral template builder](../../../scripts/build_vppp_spectral_templates.py)
then calculates per-residue finite-width resonance and conformal absorptive
templates from explicit inputs, following [@Hoferichter:2019mqg]. The
[profile calculator](../../../scripts/profile_vppp_dispersive_basis.py)
accepts those templates, conformal spectral
templates, quadrature weights, the previously derived subtraction/residue
maps, complex evaluation points, a reference-prediction residual, a local
observation-response Jacobian, and its covariance. It constructs the
once-subtracted basis and enforces the linearized form of the
[@Hoferichter:2025lcz] spectral sum rule exactly. It profiles free conformal
directions by covariance whitening and SVD projection, then reports the
allowed coefficient subspace, rank, singular values, minimum-norm local
shift, and quadratic profile curvatures.

The code's result is a local linear profile around the caller's reference
prediction, not a nonlinear likelihood fit. For cross-section data,
`observation_response` must be the Jacobian of a complete bin-integrated
Khuri--Treiman prediction with respect to the real and imaginary parts of
$a(q^2)$ at that reference point. The data themselves do not supply this
Jacobian. The nuisance-profiled rank is not reported as a data rank when a
direction is fixed only by exact basis constraints.

# Executed zero-width limit

SymPy integrates a delta-function spectral pole and reproduces the exact
once-subtracted pole in the
[conditional matching record](isoscalar-vppp-dispersive-matching.md). This
checks the transform convention and subtraction at $q^2=0$; it does not
validate any finite-width spectral model.

# Inputs still required for a physical profile

- Rescattered full KT scalar grids, on-pole $3\pi$ partial widths, and the
  source-backed two-body on-pole partial widths needed by the width builders.
  No physical $\Gamma_V(s)$ grid is currently stored. High-energy tails,
  thresholds, source overlaps, and the $\rho$--$\omega$ correction must match
  the selected fit prescription.
- A complete complex Khuri--Treiman solution and its bin-integrated
  $e^+e^-\to3\pi$ prediction, including the declared final-state-radiation
  convention, to construct `observation_response` and the reference
  residual.
- A joint likelihood input. The extracted
  [BaBar 2021 cross section](babar-2021-three-pion-fit-input.md) supplies
  114 low-energy bins with two covariance matrices and 72 high-energy bins
  with a separate error schema. The SND and CMD-2 datasets and any
  cross-experiment covariance/nuisance treatment remain to be connected.
- The nonlocal source and resonance-residue remainder, or an explicit
  external constraint for it. Otherwise total fitted pole residues absorb
  the local shifts.

The calculator has algebraic fixtures for constraint elimination,
correlated-covariance whitening, nuisance profiling, and exact degeneracy.
Those fixtures are software checks only and are not physical data or spectral
templates. Until the listed inputs are supplied, there are no fitted values,
physical singular values, or parameter intervals.

Run the symbolic limit with:

    .venv/bin/python scripts/profile_vppp_dispersive_basis.py --derive-zero-width

The full calculator requires an explicit NPZ file matching the schema in the
script docstring; it refuses to run without that file.
