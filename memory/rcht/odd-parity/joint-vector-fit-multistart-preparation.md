---
title: Joint vector fit deterministic starts and scaling
status: numerical
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Joint vector fit deterministic starts and scaling

## Inputs, conventions, and assumptions

The [joint fit baseline](babar-kkpi-global-vector-sector-fit-baseline.md) specifies the 32 coordinates, 75 residuals, units, source conventions, raw measurements, and physical assumptions [@Alda:2026nrd; @Ruiz-Femenia:2003jdx; @ParticleDataGroup:2024cfk; @BaBar:2007ceh; @Qin:2024ulb]. Masses and widths are in GeV and cross sections in nb. The shared even-sector constraints are described in the [audit](shared-even-sector-gv-constraint-audit.md); the efficient residual and derivative calculations are described in the [channel-specific Jacobian record](kkpi-global-fit-channel-jacobian.md).

## Method

The executable [`scripts/fit_joint_vector_multistart.py`](../../../scripts/fit_joint_vector_multistart.py) reads the stored nonconverged baseline as one starting point. It also uses the preliminary starting point and two variants of the stored point, with the sign of $c_4$ or of all twelve odd coordinates reversed. These choices are deterministic search inputs, not additional constraints. Each must satisfy the original coordinate bounds and the likelihood's checks of the inverse-mass and source branches, positive widths and $G_V$, and residue domain. The bounds are returned by `parameter_bounds` in the original fit script.

The initial numerical scale of each optimizer coordinate is the reciprocal of its residual-Jacobian column norm at the stored point, evaluated with the chosen quadrature; that vector is held fixed for all starts. The even physical coordinates continue to use the six `EVEN_SCALES` factors of the baseline fit. The scale choice affects the numerical search and does not alter the likelihood.

The executable records starts, scales, bounds, physical branch checks, quadrature and optimizer settings, direct input and code hashes, and Python/library versions in a separate JSON file. If requested to run the fits, it records each termination and preserves distinct and unsuccessful endpoints. No endpoint is marked accepted by this preparation procedure: acceptance requires the subsequent numerical refinement.

## Verification

The focused unit tests exercise reproducible starts, physical and rectangular-domain rejection, original bounds, deterministic scales, and preservation of distinct successful and unsuccessful simulated optimizer outcomes:

```text
.venv/bin/python -m unittest scripts.test_fit_joint_vector_multistart
```

The preparation command below checks serialization at a small quadrature order. This order is used only for the preparation check; the [high-order quadrature certification](kkpi-high-order-quadrature-certification.md) remains the reference for candidate refinement. No fit campaign or accepted minimum is reported here.

```text
.venv/bin/python scripts/fit_joint_vector_multistart.py --prepare-only --energy-order 1 --pair-order 3 --angle-order 4 --output /tmp/joint-vector-starts.json
```
