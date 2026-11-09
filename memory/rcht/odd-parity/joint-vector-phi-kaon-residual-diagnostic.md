---
title: Joint vector fit phi to kaon residual diagnostic
status: numerical
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# $\phi\to K\bar K$ residual diagnostic

## Inputs and assumptions

The accepted `basin_probes/even_shared_widths` endpoint and its residuals at quadrature $(4,128,128)$ are in the [goodness-of-fit record](joint-vector-fit-goodness-of-fit.md). The [likelihood baseline](babar-kkpi-global-vector-sector-fit-baseline.md) gives the fit assumptions and residual ordering. The [shared even-sector audit](shared-even-sector-gv-constraint-audit.md) identifies the two branching fractions and their correlation $-0.72$ [@ParticleDataGroup:2024cfk]. The calculation reads the residuals from `data/vector_sector/joint_vector_fit_adequacy.json` and the measured fractions and covariance from `data/vector_sector/gv_measurements.json`. All fractions and residuals are dimensionless. The neutral measurement is $K_L K_S$, identified with the neutral-kaon channel under the baseline conventions.

The positive symmetric square root of the covariance is used to choose one rotation of the whitened residuals. This choice treats both channel labels symmetrically. With correlated measurements, no unique decomposition of $\chi^2$ into channel contributions exists; the rotated entries each mix the physical branching-fraction differences. The ratio uncertainty below uses first-order propagation at the measured central fractions and is an approximation.

## Method

1. In the order $K^+K^-$, $K^0\bar K^0$, recover prediction minus measurement by multiplying the stored Cholesky-whitened vector by the lower Cholesky factor of the measured covariance.
2. Multiply this difference by the positive symmetric inverse square root of the covariance. The product of that inverse square root and the Cholesky factor is an orthogonal matrix. Thus the new vector has the same squared length as the stored vector. Its entries follow the channel ordering, but each combines both physical differences.
3. Independently calculate the residual in the sum and the change in the charged-to-neutral ratio. Divide the sum by its measured standard deviation and the ratio change by the linearly propagated measured standard deviation. These two standardized quantities are correlated and their squares are not a partition of the block $\chi^2$.

## Result

At the stated endpoint, the measured branching fractions are $(0.491,0.339)$ and the predictions are $(0.5084953,0.3355123)$. Prediction minus measurement is $(+0.0174953,-0.0034877)$. The stored Cholesky-whitened residuals $(3.499051,2.373849)$ rotate to $(4.105707,1.010785)$ in the symmetric convention, with squared length $17.878514$, matching the recorded block $\chi^2$.

The sum difference is $+4.01035$ measured standard deviations. The predicted minus observed charged-to-neutral ratio is $+2.27495$ standard deviations using the linearized ratio uncertainty. These comparisons describe the fixed endpoint; a refit is required to assess an alternative physical model.

## Reproduction and verification

```text
.venv/bin/python scripts/diagnose_joint_phi_kaon_residuals.py
.venv/bin/python -m unittest scripts.test_diagnose_joint_phi_kaon_residuals
```

The [script](../../../scripts/diagnose_joint_phi_kaon_residuals.py) verifies agreement with the stored block $\chi^2$. The unit tests independently check the covariance-weighted quadratic form, orthogonality, channel interchange, and the accepted endpoint. The total width and fit parameters are not refitted here.
