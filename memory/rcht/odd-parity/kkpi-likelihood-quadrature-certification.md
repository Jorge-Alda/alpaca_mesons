---
title: KKπ likelihood quadrature certification at representative fitted-region points
status: numerical
created: 2026-11-04
updated: 2026-11-04
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# KKπ likelihood quadrature certification at representative fitted-region points

## Inputs and provenance

The likelihood, parameter conventions, and saved reference candidate are
described in [the joint vector-sector baseline](babar-kkpi-global-vector-sector-fit-baseline.md)
and implemented in
[`scripts/fit_babar_kkpi_global_vector_sector.py`](../../../scripts/fit_babar_kkpi_global_vector_sector.py).
The candidate is a nonconverged fit snapshot, not an established minimum. The
model and input provenance follow `Alda:2026nrd`, `Ruiz-Femenia:2003jdx`,
`ParticleDataGroup:2024cfk`, `BaBar:2007ceh`, and `Qin:2024ulb`, as recorded by
the baseline memory entry and its input file.

The reproducible implementation is
[`scripts/certify_kkpi_quadrature.py`](../../../scripts/certify_kkpi_quadrature.py).
Its machine-readable result is
[`data/vector_sector/kkpi_quadrature_certification.json`](../../../data/vector_sector/kkpi_quadrature_certification.json),
which stores the complete point vectors, source labels, hashes, invocation,
orders, tolerances, derivative settings, comparisons, and outcomes.

The selected points are the saved reference candidate and signed displacements
along the strongest and weakest singular directions of its column-normalized
residual Jacobian. The displacement radius is $0.25$ in the normalized
coordinates and is halved until the model accepts the point. These displacements
are numerical coverage points, not confidence intervals. No saved competing
minimum with the same likelihood coordinates and input hashes was available.

## Method

The saved fit grid starts at energy order $2$, pair order $16$ per pole-map
segment, and angular order $16$ per half. The certification compares these
orders independently with $(4,16,16)$, $(2,32,16)$, and $(2,16,32)$, then
compares the combined $(4,32,32)$ grid. The physical Dalitz domain and existing
$\phi$ and $K^*$ pole maps are unchanged.

For each point, the tool compares the largest prediction change in units of the
experimental uncertainty, the largest standardized-residual change, every
scaled residual-Jacobian column, and total $\chi^2$. The tolerances are
$0.01$, $0.01$, $10^{-3}$ relative per non-negligible Jacobian column, and
$0.01$, respectively. A negligible Jacobian column uses the recorded absolute
tolerance. Analytic odd-coefficient columns and forward differences for other
parameters use the same finite-difference steps at every order.

To reproduce the stored run, use:

```text
.venv/bin/python scripts/certify_kkpi_quadrature.py --output data/vector_sector/kkpi_quadrature_certification.json --energy-orders 2 4 --pair-orders 16 32 --angle-orders 16 32
```

The saved invocation was replayed to `/tmp/kkpi_quadrature_certification_replay.json`.
The replay reproduced all five point vectors, every comparison and outcome,
the counts, and the input hashes.

## Result

All five points are explicitly `uncertified` at the maximum orders
$(4,32,32)$; none is unevaluated. For the pair refinement, all points exceed
the prediction, residual, Jacobian, and $\chi^2$ tolerances. Across the five
points, the largest prediction changes range from $0.01829$ to $0.01905$ in
experimental standard deviations, standardized-residual changes from
$0.01547$ to $0.01610$, largest relative Jacobian-column changes from
$0.00392$ to $0.00399$, and $\chi^2$ changes from $0.21988$ to $0.22753$.

For the angular refinement, predictions, residuals, and Jacobian columns pass;
the $\chi^2$ changes range from $0.01385$ to $0.01552$ and fail the $0.01$
tolerance. The combined refinement also fails: prediction changes range from
$0.01749$ to $0.01822$, residual changes from $0.01480$ to $0.01539$, largest
relative Jacobian-column changes from $0.00382$ to $0.00388$, and $\chi^2$
changes from $0.20546$ to $0.21470$.

The result establishes only that the requested tolerances are not met at the
listed points and maximum orders. It is neither a quadrature pass nor a fit
quality or optimizer-convergence statement. The JSON retains the limiting
comparisons for each point so a higher-order run can continue the audit.

## Verification

The focused suites pass with 16 tests:

```text
.venv/bin/python -m unittest scripts.test_kkpi_quadrature_certification scripts.test_certify_kkpi_quadrature
```

The pole-map tests independently integrate the physical Dalitz variables and
compare the $\phi$, $K^*_t$, and $K^*_u$ pole factors against the mapped
quadrature. The certification replay also verifies that stored point selection
and numerical outcomes are reproducible from their recorded invocation and
software/input hashes.
