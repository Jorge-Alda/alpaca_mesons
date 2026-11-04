---
title: High-order KKπ likelihood quadrature certification
status: numerical
created: 2026-11-04
updated: 2026-11-04
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# High-order KKπ likelihood quadrature certification

## Inputs and provenance

The model and saved reference candidate are documented in
[the global vector-sector baseline](babar-kkpi-global-vector-sector-fit-baseline.md).
The earlier order-32 result and the deterministic five-point selection are in
[the initial quadrature certification](kkpi-likelihood-quadrature-certification.md).
The new run reused those exact point vectors and selection labels; its source
hash is stored in the high-order result. No same-input competing minimum was
available. The reference fit itself remains a nonconverged snapshot, and a
quadrature pass makes no optimizer or fit-quality claim.

The executable is
[`scripts/certify_kkpi_quadrature.py`](../../../scripts/certify_kkpi_quadrature.py),
with bounded-memory energy-bin evaluation in
[`scripts/fit_babar_kkpi_global_vector_sector.py`](../../../scripts/fit_babar_kkpi_global_vector_sector.py).
The full result is
[`data/vector_sector/kkpi_quadrature_certification_high_order.json`](../../../data/vector_sector/kkpi_quadrature_certification_high_order.json).
Its input and software hashes, comparison details, outcomes, exact command, and
historical-point source are authoritative.

## Method

The new qualification starts at $(E,P,A)=(2,64,64)$ and ends at $(4,128,128)$,
where $E$ is energy order, $P$ is pair order per pole-map segment, and $A$ is
angular order per half. It reports four comparisons: energy $2\to4$ with
$(P,A)=(64,64)$ fixed; pair $64\to128$ with $(E,A)=(2,64)$ fixed; angle
$64\to128$ with $(E,P)=(2,64)$ fixed; and the combined
$(2,64,64)\to(4,128,128)$ comparison. The earlier $16\to32$ failures remain in
the initial result as historical diagnostics and are not mixed into this
high-order qualification's status.

The calculation uses deterministic energy-bin blocks of four bins and retains
the existing physical limits and $\phi$ and $K^*$ pole maps. Blocked predictions,
all 75 residuals, and all 32 residual-Jacobian columns match the full-grid
implementation bit-for-bit at supported low and moderate orders. The tests
also construct the order-128 pole maps one bin at a time and check their
physical bounds.

The invocation recorded in the JSON is:

```text
.venv/bin/python scripts/certify_kkpi_quadrature.py --points-from data/vector_sector/kkpi_quadrature_certification.json --energy-orders 2 4 --pair-orders 64 128 --angle-orders 64 128 --bin-block-size 4 --output data/vector_sector/kkpi_quadrature_certification_high_order.json
```

## Result

All five points—reference candidate, strong-minus, strong-plus, weak-minus,
and weak-plus—are `certified`; none is unevaluated. The unchanged tolerances
are $0.01$ for prediction changes in experimental standard deviations, $0.01$
for standardized residuals, $10^{-3}$ relative for each non-negligible
scaled-Jacobian column, and $0.01$ for total $\chi^2$.

The largest changes over every point in each comparison family are:

| Comparison | Prediction change in $\sigma$ | Standardized residual change | Relative Jacobian-column change | $\chi^2$ change |
| --- | ---: | ---: | ---: | ---: |
| Energy $2\to4$ | $1.6665\times10^{-3}$ | $1.4096\times10^{-3}$ | $3.3189\times10^{-5}$ | $1.1446\times10^{-3}$ |
| Pair $64\to128$ | $1.0372\times10^{-6}$ | $8.7729\times10^{-7}$ | $8.2261\times10^{-7}$ | $2.3755\times10^{-7}$ |
| Angle $64\to128$ | $8.7117\times10^{-7}$ | $7.3689\times10^{-7}$ | $7.4437\times10^{-7}$ | $4.3369\times10^{-7}$ |
| Combined | $1.6665\times10^{-3}$ | $1.4096\times10^{-3}$ | $3.3189\times10^{-5}$ | $1.1453\times10^{-3}$ |

Every reported comparison is below its tolerance. This resolves the quadrature
accuracy acceptance criterion for the recorded point set and orders only; it
does not change the nonconverged status of the reference fit.

## Verification

The focused suites pass with 21 tests:

```text
.venv/bin/python -m unittest scripts.test_kkpi_quadrature_blocks scripts.test_kkpi_quadrature_certification scripts.test_certify_kkpi_quadrature
```

The recorded invocation was replayed with:

```text
.venv/bin/python scripts/certify_kkpi_quadrature.py --replay data/vector_sector/kkpi_quadrature_certification_high_order.json --output /tmp/kkpi_quadrature_certification_high_order_replay.json
```

The replay reproduced the selected point vectors, all comparisons and
outcomes, counts, orders, input hashes, and software source hashes exactly.
