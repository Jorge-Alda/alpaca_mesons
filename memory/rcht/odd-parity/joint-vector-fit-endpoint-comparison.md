---
title: Joint vector fit common-order endpoint comparison
status: numerical
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Joint vector fit common-order endpoint comparison

## Inputs and method

The comparison uses the original four-start ledger, the default-order continuation, the $(2,64,64)$ and $(4,128,128)$ continuation ledgers, the high-order fixed-candidate audit, and the three-start basin-probe ledger. Their inputs and assumptions are described in the [joint fit baseline](babar-kkpi-global-vector-sector-fit-baseline.md), the [continuation record](joint-vector-fit-continuation.md), and the [basin-probe record](joint-vector-fit-basin-probes.md).

The executable [`scripts/compare_joint_vector_endpoints.py`](../../../scripts/compare_joint_vector_endpoints.py) validates ledger ancestry, input hashes, bounds, branch, scales, endpoint validity, and optimizer diagnostics. It evaluates all ten saved endpoints at common quadrature $(4,128,128)$ in four-bin blocks, ranks their common-order objectives, and computes all 45 pairwise Euclidean distances in the original optimizer scales. Presentation deduplication tolerances are 0.01 for scaled distance and 0.01 for absolute $\Delta\chi^2$. An endpoint is accepted in this comparison only if its optimizer succeeded, its domain check passed, and its own certified quadrature audit passed.

The complete ranked output, pairwise distances, sources and hashes are in [`joint_vector_endpoint_comparison.json`](../../../data/vector_sector/joint_vector_endpoint_comparison.json).

## Results

The best common-order objective is `basin_probes/even_shared_widths`, with $\chi^2=225.72669915915517$. Four endpoints pass all acceptance gates:

| Endpoint | Common-order $\chi^2$ | Optimizer optimality | Accepted |
| --- | ---: | ---: | --- |
| `basin_probes/even_shared_widths` | 225.72669915915517 | 0.8669837925352333 | yes |
| `basin_probes/babar` | 225.72708295708088 | 0.7693212435860524 | yes |
| `continuation_4128128/stored_negative_c4` | 225.72719113390542 | 1.3655702197555257 | yes |
| `basin_probes/odd` | 225.72720219590212 | 0.7095565468768541 | yes |

No pair of endpoints lies within both duplicate tolerances. The probe endpoints are close in objective at common quadrature, but their scaled distances exceed 0.01, so they remain separately recorded candidates. Their optimizer exits are successful by `ftol`; the nonzero optimality values are retained as fit diagnostics.

The original four starts remain unsuccessful and are not accepted. This finite search identifies several reproducible local candidates but does not establish global optimality.

## Verification

Three focused tests passed for common-order objective ranking, pairwise scaled distances, provenance validation, and the independent acceptance gates:

```text
.venv/bin/python -m unittest scripts.test_compare_joint_vector_endpoints
```

The common-order comparison was executed with:

```text
.venv/bin/python scripts/compare_joint_vector_endpoints.py --output data/vector_sector/joint_vector_endpoint_comparison.json
```

The final focused integration set passed, 52 tests total:

```text
.venv/bin/python -m unittest scripts.test_fit_joint_vector_multistart scripts.test_select_joint_vector_continuation scripts.test_fit_joint_vector_basin_probes scripts.test_compare_joint_vector_endpoints scripts.test_global_kkpi_running_propagator scripts.test_kkpi_quadrature_blocks scripts.test_kkpi_quadrature_certification scripts.test_certify_kkpi_quadrature
```

This includes the restart and basin-probe tests, the propagator/Jacobian and bounded-block regressions, and both quadrature-certification suites. The historical baseline and optimized fit JSON files are unchanged.
