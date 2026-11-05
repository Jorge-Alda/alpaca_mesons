---
title: Joint vector fit deterministic basin probes
status: numerical
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Joint vector fit deterministic basin probes

## Inputs and conventions

The [high-order continuation](joint-vector-fit-continuation.md) supplies the seed endpoint at $(4,128,128)$ and its fixed optimizer scale vector. The [joint likelihood baseline](babar-kkpi-global-vector-sector-fit-baseline.md) defines the raw 73-summary likelihood, bounds, physical branch, and data assumptions. The [high-order quadrature record](kkpi-high-order-quadrature-certification.md) defines the fixed-candidate comparison tolerances.

## Method

The executable [`scripts/fit_joint_vector_basin_probes.py`](../../../scripts/fit_joint_vector_basin_probes.py) makes three deterministic, named starts by perturbing separate coordinate blocks around the high-order continuation point: even/shared-width coordinates (0–11), odd coordinates (12–23), and BaBar tensor/residue/normalization coordinates (24–31). Each start has scaled Euclidean radius 0.25 using the saved scale vector, alternating signs within its block. This radius is a search choice, not an uncertainty estimate. A failed domain check halves the displacement until the existing bounds and branch checks pass; no halving was needed here.

Each start was fitted at the default $(2,16,16)$ quadrature in four-bin blocks with trust-region reflective least squares, `max_nfev=1500`, and the established tolerances. Every successful, physically valid endpoint then received the certified fixed-candidate audit from $(2,64,64)$ through $(4,128,128)$. The exact starts, scale vector, source hash, endpoint diagnostics, optimizer outputs, and refinements are in [`joint_vector_basin_probes.json`](../../../data/vector_sector/joint_vector_basin_probes.json).

## Results

All three starts terminate successfully by `ftol`, pass the physical-domain and finite-residual checks, and have no active bounds. Their default-quadrature objectives are:

| Start | $\chi^2$ | `nfev` | `njev` | Optimality |
| --- | ---: | ---: | ---: | ---: |
| `even_shared_widths` | 225.6673947799526 | 252 | 243 | 0.8669837925352333 |
| `odd` | 225.66812358493632 | 258 | 246 | 0.7095565468768541 |
| `babar` | 225.6680429491065 | 288 | 275 | 0.7693212435860524 |

Each quadrature audit reports `certified`. For the combined $(2,64,64)\to(4,128,128)$ comparison, the largest prediction shifts are $0.000663$, $0.000661$, and $0.000661$ experimental standard deviations, respectively. The largest standardized-residual shifts are $0.000661$, $0.000659$, and $0.000659$; the absolute $\chi^2$ shifts are $0.003298$, $0.003320$, and $0.003319$. Every scaled-Jacobian column passes the $10^{-3}$ relative tolerance, and every observable comparison passes its stated tolerance.

The probe objectives are slightly below the default-quadrature continuation objective $225.6709550955029$. This comparison is at the same quadrature. The high-order continuation objective and default-order probe objectives are not compared as if they came from the same discretized likelihood. Optimality values remain nonzero despite successful `ftol` termination and are retained for the final endpoint comparison. These three starts test nearby basins only; they do not establish global optimality.

## Verification

Three focused tests passed:

```text
.venv/bin/python -m unittest scripts.test_fit_joint_vector_basin_probes
```

They check deterministic, block-confined starts, physical-domain validation, deterministic halving, scale and source provenance, and preparation serialization. The actual fit command was:

```text
.venv/bin/python scripts/fit_joint_vector_basin_probes.py --source data/vector_sector/joint_vector_continuation_e4p128a128.json --output data/vector_sector/joint_vector_basin_probes.json --max-nfev 1500 --refine-successful
```
