---
title: Joint vector fit continuation from the best multi-start endpoint
status: numerical
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Joint vector fit continuation from the best multi-start endpoint

## Inputs and conventions

The [multi-start result](joint-vector-fit-multistart-results.md) documents the four original endpoints, likelihood inputs, scales, and covariance assumptions. The [seed selection record](joint-vector-continuation-seed-selection.md) identifies `stored_negative_c4` as the lowest-objective physically valid endpoint, while retaining its unsuccessful optimizer status. The fit preserves the same raw 73-summary likelihood, bounds, physical branch, fixed coordinate scales, and runtime-improved bounded-bin evaluator.

## Method

The continuation starts from the exact 32-coordinate endpoint saved under `stored_negative_c4`. It uses energy order 2, pair order 16, angle order 16, four-bin blocks, trust-region reflective least squares, the recorded scale vector, and `max_nfev=5000`. This is a new optimizer call; the attempt records its own `nfev` and the cumulative count including the prior run separately. Iteration checkpoints preserve coordinates, objective, iteration, and evaluation counts. SciPy 1.17.1 does not expose optimality in this callback, so checkpoint `optimality` is null; the final optimizer result supplies it.

The source and continuation ledgers are [`joint_vector_multistart.json`](../../../data/vector_sector/joint_vector_multistart.json) and [`joint_vector_continuation.json`](../../../data/vector_sector/joint_vector_continuation.json). The separate progress record is [`joint_vector_continuation.progress.json`](../../../data/vector_sector/joint_vector_continuation.progress.json). Input and source hashes, including the prior-ledger hash and both continuation-driver hashes, are stored in the continuation ledger.

## Result

The continuation terminates successfully with SciPy status 2 (`ftol` condition), after $n_{\rm fev}=4548$ new function evaluations and $n_{\rm jev}=4313$ Jacobian evaluations. Including the previous 1500 function evaluations, the recorded cumulative count is 6048. The objective decreases from $226.19481021926364$ to $225.6709550955029$, a reduction of $0.52385512376074$. The endpoint passes the physical-domain and finite-residual check, and no coordinate is reported as active at a bound.

The final optimality is $0.6818041314705354$, so successful termination by the relative cost criterion alone does not establish a strongly stationary point. The block contributions are $35.128275$ (odd), $12.929417$ (masses and leptonic widths), $106.835220$ (shared widths), $17.883384$ ($\phi\to K\bar K$), $52.831505$ (BaBar), and $0.063154$ (normalization constraints).

## Certified-order refits

The same endpoint was refit sequentially with four-bin blocks at $(2,64,64)$ and $(4,128,128)$, keeping the original bounds and scale vector. Both fits terminate successfully by `ftol`, pass the independent domain check, and have no active bounds. The $(2,64,64)$ fit has $\chi^2=225.72462151174813$, `nfev=439`, `njev=408`, cumulative `nfev=6487`, and optimality $1.1855750576651687$. The $(4,128,128)$ fit has $\chi^2=225.72719113390542$, `nfev=37`, `njev=29`, cumulative `nfev=6524`, and optimality $1.3655702197555257$. Their separate ledgers are [`joint_vector_continuation_e26464.json`](../../../data/vector_sector/joint_vector_continuation_e26464.json) and [`joint_vector_continuation_e4p128a128.json`](../../../data/vector_sector/joint_vector_continuation_e4p128a128.json).

At the $(4,128,128)$ endpoint's fixed coordinates, the certified audit comparing energy $2\to4$, pair $64\to128$, angle $64\to128$, and the combined $(2,64,64)\to(4,128,128)$ orders passed every recorded tolerance. In the combined comparison, the largest prediction change is $0.0006607172$ experimental standard deviations, the largest residual change is $0.0006590298$, and the absolute $\chi^2$ change is $0.0033267794$. All scaled-Jacobian columns pass. The full audit is in [`joint_vector_continuation_e4p128a128_refinement.json`](../../../data/vector_sector/joint_vector_continuation_e4p128a128_refinement.json).

These results satisfy the optimizer-success, domain, and fixed-candidate quadrature-comparison gates. The non-small gradient-based optimality values remain explicit diagnostics; the finite additional-start study is still needed before the final candidate summary.

## Verification

The default-order continuation command was:

```text
.venv/bin/python scripts/fit_joint_vector_multistart.py --resume-from data/vector_sector/joint_vector_multistart.json --resume-start stored_negative_c4 --max-nfev 5000 --output data/vector_sector/joint_vector_continuation.json
```

The restart-preparation tests passed, and a preparation-only invocation confirmed exact seed-vector and scale preservation, source hash capture, and the 1500-to-5000 evaluation settings. The full continuation output records optimizer status, final coordinates and residuals, endpoint validity, and progress checkpoints. The higher-order refit commands and audit inputs, statuses, and outcomes are recorded in the linked JSON ledgers.
