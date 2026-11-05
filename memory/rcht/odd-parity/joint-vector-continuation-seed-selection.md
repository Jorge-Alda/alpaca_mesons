---
title: Joint vector fit continuation seed selection
status: numerical
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Joint vector fit continuation seed selection

## Inputs and conventions

The [four-start result](joint-vector-fit-multistart-results.md) defines the candidate endpoints and their optimizer outcomes. The likelihood, coordinate meanings, physical branch, and data assumptions are documented in the [joint fit baseline](babar-kkpi-global-vector-sector-fit-baseline.md). The recorded optimizer-coordinate scales are used for comparing endpoint separations; they do not change the likelihood.

## Method

The executable [`scripts/select_joint_vector_continuation.py`](../../../scripts/select_joint_vector_continuation.py) reads the saved [multi-start ledger](../../../data/vector_sector/joint_vector_multistart.json). It requires exactly the four expected starts, 32 finite coordinates, finite objectives, positive finite scales, and a valid physical-domain flag for every endpoint. It ranks endpoints by $\chi^2$ and calculates all six pairwise Euclidean distances after dividing each coordinate difference by the recorded scale vector. The selected endpoint is copied without numerical changes, and the source ledger SHA-256 is stored in the separate [selection record](../../../data/vector_sector/joint_vector_continuation_selection.json).

## Result

`stored_negative_c4` has the lowest objective, $\chi^2=226.19481021926364$, among the four physically valid endpoints. Its optimizer status remains unsuccessful, and the selection record explicitly keeps `accepted=false`. It is selected only as the initial vector for continuation, not as a converged minimum.

The pairwise scaled distances are approximately 203.238 between `preliminary` and `stored`, 266.592 between `preliminary` and `stored_negative_c4`, 484.376 between `preliminary` and `stored_negative_odd`, 149.075 between `stored` and `stored_negative_c4`, 484.505 between `stored` and `stored_negative_odd`, and 552.404 between `stored_negative_c4` and `stored_negative_odd`. These distances describe separation in the chosen numerical optimizer coordinates; they do not establish basin boundaries or global coverage.

## Verification

The recorded selection command is:

```text
.venv/bin/python scripts/select_joint_vector_continuation.py --ledger data/vector_sector/joint_vector_multistart.json --output data/vector_sector/joint_vector_continuation_selection.json
```

Three focused tests passed:

```text
.venv/bin/python -m unittest scripts.test_select_joint_vector_continuation
```

They verify endpoint ordering, all six scaled distances, physical-validity rejection, exact seed-vector preservation, unsuccessful status retention, deterministic output, and source-ledger immutability.
