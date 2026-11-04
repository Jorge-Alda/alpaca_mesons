# Design

## Context

The global artifact is `data/vector_sector/babar_kkpi_global_vector_sector_fit.json`, produced by `scripts/fit_babar_kkpi_global_vector_sector.py`. Its stored likelihood has 73 measurement summaries and two normalization constraints, 32 fitted parameters, and 43 degrees of freedom. It is separate from the ten-combination inclusive BaBar fit documented in `babar-kkpi-unified-vector-coefficient-fit.md`.

## Goals / Non-Goals

**Goals:** Create one auditable memory record for the joint likelihood and link it through the odd-parity and parent topic indexes. Preserve the stored output as the authority for the reported numerical state.

**Non-Goals:** Rerun or alter the optimizer, revise data, change fit conventions, or claim a converged best fit or calibrated coefficient intervals.

## Decisions

- Treat the global JSON and its named script as the authoritative artifact pair for this combined likelihood. The memory record will state its input checksum metadata and reproduction command, so any regenerated result can be distinguished from the currently stored one.
- Inventory rows by the output's six reported blocks: 12 odd summaries, 7 mass/leptonic summaries, 7 shared total widths, 2 phi-to-kaon branching summaries, 45 inclusive BaBar bins, and 2 normalization constraints. Identify the 32 fitted coordinates, including five shared-width pulls, the conditional $G_V$ coordinate, twelve odd combinations, even-sector coordinates, tensor products, excited-pole residues, and normalization pulls.
- Describe source and mass conventions and covariance approximations directly from the script and JSON. Explicitly record that shared widths are counted once and that the BaBar Table 7 isobar component is excluded as a duplicate of events represented by inclusive BaBar rates.
- Record `success=false`, 1500 function evaluations, the maximum-evaluations message, $\chi^2=307.9848522$, 43 dof, and the optimizer p-value/model-adequacy fields. Keep the older BaBar-only value in a separate comparison callout, with an explicit statement that it is not this likelihood.

## Risks / Trade-offs

- [The JSON may be regenerated and its status may change] → Include the fit output path, current stored checksum/provenance fields, and date of record; label all reported results as the current stored snapshot.
- [Summary-level row inventory could hide input provenance] → Link the underlying JSON, script, and existing indexed records for preliminary and conditional $G_V$ inputs; use exact block counts and avoid claiming unavailable covariance matrices.

## Migration Plan

Add the atomic memory record, then update `memory/rcht/odd-parity/INDEX.md` and `memory/rcht/INDEX.md` as needed. No data migration or rollback is required; remove the new record and its index links to revert.
