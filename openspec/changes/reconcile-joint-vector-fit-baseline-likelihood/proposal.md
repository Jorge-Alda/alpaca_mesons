# Proposal

## Why

The repository contains multiple vector-sector fits with different data and conditioning, and their results can be confused. The global preliminary + conditional $G_V$ + BaBar likelihood currently terminates unsuccessfully at $\chi^2=307.985$ for 43 degrees of freedom, while the older BaBar-only inclusive fit reports $\chi^2=35.7164$ for 29 degrees of freedom; these are distinct analyses and must not be presented as competing runs of the same fit.

## What Changes

- Add an indexed numerical research record identifying the authoritative global likelihood JSON and executable script.
- Inventory its measurement rows/blocks, shared nuisances, coefficient coordinates, mass/source conventions, covariance assumptions, and exclusions of duplicate inputs.
- Record the optimizer's current nonconvergence, evaluation limit, and poor fit status separately from the earlier BaBar-only result.
- Update the relevant memory indexes so the record is discoverable from the repository's hierarchy.

## Capabilities

### New Capabilities

None.

### Modified Capabilities

None. This is a research-record and indexing change; it does not alter system behavior.

## Impact

Research memory under `memory/rcht/odd-parity/`, its parent indexes, and the existing fit output and script under `data/vector_sector/` and `scripts/`. No fit implementation or numerical result is changed by this documentation work.
