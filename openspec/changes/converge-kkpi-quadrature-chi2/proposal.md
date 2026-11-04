# Proposal

## Why

The stored certification evaluated five representative points only through pair and angular order 32. Every point remained uncertified: the pair refinement changed $\chi^2$ by $0.220$--$0.228$, and the angular refinement changed it by $0.0139$--$0.0155$, both above the unchanged $0.01$ tolerance. Higher orders are needed before `alp-uzh.2` can meet its acceptance criteria.

## What Changes

- Run a new high-order qualification with pair and angular orders 64→128, retaining independent refinements and the combined sequence. Keep the failed 16→32 results in the existing certification file as historical diagnostics, separate from the new pass/fail comparisons.
- Keep the existing uncertainty, residual, Jacobian, and $\chi^2$ tolerances fixed. Use bounded-memory energy-bin blocks for the higher-order grids without changing the physical domain, pole mappings, or likelihood.
- Store and replay the resulting comparison for every saved point, including explicit `uncertified` outcomes and limiting differences if any tolerance still fails.
- Update `alp-uzh.2` only after the saved points pass every declared comparison.

## Capabilities

### New Capabilities

None. The existing `kkpi-likelihood-quadrature` capability already defines the required comparisons and pass/failure behavior.

### Modified Capabilities

None. This follow-up changes the maximum orders used for evidence and may improve memory use; it does not change the certification requirements or tolerances. The spec delta is intentionally skipped.

## Impact

The work may extend `scripts/certify_kkpi_quadrature.py` and the grid evaluation in `scripts/fit_babar_kkpi_global_vector_sector.py`, add focused high-order and bounded-memory tests, and replace or supplement `data/vector_sector/kkpi_quadrature_certification.json`. The five selected points and their provenance remain those recorded by the existing certification. No physics model, data input, optimizer, or likelihood definition changes.
