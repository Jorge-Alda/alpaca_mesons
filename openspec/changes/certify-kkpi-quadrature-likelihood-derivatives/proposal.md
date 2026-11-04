# Proposal

## Why

The joint KKπ fit currently records a single integration order, and its optional refinement audit compares only a limited set of BaBar predictions and χ² values. That does not establish that the energy-bin and two-dimensional Dalitz integrations, especially around the narrow ϕ and K*(892) structures, are accurate enough for residuals, likelihood derivatives, or fit interpretation.

## What Changes

- Add a reproducible quadrature certification workflow over representative parameter points spanning the fitted region.
- Compare independently refined energy, pair-invariant, and angular integrations for predictions, uncertainty-normalized residuals, derivatives, and χ²; report the orders, parameter points, tolerances, and invocation.
- Keep integration certification separate from optimizer convergence and do not label a fit accepted solely because it terminates successfully.

## Capabilities

### New Capabilities
- `kkpi-likelihood-quadrature`: Reproducible certification of energy-bin and Dalitz integration accuracy for the KKπ likelihood and its derivatives.

### Modified Capabilities
- None.

## Impact

The certification will extend the joint fit tooling and its machine-readable results under `scripts/fit_babar_kkpi_global_vector_sector.py` and `data/vector_sector/`. It will use the existing pole-aware mappings and likelihood implementation, preserve the current data and physics model, and add no runtime dependency or change to the likelihood definition.
