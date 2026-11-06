# Proposal

## Why

The joint vector-sector likelihood now has four accepted, quadrature-qualified local minima, the best being `basin_probes/even_shared_widths` with $\chi^2=225.727$ for 43 degrees of freedom. A $\chi^2$ per degree of freedom near 5.2 indicates that the retained model may not describe the data, so parameter estimates and the identifiability study (alp-uzh.7) cannot be interpreted until the goodness of fit is assessed and the result is shown to be robust, or the specific model failure is identified.

## What Changes

- Add a goodness-of-fit assessment at each accepted minimum: total and per-block $\chi^2$, nominal $p$-values with an explicit statement of the degrees-of-freedom assumption, standardized residual listings, the largest pulls, and the energy dependence of the BaBar $KK\pi$ residuals in each channel.
- Add block-exclusion refits (one data block removed at a time, warm-started from the best accepted minimum) to locate which blocks are in mutual tension.
- Add a predefined sensitivity study with one refit per justified variation, each reproducing the baseline likelihood exactly when the variation is switched off:
  - BaBar fit window (upper energy cut lowered from 1.78 GeV to 1.70 and 1.60 GeV);
  - excited-resonance inputs (excited-pole masses and widths varied within PDG ranges; regulator 2.5 GeV varied);
  - tensor contribution (tensor products fixed to zero; nearby tensor pole parameters varied);
  - source hypothesis (universal, single-trace positive, single-plus-double-trace, in addition to the baseline negative single-trace);
  - mass prescription (the second inverse-mass branch);
  - covariance treatment (with and without the $\phi\to K\bar K$ correlation, the alternative $\phi$ electronic-width summary, and an added fractional systematic on BaBar bins).
- Record a verdict: either the model is adequate (possibly within a stated variation), or a concrete model failure is identified with the responsible blocks and the physics extension it requires, filed as a follow-up bead.

## Capabilities

### New Capabilities

- `joint-vector-fit-adequacy`: Reproducible goodness-of-fit assessment and predefined sensitivity study of accepted joint vector-sector minima, ending in an adequacy verdict or a documented model failure.

### Modified Capabilities

None.

## Impact

- New option flags in `scripts/fit_babar_kkpi_global_vector_sector.py` for each variation, defaulting to the current likelihood; a new driver script for block exclusion and sensitivity refits, with focused tests. Script changes are delegated to the `coder` agent.
- New ledgers under `data/vector_sector/` and indexed memory records under `memory/rcht/odd-parity/`.
- Existing accepted endpoints, historical fit outputs, and the baseline likelihood are unchanged.
- May create a follow-up bead for a model extension; unblocks alp-uzh.7 only if the verdict permits parameter interpretation.
