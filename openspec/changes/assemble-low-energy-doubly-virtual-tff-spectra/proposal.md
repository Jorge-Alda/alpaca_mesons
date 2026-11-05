# Proposal

## Why

The repository has derivations for the principal low-energy TFF channels, but no assembled source-resolved spectrum with compatible phases, uncertainty provenance, and a prescription for resonance strength shared by poles and continua. A coherent low-energy representation is needed to make the TFF usable while keeping the currently nonconverged vector fit out of validated numerical inputs.

## What Changes

- Define a common low-energy TFF representation for two-pion, three-pion, kaon, narrow $ω/φ$, and mixed $ρ$–isoscalar contributions.
- Require channel inputs, current and production conventions, complex phases, cut normalization, uncertainty provenance, and the supported virtuality domain to be recorded explicitly.
- Specify single-counting rules where dressed $ρ$ strength overlaps the two-pion continuum and where narrow $ω/φ$ poles overlap their multiparticle cuts.
- Record missing amplitudes or correlated variations as concrete blockers, and bound or explicitly omit higher isoscalar strength.
- Keep numerical resonance calibration and fit-driven uncertainty propagation conditional on completion of `alp-uzh.8`.

## Capabilities

### New Capabilities
- `low-energy-tff-spectra`: Source-resolved low-energy TFF contributions with consistent channel normalization, overlap treatment, uncertainty provenance, and declared validity domain.

### Modified Capabilities
- None.

## Impact

The change will consolidate and extend the existing TFF derivations and input records under `memory/alp/lepton-pair/`, with reproducible numerical assembly in `scripts/` and an updated index. It depends on source-specific pion, three-pion, and kaon production/current inputs; complex three-pion solutions and phase matching; and a common treatment of the $φ$ region. No validated vector-fit coefficients are available until `alp-uzh.8` is resolved.
