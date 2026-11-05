# Proposal

## Why

The repository has derivations for the principal low-energy TFF channels, but the rho line shape and possible overlap between the two-pion cut and the mixed rho-isoscalar term must be settled before numerical inputs are combined. A staged source-resolved construction is needed to avoid counting the same resonance strength twice and to keep the currently nonconverged vector fit out of validated inputs.

## What Changes

- Establish an initial representation with the broad isovector $\rho$ in the source-specific two-pion spectrum and pure isoscalar $\omega/\phi$ terms in the narrow-pole approximation; treat three-pion and kaon continua as later extensions.
- Before numerical assembly, derive where the energy-dependent $\rho$ response belongs and audit possible double counting between the $a\gamma^*\to\pi\pi$ production wave and the mixed $F^{(\rho I)}$ term.
- Decide whether the mixed $\rho I$ VVP term can be included consistently in the first approximation or must be deferred until its overlap is resolved.
- Require channel inputs, current and production conventions, complex phases, cut normalization, uncertainty provenance, and the supported virtuality domain to be recorded explicitly.
- Specify single-counting rules where rho strength overlaps the two-pion continuum and where narrow $\omega/\phi$ poles overlap their multiparticle cuts.
- Record missing amplitudes or correlated variations as concrete blockers, and bound or explicitly omit higher isoscalar strength.
- Keep numerical resonance calibration and fit-driven uncertainty propagation conditional on completion of `alp-uzh.8`.

## Capabilities

### New Capabilities
- `low-energy-tff-spectra`: Source-resolved low-energy TFF contributions with consistent channel normalization, overlap treatment, uncertainty provenance, and declared validity domain.

### Modified Capabilities
- None.

## Impact

The change will first resolve analytic component and propagator conventions, then consolidate the existing derivations and inputs under `memory/alp/lepton-pair/`, with reproducible calculations in `scripts/` and an updated index. The initial numerical model uses the two-pion isovector spectrum and narrow $\omega/\phi$ poles; mixed rho-isoscalar inclusion is gated by the overlap decision. Three-pion and kaon continua require later compatible complex inputs and a common treatment of the $\phi$ region. No validated vector-fit coefficients are available until `alp-uzh.8` is resolved.
