# Tasks

## 1. Audit inputs and declare coverage

- [x] 1.1 Create a channel input manifest covering the two-pion, three-pion, charged and neutral kaon, narrow omega/phi, and mixed rho-isoscalar terms; for each, record source basis, thresholds, conventions, phases, normalization, uncertainty provenance, and input status, then verify every entry links to an existing indexed derivation or names a missing input.
- [ ] 1.2 Audit available complex production/current amplitudes, datasets, and correlated variations for the three-pion and kaon channels; verify that every unavailable item is recorded as a concrete blocker with the needed convention and uncertainty information.
- [ ] 1.3 Declare the low-energy virtuality domain and survey higher-isoscalar strength; verify either a provenance-backed bound is included or the omission and its effect on coverage are stated without an unsupported uncertainty band.

## 2. Complete channel inputs and matching prescriptions

- [ ] 2.1 Obtain or derive compatible source-resolved two-pion, three-pion, and kaon production/current inputs, including complex phases and correlated variations; verify each input against its stated convention and preserve unavailable components as blockers.
- [ ] 2.2 Derive and document the rho-to-two-pion/mixed-rho-isoscalar matching prescription, including the dressed rho response and subtraction convention; verify the shared rho strength is represented once or the unresolved numerical match is explicitly gated.
- [ ] 2.3 Establish the omega/phi pole replacement or matching across the three-pion and kaon cuts, especially in the phi region; verify that shared resonance strength is counted once and that unsupported finite-width input is marked conditional.

## 3. Build source-resolved numerical assembly

- [ ] 3.1 Implement deterministic channel-response assembly using the declared no-mixing source coefficients and the manifest, preserving channel labels, complex phases, and completeness status; verify source linearity and prevent a blocked or incomplete channel from being reported as a complete TFF.
- [ ] 3.2 Add numerical checks for cut normalization, thresholds, Bose symmetry, spacelike and timelike boundary prescriptions, and pole/continuum replacement; verify them with executable tests against the existing symbolic derivations.
- [ ] 3.3 Document the assembly inputs, units, conventions, invocation, and supported virtuality domain in an indexed memory record; verify the record links the executable derivation, cites its inputs, and is linked from the lepton-pair index.

## 4. Propagate uncertainties and validate the deliverable

- [ ] 4.1 Propagate available correlated input variations and separately report cutoff, subtraction, matching, and omitted-channel model variations; verify the reported uncertainty components retain their provenance and no unsupported covariance is assumed.
- [ ] 4.2 Keep resonance calibration and fit-driven uncertainty propagation disabled until `alp-uzh.8` has a converged, identifiable fit; verify no coefficient from the current nonconverged likelihood enters the accepted input manifest.
- [ ] 4.3 Run the executable channel and integration checks and validate the declared scope; verify status output distinguishes the assembled low-energy contribution from a complete TFF over the unbounded lepton-loop virtualities.
