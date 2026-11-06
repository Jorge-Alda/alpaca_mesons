# Tasks

## 1. Audit inputs and declare coverage

- [x] 1.1 Create a channel input manifest covering the two-pion, three-pion, charged and neutral kaon, narrow omega/phi, and mixed rho-isoscalar terms; for each, record source basis, thresholds, conventions, phases, normalization, uncertainty provenance, and input status, then verify every entry links to an existing indexed derivation or names a missing input.
- [x] 1.2 Audit available complex production/current amplitudes, datasets, and correlated variations for the three-pion and kaon channels; verify that every unavailable item is recorded as a concrete blocker with the needed convention and uncertainty information.
- [x] 1.3 Declare the low-energy virtuality domain and survey higher-isoscalar strength; verify either a provenance-backed bound is included or the omission and its effect on coverage are stated without an unsupported uncertainty band.

## 2. Resolve rho analytic content before numerical inputs

- [x] 2.1 Map rho occurrences in $F_\pi^V$, the dipion-channel $a\gamma^*\to\pi\pi$ partial wave, its crossed $t/u$ terms, the photon response $W_\gamma(y)$, and the mixed $F^{(\rho I)}$ term; record which cuts and phases each contains, then verify the inventory against the indexed equations.
- [x] 2.2 Derive a rho-specific energy-dependent $\pi\pi$ width and self-energy convention, including the P-wave threshold behavior, $\operatorname{Im}\Delta_\rho^{\mathrm{self}}=+\sqrt{s}\,\Gamma_\rho(s)$, real-part subtractions, and mass/residue conditions; verify threshold behavior, positivity, and the stable-pole limit with an executable derivation and a paired indexed memory record.
- [x] 2.3 Derive the physical two-pion production amplitude’s isoscalar pole residue, phase, and subtraction terms, and the required continuation and compatible complex production/decay amplitudes.
- [x] 2.4 Compare the single and double discontinuities and subtraction terms from the two-pion cut with the mixed $F^{(\rho I)}$ contribution, including both photon orderings; derive a matching/subtraction prescription or establish disjointness, then verify the result with a symbolic check and a linked memory record.
- [x] 2.5 Decide whether the mixed $F^{(\rho I)}$ term can be retained consistently in the first approximation or must be deferred; base the decision on tasks 2.1–2.3, update the channel manifest with the component list and rationale, and verify that no numerical first-pass sum includes unresolved shared rho strength.

## 3. Complete inputs for the selected coverage

- [ ] 3.1 Obtain or derive compatible source-resolved two-pion current and production inputs, including complex phases, cutoff convention, and correlated variations; verify conventions against the rho treatment from task 2 and record any unavailable inputs as blockers.
- [ ] 3.2 Obtain or derive compatible complex three-pion current and ALP production amplitudes, including the physical relative phase and correlated variations; verify the common current/production conventions or document concrete blockers.
- [ ] 3.3 Obtain or derive complex charged and neutral kaon current form factors and source-resolved production amplitudes with phases and correlated variations; verify separate thresholds and charge conventions or record concrete blockers.
- [ ] 3.4 Establish the $\omega/\phi$ pole replacement or matching across the three-pion and kaon cuts, especially in the $\phi$ region; verify shared strength is represented once and unsupported finite-width input remains conditional.

## 4. Build and document source-resolved assembly

- [ ] 4.1 Implement deterministic assembly of the selected baseline and any enabled extensions using the no-mixing sources and input manifest; verify source linearity, preserve channel labels and completeness status, and prevent blocked inputs from producing a complete-TFF result.
- [ ] 4.2 Add executable checks for cut normalization, thresholds, Bose symmetry, spacelike and timelike boundary prescriptions, rho line-shape limits, and pole/continuum replacement; verify each against the symbolic derivations.
- [ ] 4.3 Document the assembled components, units, conventions, invocation, and supported virtuality domain in an indexed memory record; verify it links the executable derivations and inputs and is linked from the lepton-pair index.

## 5. Propagate uncertainties and validate the deliverable

- [ ] 5.1 Propagate available correlated input variations and separately report cutoff, subtraction, line-shape, matching, and omitted-channel model variations; verify provenance is preserved and unsupported covariances are not assumed.
- [ ] 5.2 Keep resonance calibration and fit-driven uncertainty propagation disabled until `alp-uzh.8` has a converged, identifiable fit; verify no coefficient from the current nonconverged likelihood enters the accepted input manifest.
- [ ] 5.3 Run the executable channel and integration checks and validate the declared scope; verify status output distinguishes the selected low-energy approximation from a complete TFF over unbounded lepton-loop virtualities.
