# Spec Delta

## Purpose

This capability defines the ALP-specific short-distance behavior of the doubly virtual TFF and its matching to the low-energy representation. It makes sources, scheme, thresholds, validity regimes, and uncertainties explicit so that the lepton-loop integral over unbounded virtualities is controlled.

## ADDED Requirements

### Requirement: Derive asymptotics from ALP matrix elements
The high-energy TFF MUST be expressed through renormalized ALP matrix elements for each independent no-mixing source, and MUST NOT use an asymptotic coefficient taken from $\eta$ or $\eta'$.

#### Scenario: Source-dependent coefficient
- **WHEN** the asymptotic coefficient is evaluated for supplied ALP sources
- **THEN** it is given as a linear function of the source coefficients, with each matrix element, its operator normalization, and its provenance named

### Requirement: Declare scale, scheme, and thresholds
Every short-distance result MUST state the renormalization scale and scheme, the number of active flavors, the treatment of the singlet axial current and anomaly, and the heavy-quark threshold matching.

#### Scenario: Scale crosses a flavor threshold
- **WHEN** the virtuality or matching scale crosses a heavy-quark threshold
- **THEN** the number of active flavors and the matching of the matrix elements across the threshold are stated, and the result is not silently evaluated with three massless flavors

### Requirement: State the virtuality regime of each result
Each asymptotic expression MUST identify whether it holds for fixed-ratio doubly virtual, singly virtual, or nearly symmetric kinematics, and unsupported regimes MUST be marked as such.

#### Scenario: Regime outside the derivation
- **WHEN** a TFF evaluation requests virtualities outside every derived regime
- **THEN** the result is flagged as outside the supported short-distance domain rather than extrapolated

### Requirement: Retain finite-$m_a$ corrections
The derivation MUST give the leading corrections in $m_a^2/Q^2$ or state that they are beyond the declared accuracy, and MUST remain valid for $m_a>0$ without a massless-ALP singularity.

#### Scenario: Matching scale not far above $m_a$
- **WHEN** the matching scale is not parametrically above $m_a$
- **THEN** the finite-$m_a$ terms are included or their omission is quantified

### Requirement: Match to the low-energy form without double counting
A matching prescription MUST join the asymptotic form to the low-energy representation at a declared scale, preserve the low-energy normalization and assigned low-energy discontinuities, identify added high-energy discontinuities, and count spectral strength once.

#### Scenario: Physical overlap is unknown
- **WHEN** the overlap between low- and high-energy strength is not established
- **THEN** state the required subtraction or conditional separation without claiming the overlap is resolved

#### Scenario: Euclidean diagnostic is used
- **WHEN** a sharp join or interpolation has no established analytic continuation
- **THEN** do not present it as an analytic form factor or a test of physical discontinuities

#### Scenario: Matching-scale variation
- **WHEN** the matching scale is varied over the declared range
- **THEN** the resulting TFF variation is reported separately from hadronic-input and scheme variations, and alternative completions are compared

### Requirement: Keep fit-dependent matching conditional
Matching to a numerically calibrated low-energy form MUST NOT use coefficients from the nonconverged vector fit.

#### Scenario: Calibrated low-energy input unavailable
- **WHEN** no converged low-energy calibration exists
- **THEN** the matching is tested on declared symbolic or benchmark low-energy forms and labelled conditional
