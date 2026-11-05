# Spec Delta

## Purpose

This capability defines an auditable, source-resolved low-energy representation of the doubly virtual ALP transition form factor. It makes channel inputs, overlap treatment, uncertainties, and the domain of validity explicit before numerical use.

## ADDED Requirements

### Requirement: Preserve source-resolved channel contributions
The TFF representation MUST retain each independent no-mixing ALP source and each hadronic channel as separately identifiable contributions until their amplitudes are combined with the declared source coefficients.

#### Scenario: Source-resolved evaluation
- **WHEN** a TFF contribution is evaluated for supplied ALP sources
- **THEN** the result identifies the source coefficients and the channel response used for each term

### Requirement: Declare channel conventions and normalization
Each spectral contribution MUST identify its current and production amplitude conventions, complex phase convention, cut measure and normalization, thresholds, and numerical input provenance.

#### Scenario: Channel input is incomplete
- **WHEN** a required complex amplitude, phase match, cut normalization, or uncertainty input is unavailable
- **THEN** the channel is marked incomplete with the missing input named, and no numerical prediction is presented as complete

### Requirement: Count shared resonance strength once
The representation MUST specify how resonance strength shared between explicit poles and hadronic continua is matched or replaced, so the same strength is not included twice.

#### Scenario: Dressed rho overlaps two-pion strength
- **WHEN** the dressed mixed rho–isoscalar contribution and the two-pion continuum are assembled
- **THEN** their shared rho strength is matched or allocated once, with the prescription and any unresolved matching input stated

#### Scenario: Omega or phi overlaps a multiparticle cut
- **WHEN** a narrow omega or phi pole is combined with a continuum containing that resonance
- **THEN** one representation replaces or is matched to the shared pole strength, including common three-pion and kaon strength in the phi region

### Requirement: Record uncertainty provenance and omitted strength
Every numerical channel result MUST distinguish input uncertainty from model variation and identify omitted higher isoscalar strength, with a bound when one is supported.

#### Scenario: Higher isoscalar channels are not constrained
- **WHEN** no supported bound or numerical input exists for omitted higher isoscalar strength
- **THEN** the result states that omission as an unresolved limitation and does not assign an unsupported uncertainty band

### Requirement: Declare the supported virtuality domain
Every assembled low-energy TFF result MUST state its virtuality domain and identify the continuation or matching needed outside it; a low-energy representation MUST NOT be presented as a complete form factor over unbounded loop virtualities without that input.

#### Scenario: Low-energy result used in the full lepton loop
- **WHEN** the requested loop probes virtualities beyond the declared low-energy domain
- **THEN** the missing higher-energy continuation and its matching requirements are reported explicitly

### Requirement: Keep fit-dependent calibration conditional
Numerical resonance calibration and fit-driven uncertainty propagation MUST use only a validated vector-sector fit, and MUST remain unassigned while the required fit is nonconverged or nonidentifiable.

#### Scenario: Vector-sector fit is unresolved
- **WHEN** the vector-sector fit has not passed convergence and identifiability requirements
- **THEN** its coefficients are not used as validated inputs and fit-driven uncertainties are reported as pending
