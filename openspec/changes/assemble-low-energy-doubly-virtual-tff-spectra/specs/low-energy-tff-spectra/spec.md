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

### Requirement: Declare rho line shape and analytic placement
Every rho-containing term MUST identify whether its rho dynamics enter through the pion vector form factor, a dressed propagator, or a stable-pole approximation, and state the relevant invariant and boundary prescription.

#### Scenario: Energy-dependent rho response
- **WHEN** a rho propagator is dressed with an energy-dependent self-energy
- **THEN** the width, open-channel thresholds, real-part subtraction conditions, and mass/residue convention are declared, and existing pion rescattering is not counted again

### Requirement: Count shared resonance strength once
The representation MUST distinguish production amplitudes used inside unitarity cuts from separate TFF pole or continuum terms, and MUST include shared resonance strength only once.

#### Scenario: Two-pion production and mixed rho-isoscalar term
- **WHEN** the $a\gamma^*\to\pi\pi$ production wave feeds a two-pion discontinuity and a mixed $F^{(\rho I)}$ term is considered
- **THEN** their overlap is checked at the single and double discontinuities and subtraction terms before addition, and an unresolved overlap blocks numerical inclusion of the mixed term

#### Scenario: Stable rho pole and two-pion spectrum
- **WHEN** the two-pion dispersive component already represents rho strength
- **THEN** a separate stable rho pole for that same strength is excluded or removed by a documented matching prescription

#### Scenario: Omega or phi overlaps a multiparticle cut
- **WHEN** a narrow omega or phi pole is combined with a continuum containing that resonance
- **THEN** one representation replaces or is matched to the shared pole strength, including common three-pion and kaon strength in the phi region

### Requirement: Declare approximation content
Every low-energy TFF result MUST list its included, conditional, and omitted components so a staged approximation is not mistaken for a complete spectrum.

#### Scenario: First approximation is evaluated
- **WHEN** the two-pion spectrum and narrow omega/phi poles are assembled
- **THEN** the result states whether the mixed rho-isoscalar term is included and identifies three-pion, kaon, higher-state, and short-distance contributions that remain outside the approximation

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
