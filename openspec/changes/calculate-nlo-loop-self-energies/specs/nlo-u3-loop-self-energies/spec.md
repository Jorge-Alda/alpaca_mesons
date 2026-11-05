# Spec Delta

## Purpose

Defines the renormalized one-loop meson and ALP–meson two-point functions of the $U(3)$ NLO theory with explicit resonances, with their divergence cancellation, threshold structure and pole-residue contributions, so that pole and LSZ calculations use one set of loop inputs.

## ADDED Requirements

### Requirement: Loop content restricted to declared classes
The two-point functions SHALL contain exactly the loop classes G-LO and G-ALP of alp-gii.1, at first order in $1/f_a$ for the mixed function. Loops with an internal resonance line SHALL NOT be included, and the records SHALL point to the truncation record for their size.

#### Scenario: Resonance loop requested
- **WHEN** a downstream issue asks for a loop with an internal resonance line in a two-point function
- **THEN** the records state that it is outside the calculated order and refer to the truncation uncertainty of alp-gii.1

#### Scenario: Classification of surviving topologies
- **GIVEN** the LO vertices of the retained sources
- **WHEN** the cubic and quartic vertices are enumerated
- **THEN** the records list which of tadpole and bubble topologies exist for $\Sigma_{ij}$ and $\Sigma_{i\hat a}$, and the enumeration is reproduced by the script

### Requirement: Explicit meson self-energy matrix
The records SHALL give $\Sigma_{ij}(s)$ for the full $U(3)$ nonet as an explicit function of $s$, the LO masses and mixing, and $\mu_\chi$, using the notation of Alda:2026nrd and the sign convention $Q(s)=sK-M^2+\Sigma(s)$.

#### Scenario: Neutral sector with mixing
- **WHEN** the $\pi^0$, $\eta$, $\eta'$ block is requested
- **THEN** $\Sigma_{ij}(s)$ is given in the $(3,8,0)$ basis, including the $\eta'$ loop functions with the LO mass of the neutral-meson record

#### Scenario: Degenerate-mass limit
- **WHEN** all LO meson masses are equal and the singlet is decoupled
- **THEN** the result reduces to the known $SU(3)$ chiral-logarithm form of Gasser:1984gg

### Requirement: Explicit mixed ALP–meson two-point function
The records SHALL give $\Sigma_{i\hat a}(s)$ at first order in $1/f_a$ for the ALP-source vertices of Alda:2026nrd, with flavour matrices kept non-commuting.

#### Scenario: Axial-coupling contribution
- **WHEN** only $c_q^A$ is non-zero
- **THEN** $\Sigma_{i\hat a}(s)$ is expressed through the same loop functions as $\Sigma_{ij}(s)$, and the records state whether it vanishes or is proportional to the kinetic-alignment combination of the existing records

#### Scenario: CP-odd source
- **WHEN** a CP-odd source or $\theta\neq0$ is switched on
- **THEN** the cubic vertices it generates are listed and the corresponding bubbles are either given or marked as an excluded conditional extension

### Requirement: Cancellation of divergences and scale dependence
The poles and the $\mu_\chi$ dependence of each loop function SHALL cancel against the local NLO terms built from $L_i^r(\mu)$ and $B^r(\mu)$ with the $\Gamma_i$ of alp-gii.1. The cancellation SHALL be checked symbolically, entry by entry.

#### Scenario: Scale independence of the renormalized inverse propagator
- **WHEN** $\mu_\chi$ is changed from $0.77\,$GeV to $0.5$ and $1.0\,$GeV with the $L_i^r$ shifted by their one-loop running
- **THEN** every entry of the renormalized $Q(s)$ is unchanged at one-loop order

#### Scenario: Coefficient mismatch
- **WHEN** a pole coefficient of the loop disagrees with the $\Gamma_i$ table of alp-gii.1
- **THEN** the script fails and the discrepancy is reported in the records

### Requirement: Threshold and absorptive structure
The records SHALL state, for each loop function, its non-analytic structure: the thresholds that exist, the branch cuts, and the imaginary part for $s$ above each threshold. If no bubble topology survives, the records SHALL state that the functions are real and analytic for all $s$ below the first resonance-induced threshold.

#### Scenario: Function with a cut
- **WHEN** a two-meson cut exists
- **THEN** the threshold is given as $(m_1+m_2)^2$ in GeV$^2$ with the sign of $\operatorname{Im}\Sigma$ fixed by unitarity

#### Scenario: Tadpole-only function
- **WHEN** only tadpoles contribute
- **THEN** the record states that the dependence on $s$ is polynomial and that no unitarity cut is generated at this order

### Requirement: Derivatives and pole-residue contributions
The records SHALL give $\partial_s\Sigma_{ij}$ and $\partial_s\Sigma_{i\hat a}$ and their effect on the LSZ meson components $w^{(a)}$, the physical-state residues and the perturbative NLO pole shifts, in the notation of the LSZ records.

#### Scenario: Strict NLO expansion
- **WHEN** the pole shift is expanded to first order in the one-loop terms
- **THEN** it equals the formula of the LSZ loop-correction record with $\delta Q$ replaced by the explicit renormalized $\Sigma$ and its derivative, evaluated at the LO pole

#### Scenario: Cross-check against the local limit
- **WHEN** the loop terms are set to zero
- **THEN** the result equals the tree-level LSZ solution of the existing records

### Requirement: Perturbative versus resummed treatment
The records SHALL separate the strict NLO result from any resummation. Near a threshold or a resonance pole where $\Sigma/(\text{LO inverse propagator})$ is not small, the prescription SHALL name the resummed quantity, label it as outside the strict order, and require downstream results to report both.

#### Scenario: Pole near a threshold
- **WHEN** the pole lies within the range where the perturbative correction exceeds the truncation estimate of alp-gii.1
- **THEN** the records give the criterion that triggers the resummed treatment and the form of the resummed propagator

### Requirement: Reproducible numerical evaluation
An executable script SHALL tabulate the renormalized functions and their derivatives at stated GeV inputs and SHALL reproduce the analytic expressions of the records, including the pole cancellation and the scale-independence checks.

#### Scenario: Script run with default inputs
- **WHEN** the script is run with its documented command
- **THEN** it exits with code 0, prints the loop functions at the test points, and asserts the cancellation and limit checks
