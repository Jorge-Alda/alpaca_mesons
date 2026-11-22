# Spec Delta

## Purpose

Defines the analytical, redefinition-invariant power-counting estimate of the value of $f_\pi/f_a$ at which the one-loop size of an ALP-sourced amplitude becomes comparable to its exact tree-level value in the finite-mixing LO theory, including the justification of the observable class. Conditions that involve tree-level quantities only, and the evaluation of loop integrals, are not part of the requirements.

## ADDED Requirements

### Requirement: Justified choice of observable class
The record SHALL compare at least the meson pole mass, the ALP-meson mixing entry, the ALP pole residue and an on-shell ALP-sourced amplitude, and SHALL state the class used for the counting. For each candidate it SHALL state its frame invariance, the power of $1/f_a$ at tree level, and the implication for the $\epsilon$ dependence of the ratio.

#### Scenario: Candidate rejected for lack of invariance
- **WHEN** a candidate depends on an off-shell entry that changes under a constant axial redefinition
- **THEN** the record states the change explicitly and does not use the candidate for the estimate

#### Scenario: Candidate with epsilon-independent ratio
- **WHEN** the tree and the loop term of a candidate scale with the same power of $\epsilon$ at leading order
- **THEN** the record states that this candidate cannot locate a value of $f_\pi/f_a$ and says what it can be used for

### Requirement: Exact tree quantity
The tree-level quantity SHALL be taken from the exact finite-mixing pole solution, written with $\mathbb M_{m_a^2}^{-1}$ and the invariants $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$ and $\mathcal C$, in matrix notation, and SHALL NOT be replaced by its expansion in $\epsilon$. The record SHALL state where an expansion is used and which order is kept.

#### Scenario: Exact against expansion
- **WHEN** the tree quantity is evaluated with exact rational sources at a finite $\epsilon$
- **THEN** it equals the closed formula, and the difference from its leading expansion is stated

### Requirement: Loop-to-tree size from power counting
The record SHALL give, for each loop class of the dressed-frame vertices, the power of $\epsilon$ per closed ALP line, the loop factor, and the dimensionless contraction of the invariants that multiplies them, in matrix notation without flavour indices. It SHALL NOT evaluate loop integrals or an amplitude.

#### Scenario: Power of epsilon of each class
- **WHEN** a loop class contains $n_a$ ALP lines
- **THEN** the record states its power of $\epsilon$ relative to the tree quantity, the loop factor and the invariant contraction

#### Scenario: Meson-only loops
- **WHEN** a loop contains no ALP line
- **THEN** the record states that its relative size does not depend on $\epsilon$ and treats it as the usual chiral expansion, not as a crossing

### Requirement: Naturalness assumption
The record SHALL state that the neglected next-order term is taken to have the size of the loop, and SHALL list the accidental cancellations (for example among classes at the same power of $\epsilon$, or by relations among the invariants) for which the estimate overestimates it.

#### Scenario: Cancellation among classes
- **WHEN** two loop classes have the same power of $\epsilon$ and their invariant contractions can cancel
- **THEN** the record states the cancellation condition and that the estimate does not apply there

### Requirement: Frame invariance of the estimate
The loop-to-tree ratio and the resulting estimate SHALL be shown, with an executable symbolic or exact check, to be unchanged under a constant axial redefinition of the chiral frame, and SHALL be expressed only through frame-invariant quantities. No value, representative or relation among the ALP coefficients, and no relation among the quark masses (such as $m_u=m_d$) or among the invariants (such as $\langle\lambda_3\widehat P_1\rangle=0$), SHALL be assumed.

#### Scenario: Constant axial redefinition
- **WHEN** the sources are transformed by a constant axial redefinition
- **THEN** the loop-to-tree ratio and the crossing formula take the same value in both frames

### Requirement: Closed formula for the crossing
The record SHALL give a closed formula for the $f_\pi/f_a$ at which the ratio of the leading ALP-line class to the tree quantity reaches the order-one criterion $R_{\max}$, in terms of the invariants, the unmixed masses, $m_a$, $f_\pi$, $B_0$, $m_0^2$ and the loop factor, valid for a general $R_{\max}$. No scan over parameters SHALL be used to obtain or present it.

#### Scenario: Formula against vertex counting
- **WHEN** the vertex scalings are evaluated with exact rational sources
- **THEN** the closed formula reproduces the stated powers of $\epsilon$ and the invariant contraction

### Requirement: Validity and limits
The record SHALL state the conditions of validity (non-degenerate poles, one-loop truncation), the behaviour when the ALP pole approaches a meson pole, and the limit $f_a\to\infty$.

#### Scenario: Near degeneracy
- **WHEN** the ALP pole is within the mixing entry of a meson pole
- **THEN** the record states that the closed formula fails and which exact expression replaces it

### Requirement: Exclusion of tree-only conditions
The estimate SHALL NOT be defined by a condition on tree-level quantities only, such as a large mixing angle or partial-wave unitarity.

#### Scenario: Review of the estimate
- **WHEN** the definition of the crossing is inspected
- **THEN** it is a comparison of the one-loop size and the tree-level value of the same amplitude
