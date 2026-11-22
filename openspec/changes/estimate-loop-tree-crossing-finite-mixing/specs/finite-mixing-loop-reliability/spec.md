# Spec Delta

## Purpose

Defines the analytical, redefinition-invariant estimate of the value of $f_\pi/f_a$ at which one-loop corrections to an ALP-sourced observable become comparable to its exact tree-level value in the finite-mixing LO theory, including the justification of the chosen observable. Conditions that involve tree-level quantities only are not part of the requirements.

## ADDED Requirements

### Requirement: Justified choice of observable
The record SHALL compare at least the meson pole mass, the ALP-meson mixing entry, the ALP pole residue and an on-shell ALP-sourced amplitude, and SHALL select one. The comparison SHALL state, for each candidate, whether it is invariant under constant chiral redefinitions, how the ALP enters it at tree level (the power of $1/f_a$), whether it depends on a fitted number or on the subtraction scheme, and what that implies for the $\epsilon$ dependence of the loop-to-tree ratio.

#### Scenario: Candidate rejected for lack of invariance
- **WHEN** a candidate depends on an off-shell entry that changes under a constant axial redefinition
- **THEN** the record states the change explicitly and does not use the candidate for the estimate

#### Scenario: Candidate with epsilon-independent ratio
- **WHEN** the tree and the loop term of a candidate scale with the same power of $\epsilon$ at leading order
- **THEN** the record states that this candidate cannot locate a value of $f_\pi/f_a$ and says what it can be used for

### Requirement: Exact tree quantity
The tree-level quantity of the selected observable SHALL be taken from the exact finite-mixing pole solution, written with $\mathbb M_{m_a^2}^{-1}$ and the invariants $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$ and $\mathcal C$, and SHALL NOT be replaced by its expansion in $\epsilon$. The record SHALL state where an expansion is used and which order is kept.

#### Scenario: Exact against expansion
- **WHEN** the tree quantity is evaluated with exact rational sources at a finite $\epsilon$
- **THEN** it equals the closed formula, and the difference from its leading expansion is stated

### Requirement: Loop quantity split by loop function and ALP lines
The one-loop correction to the selected observable SHALL be written as a sum of coefficients of independent loop functions, each coefficient tagged by the number of ALP lines in the loop and by its power of $\epsilon$. Notation for fields and sources SHALL follow the existing records (the multiplet $\Phi_I$, matrices without flavour indices).

#### Scenario: Power of epsilon of each class
- **WHEN** a loop class contains $n_a$ ALP lines
- **THEN** the record states its power of $\epsilon$ relative to the tree quantity and the loop factor that accompanies it

#### Scenario: Meson-only loops
- **WHEN** a loop contains no ALP line
- **THEN** the record states that its relative size does not depend on $\epsilon$ and treats it as the usual chiral expansion, not as a crossing

### Requirement: Absorbable and non-absorbable contributions
The record SHALL state which parts of the one-loop correction are absorbed by local NLO operators or by the renormalization of fitted parameters, and SHALL base the estimate only on the part that is not absorbed. It SHALL state the dependence of the estimate on the subtraction scale.

#### Scenario: Scale dependence
- **WHEN** the subtraction scale $\mu$ is changed
- **THEN** the record states which coefficient changes, by how much, and that the estimate is reported with this dependence explicit

### Requirement: Frame invariance of the estimate
The loop-to-tree ratio and the resulting estimate SHALL be shown, with an executable symbolic or exact check, to be unchanged under a constant axial redefinition of the chiral frame, and SHALL be expressed only through frame-invariant quantities. No value, representative or relation among the ALP coefficients SHALL be assumed.

#### Scenario: Constant axial redefinition
- **WHEN** the sources are transformed by a constant axial redefinition
- **THEN** the loop-to-tree ratio and the crossing formula take the same value in both frames

### Requirement: Closed formula for the crossing
The record SHALL give a closed formula for the $f_\pi/f_a$ at which the loop-to-tree ratio of the selected observable reaches the stated order-one criterion, in terms of the invariants, the unmixed masses, $m_a$, $f_\pi$, $B_0$, $m_0^2$ and the loop functions. It SHALL state the criterion value used and how the formula changes if another value is chosen. No scan over parameters SHALL be used to obtain or present it.

#### Scenario: Formula against exact ratio
- **WHEN** the exact loop-to-tree ratio is evaluated with exact rational sources at the crossing value
- **THEN** it equals the criterion value within the error bound stated for the closed formula

### Requirement: Validity and limits
The record SHALL state the conditions of validity (non-degenerate poles, one-loop truncation, absence of NLO counterterms in the ALP-line classes), the behaviour when the ALP pole approaches a meson pole, and the limit $f_a\to\infty$.

#### Scenario: Near degeneracy
- **WHEN** the ALP pole is within the mixing entry of a meson pole
- **THEN** the record states that the closed formula fails and which exact expression replaces it

### Requirement: Exclusion of tree-only conditions
The estimate SHALL NOT be defined by a condition on tree-level quantities only, such as a large mixing angle or partial-wave unitarity.

#### Scenario: Review of the estimate
- **WHEN** the definition of the crossing is inspected
- **THEN** it is a comparison of one-loop and tree-level terms of the same observable
