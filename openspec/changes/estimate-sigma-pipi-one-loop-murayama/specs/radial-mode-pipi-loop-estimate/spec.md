# Spec Delta

## Purpose

Defines the analytical size estimate of the one-loop correction to the radial-mode decay $\sigma\to\pi^+\pi^-$ in the generalized Murayama model, in terms of generic source coefficients $\widehat S_{r,s}$ and $\widehat P_{r,s}$, with the loop integrals replaced by their size and the source dependence kept exact.

## ADDED Requirements

### Requirement: Exact dependence on the invariants
The estimate SHALL be written in the generic coefficients $\widehat S_{r,s}$ and $\widehat P_{r,s}$. Their flavour traces, together with the external pion matrices, the pole vectors, $N(s)$ and the powers of $\epsilon$, SHALL be exact and symbolic.

#### Scenario: Generic sources
- **WHEN** the estimate is inspected
- **THEN** it contains the coefficients $\widehat S_{r,s}$ and $\widehat P_{r,s}$ for arbitrary $r,s$, with no up-only substitution and no expansion in $\epsilon$

#### Scenario: Internal meson line
- **WHEN** a graph has an internal meson line
- **THEN** the generator sum is carried out with the completeness relation, and the result is a combination of compact traces of the sources and external matrices (no component notation)

### Requirement: No small-coupling truncation
No step SHALL assume $\epsilon\ll1$. The selection of graphs SHALL be justified for arbitrary $\epsilon$, with large $\epsilon$ as the case of interest.

#### Scenario: Graph selection
- **WHEN** a graph is kept or dropped
- **THEN** the reason is its behaviour at large $\epsilon$ or an exact vanishing, and not a small-$\epsilon$ expansion

### Requirement: Compact trace notation
Traces SHALL be written with matrix products, commutators, anticommutators and sums over permutations or cyclic permutations only. Matrix indices, generator components and the basis $E_{ij}$ SHALL NOT be used.

#### Scenario: Inspection of records
- **WHEN** any record or script output of this change is inspected
- **THEN** it contains no matrix indices, no $E_{ij}$ and no explicit generator components

### Requirement: Analytical verification
All checks SHALL be exact symbolic identities.

#### Scenario: Test inputs
- **WHEN** a test is inspected
- **THEN** it uses no random or numerical matrices

### Requirement: Loop integrals replaced by the scalar-moment size
Each loop SHALL contribute its scalar moments with $r$ denominators replaced by $(m^2)^{2-r}/(16\pi^2)$, with $m$ the loop scale and the polynomial coefficient of the moment kept. The estimate SHALL NOT evaluate any loop integral, scalar loop function, loop coefficient, renormalization-scale derivative or loop amplitude.

#### Scenario: Bubble graph
- **WHEN** a graph has two internal lines
- **THEN** its factor is $(m^2)^{0}/(16\pi^2)$ times its source traces and symmetry factor, with no function of the external momenta

#### Scenario: Loop value required
- **WHEN** a step would need the value of a loop integral
- **THEN** the step is replaced by the size prescription and no value is computed

### Requirement: Loop scale per term and routing bound
Each term of a graph SHALL use as loop scale the largest pole scale among its lines: $m_\sigma$ for a radial line, $m_a$ for an ALP-like line, $m_P$ for the bare generator part of a meson line and $m_a$ for the rest of the meson line. For each independent product of traces the routing of the loop momentum with the largest magnitude of the coefficient SHALL be used, and the record SHALL state that the result is an upper bound per graph on the routing ambiguity, not on the loop.

#### Scenario: Meson-only graph
- **WHEN** a graph has only meson lines
- **THEN** its terms with every line in the bare generator part use $m_P$, and no scale is tied to $m_\sigma$

#### Scenario: Routing-dependent graph
- **WHEN** a coefficient vanishes in some routing
- **THEN** the largest-magnitude routing is used and the spread is reported as a diagnostic

### Requirement: Complete graph list for one loop
The record SHALL list every one-loop graph with the external legs $\sigma$, $\pi^+$, $\pi^-$ built from the vertices of the model: tadpoles, bubbles and triangles, containing at least one ALP-like or radial internal line, plus the corrections of the external legs. For each graph it SHALL state the vertices, symmetry factor and number of internal lines.

#### Scenario: Excluded topology
- **WHEN** a topology cannot be built with three external legs, such as a box with four cubic vertices
- **THEN** the record states why it is absent

#### Scenario: Meson-only loops
- **WHEN** a graph has only meson internal lines
- **THEN** only its term with the bare $U(3)$ generator sum on every line and a single $\sigma$ end is excluded, because it contains no letter of the ALP sector and its ratio to the tree depends on $\epsilon$ only through $p_\sigma^2=m_\sigma^2$; the record states it as a chiral background that is not computed, and the remaining terms (ALP-field entry and ALP-like pole parts of the meson lines) are kept

### Requirement: Comparison with the exact tree
The record SHALL define $R=|M_1^{\rm est}|/|M_0|$, with $M_0$ the exact tree $\sigma\to\pi^+\pi^-$ amplitude of the generalized model, and SHALL state that $R$ is a size ratio and not a computed correction.

#### Scenario: Vanishing tree
- **WHEN** the tree amplitude vanishes for some sources
- **THEN** the record states that $R$ is undefined there and reports only $M_1^{\rm est}$

### Requirement: Frame invariance and limits
The estimate SHALL be unchanged under a constant axial redefinition of the sources, checked on the traces. The record SHALL state the limits: absorptive parts and thresholds, cancellations between graphs, and the naturalness assumption for the neglected terms with no analysis of counterterms.

#### Scenario: Axial redefinition
- **WHEN** the sources are redefined by a constant axial transformation
- **THEN** the estimate is unchanged for the same physical legs

#### Scenario: Open decay channel in the loop
- **WHEN** an internal pair is kinematically open at $s=m_\sigma^2$
- **THEN** the record states that the estimate gives a magnitude only and carries no absorptive part
