# Spec Delta

## Purpose

Defines an independent analytical size estimate of the one-loop correction to $\sigma\to\pi^+\pi^-$ in the generalized Murayama model with the Cartesian fields $\varphi=f_a+s+ip$, written in the invariants $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$, under the same assumptions and approximations as the polar estimate.

## ADDED Requirements

### Requirement: Same assumptions, independent derivation
The estimate SHALL use the size rule, per-term loop scales, routing rule, exclusion of the meson-only chiral background, CP and flavour-diagonality assumptions, residue prescription and large-$\epsilon$ regime of `radial-mode-pipi-loop-estimate`, restated in the Cartesian invariants. The vertices, the graph list and the traces SHALL be derived from the Cartesian Lagrangian and SHALL NOT be obtained by translating polar records or through the map $(\sigma,a)\to(s,p)$.

#### Scenario: Assumption audit
- **WHEN** the records are inspected
- **THEN** each assumption is stated in terms of $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$ with a link to the polar assumption it corresponds to

#### Scenario: Independence
- **WHEN** a vertex, graph or trace record is inspected
- **THEN** it is derived from the Cartesian Lagrangian and cites no polar vertex or polar graph as its input

### Requirement: Cartesian vertices in the invariants
Vertices with $s$, $p$ and meson legs SHALL be derived for generic $\widehat S^C_{r,s}$ and $\widehat P^C_{r,s}$ up to five legs from the Cartesian Lagrangian, whose kinetic term for $(s,p)$ gives no vertex. The meson kinetic vertices are those of the chiral Lagrangian.

#### Scenario: Internal checks
- **WHEN** the vertices with $n=3,4,5$ legs are verified
- **THEN** the checks are relabelling symmetry, dimensional scaling and the up-only limit, as exact symbolic identities, with no polar vertex used

#### Scenario: Tree equivalence
- **WHEN** the on-shell cubic vertex $V(s;\pi^+,\pi^-)$ is computed
- **THEN** it is compared with the polar tree $M_0$ of the existing notes using the relation between the $\widehat S^C,\widehat P^C$ and the polar coefficients, and this tree equivalence is the only comparison with the polar results used in the derivation

### Requirement: Exact dependence on the invariants and compact traces
The estimate SHALL be written in generic coefficients with exact flavour traces, pole vectors, $N$ and powers of $\epsilon$, with no expansion in $\epsilon$. Traces SHALL use only matrix products, commutators, anticommutators and sums over permutations; no indices, generator components or $E_{ij}$. All checks SHALL be exact symbolic identities with no random or numerical matrices.

#### Scenario: Test inputs
- **WHEN** a test is inspected
- **THEN** it uses no random or numerical matrices

### Requirement: Loop integrals replaced by the scalar-moment size
Each scalar moment with $r$ denominators SHALL be replaced by $(m^2)^{2-r}/(16\pi^2)$ with the polynomial coefficient kept. No loop integral, scalar loop function, loop coefficient or amplitude SHALL be evaluated.

#### Scenario: Loop value required
- **WHEN** a step would need the value of a loop integral
- **THEN** it is replaced by the size prescription

### Requirement: Estimate, ratio and large-epsilon powers
The record SHALL give $M_1^{\rm est}$ and $R=|M_1^{\rm est}|/|M_0|$ in the Cartesian invariants, the leading powers of $\epsilon$ of each graph at fixed sources and fixed $\lambda$, and the up-only limit with its region boundary.

#### Scenario: Vanishing tree
- **WHEN** the tree amplitude vanishes
- **THEN** $R$ is stated undefined and only $M_1^{\rm est}$ is reported

### Requirement: Order-of-magnitude status and diagnostic comparison
The record SHALL state that each graph class is estimated only to order of magnitude, that cancellations between graphs with the same $\epsilon$ scaling can be missed and can differ between parametrizations, that agreement with the polar estimate (including the polar up-only prefactor $-21/(64\pi^2)$) is not expected, and that neither estimate is authoritative. A comparison with the polar estimate, graph class by graph class, SHALL be given as a diagnostic after the derivation is complete.

#### Scenario: Differing results
- **WHEN** the Cartesian and polar estimates differ in power or coefficient
- **THEN** the record reports the difference without attributing it to either parametrization being correct

#### Scenario: Possible cancellations
- **WHEN** Cartesian graphs share the same leading $\epsilon$ power
- **THEN** the record states whether their leading coefficients could cancel in the full loop, without asserting that they do

### Requirement: Frame invariance and limits
The estimate SHALL be unchanged under a constant axial redefinition of the sources, checked on the traces. The limits of the polar record (absorptive parts, thresholds, cancellations, naturalness, validity of $R\ll1$) SHALL be restated.

#### Scenario: Axial redefinition
- **WHEN** the sources are redefined by a constant axial transformation
- **THEN** the estimate is unchanged for the same physical legs
