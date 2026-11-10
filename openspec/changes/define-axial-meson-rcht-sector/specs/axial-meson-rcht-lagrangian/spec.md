# Spec Delta

## Purpose

Defines the even-parity resonance chiral Lagrangian of the light axial-vector nonet, so that every ALP application uses one set of fields, conventions and independent coefficients.

## ADDED Requirements

### Requirement: Axial field and flavor conventions
The project SHALL document the axial nonet as an antisymmetric tensor $A_{\mu\nu}$ transforming as $A\to hAh^\dagger$, with generator normalization, flavor matrix, phases and metric/Levi-Civita conventions consistent with Alda:2026nrd and the vector nonet record.

#### Scenario: Matrix content is explicit
- **WHEN** a downstream task needs the field content
- **THEN** the record gives the explicit $3\times3$ nonet matrix, the octet–singlet and ideal-basis definitions, and the sign of the strange phase

#### Scenario: Alternative representation is justified
- **WHEN** the tensor representation is chosen over the Proca (vector-field) one
- **THEN** the record states the reason and the map between them

### Requirement: Quadratic terms and nonideal mixing
The Lagrangian SHALL contain canonical kinetic terms, flavor-breaking and singlet mass terms, and a nonideal mixing angle for the neutral isoscalar states, with the mass matrix and rotation derived by a script.

#### Scenario: Mass matrix is diagonalized
- **WHEN** the neutral isoscalar mass matrix is diagonalized
- **THEN** the mixing angle, physical masses and generator rotation are given in terms of the Lagrangian parameters, with the ideal limit recovered

#### Scenario: Neglected breaking is listed
- **WHEN** a mass operator is omitted
- **THEN** the record lists it explicitly as an assumption

### Requirement: Complete even-parity operator basis
The record SHALL list every even-parity operator with one axial field needed for the ALP and strong/electromagnetic decay channels at the chosen chiral order, each with a stated normalization, and justify completeness.

#### Scenario: Completeness argument
- **WHEN** the operator list is presented
- **THEN** it states the chiral order, the building blocks used ($u_\mu$, $f_\pm^{\mu\nu}$, $\chi_\pm$, $\nabla_\mu$), and the identities (Schouten, equations of motion, Bianchi) used to remove redundant operators

### Requirement: Independent coefficients versus symmetry relations
The record SHALL separate coefficients fixed by symmetry from independent couplings and name each independent coefficient with a descriptive notation.

#### Scenario: Coefficient table
- **WHEN** the table is consulted
- **THEN** each coefficient is marked independent or fixed, with the relation and its origin (chiral symmetry, charge conjugation, isospin, SU(3) flavor limit)

#### Scenario: New notation
- **WHEN** a new symbol is introduced for a physical quantity
- **THEN** it is submitted to the user for approval before being relied on

### Requirement: Power counting and width treatment
The record SHALL state the power counting and the treatment of the axial widths, including when a state is treated as stable, when a finite-width line shape is used in exchange, and the validity limits.

#### Scenario: Broad resonance as external state
- **WHEN** an axial state is a final state
- **THEN** the record states the narrow-width condition under which it is treated as an asymptotic state, and otherwise treats it as a decay chain

#### Scenario: Width in exchange
- **WHEN** an axial resonance is exchanged
- **THEN** the propagator and width prescription are specified and gauge/Ward-identity preservation is noted
