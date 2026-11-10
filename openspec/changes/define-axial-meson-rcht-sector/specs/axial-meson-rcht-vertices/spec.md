# Spec Delta

## Purpose

Provides reproducible, script-verified Feynman rules for the even-parity axial-meson interactions with pseudoscalars, vectors, photons and ALP sources, as the single input for later ALP decay amplitudes.

## ADDED Requirements

### Requirement: Executable vertex derivations
Each vertex ($A$–$P$ mixing, $A$–$V$ mixing, $AVP$, $APP$ where allowed, $A\gamma P$, $A V\gamma$, and the ALP-source vertices) SHALL be derived by a deterministic `sympy` script with a documented invocation and a paired step-by-step memory record.

#### Scenario: Script reproduces the record
- **WHEN** the script is run as documented
- **THEN** it prints the vertex expressions listed in the memory record and exits with status zero

#### Scenario: Inputs are traceable
- **WHEN** a derivation record is read
- **THEN** every input is linked to a memory record or cited by BibTeX key

### Requirement: Channel coverage
The vertices SHALL cover every channel in the ALP axial process inventory (axial external states and axial exchange), and any channel not covered SHALL be recorded as excluded with a reason.

#### Scenario: Inventory cross-check
- **WHEN** the inventory of ALP processes with axial states is compared with the vertex list
- **THEN** each process maps to vertices or to a recorded exclusion

### Requirement: Quadratic mixing is resolved
The derivation SHALL remove or diagonalize any $A$–$P$ and $A$–$V$ quadratic mixing and give the resulting field redefinitions and normalization.

#### Scenario: No residual mixing
- **WHEN** the redefinitions are applied
- **THEN** the quadratic Lagrangian is diagonal in the physical fields and the script checks this symbolically

### Requirement: Representation equivalence and induced contacts
The record SHALL give the field redefinitions linking equivalent representations and every induced local term, and a script SHALL verify that on-shell amplitudes agree.

#### Scenario: Contact terms listed
- **WHEN** a redefinition generates local operators
- **THEN** each is listed with its coefficient as a function of independent couplings

#### Scenario: Amplitude equality
- **WHEN** the same process is computed in both representations
- **THEN** the amplitudes agree after adding the induced contacts

### Requirement: ALP source conventions and Ward identities
The ALP vertices SHALL use the no-ALP-meson-mixing source convention of Alda:2026nrd, and the vertices SHALL obey the chiral Ward identities and crossing.

#### Scenario: Ward check
- **WHEN** a source vertex is contracted with its momentum
- **THEN** the result matches the Ward identity of the Lagrangian, verified by a script

### Requirement: Coefficient and contact handoff table
The change SHALL deliver one table of all independent coefficients and induced contacts, with their vertex dependence, for use by alp-ct8.3 to alp-ct8.6, without importing numerical values from the conditional vector fits.

#### Scenario: Handoff completeness
- **WHEN** alp-ct8.3 starts
- **THEN** the table lists each coefficient, the vertices and channels where it appears, and marks no numerical input as established
