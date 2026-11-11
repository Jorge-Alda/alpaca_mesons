# Spec Delta

## Purpose

Defines how derivative and mass mixing between the ALP and neutral mesons are treated together as one inverse propagator, and how field redefinitions reorganize that mixing without changing physical poles, residues or amplitudes.

## ADDED Requirements

### Requirement: Combined inverse propagator is the physical object
The derivation SHALL write the quadratic theory in the field space $\Phi_I=(P_i,a)$ as $D(s)=sK-M^2$ and SHALL NOT define the task as separately removing kinetic mixing or mass mixing.

#### Scenario: Field redefinition $\Phi=T\Phi'$
- **WHEN** a constant invertible $T$ is applied
- **THEN** $D'(s)=T^{\mathsf T}D(s)T$, so $K'$ and $M^{2\prime}$ change individually while $\det D(s)$ zeros, null vectors and residues map covariantly

### Requirement: Invariance of poles, residues and amplitudes
The derivation SHALL show that pole positions are unchanged, eigenvectors map as $\Phi=T\Phi'$, the residues of the full propagator and the LSZ-normalized one-ALP amplitude are unchanged, and SHALL verify this symbolically.

#### Scenario: Random invertible $T$
- **WHEN** $T$ is a generic invertible matrix acting on $(P_i,a)$ with $N_f=3$ neutral mesons
- **THEN** the zeros of $\det D$, the residue of $D^{-1}$ at each pole and the amplitude $\mathcal M(a_{\rm phys},X)$ agree before and after to symbolic or exact-rational precision

### Requirement: Pole-aligned coordinates for one eigenstate
The derivation SHALL construct the coordinate choice that removes the mixing of $a'$ with all $P'_i$ at one selected pole $s=m_k^2$, state that it is not required to remove kinetic or mass mixing off that pole, and state when the pole vector $w^{(k)}(s)$ is regular.

#### Scenario: ALP-like pole
- **WHEN** the pole closest to the ALP is selected
- **THEN** the shift $P_i=P_i'+w_i^{(k)}(m_k^2)a'$ gives $m_k^2K'_{i\hat a}-M^{2\prime}_{i\hat a}=0$ and reproduces the finite-mixing record

#### Scenario: Degenerate or exceptional pole
- **WHEN** $m_k^2K-M^2$ is singular in the meson block
- **THEN** the derivation states that the shift is not defined and gives the generalized-eigenproblem alternative

### Requirement: Alignment of several poles
The derivation SHALL state necessary and sufficient conditions for a single transformation to align more than one pole, separately for transformations of the form $P_i\to P_i+w_ia$ and for general $T$, and SHALL show that the full eigenbasis alignment diagonalizes $D(s)$ at all $s$ only for a general $T$.

#### Scenario: Two poles with a shift-type map
- **WHEN** $w^{(1)}(m_1^2)\neq w^{(2)}(m_2^2)$
- **THEN** no shift of the form $P_i\to P_i+w_ia$ aligns both poles

### Requirement: Transformed vertices for multi-meson channels
The derivation SHALL give the transformed amputated vertices including nonlinear terms, and SHALL show that the physical amplitude for a multi-meson channel is independent of the coordinates.

#### Scenario: One ALP and two mesons
- **WHEN** the amplitude is computed in original coordinates with a coherent sum over mixed components and in pole-aligned coordinates with the complete transformed vertex
- **THEN** both agree; dropping the contact terms of the nonlinear map breaks the agreement, and the script reports the difference
