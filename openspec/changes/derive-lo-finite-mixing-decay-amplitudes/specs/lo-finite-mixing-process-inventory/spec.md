# Spec Delta

## Purpose

Defines the inventory of ALP decay and production processes at LO with a physical mixed pole as the parent, with quantum numbers, exact selection rules, thresholds and the rules for channels that appear or disappear at finite mixing.

## ADDED Requirements

### Requirement: Process inventory with LO status
The inventory SHALL list the channels $\gamma\gamma$; $P$, $PP$, $PPP$; $\gamma P$-type and $\gamma PP$; baryon-antibaryon; vector-meson channels; $\ell^+\ell^-$; and resonance final states. Each entry SHALL carry one status: derived in this change, derived in an existing record (linked), or deferred beyond LO with the missing ingredient named (resonance exchange, loop, NLO term, width).

#### Scenario: Deferred channel
- **WHEN** a channel needs vector-meson exchange, a loop or a transition form factor
- **THEN** it is marked deferred with that reason and no LO amplitude is claimed for it

#### Scenario: Existing record
- **WHEN** a finite-mixing amplitude already exists in the repository ($a\to3\pi$, two-meson vertex, baryon vertex)
- **THEN** the inventory links it and the new construction reproduces it as a check instead of replacing it

### Requirement: Generic pseudoscalar channels and excluded $4P$
Pseudoscalar channels SHALL be listed with generic indices ($a\to P_iP_j$, $a\to P_iP_jP_k$), not as individual charge or flavour assignments. The $a\to4P$ channels are outside this change and SHALL appear only as deferred.

#### Scenario: Generic listing
- **WHEN** the three-pseudoscalar channels are listed
- **THEN** one entry $a\to P_iP_jP_k$ covers $\pi^0\pi^+\pi^-$, $\eta\pi^+\pi^-$, $\eta K^+K^-$ and the others

#### Scenario: Four pseudoscalars
- **WHEN** the inventory lists $a\to4P$
- **THEN** it is marked deferred (outside this change) and no amplitude is claimed

### Requirement: Exact selection rules from the pole vectors
For each LO channel the records SHALL state the selection rules that hold exactly at finite mixing: conserved charge, CP and flavour content of each physical pole, the exact $K^0$ pole with zero ALP component, and which components of $v^{(n)}$ must be nonzero for the channel to be open. Rules that hold only at first order in $f_\pi/f_a$ SHALL be labelled as such.

#### Scenario: CP-conserving sources
- **WHEN** the sources are CP conserving
- **THEN** CP-even and CP-odd poles do not mix, and the allowed final states follow from the CP of each pole

#### Scenario: Channel that opens only at finite mixing
- **WHEN** a final state has zero amplitude for an unmixed ALP but nonzero overlap with a meson component of the physical pole
- **THEN** the record states the condition and the vertex or pole component that opens it

### Requirement: Thresholds in physical masses
Open-channel conditions SHALL be given in terms of the physical pole masses $m_n$ returned by the spectrum solver, and the inventory SHALL state analytically, as inequalities among those masses, when each channel opens or closes. No numerical parameter values are used.

#### Scenario: Parent below threshold
- **WHEN** the threshold inequality $m_{\rm parent}<\sum m_{\rm daughters}$ holds
- **THEN** the channel is closed and no amplitude on shell exists

### Requirement: Stable and broad daughters
The inventory SHALL classify each final state as a stable pole at LO or as a broad resonance, and SHALL state the rule for rates in terms of physical daughters (no daughter width, no spectral function at LO).

#### Scenario: Broad final state
- **WHEN** a final-state meson is a broad resonance (for example in $\pi\pi$ or $K\pi$ channels)
- **THEN** the record states that the LO rate treats it as a stable daughter and that a spectral-function treatment is outside this change
