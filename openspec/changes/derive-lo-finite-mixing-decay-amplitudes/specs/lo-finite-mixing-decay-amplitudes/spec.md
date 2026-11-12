# Spec Delta

## Purpose

Defines how LO tree amplitudes for a physical mixed pole are built: pole projection of every external leg, complete transformed vertices including nonlinear contact terms, exact internal propagators, and the symmetry checks the amplitudes must pass.

## ADDED Requirements

### Requirement: Pole projection on every external leg
The amplitude SHALL contract the amputated vertex with the $K$-orthonormal pole vector $v^{(n)}$ of the parent and of every outgoing meson that is a pole of the same $(P_i,a)$ system, using the spectrum from `lo-finite-mixing-spectrum`. The amplitude SHALL NOT expand $v^{(n)}$ in $f_\pi/f_a$, and the ALP component of an outgoing meson pole SHALL be included.

#### Scenario: Meson daughter with ALP admixture
- **WHEN** an outgoing $\pi^0$-like pole has a nonzero $a$ component
- **THEN** vertices with that leg replaced by $a$ contribute coherently to the amplitude

### Requirement: Generic pseudoscalar amplitudes
Amplitudes for pseudoscalar final states SHALL be derived for generic pseudoscalar indices, as functions of the nonet indices $i,j,k$, the flavour matrices of the generators, the pole vectors and the momenta, and SHALL NOT be derived separately for individual particle assignments. Specific channels SHALL be obtained by substituting the indices, and the substitution SHALL be checked against any existing channel-specific record. The $a\to4P$ channels are not part of this capability.

#### Scenario: Specific channel from the generic result
- **WHEN** the generic $a\to P_iP_jP_k$ amplitude is evaluated for $\pi^0\pi^+\pi^-$
- **THEN** it equals the existing finite-mixing three-pion amplitude, and no channel-specific derivation is needed for $\eta\pi^+\pi^-$ or $\eta K^+K^-$

### Requirement: Joint vertex generation
Vertices SHALL come from the expansion of the LO Lagrangian jointly in $(P_i,a)$ up to the order of the channel, including the exact meson exponential, with the ALP dependence entering through the full functions $M(x)$, $C_{R,L}(x)$ and $\Theta(x)$ ($x=a/f_a$) and NOT through an assumed linear dependence on $a/f_a$, in the conventions of Alda:2026nrd. Contact vertices of the form $a^kP^m$ with $k\ge1$ SHALL be kept, with coefficients given by the Taylor coefficients of those functions at $x=0$.

#### Scenario: Several ALP-component legs
- **WHEN** two or more legs carry an $a$ component
- **THEN** the terms with $a^2$ or higher powers in the Lagrangian expansion are included

#### Scenario: Nonlinear source functions
- **WHEN** $M(x)$, $C_{R,L}(x)$ or $\Theta(x)$ have nonzero second or higher Taylor coefficients
- **THEN** the vertices with $k\ge2$ ALP legs depend on those coefficients, and the linear-spurion result is recovered when they vanish

### Requirement: Exact internal propagators
Exchange diagrams SHALL use the exact $D^{-1}(s)$ of the combined system at the exchanged momentum, so that mixed propagation between meson and ALP components is kept off the pole.

#### Scenario: Pole versus full propagator
- **WHEN** the exchanged momentum is far from every pole
- **THEN** the amplitude uses $D^{-1}(s)$ and not only the pole sum, and the two agree on the pole to the stated tolerance

### Requirement: Complete transformed vertices
For each derived channel the amplitude SHALL also be obtained in the pole-aligned coordinates, with every nonlinear term of the exponential chiral map included, and SHALL equal the pole-projection result. The linear shift alone SHALL NOT be accepted for multi-meson channels.

#### Scenario: Agreement of the two routes
- **WHEN** a channel with at least two outgoing mesons is evaluated by both routes
- **THEN** the amplitudes agree at random kinematic points, and the contribution of the nonlinear contact terms is reported separately

### Requirement: Symmetry checks
Each channel SHALL be checked for Bose symmetry and the symmetry factor of identical daughters, crossing between decay and production, Ward identities for photon channels, invariance under constant invertible field redefinitions, and the weak-mixing limit against Alda:2026nrd.

#### Scenario: Photon channel
- **WHEN** an amplitude has an external photon
- **THEN** replacing the polarization by the photon momentum gives zero for the full amplitude at finite mixing

#### Scenario: Weak mixing
- **WHEN** $f_a/f_\pi$ is large
- **THEN** the amplitude equals the first-order result of Alda:2026nrd with a residual that scales as the next order

### Requirement: Interference and coherence
The amplitude SHALL be the coherent sum of the direct ALP-component and meson-component contributions, and the record SHALL give the relative contribution of each piece so that cancellations are visible.

#### Scenario: Cancellation
- **WHEN** direct and meson-component terms nearly cancel
- **THEN** the report shows each term and the sum with the cancellation factor

### Requirement: Limits of validity
Each amplitude record SHALL state that it is tree level at LO, that the poles are those of `lo-finite-mixing-spectrum`, and SHALL flag results where the spectrum diagnostics (near-degeneracy, $O(1)$ admixture) mark the LO model as unreliable.

#### Scenario: Near-degenerate parent
- **WHEN** the parent pole is flagged near-degenerate with a meson pole
- **THEN** the amplitude is reported for the degenerate subspace projector and flagged outside the controlled domain
