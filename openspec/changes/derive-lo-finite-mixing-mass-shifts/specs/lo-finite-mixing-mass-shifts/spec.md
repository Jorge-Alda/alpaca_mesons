# Spec Delta

## Purpose

Defines the analytical LO mass shifts that finite ALP-meson mixing induces on the pseudoscalar mesons and on the ALP, the exact relations they satisfy, and the invariance and limiting checks they must pass.

## ADDED Requirements

### Requirement: Exact relations without expansion
The record SHALL give statements valid at all orders in $f_\pi/f_a$: the charged mesons have zero shift when $c_q^A$, $c_q^P$, $\chi_0$ commute with the quark charge; the neutral-kaon combination that decouples from the ALP has its pole exactly at $m_{K^0,0}^2$; and the poles satisfy $\sum_nm_n^2=\operatorname{tr}(K^{-1}M^2)$ and $\prod_nm_n^2=\det M^2/\det K$. The trace and product SHALL be written explicitly in terms of the sources using the analytic $K^{-1}$.

#### Scenario: Charged meson
- **WHEN** the sources commute with the quark charge
- **THEN** $m_{\pi^\pm}$ and $m_{K^\pm}$ equal their unmixed values for every $f_a$

#### Scenario: Sum rule against the solver
- **WHEN** the solver spectrum is computed for exact rational sources
- **THEN** the sum and product of its poles equal the analytic expressions exactly

### Requirement: Explicit neutral-meson shifts
The record SHALL give closed analytical expressions for the shifts of the $\pi^0$, $\eta$, $\eta'$ and $K^0$ squared masses at $O(\epsilon^2)$ and $O(\epsilon^4)$, $\epsilon=f_\pi/f_a$, in terms of the unmixed masses, $m_{a,0}^2$ and the invariant combinations $c_q^{AP}$ and $c_{GA}$ rather than the individual $c_q^A$, $c_q^P$, $c_G$. The expressions SHALL state the isospin-limit form (explicit $\eta$-$\eta'$ rotation) and the general-$m_q$ form through the eigenvectors of $M^2_{PP}$.

#### Scenario: Dependence only on invariants
- **WHEN** a constant chiral redefinition is applied to the sources
- **THEN** every shift is unchanged, and it depends on the sources only through the stated invariants

#### Scenario: Weak-mixing agreement
- **WHEN** the exact LO poles are expanded in $\epsilon$
- **THEN** the closed forms reproduce the $\epsilon^2$ and $\epsilon^4$ coefficients of the weak-mixing record

### Requirement: ALP mass shift
The record SHALL give $m_a^2-m_{a,0}^2$ analytically at $O(\epsilon^2)$ and $O(\epsilon^4)$, separated into the kinetic-normalization, anomaly and meson-exchange contributions, together with the exact implicit equation for $m_a^2$ in the frame aligned with the ALP pole.

#### Scenario: No mesons
- **WHEN** all couplings that mix the ALP with mesons vanish
- **THEN** only the anomaly contribution to $m_a^2$ remains, and it vanishes for $c_G=0$

#### Scenario: Aligned frame
- **WHEN** the single-pole aligned coordinates are used
- **THEN** the ALP pole is the zero of one scalar function, and its expansion agrees with the unaligned result

### Requirement: Limits and validity domain
Every shift SHALL reduce to zero as $f_a\to\infty$ and to the unmixed chiral masses at $\epsilon=0$. The record SHALL state the level-repulsion sign rule, the condition $|m_{j,0}^2-m_{a,0}^2|\gg|D_{ja}|$ under which the expansion holds, and SHALL give the exact two-level formula for the near-degenerate case.

#### Scenario: Near degeneracy
- **WHEN** $m_{a,0}^2$ approaches an unmixed meson mass
- **THEN** the perturbative expressions are declared invalid and the two-level formula applies

### Requirement: Scope boundaries
The record SHALL state that decay constants, residues, mixing angles, other meson multiplets, NLO and loop effects and uncertainties are not covered, and SHALL NOT use phenomenological numerical inputs.

#### Scenario: Open items
- **WHEN** alp-sjp.6 is reviewed
- **THEN** the remaining items are listed and the issue stays open
