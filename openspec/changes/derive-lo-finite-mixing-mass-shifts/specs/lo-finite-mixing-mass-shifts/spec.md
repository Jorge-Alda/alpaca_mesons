# Spec Delta

## Purpose

Defines the exact LO pole equations that finite ALP-meson mixing gives for the pseudoscalar mesons and the ALP, the exact relations they satisfy, the closed approximations for isolated poles with their error control, and the invariance and limiting checks they must pass. No expansion in $f_\pi/f_a$ is part of the requirements.

## ADDED Requirements

### Requirement: Exact relations without expansion
The record SHALL give statements valid at all orders in $f_\pi/f_a$: the charged mesons have zero shift when $c_q^A$, $c_q^P$, $\chi_0$ commute with the quark charge; the neutral-kaon combination that decouples from the ALP has its pole exactly at $m_{K^0,0}^2$; and the poles satisfy $\sum_nm_n^2=\operatorname{tr}(K^{-1}M^2)$ and $\prod_nm_n^2=\det M^2/\det K$. The trace and product SHALL be written explicitly in terms of the sources using the analytic $K^{-1}$.

#### Scenario: Charged meson
- **WHEN** the sources commute with the quark charge
- **THEN** $m_{\pi^\pm}$ and $m_{K^\pm}$ equal their unmixed values for every $f_a$

#### Scenario: Sum rule against the solver
- **WHEN** the solver spectrum is computed for exact rational sources
- **THEN** the sum and product of its poles equal the analytic expressions exactly

### Requirement: Exact neutral-meson and ALP poles
The record SHALL give the exact pole equations of the neutral mesons and the ALP at finite mixing, written in the canonical frame through the invariant combinations $c_q^{AP}$, $c_{GA}$, $\chi_0$, $m_0^2$, $m_{a,0}^2$ and $\operatorname{tr}\widehat S_2$, and SHALL state that the poles are the eigenvalues of the canonical mass matrix. The record SHALL give the explicit secular polynomial and the per-meson residues for diagonal $\chi_0$, including $\pi^0$, $\eta$, $\eta'$ and the $K^0$ doublet.

#### Scenario: Dependence only on invariants
- **WHEN** a constant chiral redefinition is applied to the sources
- **THEN** every pole is unchanged, and the spectrum depends on the sources only through the stated invariants

#### Scenario: Explicit polynomial against the solver
- **WHEN** the secular polynomial is evaluated for exact rational sources
- **THEN** its roots equal the solver's poles to the working precision

### Requirement: Interlacing bounds
The record SHALL state that the poles of the canonical mass matrix interlace with the unmixed meson masses, and SHALL state the bound under which each interval holds one pole.

#### Scenario: Interlacing
- **WHEN** the unmixed meson masses and the canonical mass matrix are known
- **THEN** each unmixed mass lies between two consecutive poles, and the poles lie in the intervals stated by the record

### Requirement: Closed approximations for isolated poles
The record SHALL give a closed form for an isolated meson pole, with an exact error bound, and a closed form and a convergent iteration for an isolated ALP pole. Both SHALL be labelled as approximations under the far-apart assumption, and SHALL state the condition under which they hold.

#### Scenario: Isolated meson
- **WHEN** the meson is far from the other unmixed masses in the sense stated by the record
- **THEN** the two-level closed form approximates the pole and the error bound holds

#### Scenario: Isolated ALP
- **WHEN** the ALP entry is far from the meson masses and the convergence condition holds
- **THEN** the fixed-point iteration converges to the ALP pole

#### Scenario: Near degeneracy
- **WHEN** two unmixed masses are closer than the stated bound allows
- **THEN** the closed approximations are declared invalid and the exact eigenvalues of the canonical mass matrix apply

### Requirement: Limits and validity domain
Every pole SHALL reduce to the unmixed chiral mass in the limit $f_a\to\infty$. The record SHALL state the level-repulsion sign rule and the regimes: in the small-mixing regime no shift at order $f_\pi/f_a$ is kept, and in the large-mixing regime the finite shifts are the exact poles or their closed approximations.

#### Scenario: Decoupled ALP
- **WHEN** all couplings $\operatorname{tr}(\lambda_n\mathcal C)$ vanish
- **THEN** the ALP pole is the canonical corner entry and the meson poles are the unmixed masses

### Requirement: Scope boundaries
The record SHALL state that decay constants, residues, mixing angles, other meson multiplets, NLO and loop effects and uncertainties are not covered, and SHALL NOT use phenomenological numerical inputs.

#### Scenario: Open items
- **WHEN** alp-sjp.6 is reviewed
- **THEN** the remaining items are listed and the issue stays open
