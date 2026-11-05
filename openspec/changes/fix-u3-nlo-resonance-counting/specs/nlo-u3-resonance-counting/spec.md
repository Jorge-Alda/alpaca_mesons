# Spec Delta

## Purpose

Fixes, in a single set of auditable memory records with an executable check, the expansion, field content, loop content, subtraction scheme and validity range of the combined $U(3)$ NLO plus explicit-resonance theory, so that every alp-gii calculation uses the same order and renormalization.

## ADDED Requirements

### Requirement: Declared field content
The memory SHALL state the fields of the theory: the pseudoscalar nonet with $\eta_0$, the ALP entering through the sources of Alda:2026nrd, and explicit scalar, pseudoscalar, vector and axial-vector resonance nonets. For each nonet it SHALL state the transformation under the chiral group and the representation of spin-one fields.

#### Scenario: Downstream issue looks up a resonance field
- **WHEN** alp-gii.2 or alp-gii.3 needs the definition of a resonance nonet
- **THEN** the counting record gives its field normalization, its chiral transformation and the memory record or reference that fixes its LO Lagrangian

### Requirement: Combined power counting
The memory SHALL declare the $\delta$ counting of Kaiser:1998ds, $p^2\sim m_q\sim 1/N_c\sim\delta$. It SHALL state the $N_c$ scaling of resonance masses and couplings, and that resonance propagators are kept unexpanded in $p^2/M_R^2$. The calculated order SHALL be $\mathcal O(\delta^2)$ in the Goldstone/ALP sector.

#### Scenario: Assigning an order to an operator
- **GIVEN** any operator of the retained basis in memory/nlo/lagrangian.md, or a resonance-coupling operator
- **WHEN** its $\delta$ order is requested
- **THEN** the record and the script give the same order, with a derivation from its powers of $p$, $m_q$ and $N_c$

#### Scenario: Assigning an order to a resonance-exchange graph
- **WHEN** a tree graph with one resonance propagator between two Goldstone/ALP vertices is classified
- **THEN** the record gives the order of its low-energy limit, and states that the full propagator is kept at the declared order

### Requirement: Included and excluded loop classes
The memory SHALL list the loop classes in the declared order: one-loop graphs with only nonet-Goldstone internal lines and LO $U(3)$ vertices, including LO ALP-source vertices. Loops with any internal resonance line SHALL be listed as excluded, together with the order they would enter at.

#### Scenario: Classifying a loop graph
- **WHEN** a one-loop graph contributing to a meson or ALP–meson two-point function or to an ALP decay vertex is considered
- **THEN** the classification states whether it is included, and gives the reason in terms of the declared counting

#### Scenario: Goldstone loop in a resonance self-energy
- **WHEN** a Goldstone loop dresses a resonance propagator, which gives the resonance width
- **THEN** the record states whether this enters the strict order or only a declared line-shape prescription near the resonance pole, and links the existing width record

### Requirement: Subtraction scheme and scales
The memory SHALL fix the dimensional-regularization subtraction $\lambda(\mu_\chi)$ of Gasser:1983yg in the $U(3)$ form used by Kaiser:1998ds, and state the default value of $\mu_\chi$. It SHALL separately state the QCD-scale dependence of the singlet constants through $Z_A$.

#### Scenario: Changing the chiral scale
- **WHEN** $\mu_\chi$ is changed from $\mu_1$ to $\mu_2$
- **THEN** the executable script runs every renormalized constant with its $\Gamma_i$, and confirms that each combination labelled scale independent in the record is scale independent

#### Scenario: Singlet constant under the QCD scale
- **WHEN** $\Lambda_1$, $\Lambda_2$, $2L_5+3L_{18}$ or $2L_8-3L_{25}$ is used
- **THEN** the record states its $Z_A$ transformation and the QCD scale at which numerical inputs are quoted

### Requirement: Renormalized counterterm basis
The memory SHALL tabulate every local constant entering the quadratic sector or the ALP-source vertices: $L_4$–$L_8$, $L_{18}$, $L_{25}$, $\Lambda_1$, $\Lambda_2$, $L_{47}$ and any further ALP-source operator. For each it SHALL give the $\delta$ order, the chiral counterterm coefficient in the declared scheme, the QCD-scale behaviour and the resonance multiplets that contribute at tree level.

#### Scenario: Divergence of an included loop graph
- **WHEN** an included loop class yields a UV divergence in a two-point function
- **THEN** the table names the counterterms that absorb it, with coefficients matching Kaiser:1998ds and Herrera-Siklody:1996tqr

#### Scenario: Constant without a chiral pole
- **WHEN** a constant has zero chiral counterterm at the declared order, for example $L_{25}$, $\Lambda_{1,2}$ or $L_{47}$ under the truncation of Herrera-Siklody:1996tqr
- **THEN** the record states the truncation under which the zero holds and the order at which it fails

### Requirement: Separation of resonance exchange and residual constants
The memory SHALL define the residual local constants left after explicit resonance exchange is removed, and the scale at which that separation is made. It SHALL state that a contribution kept as an explicit resonance propagator is not also included in a local constant.

#### Scenario: Using an explicit scalar resonance with $L_5$
- **WHEN** the scalar octet is kept explicit and $L_5$ enters a two-point function
- **THEN** only the residual $L_5$ is used, as defined in the record, and the remaining matching is referred to alp-gii.3

### Requirement: Validity range and truncation uncertainty
The memory SHALL state that the framework applies to external momenta and ALP masses up to about $1\,$GeV. It SHALL give a reproducible prescription for the truncation uncertainty that covers the omitted $\mathcal O(\delta^3)$ terms and the excluded resonance loops.

#### Scenario: Prediction at a given ALP mass
- **WHEN** a downstream calculation quotes a result at an ALP mass $m_a\le 1\,$GeV
- **THEN** the truncation prescription gives an uncertainty as a function of $m_a$ and of the meson masses, with each input identified

#### Scenario: Outside the range
- **WHEN** $m_a$ exceeds about $1\,$GeV or approaches a resonance pole where the narrow-width tree propagator fails
- **THEN** the record states that the declared counting does not apply there, and which prescription, if any, replaces it

### Requirement: Executable and auditable record
Each statement of order, counterterm coefficient and scale behaviour SHALL be reproduced by a deterministic script under `scripts/`. Every memory record SHALL carry the required YAML front matter and BibTeX keys, and SHALL label assumptions explicitly.

#### Scenario: Running the check
- **WHEN** the script is run with the documented invocation
- **THEN** it prints the order table and counterterm table, and exits successfully only if all consistency checks pass
