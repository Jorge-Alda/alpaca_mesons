# Spec Delta

## Purpose

Defines how the inadequacy of the joint vector-sector fit is diagnosed and how switchable model extensions are tested against the adequacy criterion.

## ADDED Requirements

### Requirement: Diagnostics before new operators
The work MUST run the diagnostics that need no new operator before any extension is fitted, and MUST record their results.

#### Scenario: Width diagnostics
- **WHEN** the diagnostics are run at the accepted endpoint `basin_probes/even_shared_widths`
- **THEN** the record gives the on-shell $\Gamma_{K^*}$ against the tree-level RChT width at physical masses, and the minimum $\chi^2$ of the two $K^*$ width rows over any common width

#### Scenario: Channel and block diagnostics
- **WHEN** the $\phi\to K\bar K$, $\eta'$ and mass-block diagnostics are run
- **THEN** the record gives the whitened $\phi\to K\bar K$ residuals rotated to the $K^+K^-$ and $K^0\bar K^0$ channels, the refit without the $\eta'$ radiative rows with the resulting $\Gamma_{\eta'}$ and odd-block $p$-values, and the mass rows that move in the refit without the mass block

### Requirement: Switchable extensions that reduce to the baseline
Each extension MUST be a likelihood option whose default reproduces the baseline residuals to machine precision, MUST state its added coordinates with bounds, and MUST recount the degrees of freedom.

#### Scenario: Extension off
- **WHEN** the likelihood is evaluated at an accepted endpoint with every extension option at its default
- **THEN** the residual vector equals the baseline residual vector to machine precision

#### Scenario: Extension fitted
- **WHEN** an extension is refit
- **THEN** the ledger keeps its definition, added coordinates, optimizer and quadrature gate outcomes, $\chi^2$, degrees of freedom, per-block $\chi^2$ and nominal $p$-values, including unsuccessful refits

### Requirement: Operator contributions are derived before use
An extension that adds an operator contribution MUST have its amplitude derived by an executable SymPy script and recorded in a memory file before it enters the likelihood.

#### Scenario: Double-trace source
- **WHEN** the double-trace source term $y\,\mathrm{tr}(\chi_+Q)\mathbb{1}$ is considered
- **THEN** a derivation record states whether it enters the vector mass spectrum and the $KK\pi$ current, and the extension is run only if it does

### Requirement: Adequacy outcome
The work MUST re-apply the adequacy criterion of `assess-joint-vector-fit-adequacy` to each extended model and MUST end with either an adequate model or a documented reason why the model cannot meet the criterion.

#### Scenario: Criterion met
- **WHEN** an extended model has total nominal $p>0.01$ and no block with $p<10^{-3}$
- **THEN** the record names the extensions required and the accepted endpoint, and alp-uzh.7 is marked as unblocked

#### Scenario: Criterion not met
- **WHEN** all planned extensions have been tried without meeting the criterion
- **THEN** the record states the rows or blocks that remain, the irreducible $\chi^2$ they imply, and the data-level or model-level reason, with the explicit list of tried extensions
