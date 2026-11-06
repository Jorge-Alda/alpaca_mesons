# Spec Delta

## Purpose

Defines how the adequacy of accepted joint vector-sector minima is assessed and how their sensitivity to justified modelling and data-treatment choices is recorded, so that parameter estimates are interpreted only when the model describes the data.

## ADDED Requirements

### Requirement: Report goodness of fit at accepted minima
For every accepted minimum, the record MUST report the total $\chi^2$, the number of degrees of freedom and how it was counted, the nominal $p$-value, the $\chi^2$ and row count of every data block, the full standardized residual vector, and the rows with the largest absolute pulls.

#### Scenario: Accepted minimum is assessed
- **WHEN** the assessment is run on an accepted minimum
- **THEN** the output lists total and per-block $\chi^2$, degrees of freedom, nominal $p$-values, all standardized residuals, and the ranked largest pulls, together with the endpoint identifier and input hashes

#### Scenario: Endpoint is not accepted
- **WHEN** the requested endpoint is not marked accepted in the endpoint comparison ledger
- **THEN** the assessment refuses to report it as a goodness-of-fit result

### Requirement: Locate tensions between data blocks
The assessment MUST refit the likelihood with each data block removed in turn and report the change in $\chi^2$ of the remaining blocks and the prediction for the removed block.

#### Scenario: Block is removed
- **WHEN** one data block is excluded and the remaining likelihood is refit
- **THEN** the record keeps the refit's optimizer status, its $\chi^2$ per remaining block, and the pulls of the excluded block evaluated at the new endpoint

### Requirement: Predefined sensitivity variations
Each variation of fit window, resonance inputs, tensor contribution, source hypothesis, mass prescription, and covariance treatment MUST be declared before running, MUST reproduce the baseline residuals exactly when switched off, and MUST be refit with the same acceptance gates as the baseline minimum.

#### Scenario: Variation switched off
- **WHEN** the likelihood is evaluated at an accepted endpoint with all variation options at their defaults
- **THEN** the residual vector equals the baseline residual vector to machine precision

#### Scenario: Variation refit
- **WHEN** a declared variation is refit
- **THEN** the record keeps its definition and justification, optimizer and acceptance-gate outcomes, $\chi^2$, degrees of freedom, per-block $\chi^2$, and parameter shifts in units of the baseline local errors, including unsuccessful refits

### Requirement: Adequacy verdict
The record MUST state whether the retained model is adequate, adequate only under stated variations, or inadequate; an inadequate verdict MUST name the responsible blocks or rows and the model extension needed, and MUST create a follow-up issue.

#### Scenario: Model is inadequate under all variations
- **WHEN** no declared variation gives an acceptable nominal $p$-value
- **THEN** the verdict is inadequate, parameter estimates are marked as not interpretable, and a follow-up issue describes the required extension

#### Scenario: Model is adequate
- **WHEN** the baseline or a justified variation gives an acceptable nominal $p$-value
- **THEN** the verdict states the conditions under which it holds and the parameter estimates may be passed to the identifiability study
