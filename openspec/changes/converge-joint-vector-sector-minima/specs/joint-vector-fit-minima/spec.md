# Spec Delta

## Purpose

Defines reproducible evidence for minima of the joint vector-sector likelihood, including competing starts, physical constraints, numerical refinement, and optimizer outcomes.

## ADDED Requirements

### Requirement: Reproduce and retain multi-start fit outcomes
The fit record MUST identify the likelihood inputs, numerical settings, parameter scaling, bounds, permitted physical branches, and each start and its optimizer outcome so that attempted minima can be reproduced and compared.

#### Scenario: Multiple permitted starts are evaluated
- **WHEN** the joint likelihood is minimized from more than one physically permitted initial point
- **THEN** the result records each initial point, branch, bounds, scaling, quadrature settings, final coordinates, residual diagnostics, objective, and optimizer termination details, including unsuccessful attempts

#### Scenario: Competing endpoints are found
- **WHEN** starts terminate at distinct endpoint coordinates or objective values
- **THEN** the result preserves each endpoint and its diagnostics and does not silently replace competing solutions with one selected result

### Requirement: Qualify an accepted minimum
A minimum MUST NOT be reported as accepted solely because an optimizer returned parameter values; acceptance requires successful termination, a valid physical-domain point, finite residuals, and demonstrated stability under the documented numerical refinement.

#### Scenario: Optimizer does not converge
- **WHEN** an attempt reaches an evaluation limit or otherwise reports unsuccessful termination
- **THEN** its result is retained as unsuccessful and is excluded from accepted minima

#### Scenario: Candidate passes convergence and refinement checks
- **WHEN** an attempt terminates successfully at a point satisfying the recorded physical bounds and branch conditions and the objective, predictions, and residual diagnostics remain within the stated tolerances under quadrature refinement
- **THEN** the result identifies it as an accepted minimum and records the checks and tolerances supporting that status

### Requirement: Preserve the certified likelihood and practical evaluation path
The multi-start fit MUST evaluate the same raw likelihood and shared constraints, and its results MUST remain numerically compatible with the certified quadrature and validated runtime-improved evaluation path.

#### Scenario: Runtime improvements are used for repeated fits
- **WHEN** the multi-start procedure evaluates residuals and Jacobians, including at high quadrature order
- **THEN** it retains the established propagator reuse, bounded energy-bin evaluation, and channel-specific calculation paths while preserving residual, prediction, and Jacobian agreement within documented tolerances

#### Scenario: Shared even-sector constraints are respected
- **WHEN** a start or candidate is evaluated
- **THEN** the fit preserves the audited shared $G_V$ and even-sector definitions, duplicate-measurement exclusions, covariance assumptions, and physically permitted branch conditions
