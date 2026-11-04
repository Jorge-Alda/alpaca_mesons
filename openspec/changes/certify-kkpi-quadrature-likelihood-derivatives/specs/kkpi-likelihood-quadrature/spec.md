# Spec Delta

## Purpose

This capability certifies that numerical integration of the KKπ likelihood is accurate for its predictions, residuals, derivatives, and χ² across representative fitted parameter points. It makes the numerical orders, comparison tolerances, sampled points, and invocation reproducible.

## ADDED Requirements

### Requirement: Certify the KKπ energy-bin and Dalitz integrations
The certification workflow MUST refine energy-bin, pair-invariant, and angular quadrature independently while retaining the physical integration domain and the pole-aware mappings for the narrow ϕ and K*(892) structures.

#### Scenario: Independent quadrature refinement
- **WHEN** certification is run for a representative parameter point
- **THEN** the result compares predictions using at least two increasing orders for each integration dimension, with the other dimensions held fixed during each comparison

#### Scenario: Pole structures are resolved
- **WHEN** a representative point includes the physical ϕ and K*(892) contributions
- **THEN** the certification uses the existing pole-aware integration mappings and records the pole mapping and quadrature orders used

### Requirement: Certify likelihood outputs and derivatives
For every representative point, certification MUST compare the per-bin predictions, uncertainty-normalized residuals, residual Jacobian, and χ² between successive refined quadrature evaluations. Certification MUST pass only when each reported comparison is within its declared tolerance.

#### Scenario: All comparisons converge
- **WHEN** successive refinements satisfy the declared prediction, residual, Jacobian, and χ² tolerances
- **THEN** the point is reported as certified together with the comparisons and tolerances

#### Scenario: One comparison does not converge
- **WHEN** any output exceeds its declared tolerance at the maximum configured order
- **THEN** the point is reported as uncertified with the failing output, refinement orders, and observed differences

### Requirement: Cover representative fitted-region points reproducibly
The certification MUST evaluate the reference fit point and additional documented points that span the fitted region, and MUST identify each point by its source or complete parameter coordinates. Results MUST record the exact invocation, integration orders, declared tolerances, and fit/input provenance needed to reproduce the comparisons.

#### Scenario: Reproduce a certification run
- **WHEN** a downstream user reads a stored certification result
- **THEN** the user can identify all sampled points, their fit-region rationale, the quadrature settings, the tolerances, the input provenance, and the command used

#### Scenario: A sampled point is invalid
- **WHEN** a requested point cannot be evaluated within the model's physical parameter domain
- **THEN** the run records that point as unevaluated with its reason and does not count it as certified
