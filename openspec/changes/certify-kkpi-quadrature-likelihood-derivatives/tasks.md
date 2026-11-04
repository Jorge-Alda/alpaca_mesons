# Tasks

## 1. Independent refinement and output comparisons

- [x] 1.1 Have the coder extend the joint KKπ fit tooling to refine energy, pair, and angular orders independently and compare predictions, standardized residuals, Jacobian columns, and χ² against the declared tolerances; verify with focused tests that each dimension is varied independently and each comparison reports pass or failure.
- [x] 1.2 Have the coder add tests for pole-aware ϕ/K* integration, failed maximum-order convergence, and invalid sampled points; verify the focused quadrature certification test suite passes and preserves failure diagnostics.

## 2. Fitted-region sampling and reproducible results

- [ ] 2.1 Have the coder add deterministic selection of the reference point, available competing minima, and representative valid fitted-region displacements; verify repeated runs produce the same complete point vectors and provenance labels.
- [ ] 2.2 Have the coder write the machine-readable certification result with data hashes, command, tolerances, orders, mappings, derivative settings, comparisons, and outcomes; verify a saved result can be reproduced from its recorded invocation and reports every sampled point.
- [ ] 2.3 Run the certification on the documented representative point set at sufficient maximum orders and store the resulting evidence; verify all points meet the declared tolerances or retain explicit uncertified status with the limiting comparisons.
