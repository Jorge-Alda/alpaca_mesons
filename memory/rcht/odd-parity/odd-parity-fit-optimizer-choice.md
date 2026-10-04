---
title: Optimizer choice for the odd-parity coupled fit
status: derived
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: []
provenance: assumption
---

# Optimizer Choice for the Odd-Parity Coupled Fit

## Decision

The odd-parity coupling fit uses **pure `scipy`** — `scipy.linalg` for the fixed
Cholesky whitening and `scipy.optimize.least_squares` (trust-region reflective,
ordinary linear loss) for the minimization. No probabilistic-programming library
(PyMC, NumPyro, emcee, Stan) is used, and none is installed. This is an explicit
methodological choice, not a claim that sampling is never appropriate in this
project.

## Reasons

1. **The objective is already quadratic.** The likelihood is a fixed-covariance
   Gaussian assembled from analytic amplitudes, so a whitened nonlinear
   least-squares problem is the exact statement of the inference. No
   log-density, prior hierarchy, or intractable normalization is present.
2. **The binding constraint is rank, not curvature.** The connected first-release
   problem is rank deficient (one data block plus five OPE rows in an
   eleven-coupling space), as recorded in
   [the infrastructure record](correlated-odd-parity-fit-infrastructure.md).
   Sampler diagnostics such as $\hat r$ or effective sample size say nothing about
   a flat direction: NUTS would wander the six-dimensional null space and
   reproduce the starting point. The quantities actually needed are the singular
   values, null directions, and condition number of the Fisher information, which
   `scipy.linalg.svd` supplies directly and which the infrastructure reports.
3. **Determinism is a project requirement.** Reproducible research records
   demand a bitwise-repeatable result. A trust-region solve is deterministic; a
   gradient-based sampler introduces a seed dependence and version-dependent
   autodiff behavior, and would require $\hat r$/ESS/divergence bookkeeping that is
   overhead for a quadratic objective.
4. **Cost.** The likelihood evaluates in microseconds for eleven parameters.
   Routing it through a tensor-graph backend buys nothing and adds a heavy
   dependency tree to a project whose numerical requirements are numpy and scipy.

## What would change this decision

Sampling becomes appropriate only if the inference target stops being a
quadratic maximum-likelihood point estimate — for example a genuinely
non-Gaussian systematics model, a non-Gaussian data likelihood, or a
hierarchical prior over extracted tables. Even then, rank deficiency must be
resolved first, either with additional data or with explicit priors on the null
directions; sampling a flat posterior is not a remedy for it. Error propagation
from a rank-deficient Gaussian fit is a linear-algebra and profile-likelihood
task, not a sampling task.

## Provenance

No literature input is involved. The inputs are the structure of the likelihood
defined in [the infrastructure record](correlated-odd-parity-fit-infrastructure.md)
and the measured rank deficiency recorded in
[the data-implementation record](odd-parity-data-implementation.md).