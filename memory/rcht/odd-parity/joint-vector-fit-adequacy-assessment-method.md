---
title: Goodness-of-fit assessment method at accepted joint vector minima
status: numerical
created: 2026-11-06
updated: 2026-11-06
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Goodness-of-fit assessment method at accepted joint vector minima

## Purpose

Report, for each accepted minimum of the joint vector-sector likelihood, the total and per-block $\chi^2$, the degrees of freedom, nominal $p$-values, all standardized residuals with row labels, the largest pulls, and the energy dependence of the BaBar $KK\pi$ residuals.

## Inputs

- Accepted endpoints: [`joint_vector_endpoint_comparison.json`](../../../data/vector_sector/joint_vector_endpoint_comparison.json), described in [the endpoint comparison](joint-vector-fit-endpoint-comparison.md). An endpoint is used only if its `accepted` flag is true; any other requested endpoint is refused with an error.
- Likelihood: `GlobalLikelihood` of [`scripts/fit_babar_kkpi_global_vector_sector.py`](../../../scripts/fit_babar_kkpi_global_vector_sector.py), with the blocks and assumptions of [the baseline record](babar-kkpi-global-vector-sector-fit-baseline.md).

## Conventions and units

GeV and nb; residuals are standardized and dimensionless. The two $\phi\to K\bar K$ rows are Cholesky-whitened with the PDG correlation $-0.72$: the first component is the $K^+K^-$ pull alone, the second is the $K_LK_S$ pull conditional on the first. Quadrature $(4,128,128)$ with four-bin energy blocks.

## Method

1. Evaluate the 75 residuals at each accepted endpoint and check the total $\chi^2$ against the comparison ledger (tolerance $10^{-8}$).
2. Degrees of freedom: 75 residuals (73 observed summaries plus 2 normalization constraints) minus 32 free coordinates gives 43. **Assumption:** the per-block nominal $p$-value uses the block row count as degrees of freedom; fitted parameters are not subtracted, so these values only rank blocks. No correction for bounded or nonlinear parameters is applied.
3. Rank the 15 largest absolute pulls.
4. For each BaBar channel, order residuals in energy, record the sign sequence, the Wald–Wolfowitz runs count with its expectation and standard deviation, and the $\chi^2$, mean pull and sign counts in the a priori ranges $E<1.50$, $1.50$–$1.60$, $1.60$–$1.70$ and $1.70$–$1.78$ GeV.
5. Compute linearized local errors $\sqrt{\mathrm{diag}((J^TJ)^{-1})}$ from the $(4,128,128)$ Jacobian, as in the fit script. They are used only as units for parameter shifts in the sensitivity study.

## Output and verification

Executable: [`scripts/assess_joint_vector_fit_adequacy.py`](../../../scripts/assess_joint_vector_fit_adequacy.py); ledger: [`joint_vector_fit_adequacy.json`](../../../data/vector_sector/joint_vector_fit_adequacy.json), with input and source SHA-256 values.

```text
.venv/bin/python scripts/assess_joint_vector_fit_adequacy.py --output data/vector_sector/joint_vector_fit_adequacy.json
```

Checks: all four totals reproduce the comparison ledger exactly (difference 0). The unit tests in [`scripts/test_assess_joint_vector_fit_adequacy.py`](../../../scripts/test_assess_joint_vector_fit_adequacy.py) check refusal of non-accepted and unknown endpoints, block sums, row counts, labels, ordering of the pulls and the runs statistic.
