---
title: Goodness of fit at the four accepted joint vector-sector minima
status: numerical
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh]
provenance: numerical-computation
---

# Goodness of fit at the accepted minima

## Inputs and provenance

- Accepted endpoints: `data/vector_sector/joint_vector_endpoint_comparison.json` ([endpoint comparison](joint-vector-fit-endpoint-comparison.md)).
- Likelihood, blocks and baseline assumptions: [baseline record](babar-kkpi-global-vector-sector-fit-baseline.md).
- Procedure, degrees-of-freedom counting and the per-block $p$-value convention: [assessment method](joint-vector-fit-adequacy-assessment-method.md).
- Ledger: `data/vector_sector/joint_vector_fit_adequacy.json`, written by `scripts/assess_joint_vector_fit_adequacy.py`. The tables below are printed by `.venv/bin/python scripts/summarize_joint_vector_fit_adequacy.py`.
- Quadrature order $(4,128,128)$.

## Assumptions

- The total degrees of freedom are 75 residuals minus 32 free coordinates, which gives 43.
- The nominal block $p$-values use the row count of the block as its degrees of freedom. This is an assumption: the fitted parameters are shared between blocks.
- The $p$-values are not corrected for bounded or nonlinear parameters.

## Result

All four minima agree to within $10^{-3}$ in every block. The values at `basin_probes/even_shared_widths` are:

| block | rows | $\chi^2$ | nominal $p$ |
|---|---|---|---|
| odd | 12 | 35.181 | $4.4\times10^{-4}$ |
| masses and leptonic widths | 7 | 12.935 | $7.4\times10^{-2}$ |
| shared total widths | 7 | 106.842 | $4.1\times10^{-20}$ |
| $\phi\to K\bar K$ | 2 | 17.879 | $1.3\times10^{-4}$ |
| BaBar $KK\pi$ | 45 | 52.826 | $0.197$ |
| normalization constraints | 2 | 0.064 | $0.97$ |
| total | 75 | 225.7267 | $1.3\times10^{-26}$ (dof 43) |

The total $\chi^2$ of each endpoint equals `common_chi2` in the comparison ledger exactly, so the difference is 0, within the required $10^{-8}$.

The largest pulls, in the sign convention stored in the ledger, are:

| rank | row | pull | $\chi^2$ |
|---|---|---|---|
| 1 | $\Gamma_{K^{*+}}$ | $-8.27$ | 68.4 |
| 2 | $\Gamma_{K^{*0}}$ | $-4.39$ | 19.3 |
| 3 | $\phi\to K^+K^-$ (whitened component 1) | $+3.50$ | 12.2 |
| 4 | $\Gamma_{\eta'}$ | $-3.31$ | 10.9 |
| 5 | $\rho^0\to\eta\gamma$ | $-3.27$ | 10.7 |
| 6 | charged $K^*$ mass | $-3.20$ | 10.2 |
| 7 | BaBar $K_SK\pi$, 1.32–1.36 GeV | $+3.09$ | 9.5 |
| 8 | $\Gamma_\phi$ | $+2.83$ | 8.0 |
| 9 | $K^{*0}\to K^0\gamma$ | $+2.46$ | 6.1 |
| 10 | $\eta'\to\gamma\gamma$ | $+2.46$ | 6.0 |

The same ten rows lead at all four minima. The $\phi\to K\bar K$ rows are Cholesky-whitened components of the two correlated branching fractions, not the raw measurements. The $\rho$ and $\omega$ widths have $|{\rm pull}|<0.6$.

## BaBar residuals as a function of energy

At the best minimum, the bin pulls ordered by energy have these signs:

- $K\bar K\pi^0$ (21 bins, $\chi^2=22.3$): `--+-+--++-+---+-+-+-+`. There are 16 sign runs against 11.3 expected, $z=+2.2$, so the sign changes more often than chance and there is no coherent trend.
- $K_SK\pi$ (24 bins, $\chi^2=30.5$): `--+++----+-+-+--+---+--+`. There are 14 runs against 12.25 expected, $z=+0.8$.

| channel | below 1.50 GeV | 1.50–1.60 GeV | 1.60–1.70 GeV | 1.70–1.78 GeV |
|---|---|---|---|---|
| $K\bar K\pi^0$, $\chi^2$ (bins) | 11.1 (7) | 3.9 (5) | 3.8 (5) | 3.5 (4) |
| $K_SK\pi$, $\chi^2$ (bins) | 21.8 (10) | 2.6 (5) | 4.5 (5) | 1.6 (4) |

The largest BaBar contribution is in $K_SK\pi$ below 1.50 GeV. No growth of the residuals toward the upper window edge is seen.

## Verification

- Ledger check: `all_chi2_match_comparison: true`.
- Unit tests: `.venv/bin/python -m unittest scripts.test_assess_joint_vector_fit_adequacy` (5 tests OK) and `scripts.test_summarize_joint_vector_fit_adequacy` (5 tests OK).
