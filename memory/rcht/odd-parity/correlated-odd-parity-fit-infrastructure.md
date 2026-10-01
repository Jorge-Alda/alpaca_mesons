---
title: Fixed-prior correlated odd-parity fit infrastructure
status: mixed
created: 2026-09-30
updated: 2026-09-30
bibtex_keys: [Alda:2026bar, Ruiz-Femenia:2003jdx]
provenance: numerical-computation
---

# Fixed-Prior Correlated Odd-Parity Fit Infrastructure

## Evidence, conventions, and assumptions

The five OPE residual definitions and their ordering are literature inputs from the Ruiz--Femenia one-multiplet odd-parity RChT relations [@Ruiz-Femenia:2003jdx]. Natural units and Alda notation are used [@Alda:2026bar]. The Gaussian theory prior, structural scale prescription, fixed equicorrelation, and truncation fraction are **explicit assumptions authorized by the user**, not literature results. The user also explicitly authorizes diagonal published errors where covariance is absent, with source-wide normalization/luminosity errors represented as explicit nuisances rather than counted in the diagonal covariance. HLS is excluded.

The input conventions are those recorded in [the numerical kernel](odd-parity-rcht-numerical-kernel.md): $c_i,d_i$ and all five residuals are dimensionless, while $M_V,F_V,F$ are in a common mass unit. [`fixed_ope_prior_scales`](../../../scripts/implementation/odd_parity_rcht.py) takes one caller-supplied immutable `OPEParameters` reference point. The resulting covariance is computed once and remains fixed during minimization; fit parameters cannot alter it.

## Fixed OPE prior

For the existing residual order $(4c_3+c_1,\ c_1-c_2+c_5,\ c_5-c_6-\mathrm{SD},\ d_1+8d_2-\mathrm{SD},\ d_3-\mathrm{SD})$, define

$$
S_c=\frac{N_c}{64\pi^2}\frac{M_V}{\sqrt2F_V},\qquad
D_g=\frac{N_c}{64\pi^2}\frac{M_V^2}{F_V^2},\qquad
D_f=\frac{F^2}{8F_V^2},
$$

$$
\boldsymbol s=(S_c,S_c,S_c,D_g+2D_f,D_g+D_f),\qquad
\sigma_i=\tau s_i,
$$

and the fixed equicorrelated covariance

$$
\Sigma_{ij}=\sigma_i\sigma_j\big[(1-\rho)\delta_{ij}+\rho\big],
\qquad -\frac14<\rho<1.
$$

The baseline is $\tau=1/N_c$ and $\rho=0$. At $N_c=3$, $M_V=0.770\,\mathrm{GeV}$, $F=0.092\,\mathrm{GeV}$, and $F_V=\sqrt2F$, the dimensionless baseline standard deviations are

$$
\boldsymbol\sigma=(0.006625111362483023,\ 0.006625111362483023,\ 0.006625111362483023,\ 0.09711596828744847,\ 0.07628263495411514).
$$

Sensitivity choices are $\tau=0.2,0.5$ and $\rho=0.5$, always fixed before fitting. No covariance jitter, pseudoinverse, or parameter-dependent theory covariance is used.

## Generic fit method and API

[`correlated_fit.py`](../../../scripts/implementation/correlated_fit.py) provides frozen `ParameterSpec`, `NuisanceEffect`, `ObservationGroup`, `ParameterPrior`, `FixedOPEConstraint`, `FitProblem`, `FitResult`, and diagnostics records. Every covariance is finite, real, exactly symmetric, and positive definite; complex arrays and empty observation vectors are rejected. Each residual block is whitened with its fixed lower Cholesky factor, and `scipy.optimize.least_squares` minimizes the concatenated residual with ordinary quadratic loss. It reports separate data, parameter-prior, and OPE $\chi^2$ components; data-only and full Jacobian rank, singular values, null vectors, and condition number; and active bounds.

The stepwise numerical method is:

1. Supply parameter starts/bounds and observation values/covariances in each observable's own consistent units. Supply any parameter prior in the corresponding parameter units. For the odd-parity OPE block, obtain $A,b$ from the [numerical kernel](odd-parity-rcht-numerical-kernel.md) and the fixed covariance above; these are the only physics-specific inputs to the generic layer.
2. Resolve each serialized predictor ID through the runtime registry and evaluate its prediction in the observation ordering.
3. Form prediction-minus-observation, parameter-minus-mean, and $A\boldsymbol\theta-b$ residual blocks. Whiten each with its own precomputed Cholesky factor. If there is more than one observation, prior, or OPE covariance block in total, `FitProblem.block_independence_assumption` must explicitly state the assumption that cross-block covariances vanish. The implementation does not infer or represent cross-block covariance silently; correlated quantities must be assembled into one covariance block or their factorization must be declared.
4. For each declared standardized nuisance $\eta_s\sim\mathcal N(0,1)$, first apply multiplicative responses $\mu_i\exp[\log(1+f_s)\eta_s]$, then additive responses $a_{is}\eta_s$. Validate the nuisance parameter and standardized prior, shared source-to-parameter naming, vector lengths and represented units, and disjoint diagonal/nuisance source IDs.
5. Minimize the concatenated residual and diagnose the returned finite-difference Jacobian without changing any covariance.
6. Return fitted parameter values, residuals, convergence metadata, component $\chi^2$, data-only/full identifiability diagnostics, active bounds, and the exact serialized-problem digest.

Predictors are registered at runtime by stable string ID. Problems store only IDs and finite JSON configurations, never callables. `serialize_fit_problem` emits deterministic canonical JSON with schema `correlated-fit-problem-v3`, including nuisance metadata, source IDs, represented units, and the block-independence declaration or `null` for a single block; `fit_problem_digest` computes its SHA256 digest. No implicit diagonalization of covariance, undeclared block independence, robust loss, jitter, or pseudoinverse fallback exists.

Minimal toy invocation:

```python
registry = PredictorRegistry()
registry.register("toy.affine", lambda p, c: p["a"] + p["b"] * np.asarray(c["x"]))
problem = FitProblem(
    parameters=(ParameterSpec("a", 0.0), ParameterSpec("b", 1.0)),
    observations=(ObservationGroup([1.0, 3.0], [[1.0, 0.2], [0.2, 1.0]],
                                   "toy.affine", {"x": [0.0, 1.0]}),),
    block_independence_assumption=None,  # one covariance block only
)
result = fit_correlated(problem, registry)
```

This is an infrastructure example only; it contains no physics data or fitted $c_i,d_i$ values.

## Verification

Run:

```text
.venv/bin/python -m unittest scripts.implementation.test_odd_parity_rcht
.venv/bin/python -m unittest scripts.implementation.test_correlated_fit
.venv/bin/python -m unittest discover -s scripts/implementation -t .
```

The OPE tests reproduce the numerical standard deviations, check tight/broad scaling, enforce the strict equicorrelation SPD interval, and expose all six exact null directions of the $5\times11$ OPE system. Toy-fit tests compare a nonzero correlated residual directly with a lower-Cholesky solve, verify $J\boldsymbol n\simeq0$ for reported null vectors, and test prior and OPE components. Nuisance tests compare manual $\chi^2$, verify the $\eta=0$ and $\eta=1$ responses, positivity and the $O(f^2)$ difference from a linear response, and reject missing priors/parameters, inconsistent shared names, source double counting, and wrong additive lengths/units.

## Readiness matrix

| physical observable predictor | data ingestion | predictor | readiness |
|---|---|---|---|
| $\omega\to\pi^0\gamma$ | ready: one PDG branching fraction | ready: uneliminated Eq. (27) width and branching fraction | connected first-release observation |
| other $V\to P\gamma$ | missing | missing | not executable; flavor-resolved matching required |
| $P\to\gamma\gamma$ and transition form factors | missing | missing | not executable |
| $e^+e^-\to P\gamma$ and $e^+e^-\to VP$ | missing | missing | not executable |
| $\omega/\phi\to3\pi$ | missing | missing | not executable |

The fixed-covariance fit infrastructure is ready. The runtime-compatible width and branching predictors are documented in [the on-shell predictor record](omega-pi-gamma-rcht-predictor.md), and the one PDG branching datum is connected by the [data implementation](odd-parity-data-implementation.md). No physical fit has been run: one datum plus five fixed OPE rows leave the eleven-coupling problem rank deficient, while all other extracted tables still lack physical predictors. Other radiative channels must not be inferred from this formula: the source warns that $K^*$ and $\phi$ channels require flavor-resolved matching [@Ruiz-Femenia:2003jdx].
