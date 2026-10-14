---
title: Conditional low-energy odd-action calibration and adequacy
status: mixed
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# Scope and frozen policy

This executes stage D of the [maintained calibration plan](odd-breaking-low-energy-calibration-plan.md), following the [completed fixed-scenario measurement Jacobians](odd-low-energy-measurement-jacobians.md). It is a conditional calibration of the pinned low-energy summaries, with no ultraviolet matching constraint or coefficient prior.

The [PDG snapshot](../../extracts/ParticleDataGroup:2024cfk/odd-parity-low-energy-measurements.md) supplies twelve accepted odd observations. The four mass observations, three electronic-width observations and all seven original shared total widths are included once in each dataset control. Stage-3 fitted values supply initial points only. The missing joint covariance is explicitly replaced by diagonal quoted-summary errors for this diagnostic calculation; removing broad-rho or constrained-fit K-star summaries does not supply that covariance. Alda neutral central inputs remain fixed, with LO source-compatible vertex masses and physical phase-space masses treated separately.

All 64 mass/input/source alternatives remain separate. The raw parent filters distinguish the original eleven operators, soft leading parents, leading double-trace additions, single-trace spurion/derivative additions, their first-correction union, promoted leading parents and the full bounded action. Parent filters act before observable reduction. Saturated models test interpolation and numerical implementation; their improved in-sample residuals do not establish a controlled hierarchy.

# Inputs and assumptions

- [Bounded numerical SM kernels](odd-bounded-sm-numerical-kernels.md): normalized real-photon channel gradients, WZW offsets and original-coupling embedding.
- [Measurement Jacobians](odd-low-energy-measurement-jacobians.md): separate scenario definitions, measurement units, shared-width treatment, null directions and conditional ranks.
- [Power-counting proposal](odd-action-power-counting-proposal.md): proposed parent topology/derivative filtration, without numerical naturalness penalties.
- [Stage-3 calibration](../foundations/vector-stage3-calibration.md): mass/source action and starting values; its data are counted once.
- [Alda input conversion](alda-neutral-amplitude-inputs.md): scalar generators, spurions, field conventions and stated central-input approximation [@Alda:2026nrd].

The selected likelihood is quadratic and uses no inherited OPE Gaussian row, implicit OZI penalty, or robust-loss downweighting. The on-shell broad-vector approximation and partial-isospin treatment remain model limitations. Local covariance and nominal rank-based degrees of freedom do not establish a globally unique nonlinear/sign solution.

# Interpretation constraints

Width measurements leave discrete amplitude signs. Signed-amplitude target searches initialize nonlinear calibration, retaining the WZW offset. Representative parent coordinates are a coefficient gauge, not measurements of all action coefficients. Their physical amplitude map must be regenerated as the even/source nuisances move, and the retained-column image must still be contained in the representative image. An invisible direction at one reference point cannot be assumed invisible throughout a nuisance profile.

Retaining separate scenario labels requires preserving the selected inverse-action mass root and source-current signs throughout minimization; a starting point alone does not enforce this. Gauge-image failure invalidates physical-action inference even if the optimizer converges. Optimization success and scientific validity are separate recorded statuses.

Use changes of identifiable channel amplitudes for nested correction diagnostics. A coefficient-level topology decomposition can depend on the null-coordinate gauge. Held-out channel predictions require the training nullspace to annihilate the held-out row; saturated models generally lose that predictability when one independent row is withheld. No forecast is obtained by silently setting the new null coordinate to zero.

ALP uncertainty propagation remains a separate gate: this SM fit does not establish that unmeasured action directions vanish in a selected ALP target. The normalized target sensitivity and source/contact calculation must expose those directions before an uncertainty band is meaningful.

# Method and reproduction

Run the repository environment:

```text
.venv/bin/python -m unittest scripts/implementation/test_odd_low_energy_calibration.py
.venv/bin/python scripts/calibrate_odd_low_energy.py
.venv/bin/python scripts/verify_odd_low_energy_calibration.py
```

The [calibration script](../../../scripts/calibrate_odd_low_energy.py) precomputes multilinear flavor contractions of the existing normalized kernel and regenerates their mass/source-dependent amplitude map at every nuisance point. The fast contractions are checked against the original implementation at displaced mass/source inputs and each scenario start. Width derivatives and the exact inverse-mass-root separator are generated with SymPy from the observable and existing mass determinant polynomial; their step-by-step expressions are saved in the numerical transcript.

For each dataset and parent model, a pivot image supplies explicitly labelled representative odd coordinates. The script enumerates all signed amplitude targets at the ancillary optimum and uses deterministic distinct candidates to initialize the joint nonlinear quadratic minimization. Exhaustive fixed-point sign enumeration does not prove a global nonlinear minimum. Saturated images instead solve the signed-amplitude interpolation problem exactly at the ancillary optimum; that is an implementation/interpolation diagnostic.

The [verification script](../../../scripts/verify_odd_low_energy_calibration.py) checks input hashes, complete scenario/model/control coverage, inverse mass branches, source signs, image certificates and optimizer statuses. It saves per-model/control summaries and nested changes of measured predictions. The [tests](../../../scripts/implementation/test_odd_low_energy_calibration.py) exercise displaced-input kernel agreement, sign/phase covariance, detection of a nuisance-induced gauge-image lift and correct versus wrong inverse branches.

# Executed conditional result

All 64 scenarios, seven parent models and four dataset selections were run: **1792 calibrations**. Every optimizer reported convergence; all 1792 fitted inverse-mass-root and source-sign certificates passed. No fixed-gauge image failure or active numerical search bound was found. These certify the declared local calculation, not a global nonlinear minimum or the missing experimental covariance.

Only the original eleven-operator image and the soft-leading image remain unsaturated. Their measured images coincide at the fitted nuisance points: the largest bidirectional containment defect is $7.20\times10^{-15}$. Each of the separate leading double-trace and single-trace correction extensions already saturates the odd rows in every dataset selection; so do their union, promoted leading and the full bounded image. Their interpolation cannot select trace topology or establish numerical correction suppression.

The smallest achieved soft-leading total residuals across the retained scenarios are:

| Dataset selection | Achieved total $\chi^2$ |
|---|---:|
| Primary twelve summaries | 5.5675233 |
| Exclude broad-rho summaries | 0.3184 |
| Exclude constrained-fit K-star summaries | 5.4353 |
| Exclude both sets | $1.88\times10^{-7}$ |

The primary feasible example uses charged-rho and charged-K-star masses, the evaluated phi electronic width, the near-ideal mass root and the negative single-trace strange-source hypothesis. Its odd and ancillary contributions are 3.7113018 and 1.8562215; its local residual degrees of freedom are six, with nominal descriptive $p=0.4733$. The largest odd-channel pull is approximately 1.13. The nominal tail and conventional pull flags are descriptive diagnostics, not retrospective preferred-model rules. This is an achieved feasible scenario, not a selected unique central solution. The near-zero narrow-control residual has little residual testing power and does not establish that higher corrections vanish.

There are 512 unsaturated fits. Twelve approach a source-domain boundary and omit symmetric Gaussian uncertainty claims; the remaining **500** have full retained-coefficient/nuisance local identifiable prediction covariance and joint held-out null tests. The largest derivative step-halving change is $1.92\times10^{-8}$. Null-sensitive held-out predictions remain explicitly withheld. These are conditional local Gaussian experimental/input uncertainties without coefficient priors or residual rescaling, and do not replace nonlinear profiles or establish a unique amplitude-sign solution. All six implementation tests passed, and `git diff --check` was clean.

Artifacts: [calibration results](../../../data/vector_sector/odd_low_energy_calibration.json), [independent verification/local diagnostics](../../../data/vector_sector/odd_low_energy_calibration_verification.json), and [SymPy transcript](../../../data/vector_sector/odd_low_energy_calibration.md). Input hashes, all separate minima, pulls, nuisance values, representative gauges, conditioning, boundary checks and nested measured-prediction changes are retained.

# Remaining exit gates and next step

The conditional stage-D diagnostic run and subsequent [finite sign/profile continuation](odd-low-energy-profiles.md) are executed. The continuation improves 80 preceding minima and certifies 112 physical profiles; 26 capped searches continue to improve and 30 one-sided source endpoints remain unresolved. The full stage-D validation/propagation exit criterion remains open. A poor achieved residual alone does not prove rejection while the nonlinear minimum remains unproven. No unique truncation is selected by this snapshot.

The [finite neutral target/null audit](odd-alp-measured-null-audit.md) now projects fixed-nuisance odd-coefficient nulls into the available resonance target families in 448 cases. Correction-model nulls are target-sensitive in every scenario. Resolve remaining sign/source cases and compute the joint nuisance/target nulls, with missing charged/local/physical-predictor completion, before any physical ALP band. Neither saturated interpolation nor a representative zero-null coefficient gauge determines unmeasured target directions.
