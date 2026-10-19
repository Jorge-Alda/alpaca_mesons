---
title: Finite sign-basin and physical-amplitude profile validation
status: mixed
created: 2026-10-14
updated: 2026-10-19
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# Scope and frozen finite policy

This continues the [conditional calibration](odd-low-energy-calibration.md) and [maintained stage-D plan](odd-breaking-low-energy-calibration-plan.md). All 64 input/mass/source alternatives and four predeclared datasets remain separate. The diagonal PDG-summary covariance, fixed Alda central inputs, LO source-compatible vertex masses, physical phase-space masses and single counting of ancillary observations are inherited assumptions. No new prior or OPE constraint is introduced.

Strengthen the sign-basin search for all 512 unsaturated original/soft-leading cases. Pair only their training-observable likelihood evaluations while bidirectional image containment is certified at every proposed nuisance point. Pairing is not an assertion of equal action nullspaces, held-out predictions or ALP sensitivities. Each of the 256 scenario/dataset likelihood pairs uses sixteen distinct fixed-point sign starts and both previous optima. If the achieved residual improves by more than $10^{-4}$, run one further deterministic batch of sixteen starts; stop at that declared budget and label continuing improvement. This is a stronger finite search, not an exhaustive proof of a global nonlinear minimum.

Prioritize interior post-search cases with descriptive nominal $p\geq0.05$ and maximum odd-channel $|\mathrm{pull}|\leq3$ for profiles. This is a computational allocation rule, not an acceptance, exclusion or preferred-model criterion. Preserve the coverage/status of all other cases. Source-boundary approaches receive one-sided source-residue profiles rather than symmetric Gaussian errors.

# Physical profile coordinates

Choose two most-uncertain standardized training-amplitude contrasts from the full joint null-safe local amplitude covariance at the reference solution. Reference amplitude uncertainties are standardized using the fixed measurement errors and observable derivative. Freeze the contrast weights, standardization and reference amplitudes; regenerate the physical amplitude map as nuisances move. Require each contrast gradient to annihilate the complete retained-coefficient/nuisance training nullspace.

The contrast is linear in odd coefficients at fixed nuisances, so a nuisance-dependent constrained-coordinate elimination can enforce the same physical quantity throughout minimization. Raw coefficient SVD axes are not used to choose profile quantities. Constraints retain the exact inverse-action mass branch, source-current signs, positive total widths and open physical channels.

Each contrast uses seventeen targets from minus four to plus four reference standard deviations, in half-standard-deviation steps, with at most three deterministic starts per point. Extend a missing crossing only to six and then eight standard deviations. Bracket $\Delta\chi^2=1$ and $3.84$ as descriptive likelihood levels; unresolved crossings remain open. A newly lower reference is flagged unresolved and its interval brackets are withheld until re-anchoring. No such event occurred in this run. Failed optimizer/domain/image/constraint checks are unresolved points, not assigned large residuals.

One-sided source-residue targets use positive multiples $0.05,0.25,0.5,1,2,4,8$ of the reference magnitude, with an explicitly recorded numerical domain floor. Do not cross a hypothesis-defining sign or mass-root boundary to close an interval.

# Inputs and provenance

- [Conditional calibration](odd-low-energy-calibration.md): separate scenario optima, sign starts and frozen likelihood.
- [Measurement Jacobians](odd-low-energy-measurement-jacobians.md): retained parent filters, coefficient units and shared-width observations.
- [Bounded SM numerical kernels](odd-bounded-sm-numerical-kernels.md): normalized amplitudes and WZW offsets.
- [PDG snapshot](../../extracts/ParticleDataGroup:2024cfk/odd-parity-low-energy-measurements.md): measurements and absent joint covariance [@ParticleDataGroup:2024cfk].
- [Alda input conversion](alda-neutral-amplitude-inputs.md): fields, central spurions and source/physical mass distinction [@Alda:2026nrd].

The next structural gate is a separate retained-action null projection into available normalized ALP targets. Missing source/contact or physical-normalization kernels are explicit completion gates; this finite continuation does not derive all missing predictors or produce an ALP band.

# Executed result

The complete 256-pair run covers all 512 original/soft-leading model cases and performed **5888** finite sign-search fits. Eighty achieved minima improve the preceding calibration by more than $10^{-4}$; 26 pairs continue to improve in the capped second batch, so the finite search does not establish nonlinear global minima. All selected best fits converged and independently pass the inverse-mass-root, source-sign and training-image certificates. The best primary achieved $\chi^2=5.5675232799$ is unchanged; improved alternative/control minima are retained separately rather than used to select a preferred scenario.

Fifty-six computationally prioritized interior pairs receive two profiles each. All **1909** physical profile points satisfy the optimizer/domain/image and physical-contrast constraints; the maximum independently reconstructed normalized constraint defect is $1.39\times10^{-12}$ and the residual reproduction defect is $4.25\times10^{-16}$. No profile found a lower reference. All 112 contrasts bracket both sides of the descriptive $\Delta\chi^2=3.84$ level, with no open side in the declared grid/extensions. These finite, sign-conditional brackets do not establish coverage or a globally unique solution.

Six pairs whose preceding solution approached a source boundary receive one-sided source profiles. Twelve endpoints were initially valid and **30 initially remained unresolved**. A finite continuation ([script](../../../scripts/continue_odd_low_energy_profiles.py), [result](../../../data/vector_sector/odd_low_energy_continuation.json)) covered all 28 flagged pairs and performed 442 additional unconstrained fits. Six achieved minima improved by more than $10^{-4}$, at pairs 37, 45, 64, 90, 122 and 218; their $\chi^2$ values were 2.28602596, 0.365868788, 3.03653351, 348.344646, 439.847677 and 371.508355. The largest training-image defect was $9.20\times10^{-15}$. The source tolerance was not changed.

A direct physical-source-coordinate calculation with independent 80-digit Decimal reconstruction and a forward-error bound ([script](../../../scripts/continue_odd_low_energy_source_profiles.py), [result](../../../data/vector_sector/odd_low_energy_source_continuation.json)) certified all 29 previously unresolved endpoints across 220 candidate fits. The maximum Decimal constraint defect is $3.51\times10^{-17}$ GeV, while the maximum derived arithmetic bound is $8.90\times10^{-4}$ of target; all strict source signs are certified. The largest attempted training-image defect is $8.95\times10^{-11}$, below the inherited $10^{-8}$ limit.

Accepted source minima lowered the old references for pairs 90, 122, 218 and 250. Re-anchored unconstrained fits ([script](../../../scripts/reanchor_odd_low_energy_profiles.py), [result](../../../data/vector_sector/odd_low_energy_profile_reanchoring.json)) give $\chi^2=1308.568779$, 2533.137264, 1308.428487 and 2532.833120, respectively. A refreshed source grid for these four pairs ([script](../../../scripts/reprofile_odd_low_energy_reanchored.py), [result](../../../data/vector_sector/odd_low_energy_reanchored_profiles.json)) certified all 28 points; physical-amplitude profiles were omitted under the predeclared $\max|\mathrm{pull}|\leq3$ allocation rule. The new grids found lower references for pairs 90 and 218, so targeted unconstrained restarts and another source-grid pass were performed ([script](../../../scripts/restart_odd_low_energy_profile_references.py), [result](../../../data/vector_sector/odd_low_energy_targeted_profile_restart.json)). Final references are $\chi^2=1232.13594777$ and 1231.57510830. All 14 refreshed source points are valid and neither grid finds a decrease larger than $10^{-4}$; these sparse finite checks do not establish global minima or precise likelihood-level crossings. Pairs 122 and 250 remain flat across their declared grids. No symmetric Gaussian errors or closed source-domain interval is assigned.

# Reproduction and artifacts

```text
.venv/bin/python scripts/profile_odd_low_energy.py
.venv/bin/python scripts/continue_odd_low_energy_profiles.py
.venv/bin/python scripts/continue_odd_low_energy_source_profiles.py
.venv/bin/python scripts/reanchor_odd_low_energy_profiles.py
.venv/bin/python scripts/reprofile_odd_low_energy_reanchored.py
.venv/bin/python scripts/restart_odd_low_energy_profile_references.py
.venv/bin/python scripts/verify_odd_low_energy_profiles.py
.venv/bin/python scripts/plot_odd_low_energy_profiles.py
.venv/bin/python scripts/export_odd_profile_target_points.py
.venv/bin/python scripts/implementation/test_odd_low_energy_profiles.py
```

The executed sweep used four disjoint chunks, with `--start-pair` equal to 0, 64, 128 and 192 and `--limit-pairs 64`; `verify_odd_low_energy_profiles.py --merge-chunks` merges only matching input/script hashes and certifies exact coverage. Every completed pair has a checkpoint. The full single-process invocation above performs the same declared calculation. Four algebra tests verify covariance null rejection, coordinate rescaling, deterministic degenerate contrasts and exact constrained-coordinate elimination with moving nuisances.

- [Full profile results](../../../data/vector_sector/odd_low_energy_profiles.json) retain every separate sign start, best fit, profile and failed endpoint.
- [Independent verification](../../../data/vector_sector/odd_low_energy_profiles_verification.json) certifies hashes, coverage, branch/source domains, residuals and physical/source constraints.
- [Profile figure](../../../data/vector_sector/odd_low_energy_profiles.png) and [PDF](../../../data/vector_sector/odd_low_energy_profiles.pdf) show the best achieved profiled example in each control; they are illustrative conditional examples, not pooled model bands.
- [Primary target points](../../../data/vector_sector/odd_profile_primary_target_points.json) export improved nuisances with independently reconstructed original and soft-leading coefficient gauges. Training-image equivalence does not equate their action nullspaces or target predictions.
