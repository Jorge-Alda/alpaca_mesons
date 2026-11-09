---
title: Eta-prime radiative and mass-block exclusion refits of the joint vector fit
status: numerical
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Eta-prime radiative and mass-block exclusion refits

## Inputs and conventions

The accepted starting point `basin_probes/even_shared_widths` and the baseline pulls are in the [goodness-of-fit record](joint-vector-fit-goodness-of-fit.md). The [likelihood baseline](babar-kkpi-global-vector-sector-fit-baseline.md) specifies the measured rows, conventions and assumptions. The [refit procedure](joint-vector-block-exclusion-and-variation-refits.md) specifies the optimizer, bounds, fixed coordinate scales and quadrature orders. The previously computed [block-exclusion results](joint-vector-fit-block-exclusion-results.md) supply the mass-and-leptonic exclusion. Input measurements follow [@ParticleDataGroup:2024cfk] and the amplitudes follow [@Alda:2026nrd]. Masses and widths are in GeV, cross sections are in nb, and pulls and $\chi^2$ are dimensionless.

The fit assumes the baseline diagonal covariance for odd radiative rows and assumes zero unavailable cross-covariances. The $\eta'\to\rho^0\gamma$ measurement was already deferred in the baseline, so it is not a retained row. The two excluded $\eta'$ radiative observations are $\eta'\to\omega\gamma$ and $\eta'\to\gamma\gamma$; the $\eta'$ total-width measurement remains. The removed mass block contains four vector masses and three electronic widths. Per-block nominal $p$-values use the number of retained rows as degrees of freedom and do not subtract shared fitted parameters; they are a ranking convention, not calibrated probabilities.

## Method

1. Start at `basin_probes/even_shared_widths`. Remove the two retained $\eta'$ radiative rows and fit the remaining 73 residuals using scipy `least_squares`, quadrature $(2,16,16)$, and the settings of the [refit procedure](joint-vector-block-exclusion-and-variation-refits.md). Record the optimizer status, termination message, evaluation counts and Jacobian conditioning.
2. Re-evaluate the resulting parameters at $(4,128,128)$, including the excluded rows. Report retained odd and shared-width block $\chi^2$ and nominal $p$-values, the fitted $\Gamma_{\eta'}$, and each pull's change from the accepted endpoint. A stopped or ill-conditioned fit is labelled accordingly.
3. Read the pre-existing fit without the mass-and-leptonic block, which used the same start, settings and quadrature orders. Its ledger retains optimizer status and all seven removed-row pulls. Compare those pulls with the accepted endpoint and evaluate the $\eta'$ width from the stored parameters.

## Results

The machine-readable ledger is [`joint_vector_diagnostic_exclusion_refits.json`](../../../data/vector_sector/joint_vector_diagnostic_exclusion_refits.json); it retains all row changes.

Without the two $\eta'$ radiative rows the optimizer terminates successfully with `status=2` (`ftol`), after 424 evaluations. The fit-order Jacobian has no direction below the recorded singular-value threshold. The retained odd block has $\chi^2=9.685$ for 10 rows and nominal $p=0.469$; the shared-width block has $\chi^2=98.983$ for 7 rows and nominal $p=1.75\times10^{-18}$. The fitted $\Gamma_{\eta'}=0.2300$ MeV, matching its measured central value. Its pull changes from $-3.307$ to zero. The $\rho^0\to\eta\gamma$ pull moves from $-3.266$ to $+0.187$. Evaluated at the refitted point, the excluded $\eta'\to\gamma\gamma$ and $\eta'\to\omega\gamma$ pulls are $+21.828$ and $+5.111$, respectively, compared with $+2.455$ and $-0.723$ at the baseline endpoint. The retained total is $\chi^2=183.415$ for 41 nominal degrees of freedom.

For the stored mass-and-leptonic exclusion, the optimizer also terminated successfully with `status=2` (`ftol`), after 640 evaluations; its fit-order Jacobian has no direction below the recorded threshold. The retained odd block has $\chi^2=0.400$ for 12 rows and nominal $p=0.999999925$; the shared-width block has $\chi^2=22.225$ for 7 rows and nominal $p=0.00232$. The fitted $\Gamma_{\eta'}=0.22794$ MeV. Removed mass pulls change from $(+0.827,-3.198,+0.402,+0.100)$ to $(+84.160,+8.291,+408.407,-1158.513)$ for $\rho^+$, $K^{*+}$, $\omega$ and $\phi$, respectively. These are fixed-endpoint comparisons; neither exclusion establishes a global minimum or an acceptable fit including the excluded measurements.

## Reproduction and verification

```text
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/diagnose_joint_vector_exclusion_refits.py --output data/vector_sector/joint_vector_diagnostic_exclusion_refits.json
.venv/bin/python -m unittest scripts.test_diagnose_joint_vector_exclusion_refits
```

The [driver](../../../scripts/diagnose_joint_vector_exclusion_refits.py) records the fit and evaluation orders, input hashes, all row changes, and full optimizer records. Its checkpoint permits continuation of the $\eta'$ fit. The tests verify the excluded channels and retained $\eta'$ width, signed row changes, the stored mass-block optimizer status and a reduced refit's optimizer status and retained block counts. The mass-block ledger was produced using earlier versions of the option and refit scripts; as documented in the [block-exclusion results](joint-vector-fit-block-exclusion-results.md), the later changes preserve default residuals bitwise.
