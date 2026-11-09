---
title: Excited-pole and a2(1320) pole shifts of the joint vector likelihood and their linear screen
status: mixed
created: 2026-11-07
updated: 2026-11-07
bibtex_keys: [ParticleDataGroup:2024cfk, Qin:2024ulb, Alda:2026nrd, BaBar:2007ceh]
provenance: numerical-computation
---

# Excited-pole and $a_2(1320)$ pole shifts

Decision: user, 2026-11-07, recorded in `openspec/changes/assess-joint-vector-fit-adequacy/design.md`. These shifts test how much the joint fit depends on the fixed pole inputs. They do not replace the baseline in the adequacy verdict.

## Conventions and units

Masses, widths and shifts are in GeV; invariants are in GeV$^2$. Notation and conventions follow Alda:2026nrd and the [baseline likelihood](babar-kkpi-global-vector-sector-fit-baseline.md). The excited-multiplet propagators are those of Qin:2024ulb as implemented in `fit_babar_kkpi_g_coefficients.propagator`.

## Inputs

- PDG 2024 uncertainties: [pdg2024-excited-vector-and-tensor-poles.md](pdg2024-excited-vector-and-tensor-poles.md), read at run time from `data/vector_sector/pdg2024_excited_poles.json`. None are hard-coded in the scripts.
- Baseline (hard-coded) poles, i.e. the centres of the shifts: $\rho(1450)$ 1.465/0.400, $\rho(1700)$ 1.720/0.250, $K^*(1410)$ 1.414/0.232, $K^*(1680)$ 1.718/0.322, $\omega(1420)$ 1.410/0.290, $\omega(1650)$ 1.670/0.315, $\phi(1680)$ 1.680/0.150, $\phi(2170)$ 2.162/0.100, $a_2(1320)$ 1.3182/0.107 (mass/width in GeV). A unit test confirms that these values reproduce the propagator.
- Endpoint and refit settings: [block-exclusion and variation refits](joint-vector-block-exclusion-and-variation-refits.md).

## Method

1. Each of the nine states has its mass and its width shifted on their own by $+1\sigma$ and by $-1\sigma$ around the baseline value, giving 36 variations. $+$ uses the upper and $-$ the lower PDG uncertainty. The $\phi(2170)$ width uses $+0.024/-0.018$ GeV as quoted, without symmetrising. The $a_2$ width uses the OUR ESTIMATE $\pm0.005$ GeV.
2. Excited vector poles: the inherited BaBar amplitude is evaluated with the reference propagator replaced by $P(s)+\beta_k\,[P_k^{\rm shifted}(s)-P_k^{\rm baseline}(s)]$. Here $\beta_k$ is the excitation coefficient of the shifted slot (first or second excited multiplet) and $P_k=1/(M^2-s-iM\Gamma(s))$ keeps the running-width shape of the propagator. For $\rho$ states this is the P-wave $\pi\pi$ phase space times $\sqrt s/M$; for $K^*$ states the P-wave $K\pi$ phase space times $M^4/s^2$; for $\omega$ and $\phi$ states the width is constant. The replacement only happens while the BaBar block is computed.
3. $a_2(1320)$: in the $K_SK\pi$ amplitude the term $(t-u)[1/(M^2-s-iM\Gamma)-1/(M_0^2-s-iM_0\Gamma_0)](-2g_TC_T/F^3+\sqrt2 g_TD_T\,\mathcal S_\rho/F^3)$ is added, and the prediction and residual of that channel are recomputed. This is the exact replacement of the pole in the tensor term.
4. The fit modules `fit_babar_kkpi_g_coefficients.py` and `fit_babar_kkpi_global_vector_sector.py` are not modified. Their SHA-256 hashes are checked by the saved ledgers.
5. Linear screen (`scripts/screen_joint_vector_pole_shifts.py`): at the fixed best endpoint `basin_probes/even_shared_widths`, at quadrature (2,16,16) with four-bin blocks, with baseline residual $r_0$, shifted residual $r$ and baseline Jacobian $J$ (32 coordinates):
   $$\Delta\chi^2_{\rm fixed}=|r|^2-|r_0|^2,\qquad \Delta\chi^2_{\rm lin}=\min_d|r+Jd|^2-\min_d|r_0+Jd|^2 .$$
   **Approximation:** the screen is linear in the parameters and in the shift, and ignores the bounds. It is used only to order the refits and as a cross-check, never to drop a variation.
6. Refits: the same settings, high-order (4,128,128) re-evaluation, acceptance gate and basin-switch reruns as the other variations. Every start endpoint runs as its own task (`--parallel-runs`). The run parts are written to the same parts directory, and the final ledger `data/vector_sector/joint_vector_sensitivity_variations.json` holds all 47 variations.

## Outputs

- `data/vector_sector/joint_vector_pole_shift_screen.json`: both $\Delta\chi^2$ values per shift, the linear maximum parameter shift in baseline local errors, and the refit order. The fit-order $\chi^2$ at the endpoint is 225.667, with linear minimum 225.577.
- Largest linear-refit predictions: $\phi(1680)$ width $-1\sigma$ (+54.8), $\omega(1420)$ width $-1\sigma$ (+29.4), $\phi(1680)$ mass $+1\sigma$ (+15.0), $\phi(1680)$ mass $-1\sigma$ (+8.2). All $K^*$, $\phi(2170)$ and $a_2$ shifts give $|\Delta\chi^2_{\rm lin}|<0.15$. Fixed-parameter values are much larger, e.g. 1053.6 for the $\phi(1680)$ width $-1\sigma$, so refitting absorbs most of the change.
- Refit results: see the variation ledger and [the refit record](joint-vector-block-exclusion-and-variation-refits.md).

## Verification

- `scripts/test_joint_vector_likelihood_options.py`: defaults and explicit `pole_shift=None` bitwise equal to the baseline (residual and Jacobian at (2,16,16)); the single-pole propagator reproduces the reference propagator difference for all eight states (relative $10^{-9}$); zero shifts reproduce the baseline residual ($10^{-10}$); shifts change only BaBar rows, and the $a_2$ shift only $K_SK\pi$ rows; the $a_2$ shift has no effect when $g_TC_T=g_TD_T=0$; the shifted Jacobian is finite; the 36 declared shifts read the PDG file (e.g. $+0.024/-0.018$ for the $\phi(2170)$ width).
- `scripts/test_screen_joint_vector_pole_shifts.py`: the linear minimum equals the projected residual; the stored screen has 36 entries in decreasing order.
- `scripts/test_refit_joint_vector_sensitivity.py`: 47 unique variations; a pole-shift variation refits at cheap orders; run selection rule.

Invocation:
```
.venv/bin/python scripts/screen_joint_vector_pole_shifts.py --output data/vector_sector/joint_vector_pole_shift_screen.json
.venv/bin/python scripts/refit_joint_vector_sensitivity.py --study variations --output data/vector_sector/joint_vector_sensitivity_variations.json --workers 11 --parallel-runs [--only NAME ...]
```
