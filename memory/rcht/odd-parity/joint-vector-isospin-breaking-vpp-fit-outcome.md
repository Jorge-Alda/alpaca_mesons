---
title: Outcome of the joint vector fit with the isospin-breaking VPP coordinates (task 3.2)
status: numerical
created: 2026-11-09
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh]
provenance: numerical-computation
---

# Outcome of the joint vector fit with the isospin-breaking VPP coordinates

This record gives the result of task 3.2 of `extend-joint-vector-model-diagnostics`: the joint likelihood with the switchable isospin-breaking $VPP$ coordinates $r_1,r_2,r_3$ ([option record](joint-vector-isospin-breaking-vpp-likelihood-option.md), [s-channel factors](joint-vector-isospin-breaking-vpp-s-channel-factors.md)), fitted with the gates of the [fit plan](joint-vector-isospin-breaking-vpp-fit-plan.md). **The adequacy criterion is not met.** The extension is therefore recorded as inadequate, and the plan continues with section 4.

## Ledgers

All ledgers are in `data/vector_sector/`. Each records its command, optimizer settings, source hashes and every run, including unsuccessful ones.

| ledger | study | what it contains |
| --- | --- | --- |
| `joint_vector_isospin_breaking_vpp_fit.json` | scan point and zero control | four refits from the scan point $r=(-0.350,1.855,0.729)$ GeV$^{-2}$ with $r$ free, all hitting the 5000-evaluation cap; control with $r=0$, accepted |
| `joint_vector_isospin_breaking_vpp_restart_fit.json` | restart from the best scan run | one refit, `ftol` satisfied after 160 evaluations, accepted |
| `joint_vector_isospin_breaking_vpp_restart_tight.json` | restart with `ftol`, `xtol` $=10^{-12}$, `gtol` $=10^{-10}$ | hit the 5000-evaluation cap, not accepted |
| `joint_vector_isospin_breaking_vpp_r1_wide_profile.json` | $r_1$ profile over $[-0.45,-0.18]$ | four points, one unconverged |
| `joint_vector_isospin_breaking_vpp_r_profile.json` | $r_1$ profile over $[-0.325,-0.315]$ | five points, two unconverged |
| `joint_vector_isospin_breaking_vpp_r2_profile.json` | $r_2$ profile over $[1.40,1.80]$ | four points, all converged |
| `joint_vector_isospin_breaking_vpp_r3_profile.json` | $r_3$ profile over $[0.25,0.45]$ | four points, one unconverged |

## Best fit

The four scan-point refits and the restart agree on the same region:
- $\chi^2_{(4,128,128)}=98.10$ for 40 nominal degrees of freedom, so the total nominal $p=8.7\times10^{-7}$;
- $r=(-0.32,\,1.59,\,0.35)$ GeV$^{-2}$;
- $\chi^2$ by block at $(4,128,128)$: odd $18.8$ (12 rows), masses and leptonic widths $0.8$ (7), shared widths $10.5$ (7), $\phi\to K\bar K$ $2.5$ (2), BaBar $K\bar K\pi$ $65.4$ (45), normalization $0.2$ (2).

The zero-coefficient control reproduces the baseline endpoint: $\chi^2=225.67$ against the baseline $225.73$, within the control tolerance.

The baseline block values are $106.8$ for the shared widths and $17.9$ for $\phi\to K\bar K$. The option therefore removes most of the width and $\phi$ tension, and the BaBar block rises from $52.8$ to $65.4$.

## Profiles

Each profile fixes one $r_i$ and refits the other 34 coordinates from the restart endpoint.

| coordinate | range | $\Delta\chi^2$ over range | BaBar block | reading |
| --- | --- | --- | --- | --- |
| $r_1$ | $-0.38$ to $-0.18$ | $0.2$ to $0.6$ | $64$–$67$ | weakly constrained; rough nominal 1σ $\approx\pm0.1$ |
| $r_2$ | $1.40$ to $1.80$ | $0.05$ to $0.2$ | $65$–$66$ | weakly constrained; rough nominal 1σ $\approx\pm0.4$ |
| $r_3$ | $0.25$ to $0.40$ | $0.13$ to $0.43$ | $65$–$66$ | weakly constrained; rough nominal 1σ $\approx\pm0.13$ |

The $\Delta\chi^2$ values are relative to $98.10$. The 1σ estimates come from a quadratic reading of three or four converged points, with convergence noise of about $0.1$ to $0.2$ in $\chi^2$. They are rough, not profile intervals. Points that did not converge are excluded: $r_1=-0.45$ (cap), $r_1=-0.3225$ and $-0.3150$ (cap, but within 0.1 of the others), and $r_3=0.45$ (cap, optimality $46$). For $r_1=-0.3225$ and $-0.3150$ the cap was reached, so their values are upper bounds on the true profile.

## Adequacy

The criterion of [assess-joint-vector-fit-adequacy](joint-vector-fit-adequacy-assessment-method.md) is total nominal $p>0.01$ with no block below $10^{-3}$. Every converged point gives total $p\approx10^{-6}$. No single block fails $10^{-3}$, but the per-block test uses row counts as degrees of freedom and does not see the total. **The criterion is not met.**

## Interpretation

- The option removes most of the shared-width and $\phi\to K\bar K$ tension. These blocks reach $\chi^2$ close to their row counts.
- The remaining $\chi^2$ floor of about $98$ sits in the BaBar block (about $65$) and in the $\eta'$ total-width pull (about $-3.1$). Neither changes measurably over any of the three profiles.
- The $r_i$ are weakly constrained at the level of the profile noise, so they do not explain the BaBar residual.

## Assumptions and limits

- The optimizer terminates on `ftol`, not on a zero gradient, and several points stop at optimality $3$–$8$. The $\chi^2$ values are therefore upper bounds, with noise of order $0.1$ to $0.2$.
- The $r_i$ profile curvatures are rough nominal estimates. They are not profile likelihood intervals.
- The excited-$K^*$ scaling is `ground_only`. The `whole_propagator` alternative was not fitted.
- The option has the $s$-channel ground-state factors only; excited $\rho$ and $\omega$ $s$-channel terms are not included.
- Degrees of freedom are $75-35=40$ with the three added coordinates, as in the [fit plan](joint-vector-isospin-breaking-vpp-fit-plan.md).

## Consequence for the plan

The change plan (task 3.2) asks for the ledger entry and the criterion check. Both are done, and the criterion fails, so section 3 does not close the inadequacy. The next step in the plan is section 4, the double-trace source term, starting with its derivation in task 4.1. The choice of excited-$K^*$ scaling and the BaBar residual pattern remain open for later checks.
