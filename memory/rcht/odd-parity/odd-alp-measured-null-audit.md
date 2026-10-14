---
title: Finite neutral ALP resonance target projection of measured odd nullspaces
status: mixed
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# Inputs and domain

Use the [sign/profile continuation](odd-low-energy-profiles.md), [conditional calibration](odd-low-energy-calibration.md), [bounded gradients](odd-bounded-sm-numerical-kernels.md), [source-EOM/contact maps](odd-power-counted-channel-maps.md), and [LO source/LSZ assembly](odd-alp-source-lsz-assembly.md). Literature inputs are already cataloged in the bibliography and linked records. The photon contact uses the source-backed [Bauer coefficient normalization](../../extracts/Bauer:2020jbp/photon-operator-normalization.md) and has zero derivative with respect to fitted odd coefficients; independent hadronic anomaly/local coefficients are uncomputed.

This is a **finite neutral resonance odd-coefficient audit at fixed fitted nuisances**, not the complete joint coefficient/nuisance audit or complete physical predictors. Use primary-data points in all 64 scenarios and seven separate truncations. Improved primary nuisances replace the original/soft-leading points; other models retain preceding calibrated points. Separate coefficient gauges reproduce primary training amplitudes without equating action nullspaces or target predictions.

Probe $m_a^2=(0.01,0.16,1.44)$ GeV squared and the four unit UV coordinates $(c_u^{AP},c_d^{AP},c_s^{AP},c_{GA})$. Set $f_a=1$ GeV only as an overall rank normalization. Exclude neutral poles. Formal unit probes do not establish perturbative LSZ validity for a physical coupling, generic all-mass ranks, widths, thresholds, or unstable-vector line shapes. Daughter-vector projections are on-shell algebraic projections, including below physical decay thresholds.

# Method

The adapter evaluates the previously derived LO diagonal resolvent and physical $U,Y,Z$ sources. It checks $Y+2m_a^2U-2ZI=0$ and annihilation of **all 13 action relations** by each target row. It includes signed neutral-vector photon response, both photon assignments, canonical external-vector factors and inherited source/action phases. Five tests cover SM recovery, UV linearity/$f_a$ scaling, pole rejection, retained vacuum $3P$ zero, action relations and null-witness reconstruction.

For parent embedding $E$, let $A=M_{\rm SM}E$ and $B=M_{\rm target}E$. Column-balanced SVD gives a basis $N$ of $\ker A$, at tolerance $10^{-10}$. Compute $BN$ separately for each target family. Every nonzero image exports a coefficient witness $v$ with $Av=0$ and $Bv\ne0$. Nonzero fitted measured amplitudes ensure fixed-nuisance amplitude and width rows have the same odd-coefficient nullspace. A SymPy example derives the nullspace and target image in the executable artifact.

# Executed result

All **448** scenario/model cases pass coverage, input-hash, action-relation and null-witness checks. All 13 source/contact relations pass at every declared mass/UV/sector point; the largest normalized residual is $5.31\times10^{-16}$. The largest normalized measured-null witness defect is $5.57\times10^{-17}$. The independent verifier checks stacked-rank arithmetic. These sampled family null-image ranks are identical across the 64 scenarios:

| Retained model | $aVV$ | $aV\gamma$ | $a\gamma\gamma$ | Single odd vacuum $a3P$ |
|---|---:|---:|---:|---:|
| Original eleven | 0 | 0 | 0 | 0 |
| Soft leading | 0 | 0 | 0 | 0 |
| Leading plus double trace | 5 | 5 | 1 | 0 |
| Leading plus single-trace corrections | 3 | 4 | 2 | 0 |
| Soft first corrections | 13 | 13 | 3 | 0 |
| Promoted leading | 0 | 0 | 0 | 0 |
| Full bounded action | 28 | 21 | 4 | 0 |

The table pools declared UV/mass rows **within** each target family; it is not a rank for every fixed UV model or individual channel. Nonzero witnesses establish target-sensitive unmeasured odd directions at these probes. A representative zero-null coefficient gauge cannot remove their effects. Low-energy data therefore do not determine correction-model resonance target bands without additional information or explicit assumptions. No prior is introduced here.

Zero images for original, soft-leading and promoted-leading models mean odd-coefficient blindness is not a block **at this fixed-nuisance, neutral, retained tree projection and finite sampling accuracy**. They do not establish a nuisance-profiled band or all-orders statement. Reassess with operator-basis, off-shell/charged, loop, source/local-matching or precision changes. The $3P$ zero applies only to one odd insertion in vacuum at retained tree accuracy; it does not set the total ALP three-meson amplitude to zero.

# Reproduction and remaining gates

```text
.venv/bin/python scripts/export_odd_profile_target_points.py
.venv/bin/python scripts/audit_odd_alp_target_nulls.py --profile-points data/vector_sector/odd_profile_primary_target_points.json
.venv/bin/python scripts/verify_odd_alp_target_nulls.py
.venv/bin/python scripts/implementation/test_odd_alp_target_null.py
```

- [Audit artifact](../../../data/vector_sector/odd_alp_target_nulls.json) retains labeled probe identities, ranks, witnesses, source/contact certificates and hashes.
- [Independent verification](../../../data/vector_sector/odd_alp_target_null_verification.json) records coverage and per-model summaries.
- [Adapter](../../../scripts/implementation/odd_alp_target_null.py) uses saved normalized parents; it is not a complete physical predictor.

The finite continuation is executed. Full stage-D validation/propagation remains open: resolve continuing sign-search and source-boundary cases, then audit the joint coefficient/nuisance target tangent. Charged targets, independent hadronic local/anomaly matching, physical normalization/domain completion and the even-parity total $3P$ predictor remain explicit gates before physical ALP bands.
