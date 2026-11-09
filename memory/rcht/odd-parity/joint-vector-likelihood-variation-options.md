---
title: Switchable variations of the joint vector-sector likelihood
status: mixed
created: 2026-11-06
updated: 2026-11-07
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: assumption
---

# Switchable variations of the joint vector-sector likelihood

## Design decision

The variations are options of `VariedLikelihood` in [`scripts/joint_vector_likelihood_options.py`](../../../scripts/joint_vector_likelihood_options.py), a subclass of `GlobalLikelihood`. The baseline script is left byte-identical because the saved fit ledgers and [`compare_joint_vector_endpoints.py`](../../../scripts/compare_joint_vector_endpoints.py) verify its SHA-256. With every option at its default, the subclass calls the inherited methods unchanged; a unit test checks bitwise-equal residuals and Jacobians at the accepted endpoint `basin_probes/even_shared_widths`.

## Options (defaults reproduce the [baseline](babar-kkpi-global-vector-sector-fit-baseline.md))

| option | default | variation |
| --- | --- | --- |
| BaBar upper energy | 1.78 GeV | bins with upper edge above the cut are removed from the residual (8 bins at 1.70, 18 at 1.60) |
| regulator $\Lambda$ | 2.5 GeV | amplitude factor $e^{-E^6/\Lambda^6}$ multiplying every bin integrand squared; implemented by rescaling the bin weights by $e^{-2E^6(\Lambda^{-6}-2.5^{-6})}$ |
| source hypothesis | negative single trace | universal ($z_V=0$, domain $F_V>0$), positive single trace (light and strange residues positive) |
| inverse-mass branch | near ideal | second root; $(\Delta_V,\kappa_V)$ recentred on the stage-3 root of the PDG charged $\rho$/$K^*$ scenario |
| $\phi\to K\bar K$ correlation | $-0.72$ | 0 |
| $\phi$ electronic width | PDG OUR EVALUATION $1.27\pm0.04$ keV | PDG OUR AVERAGE $1.251\pm0.021$ keV (`pdg2024.json`, `alternatives.phi_ee_average`) |
| BaBar fractional systematic | 0 | $\sigma_i^2\to\sigma_i^2+(f\,\sigma^{\rm obs}_i)^2$, uncorrelated between bins |
| excluded blocks | none | any of the six blocks removed |
| pole shift | none | mass or width of one fixed excited vector pole or of the $a_2(1320)$ moved by a given amount; see [the pole-shift record](joint-vector-pole-shift-variations.md) |

Tensor products fixed to zero are a fixed-coordinate choice of the refit driver, not a likelihood option.

**Assumption:** the 5% bin-uncorrelated BaBar systematic has no source in the inclusive tables used here (they carry statistical errors only; normalization uncertainties of 6.5% and 5.1% are separate pulls). Its size is fixed a priori and is labelled as an assumption.

The source hypotheses and their domains follow [the stage-3 calibration](../foundations/vector-stage3-calibration.md) and [the stage-5 specification](vector-stage5-joint-fit-specification.md). For the universal source the BaBar current is the single-trace current at $z_V=0$, which is a special case of the implemented formula.

## Not implemented

- Single-plus-double-trace source: deferred by user decision on 2026-11-07. The double-trace source term $y\,\mathrm{tr}(\chi_+Q)\mathbb{1}$ has not been derived for the inclusive $KK\pi$ current, and the likelihood has no $y$ coordinate. It is a candidate model extension for the follow-up bead.

The excited-pole and tensor-pole variations, listed here as missing on 2026-11-06, are implemented since 2026-11-07 ([pole shifts](joint-vector-pole-shift-variations.md)).

## Verification

[`scripts/test_joint_vector_likelihood_options.py`](../../../scripts/test_joint_vector_likelihood_options.py): bitwise equality at defaults; each option changes only the intended rows; masks give the expected row counts; source and branch starts lie in their domains; explicit `pole_shift=None` is bitwise equal to the baseline; pole-shift checks listed in [the pole-shift record](joint-vector-pole-shift-variations.md).
