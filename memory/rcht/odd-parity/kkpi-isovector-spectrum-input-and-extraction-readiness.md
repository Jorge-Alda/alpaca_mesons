---
title: KK pi isovector spectrum input and extraction readiness
status: mixed
created: 2026-10-29
updated: 2026-10-30
bibtex_keys: [BaBar:2007ceh, GomezDumm:2010ra, Alda:2026nrd]
provenance: assumption
---

# Candidate measured input

[@BaBar:2007ceh] reports an isovector projection of $e^+e^-$ production through $K^*(892)K$, extracted from the Dalitz distribution for $K^0_S K^\pm\pi^\mp$. Nineteen Table 7 bins from $1.40$ to $1.78\,\mathrm{GeV}$ overlap the tau kinematic range, with the last published bin straddling its endpoint. We have transcribed their cross sections and statistical errors to the [CSV](../../../data/vector_sector/babar_2007_isovector_kstar_k_tau_region.csv); the paired [JSON metadata](../../../data/vector_sector/babar_2007_isovector_kstar_k_tau_region.json) records units, table and section provenance, published bin interpretation, systematics, and the endpoint qualification. The [deterministic extraction script](../../../scripts/extract_babar_kkpi_isovector_data.py) regenerates these files.

The paper assigns a 15% systematic uncertainty to the Dalitz fit, in addition to Table 1's 5.1% total cross-section systematic. It does not provide a bin-to-bin covariance for either source. The metadata records the individual Table 1 components and keeps their correlations unspecified. It also records a common multiplicative nuisance for the 15% uncertainty only as an optional fully correlated approximation; this is not a published covariance. The +4.2% correction listed in Table 1 has already been applied to the reported cross section.

# What the data can constrain

The isovector projection is a promising energy-dependent input for the $KK\pi$ vector-current coefficient map in the [rank audit](kkpi-vector-current-coefficient-rank.md). The spectrum provides substantially more leverage than an integrated branching fraction. It is a Dalitz-fit component for $K^*(892)K$, not an inclusive $KK\pi$ rate. A numerical conversion to the charged-current tau mode requires an explicit CVC/isospin normalization, charge-state conventions, and consistent treatment of the $K^*$ decay fractions.

The exact coefficient-function rank audit identified the additional VJP direction $c_1-c_2+c_5$. A direct crossed $K^*\to K\gamma^*$ quasi-two-body fit now profiles that direction and $c_4$; its coefficient results, profiles, resonance-model sensitivity, and retained nulls are recorded in [the coefficient-fit record](babar-kkpi-rcht-coefficient-fit.md). The full inclusive weak-current predictor below remains separate. BaBar reports a broad $\rho'$ and notes that a single-resonance description may be insufficient, so the coefficient profiles retain excited-pole and fit-window sensitivities.

# Remaining combined-fit gate

The direct crossed $K^*K$ coefficient fits are complete as reported in
[babar-kkpi-rcht-coefficient-fit.md](babar-kkpi-rcht-coefficient-fit.md).
The earlier proposed inclusive route remains separate: it requires a
charge-resolved current predictor, the experimental Dalitz projection, and
a cited CVC map before combining those rates with tau data. The current
quasi-two-body coefficient profiles retain the measured statistical errors,
the 15% Dalitz systematic, the 5.1% Table 1 uncertainty, and the explicit
excited-pole/window variants; these remain conditional because the published
bin covariance and a common excited-vector/RChT fit are unavailable.

The data input and extraction path were established on 2026-10-29. The subsequent coefficient fit is documented in [babar-kkpi-rcht-coefficient-fit.md](babar-kkpi-rcht-coefficient-fit.md); it constrains selected coefficient coordinates conditionally and leaves the stated coefficient null directions unresolved.
