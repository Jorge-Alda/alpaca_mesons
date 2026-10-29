---
title: KK pi isovector spectrum input and extraction readiness
status: mixed
created: 2026-10-29
updated: 2026-10-29
bibtex_keys: [BaBar:2007ceh, GomezDumm:2010ra, Alda:2026nrd]
provenance: assumption
---

# Candidate measured input

[@BaBar:2007ceh] reports an isovector projection of $e^+e^-$ production through $K^*(892)K$, extracted from the Dalitz distribution for $K^0_S K^\pm\pi^\mp$. Nineteen Table 7 bins from $1.40$ to $1.78\,\mathrm{GeV}$ overlap the tau kinematic range, with the last published bin straddling its endpoint. We have transcribed their cross sections and statistical errors to the [CSV](../../../data/vector_sector/babar_2007_isovector_kstar_k_tau_region.csv); the paired [JSON metadata](../../../data/vector_sector/babar_2007_isovector_kstar_k_tau_region.json) records units, table and section provenance, published bin interpretation, systematics, and the endpoint qualification. The [deterministic extraction script](../../../scripts/extract_babar_kkpi_isovector_data.py) regenerates these files.

The paper assigns a 15% systematic uncertainty to the Dalitz fit, in addition to Table 1's 5.1% total cross-section systematic. It does not provide a bin-to-bin covariance for either source. The metadata records the individual Table 1 components and keeps their correlations unspecified. It also records a common multiplicative nuisance for the 15% uncertainty only as an optional fully correlated approximation; this is not a published covariance. The +4.2% correction listed in Table 1 has already been applied to the reported cross section.

# What the data can constrain

The isovector projection is a promising energy-dependent input for the $KK\pi$ vector-current coefficient map in the [rank audit](kkpi-vector-current-coefficient-rank.md). The spectrum provides substantially more leverage than an integrated branching fraction. It is a Dalitz-fit component for $K^*(892)K$, not an inclusive $KK\pi$ rate. A numerical conversion to the charged-current tau mode requires an explicit CVC/isospin normalization, charge-state conventions, and consistent treatment of the $K^*$ decay fractions.

This table alone does not yield a coefficient determination. The low-energy single-multiplet form-factor expressions in [@GomezDumm:2010ra] do not include an explicit excited isovector vector contribution. BaBar describes the isovector spectrum with a broad $\rho'$ and notes the single-resonance description may be insufficient. The current repository has the exact coefficient-function rank audit but no physical, finite-width $KK\pi$ likelihood predictor, CVC implementation, or fitted energy-dependent higher-resonance term. A fit that omits these ingredients would risk assigning the $ρ'$ shape to the low-energy missing coefficient, especially the VJP combination $c_1-c_2+c_5$.

# Extraction gate

Before treating this as coefficient information, construct the charged-current or neutral-current observable from the canonical amplitudes in the [overview](alp-vector-channel-coefficient-overview.md), integrate over the measured bins, and include the $K^*$ line shape and excited-vector terms. Map the isovector normalization with a cited CVC convention. Then fit jointly with the preliminary coefficient combinations while retaining its unconstrained directions and profile the normalization, resonance, and axial contributions. Keep the Table 1 and Dalitz-fit systematics distinct until a covariance prescription is justified. The first result should be a nuisance-profiled sensitivity/rank report; a coefficient profile is meaningful only if it survives reasonable excited-vector and covariance variants.

The present step establishes a source-backed data input and a concrete path to the fit. It does not report a fitted coefficient or likelihood rank.
