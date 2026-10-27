---
title: BaBar 2021 three-pion cross section and covariance inputs
status: literature
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [BABAR:2021cde]
provenance: literature
---

# Cross-section data and uncertainties

[@BABAR:2021cde] measures $e^+e^-\to\pi^+\pi^-\pi^0$ from 0.62 to
3.5 GeV. Table 4 covers energies below 1.1 GeV and reports bin edges, the
cross section, the statistical error, a correlated systematic error, and a
residual systematic error. The paper links the statistical covariance matrix
`rec_cov_cr.dat` and residual systematic covariance matrix
`matrix_cov_syst_sum.dat` from its unfolding repository. The correlated
systematic includes luminosity, radiative correction, detection efficiency,
and unfolding; it is 1.3% in the 0.73--1.03 GeV region.

Table 5 covers energies above 1.1 GeV and reports statistical and systematic
errors. The systematic includes luminosity, radiative correction, detection
efficiency, and background subtraction. For its integrated $a_\mu$ result,
the paper conservatively treats these systematic errors as fully correlated;
that prescription should not automatically be assumed for a new cross-section
fit without documenting it.

The source and covariance files are archived under
[the source index](../../sources/BABAR:2021cde/INDEX.md). A reproducible
extraction produces [separate low/high tables](../../../data/odd_parity/babar_2021_3pi_below_1p1.csv)
and a [checksum manifest](../../../data/odd_parity/babar_2021_3pi_provenance.json)
with [the extraction script](../../../scripts/extract_babar_2021_three_pion.py).
The extraction does not combine error components or impose an unstated
cross-region covariance.
