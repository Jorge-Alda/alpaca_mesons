---
title: Extracted BaBar 2021 three-pion fit data and covariance
status: numerical
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [BABAR:2021cde]
provenance: numerical-computation
---

# Inputs and extraction

The executable [extractor](../../../scripts/extract_babar_2021_three_pion.py)
parses the source tables in [@BABAR:2021cde] and the two covariance matrices
linked by that paper. It writes separate low- and high-energy tables plus a
checksum manifest; the source-specific error schemas are not merged.

# Extracted contents

- Table 4 provides 114 bins from 0.62 to 1.10 GeV. Its CSV keeps the
  statistical error, correlated systematic-error vector, and residual
  systematic-error vector separately. The 114 by 114 statistical and
  residual-systematic covariance matrices are supplied as separate files.
- Table 5 provides 72 bins from 1.10 to 3.50 GeV, with statistical and total
  systematic errors kept separately. The source treats systematic errors as
  fully correlated for its $a_\mu$ integral over 1.1--2.0 GeV; this is not
  imposed as a fit covariance for the complete table.
- Luminosity, radiative-correction, and detection-efficiency uncertainties
  have common components across the two regions. No cross-region covariance
  is invented. The two bins from 3.0 to 3.2 GeV have the stated $J/\psi$
  contribution subtracted.

All cross sections remain in the source's Born convention. The low-energy
covariance diagonals agree with the table errors within its 0.001 nb rounding:
the maximum differences are 0.000492 nb statistical and 0.000500 nb residual
systematic. The values, source hashes, covariance hashes, and extraction
method are in
[the provenance manifest](../../../data/odd_parity/babar_2021_3pi_provenance.json).

# Fit relevance and limit

Unlike the cached 2004 BaBar table, this dataset covers the $\omega$ and
$\phi$ peaks. It is now a usable measurement input for a three-pion fit.
It does not supply the finite-width Khuri--Treiman predictor or its response
to the two VPPP combinations. No fit or coefficient profile is inferred from
the data extraction alone.

Reproduce the extraction with:

    .venv/bin/python scripts/extract_babar_2021_three_pion.py

The extraction tests check bin counts, first/last bins, covariance symmetry
and positive semidefiniteness, and agreement of covariance diagonals with the
rounded published errors.
