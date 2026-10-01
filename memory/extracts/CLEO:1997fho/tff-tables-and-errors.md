---
title: CLEO light-pseudoscalar transition-form-factor tables
status: literature
created: 2026-09-30
updated: 2026-09-30
bibtex_keys: [CLEO:1997fho]
provenance: literature
---

The arXiv source archive is cached at `memory/sources/CLEO:1997fho/arxiv-source.tar`. Tables I--X contain the measured space-like meson-photon transition form factors for $\pi^0$, $\eta$, and $\eta'$; the paper separates statistical and systematic contributions and reports fitted form-factor information [@CLEO:1997fho]. The observable for the first fit release is the tabulated $Q^2F_{P\gamma^*\gamma}(Q^2)$ point, constraining the TFF normalization and its vector-exchange slope. The source archive contains no bin-to-bin covariance matrix, so only the quoted diagonal errors may be used unless a collaboration data release is found.
The literal 63 rows across the retained tables are stored in [`data/odd_parity/cleo_tff.csv`](../../../data/odd_parity/cleo_tff.csv), with provenance and versioned likelihood metadata in the adjacent JSON sidecar. The source column is labelled $0.01\times\mathrm{GeV}$; the CSV converts central values and both errors to GeV. This table has a validated loader but no connected predictor.
