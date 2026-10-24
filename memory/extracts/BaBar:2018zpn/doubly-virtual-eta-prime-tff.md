---
title: BaBar double-tag doubly virtual eta-prime transition form factor
status: literature
created: 2026-10-24
updated: 2026-10-24
bibtex_keys: [BaBar:2018zpn]
provenance: literature
---

# Measurement and convention

BaBar:2018zpn measures $e^+e^-\to e^+e^-\eta'$ in double-tag mode and
reports the spacelike doubly virtual form factor for $2<Q_1^2,Q_2^2<60$
$\mathrm{GeV}^2$. Its amplitude convention is
$$
T=-i4\pi\alpha\,\epsilon_{\mu\nu\beta\gamma}
\varepsilon_1^\mu\varepsilon_2^\nu q_1^\beta q_2^\gamma
F_{\eta'}(Q_1^2,Q_2^2),\qquad Q_i^2=-q_i^2.
$$
The scalar is reported in $10^{-3}\,\mathrm{GeV}^{-1}$. A fit using the
repository's TFF must compare magnitudes unless its epsilon and amplitude
phase conventions have been explicitly converted. The data concern the
physical $\eta'$ state, so the predictor requires its physical flavor/mixing
matrix and source-dependent mass-spurion traces; pion-only source factors
cannot be reused [@BaBar:2018zpn].

# Published table

The paper gives efficiency-weighted average virtualities for five broad
two-dimensional bins. The statistical errors are asymmetric; each point also
has separate systematic and model errors. Two bins have low-significance
central values and published 90% confidence upper limits.

| $\overline{Q_1^2}$ | $\overline{Q_2^2}$ | $10^3F_{\eta'}$ | Statistical error | Systematic | Model | 90% upper limit |
|---:|---:|---:|---:|---:|---:|---:|
| 6.48 | 6.48 | 14.32 | $^{+1.95}_{-1.89}$ | 0.83 | 0.14 | — |
| 16.85 | 16.85 | 5.35 | $^{+1.71}_{-2.15}$ | 0.31 | 0.42 | 14.53 |
| 14.83 | 4.27 | 8.24 | $^{+1.16}_{-1.13}$ | 0.48 | 0.65 | — |
| 38.11 | 14.95 | 6.07 | $^{+1.09}_{-1.07}$ | 0.35 | 1.21 | — |
| 45.63 | 45.63 | 8.71 | $^{+3.96}_{-8.71}$ | 0.50 | 1.04 | 32.03 |

The virtualities and form-factor errors in the table use $\mathrm{GeV}^2$
and $10^{-3}\,\mathrm{GeV}^{-1}$, respectively. The form-factor value is
extracted from a two-dimensional bin-integrated cross section using the
paper's stated model procedure and then quoted at the efficiency-weighted
average virtualities. The model uncertainty reflects the extraction
dependence. Treating those points as pointwise measurements at the averages
is therefore an additional fit approximation. The source does not provide a
bin-to-bin covariance matrix; any diagonal likelihood must be labelled an
approximation. Retain the upper limits in a censoring or sensitivity analysis
rather than treating their central values as ordinary high-significance
Gaussian observations [@BaBar:2018zpn, Table I and Sec. V].

# Source and provenance

The arXiv source archive is cached at
`memory/sources/BaBar:2018zpn/arxiv-source.tar.gz`; the tabulated results are
in `journal.tex`, Table I. The source-provided INSPIRE-HEP BibTeX entry is
`BaBar:2018zpn`.
