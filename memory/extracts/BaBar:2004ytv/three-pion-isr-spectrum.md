---
title: BaBar ISR three-pion line shape
status: literature
created: 2026-09-30
updated: 2026-09-30
bibtex_keys: [BaBar:2004ytv]
provenance: literature
---

The arXiv source archive is cached at `memory/sources/BaBar:2004ytv/arxiv-source.tar`. BaBar measures the ISR Born cross section $\sigma(e^+e^-\to\pi^+\pi^-\pi^0)$ over the $\omega$, $\phi$, and higher-vector region; the mass spectrum and fitted products $\mathcal B(V\to e^+e^-)\mathcal B(V\to3\pi)$ constrain timelike vector mixing, resonance phases, and the $V\to3\pi$ amplitude [@BaBar:2004ytv]. The paper tabulates spectrum points and uncertainties but supplies no full covariance matrix in the source archive.
The 78 tabulated Born-spectrum rows are stored in [`data/odd_parity/babar_3pi_born_spectrum.csv`](../../../data/odd_parity/babar_3pi_born_spectrum.csv), retaining the $N_{corr}$, efficiency, luminosity, and cross-section error columns. The source states that the cross-section systematic errors are fully correlated across mass bins. Accordingly, only the statistical cross-section error is diagonal; the sidecar records one additive standardized nuisance with the tabulated `sigma_syst_nb` response vector in nb, plus a separate unconnected $1.2\%$ luminosity normalization nuisance [@BaBar:2004ytv].
