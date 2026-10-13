---
title: Fixed PDG 2024 vector mass and electromagnetic-residue inputs
status: literature
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [ParticleDataGroup:2024cfk]
provenance: literature
---

# Vector-sector calibration inputs

The auditable machine-readable snapshot is [pdg2024.json](../../../data/vector_sector/pdg2024.json). It pins the official 2024 listing URLs and PDF SHA-256 hashes. PDFs dated 31 May 2024 were accessed/cached on 13 October 2026 under the ignored `memory/sources/ParticleDataGroup:2024cfk/` branch. The selected quantities are listing fitted masses, not complex S-matrix pole masses. [@ParticleDataGroup:2024cfk]

| Mass category | Value (MeV) | PDF page | Stated error scale factor |
|---|---:|---:|---:|
| rho neutral only, electron-positron data | $775.26\pm0.23$ | rho 1 | none stated |
| rho charged only, tau and electron-positron data | $775.11\pm0.34$ | rho 2 | none stated |
| K-star charged only, hadroproduced | $891.67\pm0.26$ | K-star 1 | none stated |
| K-star neutral only | $895.55\pm0.20$ | K-star 2 | $1.7$ |
| omega | $782.66\pm0.13$ | omega 1 | $2.0$ |
| phi | $1019.461\pm0.016$ | phi 1 | none stated |

All are OUR AVERAGE entries. Quoted errors already include the stated PDG scale factors; do not multiply again. Charged/neutral alternatives remain separate scenarios. These categories do not establish an isospin-limit mass or an electromagnetic subtraction. [@ParticleDataGroup:2024cfk]

| Direct electron-positron partial width | Value (keV) | PDF page and category |
|---|---:|---|
| rho | $7.04\pm0.06$ | rho 14, OUR FIT and identical OUR AVERAGE |
| omega | $0.60\pm0.02$ | omega 5, OUR EVALUATION |
| phi | $1.27\pm0.04$ | phi 6, OUR EVALUATION |
| phi alternative | $1.251\pm0.021$ | phi 6, OUR AVERAGE; error includes scale $1.1$ |

The alternative phi average is a separate sensitivity input, not an additional independent datum. The phi page was visually inspected to distinguish OUR EVALUATION from OUR AVERAGE. No scale factor is separately stated for the selected direct-width evaluations. The omega electron branching fraction is $(7.38\pm0.22)10^{-5}$ (scale $1.9$), whereas its total width is $8.68\pm0.13$ MeV. Their central product is $0.640584$ keV, different from the explicitly listed direct evaluation. Accordingly this snapshot uses the direct evaluation and does not manufacture a new independently propagated product error. [@ParticleDataGroup:2024cfk]

## Covariance and modeling assumptions

**Assumption for diagnostic likelihood:** use independent Gaussian quoted errors only as a conditional adequacy test. The selected listing summaries supply no complete joint covariance across masses and partial widths. Direct partial-width summaries avoid newly treating branching and total width as independent, but their underlying measurements still carry normalization and line-shape dependencies. For example, omega page 5 notes a measurement using $B_{3\pi}=0.891\pm0.007$ and total width $8.44\pm0.09$ MeV; phi page 6 notes a branching-sum assumption for CMD-2 and a KLOE extraction using the 2004 total width. Rho mass and electron width share some contributing experiments. A precision fit requires original likelihoods/covariances or declared correlation sensitivity. [@ParticleDataGroup:2024cfk]

**Assumed numerical convention:** the script may use fixed on-shell $\alpha_{\rm em}=1/137.035999084$. This numerical input is a stated assumption here, not a newly sourced result from the particle listings. Matching tree-level masses to these fitted masses omits self-energy and electromagnetic corrections; input errors and that model limitation must remain separate.
