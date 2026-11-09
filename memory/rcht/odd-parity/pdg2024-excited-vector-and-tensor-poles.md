---
title: PDG 2024 masses and widths of excited vector resonances and a2(1320)
status: literature
created: 2026-11-07
updated: 2026-11-07
bibtex_keys: [ParticleDataGroup:2024cfk]
provenance: literature
---

# PDG 2024 excited vector and tensor inputs

Source: ParticleDataGroup:2024cfk particle listings (PDF per state under `memory/sources/ParticleDataGroup:2024cfk/`, checksums in `data/vector_sector/pdg2024_excited_poles.json`). Machine-readable copy in GeV: [pdg2024_excited_poles.json](../../../data/vector_sector/pdg2024_excited_poles.json).

These are the listing masses and widths (Breit-Wigner-type), not complex S-matrix pole positions. PDG gives one combined uncertainty for each of these headline values; no statistical/systematic split is given, so none is recorded. Values in MeV.

| State | Mass | Type | Width | Type |
|---|---|---|---|---|
| $\rho(1450)$ | $1465 \pm 25$ | OUR ESTIMATE | $400 \pm 60$ | OUR ESTIMATE |
| $\rho(1700)$ | $1720 \pm 20$ | OUR ESTIMATE | $250 \pm 100$ | OUR ESTIMATE |
| $K^*(1410)$ | $1414 \pm 15$ (scale factor 1.3) | OUR AVERAGE | $232 \pm 21$ (scale factor 1.1) | OUR AVERAGE |
| $K^*(1680)$ | $1718 \pm 18$ | OUR AVERAGE | $320 \pm 110$ (scale factor 4.2) | OUR AVERAGE |
| $\omega(1420)$ | $1410 \pm 60$ | OUR ESTIMATE | $290 \pm 190$ | OUR ESTIMATE |
| $\omega(1650)$ | $1670 \pm 30$ | OUR ESTIMATE | $315 \pm 35$ | OUR ESTIMATE |
| $\phi(1680)$ | $1680 \pm 20$ | OUR ESTIMATE | $150 \pm 50$ | OUR ESTIMATE |
| $\phi(2170)$ | $2164 \pm 6$ | OUR AVERAGE | $106\,{}^{+24}_{-18}$ (scale factor 2.0) | OUR AVERAGE |
| $a_2(1320)$ | $1318.2 \pm 0.6$ (scale factor 1.2) | OUR AVERAGE | $107 \pm 5$ | OUR ESTIMATE |

Notes from the listings.
- The $\rho(1450)$ mass and width and the $\phi(1680)$ width are labelled by PDG as an educated guess, with an error larger than that of the average of published values.
- For $a_2(1320)$ the width block prints both $107 \pm 5$ (OUR ESTIMATE) and $107.8 \pm 1.2$ (OUR AVERAGE).
- The $a_2(1320)$ listing also gives a T-matrix pole: $(1305\text{--}1321) - i\,(52\text{--}58)$ MeV (OUR ESTIMATE). The $K^*(1410)$ T-matrix pole is $(1368 \pm 38) - i\,(106\,{}^{+48}_{-59})$ MeV (OUR ESTIMATE); this is a single-analysis value.
- No PDG 2024 value of type OUR FIT occurs for these masses and widths.

## Comparison with hard-coded values in scripts

Checked in `scripts/fit_babar_kkpi_g_coefficients.py` (`propagator`) and `scripts/fit_babar_kkpi_global_vector_sector.py` (line 395, $a_2$). Values in GeV.

| Pole | Script | PDG 2024 | Status |
|---|---|---|---|
| $\rho(1450)$ | 1.465, 0.400 | 1.465, 0.400 | agree |
| $\rho(1700)$ | 1.720, 0.250 | 1.720, 0.250 | agree |
| $K^*(1410)$ | 1.414, 0.232 | 1.414, 0.232 | agree |
| $K^*(1680)$ | 1.718, 0.322 | 1.718, 0.320 | width differs by 0.002 GeV (well within the $\pm 0.110$ GeV uncertainty) |
| $\omega(1420)$ | 1.410, 0.290 | 1.410, 0.290 | agree |
| $\omega(1650)$ | 1.670, 0.315 | 1.670, 0.315 | agree |
| $\phi(1680)$ | 1.680, 0.150 | 1.680, 0.150 | agree |
| $\phi(2170)$ | 2.162, 0.100 | 2.164, 0.106 | mass differs by 0.002 GeV, width by 0.006 GeV (both within uncertainties) |
| $a_2(1320)$ | 1.3182, 0.107 | 1.3182, 0.107 | agree (width is the OUR ESTIMATE value) |

The docstring of `propagator` attributes the excited-state inputs to the PDG figures printed in the Qin:2024ulb appendix, so the two differing entries probably come from an earlier PDG edition; this was not checked. In the script, the $\omega$ and $\phi$ excited poles are labelled `omega_prime`/`phi_prime` (first) and the unnamed `second_mass, second_width` tuple (second); the $\rho$ and $K^*$ second poles are likewise unnamed. The $a_2$ pole is hard-coded in-line without a parameter name.
