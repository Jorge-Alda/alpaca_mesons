---
title: BaBar isovector KK star K spectrum
status: literature
created: 2026-10-29
updated: 2026-10-29
bibtex_keys: [BaBar:2007ceh]
provenance: literature
---

# Observable and amplitude separation

[@BaBar:2007ceh] analyzes ISR production of $K^0_S K^\pm\pi^\mp$. A Dalitz fit with charged and neutral $K^*(892)$ bands separates isoscalar and isovector $K^*(892)K$ amplitudes. It reports the two component cross sections and their relative phase as functions of total center-of-mass energy (Sec. VI, Tables 6 and 7, and Fig. 18). The charged and neutral $K^*K$ channels involve opposite interference signs between the two isospin components (Eq. (13)).

The isovector table covers $1.40$--$3.08$ GeV in variable-width bins. The overlap with the $τ$ endpoint contains the bins from $1.40$ GeV to $1.78$ GeV; its last bin reaches $1.78$ GeV and therefore straddles the physical $τ$ endpoint. Preserve that published bin and flag the partial overlap rather than shortening its central value. The table lists statistical errors. The paper additionally assigns a 15% systematic uncertainty to the Dalitz-fit procedure but does not give its bin-to-bin covariance. Any use of a fully correlated 15% nuisance is therefore a stated likelihood assumption, not a measured covariance.

# Relevance and limitations for coefficient extraction

The isovector projection is directly relevant to the vector-current $KK\pi$ coefficient audit and provides energy-dependent data rather than only a branching fraction. Its energy interval overlaps the kinematic range for $τ\to KK\pi\nu_τ$, subject to the endpoint bin caveat. Connecting the neutral-current projection to the charged-current channel still requires an explicit CVC/isospin normalization and consistent charge-state and $K^*$ branching conventions.

The spectrum is not a direct measurement of the low-energy $c_i,d_i,g_i$. The paper describes the isovector component by a broad $ρ'$ (identified with the $ρ(1450)$) and notes that a single structure may be inadequate; its fitted width may mimic a more complex shape. This excited-vector contribution is absent from the single-multiplet low-energy coefficient predictor currently present in the repository. The table can be preserved as a candidate input and used for a model-comparison or conditional fit only after adding an explicit excited-vector component and its nuisance parameters. It must not presently be counted as an extraction of the missing VJP or VPPP coefficient directions.

The paper's reported isovector cross sections and relative phase are experimental fit products; the extracted isospin components have a common Dalitz-model systematic. The relative phase is also relevant: the paper finds that the isovector shape alone can look phase-space-like, whereas the phase indicates a resonant contribution.

# Provenance

Source: [@BaBar:2007ceh], arXiv:0710.4451, Sec. VI and Sec. X.1, especially Table 7, Eqs. (13)--(14), and the discussion of the broad $ρ'$ in Sec. XI. The source archive is indexed at [memory/sources/BaBar:2007ceh](../../sources/BaBar:2007ceh/INDEX.md). A machine-readable transcription of the tau-overlap bins is stored under `data/vector_sector/` with its source table, units, and uncertainty assumptions recorded in the accompanying metadata.
