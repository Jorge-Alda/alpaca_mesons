---
title: Higher isoscalar strength and source-specific short-distance input gates
status: mixed
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv]
provenance: assumption
---

# Scope and evidence

This is the stage-5 feasibility audit of the [isoscalar TFF plan](alp-isoscalar-tff-plan.md), not a short-distance derivation or an effective-pole fit. The inputs are the retained [tree TFF](alp-doubly-virtual-tff.md) and the [lepton-pair calculation plan](alp-lepton-pair-calculation-plan.md). The latter already requires free runtime ALP mass and couplings; selecting a benchmark coupling is therefore unnecessary. The user specifies $2m_e<m_a<2m_\tau$ and the virtualities needed for the complete two-photon lepton loop. The virtuality integration is unbounded, so a finite low-energy range cannot replace the high-energy continuation and renormalization prescription.

Alda:2026nrd fixes the no-mixing source convention. Messerli:2025rnv supplies a meson dispersive strategy, not source-specific ALP spectral strengths or high-virtuality coefficients [@Alda:2026nrd; @Messerli:2025rnv]. This audit relies on these existing repository records, with both bibliography entries inspected; it introduces no additional literature result.

# Higher isoscalar strength: remaining inputs

An effective-pole fit cannot be identified from the present plan. First specify the low-energy contributions already included, their subtraction terms, and the separation from explicit vector poles. Next supply source-resolved single- or double-virtual information, or independent source-specific matching conditions, with uncertainties and correlations. A fit to physical-meson data alone does not establish the ALP response to every independent quark or gluon source.

Any future effective-pole description must state which omitted channels it represents, the allowed virtuality range, the number of independent parameters constrained by those inputs, and comparison with an alternative description. Its uncertainty remains a model uncertainty. No higher-state masses, residues, or correlations are selected in this audit.

# Short-distance matching: remaining inputs

For every independent supplied hadronic-scale ALP source, identify the relevant QCD current matrix element, operator normalization, renormalization scale and scheme, and the short-distance expansion appropriate to the virtuality limit being tested. The single-virtual limit and fixed-ratio doubly virtual limit require separate verification. Establish the matching to the no-mixing low-energy source without reconstructing an ALP gluonic or singlet coefficient from pion or physical-meson matching.

The retained tree expression is a low-energy approximation. Expanding its stable poles is an algebraic property of that approximation and does not determine the QCD coefficient. Meson asymptotic coefficients, meson decay constants, or the existing source matrices alone cannot be substituted for an independently established ALP-source matching calculation. No asymptotic coefficient is assigned here.

# Acceptance and next step

Stage 5 remains open. A subsequent [executable hard-current derivation](alp-hard-current-ope.md) now fixes the leading local parity-odd Wilson coefficient and nearly symmetric Euclidean scalar projection, while retaining explicit ALP axial-current matrix elements. The next step is their renormalized source-to-current matching, the other hard limits and low-energy pole/continuum separation. Only after these inputs exist should a matching scale be selected and varied, or effective-pole parameters fitted. Check the resulting source-specific asymptotics in the relevant virtuality limits while preserving the independently fixed low-energy normalization and discontinuities within the declared uncertainty.

The feasible independent work in the meantime is the formal source/current bookkeeping, three-pion unitarity normalization, and virtual-current kaon amplitude derivation in stages 1, 3, and 4. Their symbolic parts do not require invented numerical three-pion or kaon grids. Completion of those parts does not close the missing short-distance inputs.
