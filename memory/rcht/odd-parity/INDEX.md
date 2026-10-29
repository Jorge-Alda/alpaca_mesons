---
title: Odd-parity RChT memory index
status: mixed
created: 2026-10-01
updated: 2026-10-29
bibtex_keys: []
provenance: assumption
---

# Odd-parity RChT

Odd-parity operator extensions, coupling prescriptions, numerical kernels, calibration data, and fit planning.

- [levi-civita-convention.md](levi-civita-convention.md) — canonical upper-index orientation, action-preserving reference conversion and explicit legacy-output adapter.
All displayed physical tensors use the canonical convention in
[levi-civita-convention.md](levi-civita-convention.md). Older linked scripts,
fits and signed artifacts are explicitly reference-backend inputs; their
odd coefficients and tensor scalars require that record's adapter before
being combined with canonical photon or resonance amplitudes.
- [odd-parity-physics-draft.md](odd-parity-physics-draft.md) — physics-only manuscript narrative and precise physical wording, direct Alda:2026nrd conventions, audited photon contact and appendix/missing-work comments.
- [two-source-alp-rcht-odd-parity-extension.md](two-source-alp-rcht-odd-parity-extension.md) — two-source anomaly seed and ALP insertion.
- [wzw-gamma3pi-contribution.md](wzw-gamma3pi-contribution.md) — gauged-WZW $\gamma3\pi$ vertex and the photon--vector-mixing contribution it induces.
- [vector-three-pion-contact.md](vector-three-pion-contact.md) — symmetry-allowed local $V3\pi$ contact, its nonzero vertex and its IBP/leading-EOM reduction into the $VPPP$ basis plus an SU($N_f$) double trace.
- [isoscalar-vppp-three-pion-projection.md](isoscalar-vppp-three-pion-projection.md) — explicit five-operator isoscalar charged-three-pion projection, two coefficient combinations, rank/null directions, and consequences for off-shell TFF input.
- [isoscalar-vppp-photon-response.md](isoscalar-vppp-photon-response.md) — tree photon-source insertion for the local isoscalar VPPP vertex and its exact pole/contact degeneracy.
- [isoscalar-vppp-dispersive-matching.md](isoscalar-vppp-dispersive-matching.md) — conditional map to the $a(q^2)$ subtraction and vector residues, including rank and sum-rule mismatch.
- [isoscalar-vppp-fit-preparation.md](isoscalar-vppp-fit-preparation.md) — parameterization, likelihood, identifiability checks, data audit, and prerequisites for extracting the two local combinations.
- [kkpi-vector-current-coefficient-rank.md](kkpi-vector-current-coefficient-rank.md) — exact KK pi vector-current VJP/VVP/VPPP coefficient maps, charged/neutral ranks, cross-sector compensation and combination with the preliminary and isoscalar spans.
- [kkpi-isovector-spectrum-input-and-extraction-readiness.md](kkpi-isovector-spectrum-input-and-extraction-readiness.md) — BaBar's isovector $K^*(892)K$ bins in the tau-overlap region, covariance limits, and missing higher-resonance/CVC ingredients for a coefficient fit.
- [babar-2021-three-pion-fit-input.md](babar-2021-three-pion-fit-input.md) — extracted 114-bin low-energy cross section with covariance and 72-bin higher-energy table for the new fit.
- [vppp-dispersive-profile-input-contract.md](vppp-dispersive-profile-input-contract.md) — finite-width spectral/profile calculator, exact sum-rule elimination, zero-width check, and remaining physical inputs.
- [vppp-finite-width-spectral-templates.md](vppp-finite-width-spectral-templates.md) — SymPy-derived finite-width resonance and conformal spectral bases, upper-lip convention, and required width/tail inputs.
- [vppp-two-body-partial-width-inputs.md](vppp-two-body-partial-width-inputs.md) — conditional $k^3$ and $k^3/s$ running-width shapes for caller-supplied two-body channels; rescattered $3\pi$ widths remain external.
- [vppp-three-pion-running-width.md](vppp-three-pion-running-width.md) — phase-space integral of a supplied full KT scalar and its normalized running-width relation to the vector-current cross section.
- [u3-nonet-odd-parity-extension-options.md](u3-nonet-odd-parity-extension-options.md) — U(3) nonet and OZI options.
- [odd-parity-coupling-prescriptions-and-direct-alp-vgamma.md](odd-parity-coupling-prescriptions-and-direct-alp-vgamma.md) — $c_i,d_i$ prescriptions and direct $aV\gamma$ terms.
- [odd-prescription-matching.md](odd-prescription-matching.md) — neutral pole matching for Aloni:2018vki and Bai:2025fvl, derived even-parity $G_V$, split flat-schema inputs and the separate rank-five octet OPE family.
- [odd-parity-rcht-implementation-decisions.md](odd-parity-rcht-implementation-decisions.md) — selected implementation baseline.
- [odd-parity-rcht-numerical-kernel.md](odd-parity-rcht-numerical-kernel.md) — numerical odd-parity RChT kernel.
- [omega-pi-gamma-rcht-predictor.md](omega-pi-gamma-rcht-predictor.md) — on-shell $\omega\to\pi^0\gamma$ predictor.
- [odd-parity-coefficient-tabulation-status.md](odd-parity-coefficient-tabulation-status.md) — status of numerical coefficient lists.
- [odd-parity-data-implementation.md](odd-parity-data-implementation.md) — validated data loaders and readiness.
- [odd-parity-global-fit-dataset-plan.md](odd-parity-global-fit-dataset-plan.md) — staged global-fit dataset plan.
- [correlated-odd-parity-fit-infrastructure.md](correlated-odd-parity-fit-infrastructure.md) — fixed-prior correlated fit infrastructure.
- [odd-parity-fit-optimizer-choice.md](odd-parity-fit-optimizer-choice.md) — pure-scipy optimizer decision for the coupled fit.
- [babar-3pi-spectrum-fit-blocker.md](babar-3pi-spectrum-fit-blocker.md) — why the BaBar ISR table misses the $\omega/\phi$ peaks and needs $\omega'/\omega''$.
- [three-pion-resonance-separability.md](three-pion-resonance-separability.md) — numerical proof that no BaBar mass window is free of $\omega'/\omega''$.
- [kloe-phi-three-pion-fit-blocker.md](kloe-phi-three-pion-fit-blocker.md) — why the held KLOE $\phi\to3\pi$ table cannot constrain $d_1,d_2$ at tree level.
- [dispersive-three-pion-treatment-plan.md](dispersive-three-pion-treatment-plan.md) — scoped dispersive $V\to3\pi$ project (c) and staged alternative (c').
- [radiative-channel-saturation-probe.md](radiative-channel-saturation-probe.md) — rank-is-not-information criterion and the cheap radiative probe.
- [radiative-sm-flavor-algebra.md](radiative-sm-flavor-algebra.md) — SM flavor algebra for $V\to P\gamma$ and $P\to V\gamma$: charge and quark-mass weights, per-channel $F_Q,F_M$; every channel distinct.
- [radiative-sm-amplitude.md](radiative-sm-amplitude.md) — action-derived one-photon amplitude with linear $f_-$ and $c_4$, complete neutral-vector exchange, canonical vector conventions, and common/split-pole rank results.
- [radiative-rank-saturation.md](radiative-rank-saturation.md) — action-based common-pole rank result, split-pole sensitivity, TFF rank, and assumptions.
- [radiative-width-kinematics.md](radiative-width-kinematics.md) — per-operator kinematic factors for the direct and VVP-exchange operators, the single Lorentz invariant, the uniqueness selection rule, the assembled per-channel width, and reproduction of the published $\omega\to\pi^0\gamma$ bracket.
- [radiative-tff-kinematics.md](radiative-tff-kinematics.md) — Fourier-consistent VJP numerator and its complete two-photon placement in the TFF amplitude.
- [radiative-tff-flavor.md](radiative-tff-flavor.md) — $t_2$ flavor factors for $\pi^0,\eta,\eta'\to\gamma^*\gamma$, with stated mixing convention; full amplitude and mass-spurion factors are in the amplitude record.
- [radiative-tff-amplitude.md](radiative-tff-amplitude.md) — normalized tree-level $P\to\gamma^*\gamma$ tensor-RChT amplitude, WZW calibration, data Observation connection, and rank consequences for equal versus independent mass-spurion inputs.
- [next-steps.md](next-steps.md) — current results and remaining source-grounded meson inputs, asymptotic behavior, and matching tasks.
- [p-to-gamma-gamma.md](p-to-gamma-gamma.md) — real-photon operator audit, identical-photon width, pion normalization and OPE coupling sensitivity.
- [alda-neutral-amplitude-inputs.md](alda-neutral-amplitude-inputs.md) — reconstructed arXiv-v1 neutral generators, field conversion and explicit tree/physical kinematics.
- [vector-stage4-amplitudes.md](vector-stage4-amplitudes.md) — action-derived spectral source insertion into every radiative and two-photon vector slot.
- [vector-stage4-rank.md](vector-stage4-rank.md) — all 64 retained action scenarios, inherited-OPE sensitivity ranks and their interpretation limits.
- [vector-stage5-joint-fit-specification.md](vector-stage5-joint-fit-specification.md) — joint likelihood, covariance/source treatment, tree-level precision decision and remaining completion gates.
- [vector-stage5-preflight.md](vector-stage5-preflight.md) — executable data-only ranks and spectral single-virtual asymptotic coefficients without inherited OPE priors.
- [vector-chiral-octet-matching.md](vector-chiral-octet-matching.md) — enlarged-action chiral projection, explicit VVP OPE solution and common-pole TFF incompatibility.
- [vector-physical-tff-matching.md](vector-physical-tff-matching.md) — physical-channel hard-ray conditions and the source-breaking obstruction in all 64 retained scenarios.
- [vector-matching-mismatch-interpretation.md](vector-matching-mismatch-interpretation.md) — saved explanation of both mismatches, their possible origins, remedies and low-energy scope limits.
- [odd-breaking-low-energy-calibration-plan.md](odd-breaking-low-energy-calibration-plan.md) — selected staged one-spurion single/double-trace odd extension, low-energy fit and ALP propagation plan.
- [odd-breaking-flavor-enumeration.md](odd-breaking-flavor-enumeration.md) — exact abstract single/double traces, C/Hermitian and Bose sectors, singlet and neutral-photon projections, and literature dictionary audit.
- [odd-breaking-lorentz-reduction.md](odd-breaking-lorentz-reduction.md) — exact Lorentz/IBP and linear Bianchi reduction with independent sources, source-gradient and scalar projections.
- [odd-breaking-vertex-catalogue.md](odd-breaking-vertex-catalogue.md) — 205 bounded source-resolved coefficient tensors, dimensions, normalization, target limits and stage-B handoff.
- [odd-trace-source-eom-audit.md](odd-trace-source-eom-audit.md) — executed neutral U(3) source EOM, induced topological contacts and the remaining full-action count milestone.
- [odd-covariant-closure.md](odd-covariant-closure.md) — verified covariant identities, source-gradient remainders and the selected explicit-vector counting convention.
- [odd-goldstone-eom-multiplier-reduction.md](odd-goldstone-eom-multiplier-reduction.md) — exact covariant LO Goldstone EOM multiplier matrices and their two-spurion closure obstruction.
- [odd-leading-covariant-operator-count.md](odd-leading-covariant-operator-count.md) — complete leading covariant quotient: 40 independent operators, original rank eleven and 29 new representatives.
- [odd-covariant-source-gauge-reduction.md](odd-covariant-source-gauge-reduction.md) — independent-source frame, ordered derivatives, exact trace-jet IBP and cross-block completion.
- [odd-three-flavor-trace-gram.md](odd-three-flavor-trace-gram.md) — exact U(3) color Gram closure, including the unresolved direction beyond finite source probes.
- [odd-higher-covariant-operator-count.md](odd-higher-covariant-operator-count.md) — complete bounded action count: 228 independent operators, original rank eleven, 217 new and an explicit basis.
- [odd-trace-scalar-kernels.md](odd-trace-scalar-kernels.md) — normalized leading single/double-trace SM kernels, tests and 15 additional formal scalar EM directions.
- [alp-target-amplitude-audit.md](alp-target-amplitude-audit.md) — inventory of the four ALP target amplitudes, fresh blind-direction checks and missing-calculation decisions.
- [alp-vector-channel-coefficient-overview.md](alp-vector-channel-coefficient-overview.md) — coefficient map for even, VVP, VJP, and VPPP operators across ALP channels; exact TFF insensitivity to the common $c_2,c_5,c_6$ shift.
- [alp-two-boson-blind-direction-audit.md](alp-two-boson-blind-direction-audit.md) — exact VJP/VVP null projections for canonical vector and photon tensor legs.
- [alp-vv-completion-plan.md](alp-vv-completion-plan.md) — steps to assemble the normalized physical ALP two-vector amplitude and width.
- [alp-vgamma-completion-plan.md](alp-vgamma-completion-plan.md) — steps to complete the ALP vector-photon amplitude including exchange and source/contact terms.
- [odd-operator-alp-relevance.md](odd-operator-alp-relevance.md) — target-relevant operator families, projection limits and candidate remedies for the two matching mismatches.
- [odd-higher-two-boson-combinations.md](odd-higher-two-boson-combinations.md) — exact stripped higher two-boson coefficient maps: 66 abstract and 51 neutral combinations, with independent ALP source weights.
- [odd-action-power-counting-proposal.md](odd-action-power-counting-proposal.md) — proposed momentum/quark-mass/OZI hierarchy, first-correction truncation, resonance-domain limitations and required parent-column rank calculation.
- [odd-bounded-hilbert-series.md](odd-bounded-hilbert-series.md) — character-derived raw counts and finite graded quotient presentation, with inherited covariant syzygies explicitly identified.
- [odd-hilbert-power-counting.md](odd-hilbert-power-counting.md) — homogeneous chiral grading and executed parent-filtration dimensions 12, 101 and 228.
- [odd-resonance-derivative-counting.md](odd-resonance-derivative-counting.md) — exact hard/residual momentum split and conservative derivative-parent promotion for on-shell vectors.
- [odd-power-counted-channel-maps.md](odd-power-counted-channel-maps.md) — joint normalized neutral cubic parent maps, source-EOM/contact consistency and nested truncation ranks.
- [odd-charged-cubic-kernels.md](odd-charged-cubic-kernels.md) — all charged raw parent gradients, composite f-minus derivative momentum, Ward and LO EOM checks.
- [odd-alp-source-lsz-assembly.md](odd-alp-source-lsz-assembly.md) — general diagonal UV profile, source EOM, 964 LSZ identities and signed photon matching in Bauer's coefficient normalization.
- [odd-bounded-sm-numerical-kernels.md](odd-bounded-sm-numerical-kernels.md) — reusable neutral/charged 241-parent gradients, original-coupling embedding and independent symbolic/legacy checks.
- [odd-low-energy-measurement-jacobians.md](odd-low-energy-measurement-jacobians.md) — all 64 conditional measured-width Jacobians, nuisance projection, nested ranks and data-adequacy controls.
- [vector-dalitz-form-factor-fit-plan.md](vector-dalitz-form-factor-fit-plan.md) — proposed joint fit of $\omega\to\pi^0\ell^+\ell^-$ spectra, with finite-virtuality sensitivity and null-space rank checks.
- [vector-dalitz-conditional-fit.md](vector-dalitz-conditional-fit.md) — conditional A2/NA60 fit, one resolved old-null shape direction and the poor-fit limitation.
- [vector-dalitz-pion-tff-joint-fit.md](vector-dalitz-pion-tff-joint-fit.md) — conditional joint fit of A2/NA60 and CLEO/BaBar pion TFF data, one resolved old-null direction, and strong model tension.
- [doubly-virtual-etap-joint-fit.md](doubly-virtual-etap-joint-fit.md) — adds BaBar double-tag $\eta'$ data, reports no extra old-null direction, and exposes a large conflict with the retained tree model.
- [vector-tff-null-lifting-recommendations.md](vector-tff-null-lifting-recommendations.md) — defines short-distance-compatible TFF behavior, concrete Dalitz-model improvements, and staged null-lifting criteria.
- [improved-pion-tff-two-vector-matching.md](improved-pion-tff-two-vector-matching.md) — executable pion TFF hard limits, odd-coefficient derivatives, source-based dispersive continuum and symbolic low-energy moment matching.
- [pion-tff-spectral-input-contract.md](pion-tff-spectral-input-contract.md) — source-normalized complex spectral-grid contract and correlated pion-TFF moment integrals; physical amplitude grids remain unselected.
- [pion-tff-spectral-input-source-audit.md](pion-tff-spectral-input-source-audit.md) — audit of the Hoferichter:2018kwz supplementary table and local source files; identifies why they do not yet supply the complex spectral amplitudes.
- [pion-vector-form-factor-central-grid.md](pion-vector-form-factor-central-grid.md) — calculated central $F_\pi^V(x)$ grid from the exact-isospin three-subtracted Omnès prescription, with explicit phase continuation and numerical limitations.
- [odd-low-energy-calibration.md](odd-low-energy-calibration.md) — historical 1792 conditional joint fits with an identified overly restrictive double-trace residue-sign condition; corrected focused refit and remaining validation gates.
- [odd-preliminary-fit.md](odd-preliminary-fit.md) — central original-eleven-operator provisional fit with a minimal parameter-only JSON and conditional $G_V$ value; research archives retain the corrected double-trace refit and targeted input sensitivities.
- [preliminary-fit-lagrangian-and-null-directions.md](preliminary-fit-lagrangian-and-null-directions.md) — full retained even and original eleven-operator odd Lagrangians at the preliminary fit, including its four fixed-nuisance coefficient null directions.
- [alp-vv-amplitude.md](alp-vv-amplitude.md) — analytic no-mixing $a\to\omega\omega,\omega\phi,\phi\phi$ amplitudes with VJP photon conversion and the fitted one-trace even-vector source.
- [alp-vgamma-amplitude.md](alp-vgamma-amplitude.md) — no-mixing $a\to\omega\gamma,\phi\gamma$ amplitudes with the broken even-sector photon response and explicit null-direction projections.
- [alp-vv-width.md](alp-vv-width.md) — polarization sum, phase space, identical-particle factors and analytic narrow-vector widths for the three channels.
- [odd-low-energy-profiles.md](odd-low-energy-profiles.md) — finite sign searches, re-anchored source-boundary fits and independently certified source profiles; global minima and precise level crossings remain unproven.
- [odd-alp-measured-null-audit.md](odd-alp-measured-null-audit.md) — 448 finite neutral target/null projections, source/contact certificates and correction-model witnesses.

- [vppp-external-kt-grid.md](vppp-external-kt-grid.md) — external one-subtraction KT grid adapter, pion-unit correction, prepared Madrid phase alternatives and historical synthetic software diagnostic.
- [vppp-pipi-phase-input.md](vppp-pipi-phase-input.md) — corrected Madrid P-wave phase tables in pion units, mass assumptions and fixed KT-continuation caveat.

- [vppp-physical-kt-convergence.md](vppp-physical-kt-convergence.md) — physical-input KT iteration, external phase-space quadrature, cutoff and Madrid-phase sensitivity checks.
