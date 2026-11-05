---
title: Channel input manifest for the low-energy doubly virtual ALP TFF
status: mixed
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx, Niecknig:2012sj, Hoferichter:2023bjm, Hoferichter:2025lcz]
provenance: assumption
---

# Scope and common conventions

This manifest records the inputs required to assemble the low-energy part of
the doubly virtual ALP transition form factor (TFF). It summarizes the
linked, indexed derivations in this repository; no external source was
consulted for this inventory. Use the unnormalized scalar convention
$i\mathcal M=i e^2 R F_{a\gamma^*\gamma^*}$, the Alda:2026nrd no-mixing
source basis, natural units, metric $(+---)$, and $\epsilon^{0123}=+1$
[@Alda:2026nrd; @Messerli:2025rnv]. Electromagnetic charge factors are
stripped from one-current matrix elements where specified in the channel
records.

An entry marked *formal* has a derivation or modeling prescription but
lacks the physical numerical inputs needed for a spectral prediction.
Uncertainty provenance distinguishes known model variations from missing
input covariance; it does not imply that a numerical uncertainty has been
evaluated. All file links below point to records indexed in
[`INDEX.md`](INDEX.md).

## Channel inventory

| Contribution | Source basis and threshold/support | Current and production conventions | Phase and normalization | Uncertainty provenance | Input status and record |
|---|---|---|---|---|---|
| Two-pion isovector cut | Independent no-mixing ALP sources $g_A(\mu)$; $x\geq4M_\pi^2$ in the isospin-symmetric approximation. | $J^\mu_{\pi\pi}=(p_+-p_-)^\mu[F_\pi^V(x)]^*$ with $F_\pi^V(0)=1$; the production wave is $\mathcal F_{a\pi\pi\gamma^*}=f_{1,\mathrm{current}}=f_1/e$ in the real-epsilon crossing convention. | Single cut: $x\sigma_\pi^3/(96\pi)\,[F_\pi^V]^*\mathcal F$. Per-ordering double density: $x\sigma_\pi^3/(192\pi)\operatorname{Im}([F_\pi^V]^*\mathcal F)$, with the exchanged photon ordering included separately. The pion phase and production phase must be compatible. | Cutoff or smooth-cutoff changes are model variations. Numerical phase, rescattering, production uncertainties, and their correlations are not selected. | **Formal; numerical density incomplete.** The source-resolved production wave, numerical $F_\pi^V$, phase input, local anomalous completion, and single-spectral/subtraction audit remain required. See [two-pion spectral construction](alp-isovector-tff-spectral-construction.md), [pion current and production wave](alp-pion-vector-and-production-form-factors.md), and [source-resolved pion amplitude](alp-isovector-pion-production-amplitude.md). |
| Charged three-pion cut $\pi^+\pi^-\pi^0$ | Independent no-mixing ALP sources; $x\geq(m_++m_-+m_0)^2$. The retained state is the charged three-pion channel. | The current side is $\langle3\pi|J^\mu_{\rm em}|0\rangle$ and the ALP side is the crossed $\langle3\pi|J^\nu_{\rm em}|a\rangle$, with charge stripped and the full electromagnetic current retained on the ALP side. Born and strict vector-pole production pieces are derived; the complete complex production amplitude is not supplied. | Use $S=1+iT$, the recorded epsilon convention, relativistic phase space, and the conditional three-body discontinuity and spin-one weight. Keep $\operatorname{Disc}F$ unless Schwarz reflection conditions justify $2i\operatorname{Im}F$. The conventional strict-pole crossing factor is $\zeta_i=-1$; the physical decay normalization phase relative to production remains unresolved. | Published upper/lower curves are sensitivity bounds, not correlated replicas. A complex Dalitz grid, common phase match, and joint variations of phase input, subtraction prescription, fit parameters, and normalization are unavailable. | **Blocked for numerical spectrum.** Missing items: compatible complex $V_i\to3\pi$ grid and current amplitude, ALP-relative physical phase, and correlated variations; regular terms and contact completion also remain open. See [three-pion discontinuity](alp-three-pion-discontinuity.md), [three-pion phase match](alp-three-pion-pole-phase.md), and [numerical and phase input status](alp-three-pion-numerical-input-status.md). |
| Charged kaon cut $K^+K^-$ | Independent no-mixing ALP sources; $s\geq4M_{K^+}^2$. | $\langle K^+K^-|J^\mu_{\rm em}|0\rangle=(p_1-p_2)^\mu F_{K^+}^V(s)$, with $F_{K^+}^V(0)=1$. The retained virtual-current production amplitude is a symbolic odd-resonance tree result; the current and production factors must use the same charge-stripped convention. | Conditional distinct-particle cut coefficient is $s(1-4M_{K^+}^2/s)^{3/2}/(96\pi)$ times the conjugate current form factor and the reduced production wave. No $1/2!$ factor. The physical relative phase is not specified. | No selected complex timelike form-factor covariance, rescattering uncertainty, or correlated production variations. The retained tree truncation omits independent anomalous terms and loops. | **Formal; numerical cut incomplete.** Missing complex $F_{K^+}^V$, completed production amplitude and phase, and correlated uncertainties. See [virtual-photon kaon amplitude](alp-virtual-photon-kaon-amplitude.md). |
| Neutral kaon cut $K^0\bar K^0$ | Independent no-mixing ALP sources; $s\geq4M_{K^0}^2$, kept distinct from the charged threshold. | $\langle K^0\bar K^0|J^\mu_{\rm em}|0\rangle=(p_1-p_2)^\mu F_{K^0}^V(s)$, with $F_{K^0}^V(0)=0$. The retained virtual-current production amplitude is symbolic and source-resolved. | Same conditional distinct-particle $1/(96\pi)$ cut convention as the charged channel, with its own threshold and current form factor. No $1/2!$ factor. A three-flavor isoscalar-current half-sum requires consistent current and isospin conventions. | No complex timelike form-factor covariance or correlated production variations is selected. Strong-isospin use with physical mass splittings would be an additional approximation. | **Formal; numerical cut incomplete.** Missing complex $F_{K^0}^V$, completed production amplitude and phase, and correlated uncertainties. See [virtual-photon kaon amplitude](alp-virtual-photon-kaon-amplitude.md). |
| Pure narrow $\omega$ and $\phi$ poles | Independent no-mixing ALP sources $U_a,Y_a$; stable poles at the selected vector masses in each photon virtuality. There is no multiparticle threshold in the stable-pole term itself. | Signed neutral-vector photon residues and the ordered VJP/VVP numerators use the unnormalized scalar TFF convention. The meson-pole WZW contribution already absorbed in $\widehat c_\gamma$ is not added again. | Stable-pole denominators define the adopted approximation; their delta-function spectral weights are alternative representations of the same pole and must not be added again. The tree residues are signed, but no physical finite-width phase is supplied. | No quantitative error for the narrow-pole approximation, spectral tails, or correlated residue variations is selected. Uncertainty is model-dependent, especially near timelike peaks. | **Conditional model input.** Retain only where the same strength is not represented by a continuum; replace or match in the $\omega/\phi$ cuts, particularly the shared $\phi\to3\pi$ and $K\bar K$ strength. See [tree TFF](alp-doubly-virtual-tff.md), [mixed-vector strategy](alp-tff-mixed-vector-strategy.md), and [finite-width pole replacement](alp-isoscalar-vector-self-energy.md). |
| Mixed dressed $\rho I$ VVP, $I=\omega,\phi$ | Independent no-mixing ALP sources; the dressed $\rho$ response has the $\pi\pi$ support beginning at $s=4M_\pi^2$, while the other leg is the selected narrow $I$ pole. | Retain the ordered mixed VVP source numerator with $W_\rho^{\rm dress}(s)=R_\rho/[M_{0\rho}^2-s-\Delta_\rho^{\rm self}(s)]$ and $W_I^{\rm pole}$ on the other leg. This response-level term is not the full two-pion double-spectral TFF. | The self-energy convention is $\operatorname{Im}\Delta_\rho^{\rm self}(s+i0)=+\sqrt{s}\,\Gamma_{\rho\to\pi\pi}(s)$; the real part requires a subtracted dispersion relation. The ordered VVP numerator is retained as written in the mixed-term record; it has no independent cut prefactor. | The rho subtraction conditions, numerical self-energy, source-specific match to the two-pion continuum, and their uncertainties are not available. | **Formal; quantitative overlap unresolved.** The mixed term is generally nonzero and must be retained, but its shared rho strength must be matched to the two-pion continuum once. See [mixed rho-isoscalar VVP](alp-tff-mixed-rho-isoscalar-vvp.md), [mixed-vector strategy](alp-tff-mixed-vector-strategy.md), and [isoscalar source decomposition](alp-isoscalar-tff-source-decomposition.md). |

## Global calibration and completeness status

No channel entry above supplies a complete TFF over the unbounded virtualities
of the lepton loop. Higher isoscalar strength, subtraction terms, and the
short-distance continuation remain outside this inventory and require
separate matching or a stated omission. Numerical resonance calibration and
fit-driven uncertainty propagation remain conditional on `alp-uzh.8`; the
current nonconverged vector-sector likelihood is not accepted as a numerical
input. This manifest is an input-status record, not a numerical prediction.
