---
title: Channel input manifest for the low-energy doubly virtual ALP TFF
status: mixed
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx, Niecknig:2012sj, Hoferichter:2023bjm, Hoferichter:2025lcz, BABAR:2021cde, KLOE:2003kas, BaBar:2007ceh, Qin:2024ulb, ParticleDataGroup:2024cfk]
provenance: assumption
---

# Scope and common conventions

This manifest records the inputs required to assemble the low-energy part of
the doubly virtual ALP transition form factor (TFF). It summarizes the
linked, indexed derivations, literature extracts, and data files in this
repository; no additional external source was consulted for this audit. Use
the unnormalized scalar convention
$i\mathcal M=i e^2 R F_{a\gamma^*\gamma^*}$, the Alda:2026nrd no-mixing
source basis, natural units, metric $(+---)$, and $\epsilon^{0123}=+1$
[@Alda:2026nrd; @Messerli:2025rnv]. Electromagnetic charge factors are
stripped from one-current matrix elements where specified in the channel
records.

An entry marked *formal* has a derivation or modeling prescription but
lacks the physical numerical inputs needed for a spectral prediction.
Uncertainty provenance distinguishes known model variations from missing
input covariance; it does not imply that a numerical uncertainty has been
evaluated. All project record links below lead to records named in their
applicable topic or extract indexes.

## Channel inventory

| Contribution | Source basis and threshold/support | Current and production conventions | Phase and normalization | Uncertainty provenance | Input status and record |
|---|---|---|---|---|---|
| Two-pion isovector cut | Independent no-mixing ALP sources $g_A(\mu)$; $x\geq4M_\pi^2$ in the isospin-symmetric approximation. | $J^\mu_{\pi\pi}=(p_+-p_-)^\mu[F_\pi^V(x)]^*$ with $F_\pi^V(0)=1$; the production wave is $\mathcal F_{a\pi\pi\gamma^*}=f_{1,\mathrm{current}}=f_1/e$ in the real-epsilon crossing convention. | Single cut: $x\sigma_\pi^3/(96\pi)\,[F_\pi^V]^*\mathcal F$. Per-ordering double density: $x\sigma_\pi^3/(192\pi)\operatorname{Im}([F_\pi^V]^*\mathcal F)$, with the exchanged photon ordering included separately. The pion phase and production phase must be compatible. | Cutoff or smooth-cutoff changes are model variations. Numerical phase, rescattering, production uncertainties, and their correlations are not selected. | **Formal; numerical density incomplete.** The source-resolved production wave, numerical $F_\pi^V$, phase input, local anomalous completion, and single-spectral/subtraction audit remain required. See [two-pion spectral construction](alp-isovector-tff-spectral-construction.md), [pion current and production wave](alp-pion-vector-and-production-form-factors.md), and [source-resolved pion amplitude](alp-isovector-pion-production-amplitude.md). |
| Charged three-pion cut $\pi^+\pi^-\pi^0$ | Independent no-mixing ALP sources; $x\geq(m_++m_-+m_0)^2$. The retained state is the charged three-pion channel. | The current side has a published dispersive Khuri--Treiman parameterization and fitted normalization; the ALP side has derived Born and strict vector-pole pieces, but no complete physical complex production amplitude. The full electromagnetic current is retained on the ALP-current side. | Use $S=1+iT$, the recorded epsilon convention, relativistic phase space, and the conditional three-body discontinuity and spin-one weight. Keep $\operatorname{Disc}F$ unless Schwarz reflection conditions justify $2i\operatorname{Im}F$. The conventional strict-pole crossing factor is $\zeta_i=-1$; the physical decay normalization phase relative to production remains unresolved. | BaBar 2021 low-energy data include separate statistical and residual-systematic covariance matrices; the higher-energy table has no supplied covariance, and no cross-region covariance is available. Hoferichter:2025lcz reports marginal fit errors and model variants, not complete correlated replicas; published $V_i\to3\pi$ upper/lower curves are sensitivity bounds, not covariance replicas. | **Current data available; numerical ALP spectrum blocked.** Missing items: a numerical complex current partial-wave grid with complete correlated fit/phase variations; a complex $V_i\to3\pi$ grid and physical phase match; the completed source-dependent ALP production amplitude and its variations. See [three-pion discontinuity](alp-three-pion-discontinuity.md), [three-pion phase match](alp-three-pion-pole-phase.md), [numerical and phase input status](alp-three-pion-numerical-input-status.md), [Hoferichter:2023bjm extract](../../extracts/Hoferichter:2023bjm/three-pion-amplitude-and-normalization.md), [Hoferichter:2025lcz extract](../../extracts/Hoferichter:2025lcz/gamma-star-to-three-pion-fit.md), and [BaBar:2021cde data audit](../../extracts/BABAR:2021cde/three-pion-cross-section-and-covariance.md). |
| Charged kaon cut $K^+K^-$ | Independent no-mixing ALP sources; $s\geq4M_{K^+}^2$. | $\langle K^+K^-|J^\mu_{\rm em}|0\rangle=(p_1-p_2)^\mu F_{K^+}^V(s)$, with $F_{K^+}^V(0)=1$. The retained virtual-current production amplitude is a symbolic odd-resonance tree result; no direct complex $K^+K^-$ current grid is present in the repository inputs. | Conditional distinct-particle cut coefficient is $s(1-4M_{K^+}^2/s)^{3/2}/(96\pi)$ times the conjugate current form factor and the reduced production wave. No $1/2!$ factor. The physical relative phase is not specified. | The stored $\phi\to K^+K^-$ and $\phi\to K_LK_S$ branching fractions have a reported correlation of $-0.72$, but are on-shell decay fractions, not a complex continuum form factor or its covariance. No correlated production variations are selected. | **Formal; numerical cut incomplete.** Missing complex timelike $F_{K^+}^V$, completed production amplitude and phase, and correlated spectral variations. See [virtual-photon kaon amplitude](alp-virtual-photon-kaon-amplitude.md) and the [PDG vector-to-two-pseudoscalar input extract](../../extracts/ParticleDataGroup:2024cfk/vector-to-two-pseudoscalar-inputs.md). |
| Neutral kaon cut $K^0\bar K^0$ | Independent no-mixing ALP sources; $s\geq4M_{K^0}^2$, kept distinct from the charged threshold. | $\langle K^0\bar K^0|J^\mu_{\rm em}|0\rangle=(p_1-p_2)^\mu F_{K^0}^V(s)$, with $F_{K^0}^V(0)=0$. The retained virtual-current production amplitude is symbolic and source-resolved; no direct complex neutral-kaon current grid is present in the repository inputs. | Same conditional distinct-particle $1/(96\pi)$ cut convention as the charged channel, with its own threshold and current form factor. No $1/2!$ factor. A three-flavor isoscalar-current half-sum requires consistent current and isospin conventions. | The stored $\phi\to K^+K^-$ and $\phi\to K_LK_S$ fractions constrain on-shell branching fractions, not the neutral-kaon current phase or continuum covariance. No correlated production variations are selected. | **Formal; numerical cut incomplete.** Missing complex timelike $F_{K^0}^V$, completed production amplitude and phase, and correlated spectral variations. See [virtual-photon kaon amplitude](alp-virtual-photon-kaon-amplitude.md) and the [PDG vector-to-two-pseudoscalar input extract](../../extracts/ParticleDataGroup:2024cfk/vector-to-two-pseudoscalar-inputs.md). |
| Pure narrow $\omega$ and $\phi$ poles | Independent no-mixing ALP sources $U_a,Y_a$; stable poles at the selected vector masses in each photon virtuality. There is no multiparticle threshold in the stable-pole term itself. | Signed neutral-vector photon residues and the ordered VJP/VVP numerators use the unnormalized scalar TFF convention. The meson-pole WZW contribution already absorbed in $\widehat c_\gamma$ is not added again. | Stable-pole denominators define the adopted approximation; their delta-function spectral weights are alternative representations of the same pole and must not be added again. The tree residues are signed, but no physical finite-width phase is supplied. | No quantitative error for the narrow-pole approximation, spectral tails, or correlated residue variations is selected. Uncertainty is model-dependent, especially near timelike peaks. | **Conditional model input.** Retain only where the same strength is not represented by a continuum; replace or match in the $\omega/\phi$ cuts, particularly the shared $\phi\to3\pi$ and $K\bar K$ strength. See [tree TFF](alp-doubly-virtual-tff.md), [mixed-vector strategy](alp-tff-mixed-vector-strategy.md), and [finite-width pole replacement](alp-isoscalar-vector-self-energy.md). |
| Mixed dressed $\rho I$ VVP, $I=\omega,\phi$ | Independent no-mixing ALP sources; the dressed $\rho$ response has the $\pi\pi$ support beginning at $s=4M_\pi^2$, while the other leg is the selected narrow $I$ pole. | Retain the ordered mixed VVP source numerator with $W_\rho^{\rm dress}(s)=R_\rho/[M_{0\rho}^2-s-\Delta_\rho^{\rm self}(s)]$ and $W_I^{\rm pole}$ on the other leg. This response-level term is not the full two-pion double-spectral TFF. | The self-energy convention is $\operatorname{Im}\Delta_\rho^{\rm self}(s+i0)=+\sqrt{s}\,\Gamma_{\rho\to\pi\pi}(s)$; the real part requires a subtracted dispersion relation. The ordered VVP numerator is retained as written in the mixed-term record; it has no independent cut prefactor. | The rho subtraction conditions, numerical self-energy, source-specific match to the two-pion continuum, and their uncertainties are not available. | **Formal; quantitative overlap unresolved.** The mixed term is generally nonzero and must be retained, but its shared rho strength must be matched to the two-pion continuum once. See [mixed rho-isoscalar VVP](alp-tff-mixed-rho-isoscalar-vvp.md), [mixed-vector strategy](alp-tff-mixed-vector-strategy.md), and [isoscalar source decomposition](alp-isoscalar-tff-source-decomposition.md). |

## Three-pion and kaon data audit

### Three-pion current, decay, and ALP production amplitudes

- **Electromagnetic current:** Hoferichter:2023bjm and Hoferichter:2025lcz give the dispersive Khuri--Treiman representation and fitted normalization for $\gamma^*\to3\pi$. Their published tables and fit descriptions do not provide a machine-readable complex $f_1(s,q^2)$ grid or the complete joint fit covariance needed for correlated spectral replicas [@Hoferichter:2023bjm; @Hoferichter:2025lcz].
- **Cross-section data:** the repository contains the extracted BaBar 2021 Born cross sections from $0.62$ to $3.50$ GeV. The 114 bins below $1.10$ GeV retain separate statistical and residual-systematic covariance matrices and a separate correlated-systematic vector. The 72 bins above $1.10$ GeV retain statistical and total systematic errors but have no supplied bin covariance; no cross-region covariance is assigned. These data can constrain a reproduced current-amplitude fit, but are not themselves a complex amplitude [@BABAR:2021cde].
- **$\phi\to3\pi$ data:** the local KLOE file is a fit-parameter table, not binned Dalitz data. It has no full covariance matrix; CMD-2 tabulated Dalitz inputs are not present in the indexed local material [@KLOE:2003kas].
- **Vector decay amplitude:** Niecknig:2012sj supplies the Khuri--Treiman equations, shape plots, and upper/lower sensitivity curves, but not a machine-readable complex $V_i\to3\pi$ Dalitz grid or correlated parameter replicas. Widths and Dalitz intensities do not fix its physical phase relative to the signed ALP production residue [@Niecknig:2012sj].
- **ALP production:** the repository derives the Born current, strict $\omega/\phi$ pole residues, and local VPPP tensor terms. The complete physical source-dependent amplitude still needs the matched contact/subtraction terms, rescattering, a compatible complex decay/current input, and the relative phase and correlated variations. The required conventions are the charge-stripped electromagnetic current, the recorded epsilon tensor and crossing choice, and the common LSZ/S-matrix phase convention; missing fit and phase variations must be recorded rather than replaced by independent marginal errors.

### Charged and neutral kaon current, production, and nearby data

- **Two-body current:** the repository has the definitions and charge normalizations $F_{K^+}^V(0)=1$ and $F_{K^0}^V(0)=0$, but its local data inputs contain no complex timelike $F_{K^+}^V(s)$ or $F_{K^0}^V(s)$ grid, current phase, or spectral covariance. The required convention is $\langle K\bar K|J^\mu_{\rm em}|0\rangle=(p_1-p_2)^\mu F_K^V(s)$, with separate charged and neutral thresholds.
- **On-shell $\phi$ constraints:** the stored PDG 2024 $\phi\to K^+K^-$ and $\phi\to K_LK_S$ branching fractions have a reported $-0.72$ correlation. These constrain on-shell decay fractions and can inform a $\phi$ width input; they do not give either complex continuum form factor or its phase [@ParticleDataGroup:2024cfk].
- **Nearby $K\bar K\pi$ data:** local BaBar inputs include the inclusive $K^+K^-\pi^0$ spectrum (21 bins to $1.78$ GeV, with a $6.5\%$ common systematic and no bin covariance) and an isobar-only $K^*(892)K$ spectrum (statistical errors plus a $15\%$ Dalitz-fit systematic without bin covariance). The other inclusive charge mode has its own $5.1\%$ normalization systematic. These are three-body observables; they do not determine the two-body $K\bar K$ current form factors or the ALP $a\gamma^*\to K\bar K$ amplitude [@BaBar:2007ceh].
- **Related model input:** Qin:2024ulb supplies model form factors for $e^+e^-\to K\bar K\pi$ and fit assumptions for those channels, not charge-resolved $K\bar K$ current form factors [@Qin:2024ulb].
- **ALP production:** the virtual-current $a\gamma^*\to K^+K^-$ and $a\gamma^*\to K^0\bar K^0$ tree amplitudes are derived symbolically, including their different thresholds and the charged-channel $c_4$ term. Physical complex phases, local anomalous terms, rescattering, and correlated variations remain missing. A numerical cut requires these production amplitudes and the two current form factors in the same real-epsilon crossing convention, with the $\phi$ strength coordinated with the three-pion channel.

## Effect of this audit on the virtuality and higher-state assessment

The BaBar three-pion data extend to $3.50$ GeV and improve the possible
normalization of the electromagnetic-current side; the correlated PDG
$\phi\to K\bar K$ branching fractions constrain on-shell $\phi$ decays.
Neither supplies the source-specific ALP production amplitude or its
double-virtual dependence. The available $K\bar K\pi$ measurements are
different final states. Therefore these inputs do not fix the two-pion
cutoff $\Lambda$, bound omitted higher-isoscalar ALP TFF strength, or provide
the continuation over the unbounded lepton-loop virtualities. The task 1.3
coverage conclusion is unchanged.

## Declared virtuality coverage and higher isoscalar strength

The requested ALP mass range is $2m_e<m_a<2m_\tau$, but this does not bound
the photon virtualities in the complete two-photon lepton loop: the loop
integrates over unbounded virtualities. The present low-energy construction
therefore has no common finite numerical virtuality domain. Its coverage is
limited to the formal hadronic cuts and pole terms listed above. In
particular, the two-pion double-spectral formula requires an upper cutoff
$\Lambda$ on its spectral variables, but no value or matching prescription
for $\Lambda$ has been selected. The three-pion and kaon thresholds identify
their cut supports, not an upper virtuality limit or a validated numerical
region. These distinctions follow the [isoscalar TFF plan](alp-isoscalar-tff-plan.md)
and [higher-isoscalar input audit](alp-isoscalar-tff-input-gates.md).

Higher isoscalar strength is omitted from the present low-energy inventory;
there is no provenance-backed bound on its source-specific TFF contribution
and correlated uncertainty, so no quantitative uncertainty band is assigned.
Source-specific virtual
TFF information or matching conditions with uncertainties are still needed
before an effective-pole description could be identified. The existing
[hard-current expansion](alp-hard-current-ope.md) controls only the nearly
symmetric Euclidean hard limit at leading order and retains ALP axial-current
matrix elements that have not been matched. It does not bound omitted
strength over the singly virtual, general unequal-virtuality, or timelike
regions. Consequently the omitted strength prevents a complete isoscalar
TFF and a complete prediction over the full lepton-loop domain; it is a
coverage limitation, not a bounded residual error.

Numerical resonance calibration and fit-driven uncertainty propagation
remain conditional on `alp-uzh.8`; the current nonconverged vector-sector
likelihood is not accepted as a numerical input. This manifest is an
input-status record, not a numerical prediction.
