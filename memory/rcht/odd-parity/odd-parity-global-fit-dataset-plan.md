---
title: Minimal global-fit dataset plan for odd-parity RChT
status: mixed
created: 2026-09-30
updated: 2026-10-04
bibtex_keys: [Alda:2026bar, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, Achasov:1999wp, CLEO:1997fho, BaBar:2009rrj, KLOE:2003kas, Akhmetshin:2006sc, BaBar:2004ytv, Niecknig:2012sj]
provenance: assumption
---

# Recommended first release

This is a data-selection recommendation, not a fit result. It uses the single-trace, no-HLS odd-parity RChT baseline already recorded in `odd-parity-rcht-implementation-decisions.md`; the fitted objects are identifiable combinations of $c_i,d_i$, not each operator coefficient separately [@Alda:2026bar; @Ruiz-Femenia:2003jdx].

## Tier 0: staged target and executable first-release subset

| tier | source/data product | observable to ingest | direct coupling information | covariance/readiness |
|---|---|---|---|---|
| 0A | PDG 2024 light-meson listings | $V\to P\gamma$ branching fractions/partial widths for $\rho,\omega,\phi\to\pi,\eta,\eta'\,\gamma$; $P\to\gamma\gamma$ widths | absolute $VJP+VVP+$photon-vector combinations and anomaly normalization | only $\mathcal B(\omega\to\pi^0\gamma)$ is connected; use PDG correlations where supplied and the user-authorized diagonal approximation otherwise [@ParticleDataGroup:2024cfk] |
| 0B | CLEO, `CLEO:1997fho`, Tables I--X | space-like $Q^2F_{P\gamma^*\gamma}(Q^2)$ for $P=\pi^0,\eta,\eta'$ | normalization convention plus high-$Q^2$ vector-exchange combinations | published statistical/systematic components; no full bin covariance in the source archive; use diagonal errors only unless an external release is located [@CLEO:1997fho] |
| 0C | BaBar, `BaBar:2009rrj`, Table 11 | binned $Q^2F_{\pi^0\gamma^*\gamma}(Q^2)$, $4$--$40\,\mathrm{GeV}^2$ | $\pi^0$ TFF normalization and slope/curvature leverage | table gives statistical and systematic errors; no covariance matrix in the source archive; diagonal only [@BaBar:2009rrj] |
| 0D | KLOE/CMD-2 and BaBar $3\pi$ spectra, `KLOE:2003kas`, `Akhmetshin:2006sc`, `BaBar:2004ytv` | $\phi\to\pi^+\pi^-\pi^0$ Dalitz density and $e^+e^-\to3\pi$ Born cross section | $VVP$ three-pseudoscalar shape, direct term, $\rho$ exchange, and timelike $\omega/\phi$ propagator | KLOE paper gives fit errors and binning but not a covariance matrix; BaBar paper tabulates spectrum uncertainties but no covariance; diagonal first release [@KLOE:2003kas; @Akhmetshin:2006sc; @BaBar:2004ytv] |

Tier 0 is a staged target, not a claim that every listed table belongs in the executable first release. The current first-release builder contains only the validated PDG $\omega\to\pi^0\gamma$ branching fraction and the fixed OPE constraint. The CLEO, BaBar, and KLOE tables have checksum-validated loaders but are not observations until their physical predictors exist. No digitized plots are used and no complete world-average claim is made.

## Tier 1: precision expansion

1. Add all PDG radiative averages, retaining experiment labels (CMD-2/SND/KLOE/BESIII/BaBar) as provenance rather than treating the average as independent measurements [@ParticleDataGroup:2024cfk; @Achasov:1999wp].
2. Add low-$Q^2$ CELLO/CLEO and modern $\eta,\eta'$ TFF points. The exact machine-readable release and any bin-to-bin covariance must be verified before inclusion; do not reconstruct it from plots [@CLEO:1997fho].
3. Add the KLOE and CMD-2 $\phi\to3\pi$ Dalitz samples jointly. The dispersive comparison uses approximately $1.8\times10^3$ binned KLOE points and nearly $2\times10^6$ events, but this is not a substitute for an experimental covariance matrix [@Niecknig:2012sj].
4. Add the $e^+e^-\to3\pi$ line shape from BaBar ISR, including the quoted normalization/luminosity systematics as nuisance parameters rather than correlated point errors [@BaBar:2004ytv].

## Tier 2: optional channels

Add direct $e^+e^-\to P\gamma$ and $e^+e^-\to VP$ cross sections only after the experiment's machine-readable tables are obtained. Their observables are the Born cross section versus $s$, resonance interference phases, and channel-dependent products $\Gamma(V\to e^+e^-)\Gamma(V\to P\gamma)$ or $\Gamma(V\to VP)$. No such table/covariance is currently cached here; the first release must not fabricate it from plot pixels.

## Fit conventions and blockers

Use PDG correlations when supplied. The user authorizes a documented diagonal Gaussian approximation to published errors where covariance is absent, with the mean absolute lower/upper error for asymmetric entries. Identified normalization/luminosity sources must be explicit nuisances and excluded from the diagonal covariance. Known correlated systematics are likewise explicit nuisances: the BaBar three-pion mass-bin systematic uses its tabulated additive response vector, while only the published statistical error is diagonal. This is a data approximation, not a claim of independence. The selected fixed OPE prior is recorded in [the fit-infrastructure record](correlated-odd-parity-fit-infrastructure.md); additional isospin/OZI-breaking and line-shape uncertainties still require explicit models. HLS is not used. The [validated data implementation](odd-parity-data-implementation.md) connects only the omega branching datum. Physical predictors for $P\to\gamma\gamma$/TFF, $e^+e^-\to P\gamma/VP$, and $\omega/\phi\to3\pi$ remain missing, as do covariance matrices for several targets. Other radiative channels must not be generated by automatic flavor substitution because Ruiz--Femenia et al. warn that $K^*$ and $\phi$ require flavor-resolved matching [@Ruiz-Femenia:2003jdx]. No physical fit has been run: the one-datum-plus-OPE problem is rank deficient and the remaining loaded data are not predictor-ready.

## Structural finding: multi-body versus two-body final states

Both multi-body datasets in Tier 0 are now blocked for one shared reason, and adding an observable adds only a rank direction, not information ([the saturation criterion](radiative-channel-saturation-probe.md)):

- The BaBar $3\pi$ spectrum has no mass window in which the single-multiplet baseline can be isolated: $\omega',\omega''$ own at least 72.6% of the coherent sum at all 78 bins ([separability record](three-pion-resonance-separability.md), [blocker](babar-3pi-spectrum-fit-blocker.md)).
- The held KLOE $\phi\to3\pi$ material is a four-number amplitude decomposition that a tree-level amplitude mispredicts, because crossed-channel three-body cuts shift and attenuate the $\rho$ band ([blocker](kloe-phi-three-pion-fit-blocker.md)).
- The already-connected $\mathcal B(\omega\to\pi^0\gamma)$ adds rank one but has normalized residual $\approx0.045$, i.e. it is saturated.

Tier 0A radiative widths are therefore the only remaining route for the first release, because two-body final states are where the tree-level baseline is best controlled. Whether they are also saturated is untested and is the subject of the [saturation probe](radiative-channel-saturation-probe.md). If they fail, the multi-body dispersive treatment ([plan](dispersive-three-pion-treatment-plan.md)) becomes the only route.
