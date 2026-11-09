---
title: Synthesis of diagnostics for the joint vector-fit inadequacy
status: mixed
created: 2026-11-09
updated: 2026-11-09
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, ParticleDataGroup:2024cfk, BaBar:2007ceh, Qin:2024ulb]
provenance: numerical-computation
---

# Synthesis of diagnostics for the joint vector-fit inadequacy

## Inputs and scope

This record combines three diagnostics at the accepted endpoint
`basin_probes/even_shared_widths`. The starting fit and its block contributions
are documented in the [goodness-of-fit record](joint-vector-fit-goodness-of-fit.md).
The inputs, assumptions, calculations and machine-readable outputs are given in
the [K-star width diagnostic](joint-vector-kstar-width-diagnostic.md), the
[$\phi\to K\bar K$ channel diagnostic](joint-vector-phi-kaon-residual-diagnostic.md),
and the [$\eta'$ and mass exclusion refits](joint-vector-eta-prime-and-mass-exclusion-diagnostics.md).
Their ledgers are `data/vector_sector/joint_vector_fit_adequacy.json` and
`data/vector_sector/joint_vector_diagnostic_exclusion_refits.json`. The fit
model conventions follow [@Alda:2026nrd; @Ecker:1988te; @Ecker:1989yg], and the
width and branching-fraction measurements follow [@ParticleDataGroup:2024cfk].

The endpoint diagnostics hold the fitted coordinates fixed. Exclusion refits
are conditional fits to the retained rows and do not demonstrate that the
full likelihood is adequate. The nominal block $p$-values use the block row
counts as degrees of freedom, without subtracting shared fitted parameters.

## Results

### $K^*$ widths

At the accepted endpoint, the propagator on-shell width is $44.9454$ MeV. The
tree-level widths using the physical charged and neutral parent masses are
$44.7012$ MeV and $44.6328$ MeV. The measured widths are $51.4\pm0.8$ MeV and
$47.3\pm0.5$ MeV. If both measurements are described by one freely chosen
common width, its best value is $48.4517$ MeV and the irreducible minimum is
$\chi^2=18.8876$ for one nominal degree of freedom. The two width errors are
treated as uncorrelated because their cross-covariance is unavailable.

This is evidence for prioritizing candidate 1 in the
[candidate-causes record](joint-vector-fit-inadequacy-candidate-causes.md):
the baseline shared-width treatment cannot reproduce both central values even
after optimizing the common width. The close tree-level and propagator widths
do not support a large shift from the on-shell prescription alone. The
irreducible discrepancy is a property of the two measured rows under the
common-width hypothesis; it does not by itself establish which extension is
physically justified.

### $\phi\to K\bar K$ channels

The fixed-endpoint prediction-minus-measurement differences for
$(K^+K^-,K^0\bar K^0)$ are $(+0.0174953,-0.0034877)$. Because the measured
fractions have correlation $-0.72$, the symmetric covariance-whitened
components $(4.105707,1.010785)$ mix the two channels; they are not independent
channel pulls. The squared norm reproduces the block value
$\chi^2=17.878514$. The residual in the sum is $+4.01035$ measured standard
deviations, while the charged-to-neutral ratio difference is $+2.27495$
standard deviations using first-order error propagation.

The larger sum displacement gives some support to candidate 2, SU(3) breaking
in the vector-pseudoscalar-pseudoscalar coupling, which changes the overall
$\phi$ kaon rate. The ratio displacement also leaves candidate 4, even-sector
isospin breaking, worth testing. Since the ratio uncertainty is linearized
and correlated with the sum, these two comparisons are diagnostic directions,
not independent significance tests. A refit is required to separate the
candidate effects.

### $\eta'$ and mass exclusions

Removing the two retained $\eta'$ radiative rows,
$\eta'\to\omega\gamma$ and $\eta'\to\gamma\gamma$, gives a successful
optimizer termination (status 2, `ftol`, 424 evaluations). The fitted
$\Gamma_{\eta'}$ is $0.2300$ MeV. The retained odd-block nominal $p$-value
becomes $0.469$, while the shared-width block remains at
$1.75\times10^{-18}$. The $\rho^0\to\eta\gamma$ pull moves from $-3.266$ to
$+0.187$. This supports candidate 5 as an explanation for part of the odd
block tension, but the shared-width problem remains and the full-data model
has not been repaired by excluding observations.

Removing the mass-and-leptonic block also terminates successfully (status 2,
`ftol`, 640 evaluations). The retained odd-block nominal $p$-value is
$0.999999925$, while the shared-width block remains poor at $0.00232$. On
re-evaluation of the excluded mass rows, the $\rho^+$, $K^{*+}$, $\omega$ and
$\phi$ mass pulls are $+84.16$, $+8.29$, $+408.41$ and $-1158.51$. This
conditional result shows large movement when the mass information is removed;
it does not indicate that the mass measurements themselves should be dropped.

## Candidate assessment and next model tests

The diagnostics favour testing the common $K^*$ width/mass treatment first
(candidate 1). The $\phi$ sum residual gives a secondary reason to test the
VPP coupling ratio (candidate 2); its channel-ratio component also motivates
the later isospin-breaking check (candidate 4). Excluding $\eta'$ radiative
rows relieves the odd-block discrepancy (candidate 5), but not the shared
width tension. Removing the mass block creates extreme conditional mass
pulls and does not point to a viable removal of that block. These findings
rank tests; they do not establish an adequate extended model or a unique
physical cause.

## Reproduction and verification

The component calculations and focused unit-test commands are recorded in
each linked diagnostic record. The exclusion-refit ledger records optimizer
status, high-order block summaries and all signed pull changes. Its
reproduction command is:

```text
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/diagnose_joint_vector_exclusion_refits.py --output data/vector_sector/joint_vector_diagnostic_exclusion_refits.json
```
