---
title: Odd flavor-breaking and single/double-trace extension with low-energy calibration
status: mixed
created: 2026-10-13
updated: 2026-10-19
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Roig:2013baa, Kaiser:2000gs, ParticleDataGroup:2024cfk]
provenance: assumption
---

# Selected objective

User instruction on 2026-10-13: plan a combination of odd flavor breaking,
single/double interactions and low-energy calibration. Here “single/double”
means single- and double-**flavor-trace** interactions, including the associated
independent singlet structures; it does not mean new resonance multiplets.

Derive a declared one-vector-nonet odd-action extension, calibrate only
low-energy observables within its validated domain, then determine what it
predicts for tree-level ALP targets. Do not claim full ultraviolet matching.
Loops, dynamical mixing, additional resonance multiplets and dispersive
three-pion calibration remain separate work. HLS stays excluded.

## Indexed inputs and preserved choices

Use [mismatch interpretation](vector-matching-mismatch-interpretation.md),
[physical matching and definitions](vector-physical-tff-matching.md),
[chiral-octet matching](vector-chiral-octet-matching.md),
[U(3) extension options](u3-nonet-odd-parity-extension-options.md),
[Roig operator comparison](../../extracts/Roig:2013baa/anomalous-rcht-bases-and-double-trace.md),
[stage-4 amplitudes](vector-stage4-amplitudes.md),
[Alda input conversion](alda-neutral-amplitude-inputs.md),
[fit specification](vector-stage5-joint-fit-specification.md), and
[ALP target audit](alp-target-amplitude-audit.md).

Preserve canonical tensor fields, Alda arXiv-v1 source conventions, the
positive-strange phi phase, independent $M_V$ normalization, and WZW/current
assignment. Retain both mass branches and source-sign hypotheses. Keep
charged/neutral scenarios separate, preferring the near-ideal branch only
when adequacy and truncation considerations are equivalent. The existing
partial-isospin approximation remains labelled rather than upgraded silently.

## Stage A — define and reduce a bounded operator space

Initial execution on 2026-10-13: the
[source-resolved vertex catalogue](odd-breaking-vertex-catalogue.md),
[flavor enumeration](odd-breaking-flavor-enumeration.md) and
[Lorentz/IBP reduction](odd-breaking-lorentz-reduction.md) are saved with
executable exact SymPy results. The bounded block space has 205 coefficient
tensors before physical-source/kinematic projection; these are not 205 fit
parameters. No EOM or OPE prior was applied to that vertex catalogue.
Execution on 2026-10-14: the [complete covariant quotient](odd-higher-covariant-operator-count.md)
now establishes 228 independent operators, including all eleven original
operators, hence 217 additions. The basis, ordered derivatives, cross-block
relations and exact three-flavor certificates are saved. This completes the
operator-count milestone; physical predictor normalization remains in stage B.

**Deliverables:** an atomic operator-basis record, an executable SymPy
reduction/expansion script and a table of dimensions, trace topology,
spurion count, source dependence and target channels for every retained term.

1. Reconstruct the original seven VJP/four VVP operators and the literature
   single/double-trace alternatives before adding terms. Enumerate both
   derivative and mass/source structures. The concrete
   $i\epsilon\langle VV\rangle\langle\chi_-\rangle$ versus
   $i\epsilon\langle\{V,V\}\chi_-\rangle$ pair is an input example,
   not a complete basis or two automatically independent additions to the
   existing $d_i$ basis [@Roig:2013baa].
2. Proposed scope: parity- and charge-conjugation-allowed VJP and VVP
   one-pseudoscalar/source vertices with the baseline derivative ceiling,
   at most one quark-spurion insertion in total, and at most two flavor
   traces. Include zero-spurion independent singlet/double-trace structures
   and one-spurion single/double-trace structures. This is a selected model
   space, not a complete arbitrary-order $U(3)$ action.
3. For derivative seeds, enumerate inequivalent placements of the dynamical
   $\chi_+$ block and trace partitions. Separate its traceless flavor-breaking
   part from singlet rescalings. Do not insert a second spurion into a seed
   already containing $\chi_-$ under the one-spurion ceiling. Keep powers of
   $M_V$ explicit so all fitted coefficients have declared dimensions.
4. Reduce by cyclicity, Hermitian conjugation, parity/C, Schouten identities,
   integration by parts and valid $U(3)$ trace identities. Use leading EOM
   only with an explicit account of induced source/contact terms. Identify
   single-trace terms already represented by $c_i,d_i$ to avoid duplicate
   parameters. Recompute which combinations are independent after restricting
   to the chosen amplitudes; operator independence and data identifiability
   are separate questions.
5. Large-$N_c$/OZI counting records one extra trace as nominal suppression,
   separately from spurion counting. It is not a fixed numerical Gaussian
   prior [@Kaiser:2000gs]. Independent $X_{TC}$ coefficient functions and
   derivative-$X_{TC}$ terms remain outside this initial space and must appear
   in the ALP truncation statement.

**Exit criterion:** a reproducible, nonredundant declared action with sufficient
channel flavor freedom to test low-energy adequacy. Do not select operators
solely because they force a desired ultraviolet cancellation.

## Stage B — derive normalized SM and ALP kernels

Execution on 2026-10-14: the
[neutral independent-source EOM audit](odd-trace-source-eom-audit.md)
derives the required $X_{TC}$ contacts and three relations in its
divergence/mass subsector. The
[leading single/double-trace scalar kernels](odd-trace-scalar-kernels.md)
are normalized, implemented and pass 32 new/regression tests. They add 15
formal scalar EM directions beyond the baseline, not an action count.
The [complete bounded action quotient](odd-higher-covariant-operator-count.md)
has now been established: 217 new independent operators. Normalized SM
cubic kernels, measured combinations, calibration and finite neutral
ALP-source projections are now executed as recorded in stage D. Higher-field,
charged-target and complete physical-predictor directions remain separate
completion gates; the cubic audit does not establish them.

**Deliverables:** action-derived scripts, atomic amplitude records and tested
extensions of the existing implementation kernels.

1. Expand the full covariant blocks before replacing them by vacuum spurions.
   Differentiate independent external sources before setting them to zero.
   Determine all source/contact insertions, including those arising from a
   new $\chi_+$ or singlet structure; a vacuum flavor replacement alone does
   not determine the dynamical ALP insertion.
2. Derive V-to-P-photon, crossed P-to-V-photon, and real/two-virtual photon
   kernels with both photon assignments and the action's matrix response.
   Preserve external canonical tensor normalization and distinguish it from
   the internal pole response. Derive the small-virtuality normalization and
   slope analytically where available.
3. Check Ward and Bose identities, crossing, dimensions, phase/basis
   covariance and recovery of the existing kernels when new coefficients
   vanish. Add flavor-sensitive tests distinguishing derivative versus
   mass double traces and singlet versus octet channels.
4. Recompute sensitivity of $a\to VV,V\gamma,\gamma\gamma,3P$ to all retained
   new combinations, using the independent Alda sources. The old blind-direction
   result is not assumed to survive an extended operator space. Identify any
   new combination invisible in the SM dataset but entering an ALP target.
   Construct normalized predictors before interpreting that sensitivity as
   an observable uncertainty.
5. Re-evaluate hard expansions as labelled diagnostics. A mass-type double
   trace cannot alone supply the missing leading double-virtual constant
   cancellation; derivative structures must be assessed separately. No hard
   cancellation or old affine OPE target enters the low-energy likelihood.

**Exit criterion:** normalized SM predictors and a documented ALP sensitivity
map, with uncalibratable target directions explicitly exposed.

## Stage C — assemble a low-energy data snapshot

**Primary dataset:** mass and electronic-width inputs plus radiative decay
branching fractions/partial widths and $\pi^0,\eta,\eta'\to\gamma\gamma$
normalizations. The [PDG 2024 snapshot](../../extracts/ParticleDataGroup:2024cfk/odd-parity-low-energy-measurements.md)
now supplies twelve summary observations and seven shared total widths,
with archived sources, hashes, pages and units. Its inclusive
$\eta'\to\rho^0\gamma$ entry is deferred; the direct pion partial width
is used once without adding a derived mean-life datum. Use the actual
measured observable rather than converting a branching ratio with a fixed
total width into an allegedly independent partial width.

**Optional low-virtuality information:** only genuine near-origin TFF slopes
or binned measurements within a predeclared domain. Extract their covariance
and normalization policy, and validate the model/domain with nested truncations
and withheld data before adoption. Do not rebrand CLEO/BaBar hard tables as
low-energy input, or choose a $Q^2$ cutoff after inspecting fit residuals.
Without adequate low-virtuality data, the primary dataset stands on its own.

Audit measurement overlap and correlated total widths, normalization errors,
and mass/electronic/radiative summary correlations. Keep standardized shared
nuisances once and separate from diagonal errors. Where covariance is absent,
state conditional independence explicitly. Alda neutral-input covariance is
an input-propagation task, not an extra fit datum. Missing correlations must
be labelled; source-breaking/model uncertainty stays separate from input error.

Resolved resonance line shapes, hard spacelike TFFs and three-pion amplitude
decompositions are excluded from this first calibration. Broad-vector NWA
limitations are labelled or the affected datum is deferred until its measured
spectral observable can be predicted. No dispersive three-pion completion is
required to execute this narrowed dataset.

**Exit criterion:** frozen, nonduplicated observations and explicit covariance/
independence assumptions, with a model-domain justification for every datum.

## Stage D — calibrate only identifiable low-energy combinations

The [fixed-scenario measured-width Jacobians](odd-low-energy-measurement-jacobians.md)
are now executed for all 64 retained branches/source scenarios. Their
conditional soft ranks are $(7,12,12)$ and promoted ranks $(12,12,12)$.
The controls excluding broad-rho and constrained-fit K-star summaries reach
$(6,7,7)$ and $(7,7,7)$. The twelve-summary model is already saturated at
first corrections; no in-sample residual gain establishes power-counting
adequacy. The full action has 216 locally unmeasured directions in this
snapshot. Measurement, nuisance and visible-combination maps are saved.
The [conditional calibration](odd-low-energy-calibration.md) now executes
all 1792 scenario/model/control fits. All branch/source and image
certificates pass; 500 interior unsaturated fits have local observable
covariance and joint held-out null diagnostics; twelve source-boundary
fits omit symmetric uncertainties. The best achieved soft-leading
primary total $\chi^2$ is 5.5675 for six nominal degrees of freedom.
Separate double-trace and single-trace correction models already saturate
every odd dataset control. This is a completed conditional diagnostic
run, not selection of a validated preferred truncation or completion of
the full stage-D exit criterion. The [finite sign/profile continuation](odd-low-energy-profiles.md)
now performs 5888 initial fits and 112 physical profiles. A continuation adds
442 finite fits, validates all 29 previously unresolved source endpoints using
direct source coordinates and an independent arithmetic certificate, and
re-anchors four lower references. Targeted restarts and refreshed source grids
find no further decrease above $10^{-4}$ for pairs 90 and 218; finite searches
do not prove global minima. The [finite neutral target audit](odd-alp-measured-null-audit.md)
certifies 448 fixed-nuisance cases; correction-model measured nulls affect
the sampled two-boson resonance targets. Joint nuisance/target projections
and physical-predictor completion remain open.

1. Build the measurement-weighted Jacobian for each source/branch scenario;
   recompute rank after source/input nuisance treatment. The old 7/8 amplitude
   ranks do not determine identifiability of the enlarged action.
2. Fit nested models: original odd vertices, independent double-trace/singlet
   vertices, one-spurion single-trace breaking, and their combined extension.
   Keep coefficient topologies distinguishable; select between truncations
   using adequacy, held-out channels, conditioning and correction sizes,
   rather than rank gain alone. Fit identifiable combinations; retain null
   directions explicitly instead of fixing invisible couplings to measured
   values.
3. Jointly float the declared even mass/source and visible odd coefficients.
   Stage-3 estimates supply starts; their data are counted only once, not also
   as independent fitted-parameter priors. Freeze covariance policy before
   minimization and use the existing SciPy quadratic-loss infrastructure.
4. No inherited OPE Gaussian prior, implicit OZI/naturalness penalty or
   ultraviolet cancellation constraint is included. If the data cannot
   constrain the selected extension, reduce to data-identifiable combinations
   and report missing information; do not conceal it with an undeclared prior.
5. Report convergence, data/nuisance contributions, identifiable combinations,
   profiles or local uncertainties as justified, active bounds, residuals,
   and nested-model/input sensitivity separately for every retained scenario.
   A fit must not pool sign or mass branches into an unlabelled central value.

**Exit criterion:** adequate low-energy predictions with transparent
identifiability and uncertainty limits. If no truncation succeeds, record the
failure and reassess the action/domain before producing an ALP band.

## Stage E — propagate to tree-level ALP predictions

Propagate the calibrated combinations through the complete normalized ALP
predictors from stage B, retaining Alda's pointwise mixing-validity conditions.
Reassess which target observables depend on unconstrained new singlet/source
directions. If one does, its uncertainty remains open until an explicit
assumption or additional source-sensitive information is supplied; no unique
prediction is manufactured by zeroing that direction.

Separate experimental/input uncertainty from odd-operator truncation,
source-sign/mass-branch alternatives and the remaining partial-isospin/NWA
limitations. Publish a stated low-energy domain and label extrapolation beyond
it. Hard TFF agreement and full singlet/finite-mass OPE matching remain future
objectives. A later resonance extension must reproduce the established
low-energy amplitudes without double counting local resonance-saturated terms.

## Next executable task

The user simplified the working fit on 2026-10-19. Use the
[focused original-operator preliminary result](odd-preliminary-fit.md) with
one fixed mass/data setup, separate source hypotheses, and targeted input
and data-removal sensitivities. Reuse existing validated minima and do not
repeat the full scenario grid. Resolve only material issues exposed by this
compact result before further fit optimization. The remaining ALP predictor
and uncertainty work below is a separate task from obtaining this
conditional low-energy working fit.

The [power-counting continuation](odd-power-counted-channel-maps.md) has now
combined all leading/higher neutral cubic parents with source/Fourier/canonical
normalization, both photon placements and the retained LO topological contact.
It recovers the original scalar $C,H$ kernels and verifies every action
relation only after source-EOM compatibility. The generic neutral source-family
first-correction ranks are $(12,19,6)$ under soft counting and $(15,23,7)$
after hard derivative-parent promotion. These are not measured Jacobian ranks.
The [derivative prescription](odd-resonance-derivative-counting.md) defines
the kinematic distinction explicitly. The
[source/LSZ assembly](odd-alp-source-lsz-assembly.md),
[charged higher SM kernels](odd-charged-cubic-kernels.md) and
[signed photon/action matching](../../alp/photon/alp-photon-odd-action-phase-matching.md)
are now derived. The [measurement Jacobians](odd-low-energy-measurement-jacobians.md)
and [bounded numerical SM kernels](odd-bounded-sm-numerical-kernels.md) are
now implemented. The [conditional calibration](odd-low-energy-calibration.md)
has executed their joint low-energy fit and local uncertainty/held-out
identifiability diagnostics in all 64 scenarios. The [finite profiles](odd-low-energy-profiles.md)
and [neutral target audit](odd-alp-measured-null-audit.md) are now executed.
The later source continuation validates all 29 previously unresolved endpoints
and re-anchors four references, as recorded in the profile result. The working
fit now follows the focused task above. Full joint nuisance/target null
projections and missing charged/local/physical predictors remain prerequisites
for ALP uncertainty propagation. No unique ALP band is established. The earlier intermediate statements below refer
to their original execution stages, not completion of physical predictors.

The [power-counting proposal](odd-action-power-counting-proposal.md), saved
at the user's request on 2026-10-14, organizes leading single traces, leading
double-trace corrections and single-trace spurion/derivative corrections,
while deferring their mixed double-trace corrections. It is a proposed
truncation with unestablished numerical resonance-domain adequacy; its
parent filtration, source-family ranks and conditional measured-data ranks
are now computed. Numerical correction sizes and adequacy remain open. Apply any selected hierarchy
to parent columns before amplitude reduction, with no implicit coefficient
prior or claim of complete next-order resonance accuracy.

The [complete covariant quotient](odd-higher-covariant-operator-count.md)
contains 228 independent operators: 40 leading and 188 higher, with original
rank eleven and 217 additions. It uses the user's 2026-10-14 convention:
preserve explicit vectors and the fixed derivative ceiling; reduce chiral
identities and source-preserving LO Goldstone EOM. All induced contacts and
higher-field remainders are retained. No vector EOM or discarded two-spurion
image is used to force this count.

Map the saved basis into normalized SM and ALP predictors. Extend the existing
leading trace kernels with the one-spurion derivative and higher-derivative
mass/source blocks, including both photon assignments and fixed topological
contacts. Determine which combinations enter radiative decays, real/virtual
TFFs and each ALP target, and which remain invisible. Use those projections
to select identifiable low-energy fit combinations rather than fitting all
217 new action coefficients.

Execution on 2026-10-14: the [higher two-boson projection](odd-higher-two-boson-combinations.md)
reduces 181 derivative/mass-source representatives to 66 abstract stripped
combinations and 51 neutral combinations. The seven higher $f_-$ representatives
vanish for electromagnetic-compatible diagonal scalar sources at this cubic
order. The 51-combination map keeps independent ALP derivative/mass weights
and symbolic quark masses; it is not a normalized predictor or fit rank.
The joint normalized source/channel map now performs that combination;
fixed-scenario measured sensitivities are completed in the linked Jacobian
record. Conditional calibration diagnostics are executed; full validation
and ALP propagation remain open.

The count and basis, bounded SM kernels and conditional low-energy
calibration diagnostics are complete. Full calibration validation,
ALP-target identifiability and physical ALP propagation remain open.
