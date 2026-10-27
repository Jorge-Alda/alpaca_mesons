---
title: Plan and status audit for the low-energy isoscalar ALP TFF
status: mixed
created: 2026-10-25
updated: 2026-10-27
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Hoferichter:2023bjm, Hoferichter:2025lcz, Niecknig:2012sj, Ruiz-Femenia:2003jdx, GomezDumm:2010ra, Ecker:1988te, Ecker:1989yg]
provenance: assumption
---

# Scope and conventions

This is a staged plan and implementation-status audit for the low-energy
isoscalar contribution to the doubly virtual ALP transition form factor (TFF)
needed in the two-photon part of $a\to\ell^+\ell^-$. It complements the
[isovector construction](alp-isovector-tff-spectral-construction.md). Use
the unnormalized scalar TFF, tensor convention, and no-mixing ALP sources in
the [source-resolved tree TFF record](alp-doubly-virtual-tff.md), following
Alda:2026nrd [@Alda:2026nrd]. The dispersive decomposition is a method input
from Messerli:2025rnv, not a source of ALP spectral weights
[@Messerli:2025rnv]. Do not import the eta or eta-prime weights.

For every independent hadronic-scale ALP source, construct its response first
and combine responses only with the source coefficients in the stated
convention. Keep the explicit neutral-vector poles, hadronic cuts, subtraction
terms, and short-distance contribution distinct until their overlap is
specified. In particular, do not add a dispersive omega/phi contribution to
the tree pole expression without a matching prescription.

The specified application is $2m_e<m_a<2m_\tau$, with the full photon
virtualities required by the two-photon lepton loop. Its integration is
unbounded and crosses timelike hadronic singularities. A low-energy or
narrow-$\phi$ approximation alone therefore cannot complete this application;
the high-energy continuation, additional open cuts, and ultraviolet/local-term
matching must also be controlled. This scope is a user-specified requirement,
not an assumption about the validity range of the retained tree expression.

# Status summary

| Component | Status | Evidence and next gate |
|---|---|---|
| No-mixing ALP sources and scalar TFF convention | **Complete** for the stated LO source assembly and tensor convention. | The source/LSZ relations and contact bookkeeping are derived in the [LO source assembly](../../rcht/odd-parity/odd-alp-source-lsz-assembly.md); the full retained tree TFF is in [the tree TFF record](alp-doubly-virtual-tff.md). |
| Tree-level omega/phi responses | **Complete** within the retained stable-pole, original-operator truncation. | The source-resolved tree expression sums neutral-vector responses. This is a pole approximation, not a complete low-energy isoscalar spectral function; widths, cuts, additional operators, and short-distance matching are absent. |
| Three-pion dispersive framework | **Partly complete** as a literature method and fit-input audit; no ALP spectral predictor exists. | The Khuri--Treiman equations and fitted current normalization are documented in [the 2023 amplitude digest](../../extracts/Hoferichter:2023bjm/three-pion-amplitude-and-normalization.md) and [the 2025 fit digest](../../extracts/Hoferichter:2025lcz/gamma-star-to-three-pion-fit.md). The framework has not been implemented for the ALP source. |
| Source-specific $a\gamma^*\to3\pi$ amplitude | **The exact-isospin Born current, its spin-one overlap, and the local $VPPP$ source contribution are derived; the physical isoscalar amplitude remains incomplete.** | The [trace-based Born projection](alp-three-pion-isoscalar-born-projection.md) gives the charged $I=0$ component, and the [spin-one overlap](alp-three-pion-born-spin-one-projection.md) includes the analytic pion-denominator integral. The [contact derivation](alp-isoscalar-vector-three-pion-contact.md) derives the full local endpoint and contracts it with the unconstrained source-resolved ALP vertex, including the complementary tensor-projector term. Its independent $g_i$ remain unmatched. Rescattering and matched complex decay functions remain absent; the SM current partial wave does not determine the ALP transition amplitude. |
| $aV_i\gamma^*$ followed by $V_i\to3\pi$ | **Pole factorization, partial-wave representation, local $VPPP$ tree term, and conventional crossing/LSZ factor derived; numerical decay input and physical relative phase remain unavailable.** | The [partial-wave record](alp-three-pion-vector-partial-wave.md) gives the Khuri--Treiman relation, the [off-shell tree record](alp-isoscalar-vector-three-pion-offshell.md) derives retained sequential rho exchange, and the [contact derivation](alp-isoscalar-vector-three-pion-contact.md) adds the five-operator local endpoint and full tensor propagation. The [pole-phase record](alp-three-pion-pole-phase.md) fixes $\zeta_i=-1$ for the declared routing and decay convention; it does not infer the physical phase from decay rates. The [numerical-input audit](alp-three-pion-numerical-input-status.md) finds no machine-readable complex grid or correlated replicas in the cited records. The $g_i$ require compatible matching or independent input; rescattering and the ideal-mixing single-trace phi limit remain as stated in the linked records. |
| Numerical three-pion spectral density | **Requires additional inputs.** | Needs complex source-dependent ALP production amplitudes, a compatible complex $\gamma^*\to3\pi$ current amplitude, phase conventions, and correlated variations of fit and phase inputs. The cited 2023/2025 records do not supply a complete complex partial-wave grid or correlated replicas. |
| Three-pion contribution to the doubly virtual TFF | **Conditional unitarity and dispersion normalization derived and executed.** | The [discontinuity record](alp-three-pion-discontinuity.md) derives the charged-state phase space, scalar projector, spin-one angular weight, subtracted single/double representation, and stable-pole replacement identity. Physical spectral weights, double discontinuity, subtractions and pole/continuum separation remain open. |
| $K\bar K$ cut and phi region | **Retained virtual-current production and conditional cut normalization derived and executed.** | The [virtual-kaon record](alp-virtual-photon-kaon-amplitude.md) retains the independent virtuality terms, charged-kaon $c_4$, full photon response and crossed $K^*$ exchanges. Complex current form factors, rescattering, correlated uncertainties and a common phi-region treatment are still required. |
| Higher isoscalar strength | **Input requirements audited; numerical model not selected.** | The [input-requirement record](alp-isoscalar-tff-input-gates.md) specifies the low-energy separation and independent source-resolved constraints required before fitting omitted strength. No effective-pole parameters have been invented. |
| Short-distance matching | **Leading local hard-current coefficient derived and executed; source matching remains open.** | The [hard-current OPE record](alp-hard-current-ope.md) derives the tree-level parity-odd coefficient and nearly symmetric deep-Euclidean scalar projection with explicit ALP axial-current matrix elements. Their renormalized source matching, other hard limits, flavor thresholds and continuation are not fixed. |
| Two-photon lepton loop | **Unnormalized linear functional and ultraviolet bookkeeping derived and executed; full numerical result open.** | The [unnormalized-loop record](alp-unnormalized-lepton-loop.md) derives the coherent width, two-photon discontinuity and pointlike logarithm/local subtraction. A completed TFF, its analytic continuation and matched local finite part remain necessary. |
| Reproducibility and validation | **Executed checks now include the trace-based Born $I=0$ spin-one orientation projection.** | Its five focused tests include an independent direct-current angular quadrature and epsilon-tensor sign check. The paired scripts and tests for the other conditional derivations pass. Physical spectral integration, correlated input propagation, high-energy continuation and full-loop numerics are not yet validated. |

# Staged work plan

## Delegated continuation on 2026-10-25

Each open stage is assigned separately and its result is incorporated here
before dependent work proceeds. The current continuation preserves arbitrary
no-mixing sources and uses the supported ALP-mass and photon-virtuality
range fixed above. Existing uncommitted manuscript changes are outside this work.

| Open stage | Assignment and progress |
|---|---|
| 1: source/current decomposition | Coder completed the [retained-tree decomposition](alp-isoscalar-tff-source-decomposition.md); the executed derivation and four tests verify recombination, source linearity, Bose symmetry, Ward identities and the stored stable-pole normalization. Mixed rho--isoscalar VVP selections overlap and must be subtracted once. The total electromagnetic contact remains separate. |
| 2: complex three-pion input and relative phase | The [input/phase audit](alp-three-pion-numerical-input-status.md) finds no machine-readable complex decay grid or correlated replicas in the cited sources. The [strict-pole crossing/LSZ derivation](alp-three-pion-pole-phase.md) fixes the conventional factor $\zeta_i=-1$ for the declared routing, but not the physical overall decay phase from rate data. The assembled [production current](alp-three-pion-production-amplitude.md) now uses that convention. A numerical Khuri--Treiman solution still requires a declared phase-shift input, subtraction prescription, normalization and phase match; correlated variations remain separate. |
| 3: three-pion unitarity and dispersion | The [conditional discontinuity derivation](alp-three-pion-discontinuity.md) fixes the kinematic normalization. The [full-charge Born $I=0$ component](alp-three-pion-isoscalar-born-projection.md), its [spin-one orientation overlap](alp-three-pion-born-spin-one-projection.md), and the [local vector-contact production term](alp-isoscalar-vector-three-pion-contact.md) are derived. The latter keeps independent $g_i$ and includes the regular tensor-projector contraction. Physical coefficient matching, rescattering/current inputs, the double discontinuity, subtraction information and the pole/continuum prescription remain open. |
| 4: virtual kaon production and phi region | Coder completed the [retained virtual-kaon derivation](alp-virtual-photon-kaon-amplitude.md); derivation and four tests pass. Physical complex form factors, completed production, rescattering and common phi-region matching remain open. |
| 4a: common finite-width vector response | The [self-energy record](alp-isoscalar-vector-self-energy.md) derives the absorptive sign, positive spectral density, subtracted self-energy dispersion and stable-pole replacement. It is formal and diagonal; physical decay amplitudes, subtraction conditions, and possible omega--phi mixing remain open. |
| 5: higher strength and short distance | Scientific consultation established the [input requirements](alp-isoscalar-tff-input-gates.md). A second-iteration coder completed the [leading local OPE derivation](alp-hard-current-ope.md); derivation and four tests pass. Independent source matrix elements, other hard limits, higher strength and full matching remain open. |
| 6: lepton loop | Coder completed the [unnormalized-loop derivation](alp-unnormalized-lepton-loop.md); derivation and four tests pass. Numerical evaluation depends on the completed TFF and finite local matching. |
| Reproducibility | Each new derivation must have an executed repository-Python script, a paired atomic record, and explicit limits. |

## Framework adaptation: keep the two three-pion amplitudes distinct

The current amplitude used in the standard three-pion Khuri--Treiman
construction is
$$
\langle 0|J^\mu_{\rm em}|\pi^+\pi^-\pi^0\rangle
=-\epsilon^{\mu\nu\alpha\beta}p_{+\nu}p_{-\alpha}p_{0\beta}
\mathcal F(s,t,u;q^2),
$$
as recorded in the [Hoferichter:2025lcz digest](../../extracts/Hoferichter:2025lcz/gamma-star-to-three-pion-fit.md).
The other side of the three-pion unitarity cut is a different matrix element,
$$
\langle\pi^+\pi^-\pi^0|J^\nu_{\rm em}|a\rangle,
$$
crossed to $a\gamma^*\to3\pi$. It contains an ALP and three pseudoscalars,
has ordinary-parity vector covariants at leading chiral order, and depends on
the additional invariants involving the ALP and current momenta. The SM scalar
$\mathcal F(s,t,u;q^2)$ therefore cannot be converted into this matrix element
by replacing its normalization with ALP couplings. A Khuri--Treiman solution
can describe a selected $V\to3\pi$ subamplitude after a partial-wave
projection; it is not by itself the full five-point ALP production amplitude.

The LO Born step gauges the existing [one-ALP $a\Phi^3$ chiral
interaction](../three-pseudoscalar/alp-phi3-lo-chipt-derivation.md), retains
the charged-pion pole graphs and contact term, and checks the Ward identity.
The executed result is recorded in the [LO current-amplitude
record](alp-three-pion-current-born-amplitude.md). The resonance pole and
decay partial-wave steps are recorded separately in the [spin-one
pole-factorization record](alp-three-pion-vector-pole-factorization.md) and
the [vector partial-wave record](alp-three-pion-vector-partial-wave.md).

In exact light-quark isospin, the isoscalar electromagnetic generator is
$Q^{(0)}_\ell=\frac16\mathbf1_2$. It commutes with the charged-pion
generators, so its minimal gauged LO contribution vanishes. The executed
matrix-trace calculation also finds equal shifted source vertices and a zero
contact coefficient, making the Born current transverse. These are results
for the stated LO chiral contribution; they do not exclude a resonance
mechanism in which an ALP-vector-photon vertex produces an isoscalar vector
which then decays to three pions. The existing
[$a\to V\gamma$ record](../../rcht/odd-parity/alp-vgamma-amplitude.md) is
for an on-shell vector and a real photon; the literature $V\to3\pi$ record
also lacks the direct contact terms needed for a complete amplitude. These
records do not yet close the off-shell $a\gamma^*\to3\pi$ isoscalar
amplitude. Keep arbitrary diagonal ALP sources until the isospin projection
is explicitly applied. Selecting $V_i=\omega,\phi$ on the three-pion leg
selects the isoscalar quantum numbers on that cut; it does not project the
other electromagnetic current onto $Q^{(0)}_\ell$. Retain the full charge
matrix and all neutral intermediate-vector responses, including mixed
isovector-isoscalar trace terms. The zero of the minimal isoscalar-current
Born term alone does not establish zero overlap with the isoscalar three-pion
cut; that requires the explicit three-pion isospin and partial-wave
projection.

## 1. Freeze source basis and isoscalar current decomposition

Use the no-mixing source matrices, the unnormalized TFF, and the common
electromagnetic normalization already defined in the linked tree-TFF and
source/LSZ records. Build an independent response for each source. Explicitly
state which current and photon placements define the isoscalar piece, and how
contact contributions and mixed isovector-isoscalar vector terms are assigned.
The tree-level neutral-vector calculation is the normalization and limiting
check for this stage, not a complete spectral representation.

**Acceptance:** source linearity, current conservation, Bose symmetry under
the two photon labels, and agreement with the existing tree expression in its
stable-pole limit.

**Progress:** the [executed source/current decomposition](alp-isoscalar-tff-source-decomposition.md)
meets these checks for the retained resonance expression. It separates
full-three-flavor isoscalar/isovector photon currents from omega/phi
hadronic-pole selection and derives the mixed VVP overlap explicitly.
The supplied electromagnetic contact fixes the total current pair only;
independent projected contact weights require an external-source completion
or an explicitly declared allocation. No allocation is imposed here.
The [assembled production-current derivation](alp-three-pion-production-amplitude.md)
now combines the gauged LO Born term with the strict omega/phi pole residues
and verifies the photon Ward identity. It does not yet project the Born term
onto the isoscalar three-pion state or supply the complex rescattered decay
functions, so the physical three-pion discontinuity remains open.

## 2. Spin-one pole factorization: derived

The [pole-factorization record](alp-three-pion-vector-pole-factorization.md)
derives the full free antisymmetric-tensor response from its six-component
Hessian and contracts it with the transverse spin-one lift of the on-shell
$V_i\to3\pi$ current. The local tensor-projector term vanishes for that
declared endpoint, and the executed SymPy calculation verifies the three-pion
current transversality, photon Ward identity, vector transversality, and
projector reduction. The source-resolved production residue uses the ordered
VJP/VVP numerators with the ALP mass and the full neutral-vector photon
response.

This result is restricted to the strict $V_i$ pole residue. The relative
crossing/LSZ phase must be matched to the selected decay-current convention.
On-shell $V_i\to3\pi$ data do not determine regular off-shell terms or direct
$V3\pi$ contacts; do not extend the pole numerator into those terms without a
separate source model.

**Completed step:** the [partial-wave record](alp-three-pion-vector-partial-wave.md)
derives the Bose-symmetric Khuri--Treiman representation and identifies the
on-shell decay function multiplying the strict pole residue. The symbolic
pole derivation was executed with approval. This does not provide a numerical
complex decay grid or its correlated variations.

**Next gate:** produce or obtain the complex, normalized $V_i\to3\pi$
solution and correlated phase/fit variations, and match its S-matrix phase to
the production residue. Then derive the three-pion discontinuity using the
full electromagnetic charge matrix and mixed isovector-isoscalar traces on
the ALP-current side; use the SM $\gamma^*\to3\pi$ amplitude only on the
current side of the unitarity cut [@Hoferichter:2023bjm; @Hoferichter:2025lcz].

## 3. Derive the three-pion discontinuity and TFF terms

Apply unitarity to the compatible $\gamma^*\to3\pi$ current amplitude and
source-specific ALP production amplitude. Derive the phase-space factor,
charge sum, and identical-particle convention rather than importing a
two-pion prefactor. Construct all single- and double-spectral terms needed
for the TFF, with both photon orderings and relevant subtraction
contributions. Define the separation from explicit $\omega,\phi$ poles
before adding either representation.

**Acceptance:** the discontinuity reconstructed from the dispersive TFF
matches the unitarity result; the TFF is Bose symmetric; and the pole/continuum
prescription avoids double counting.

**Progress:** the [conditional discontinuity derivation](alp-three-pion-discontinuity.md)
establishes the relativistic three-body measure, absence of a charged-state
factorial, epsilon-tensor scalar projection, and spin-one orientation factor.
Its subtracted Cauchy representation has the expected discontinuity and
photon exchange, conditional on compatible physical spectral inputs and
convergence. It does not derive an unknown double discontinuity from the
single three-pion cut. The stable-pole distribution identity specifies how
a pole must be removed when replaced by a spectral contribution; the
physical finite-width separation remains to be chosen and validated.

## 4. Add the $K\bar K$ cut and resolve the phi region

Derive source-resolved $a\gamma^*\to K^+K^-$ and
$a\gamma^*\to K^0\bar K^0$ amplitudes, along with the corresponding
electromagnetic-current form factors. Include physical thresholds and
declared isospin-breaking effects if the precision requires them. Choose a
dispersive phi-region treatment or a narrow-pole approximation for the
specified mass and virtuality range, then match it to the $K\bar K$ continuum
without overlap. Existing $a\to\gamma K\bar K$ amplitudes supply only a
real-photon starting point.

**Acceptance:** channel normalization, phases, thresholds, and phi
pole/continuum separation are fixed and checked against available data or
source records.

**Progress:** the [virtual-kaon derivation](alp-virtual-photon-kaon-amplitude.md)
extends the original eleven-operator resonance contribution with independent
virtuality terms and the charged-kaon $c_4$ contribution. It derives the
allowed neutral/$K^*$ exchanges, Ward identity, real-photon limit, source
linearity and conditional distinct-kaon cut normalization, with separate
charged and neutral thresholds. The complex kaon form factors, complete
production amplitudes and shared three-pion/kaon phi description are still
needed; these executed tree checks do not close the physical cut.

## 5. Add higher isoscalar strength and short-distance matching

After the low-energy cuts are controlled, add higher isoscalar spectral
strength using declared data or effective poles. Independently derive the
short-distance coefficient for each ALP source. Match the low- and
high-virtuality descriptions at a stated scale and vary the matching scale,
cutoff, and higher-state parameterization to estimate model uncertainty
[@Messerli:2025rnv].

**Acceptance:** the matched TFF has the required source-specific asymptotics;
the low-energy normalization and discontinuities remain unchanged within the
stated matching uncertainty.

**Progress:** the [hard-current OPE derivation](alp-hard-current-ope.md)
calculates both free hard-quark contractions and the gamma-five trace,
fixing the leading parity-odd local Wilson coefficient at leading order in
$\alpha_s$. Its scalar projection retains the explicit vacuum-to-ALP axial
current matrix elements. The nearly symmetric Euclidean limit with fixed
soft ALP momentum is controlled at this order; the local term does not
establish arbitrary unequal-ratio or singly virtual asymptotics. Those
limits, the source-specific renormalized matrix elements, active-flavor
thresholds, higher-state information and timelike continuation remain open.

## 6. Insert the TFF in the lepton loop

Use the unnormalized full TFF in the two-photon loop so source cancellations
at zero virtuality do not produce a singular normalization. Add any local
ALP--lepton contribution only with its scale and matching prescription
stated. Identify additional open hadronic cuts and the ALP-mass range over
which the result is intended to apply, following the [lepton-pair plan](alp-lepton-pair-calculation-plan.md).

**Acceptance:** recover the pointlike-TFF limit and the known two-real-photon
absorptive part, and report separate uncertainties from source matching,
three-pion input, kaon input, pole/continuum separation, higher states, and
short-distance matching.

**Progress:** the [unnormalized-loop derivation](alp-unnormalized-lepton-loop.md)
removes the real-photon denominator before integration, retains source
linearity even when the real-photon value cancels, and derives the coherent
width with a consistently converted local coefficient. It reproduces the
two-real-photon discontinuity and derives the logarithmic ultraviolet
divergence of a constant TFF. Thus a finite pointlike limit requires a
specified renormalization and finite local matching; it is not a universal
finite number. Full spectral-loop numerics and the uncertainty budget remain
open.

# Immediate prerequisites

The next work should not be a numerical isoscalar spectral integral. The
retained source-resolved Born and strict-pole production current is now
available, but it is not yet the full isoscalar production amplitude. First
obtain or produce:

1. the missing physical regular production terms, direct $V3\pi$ inputs,
   and subtraction/left-hand-cut information. The retained Born $I=0$
   spin-one orientation overlap is derived, and the retained affine numerator's
   conditional regular quotient is recorded; neither supplies the missing
   physical terms;
2. a complex three-pion partial-wave solution with correlated fit and phase
   variations, or the inputs and approval to reproduce its fit and dispersive
   equation;
3. the matched current-amplitude normalization and phase convention needed
   for the three-pion unitarity product;
4. the user-specified scope $2m_e<m_a<2m_\tau$ with the full loop virtualities
   is fixed; establish the additional open channels and a high-energy
   continuation over that scope, rather than restricting the answer to a
   narrow phi approximation;
5. source-resolved virtual-photon kaon amplitudes and kaon-current form
   factors before claiming a $K\bar K$ spectral contribution;
6. ALP-source short-distance matrix elements before claiming a matched
   high-energy tail.

Items 1--3 remain prerequisites to the leading low-energy three-pion
discontinuity; the vector decay partial-wave representation is now derived,
but its complex numerical solution and matched phase are still missing.
Items 4--6 become prerequisites as the calculation expands to the
phi/kaon region and a full TFF for the lepton loop.

# Literature boundary

Messerli:2025rnv supplies the decomposition strategy, not ALP spectral
weights or source amplitudes [@Messerli:2025rnv]. Hoferichter:2023bjm and
Hoferichter:2025lcz supply the three-pion current dispersive framework and
fit information, but not a ready complex grid with correlated replicas
[@Hoferichter:2023bjm; @Hoferichter:2025lcz]. Alda:2026nrd supplies the
no-mixing source convention [@Alda:2026nrd]. The present records do not
establish the ALP-specific three-pion amplitude, the virtual-photon kaon
amplitudes, or the ALP-source high-energy coefficient.
