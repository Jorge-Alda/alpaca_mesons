---
title: Plan and status audit for the low-energy isoscalar ALP TFF
status: mixed
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Hoferichter:2023bjm, Hoferichter:2025lcz, Niecknig:2012sj]
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

# Status summary

| Component | Status | Evidence and next gate |
|---|---|---|
| No-mixing ALP sources and scalar TFF convention | **Complete** for the stated LO source assembly and tensor convention. | The source/LSZ relations and contact bookkeeping are derived in the [LO source assembly](../../rcht/odd-parity/odd-alp-source-lsz-assembly.md); the full retained tree TFF is in [the tree TFF record](alp-doubly-virtual-tff.md). |
| Tree-level omega/phi responses | **Complete** within the retained stable-pole, original-operator truncation. | The source-resolved tree expression sums neutral-vector responses. This is a pole approximation, not a complete low-energy isoscalar spectral function; widths, cuts, additional operators, and short-distance matching are absent. |
| Three-pion dispersive framework | **Partly complete** as a literature method and fit-input audit; no ALP spectral predictor exists. | The Khuri--Treiman equations and fitted current normalization are documented in [the 2023 amplitude digest](../../extracts/Hoferichter:2023bjm/three-pion-amplitude-and-normalization.md) and [the 2025 fit digest](../../extracts/Hoferichter:2025lcz/gamma-star-to-three-pion-fit.md). The framework has not been implemented for the ALP source. |
| Source-specific $a\gamma^*\to3\pi$ amplitude | **LO Born contribution derived and executed; full isoscalar amplitude requires additional derivation and inputs.** | The [LO current record](alp-three-pion-current-born-amplitude.md) derives the gauged $a\Phi^3$ contribution and Ward identity. Its minimal light-quark isoscalar-current part vanishes. The SM current partial wave does not determine the ALP transition amplitude. The existing [three-pion treatment plan](../../rcht/odd-parity/dispersive-three-pion-treatment-plan.md) also notes missing direct $V3\pi$ operators in the original one-pseudoscalar basis. |
| Numerical three-pion spectral density | **Requires additional inputs.** | Needs complex source-dependent ALP production amplitudes, a compatible complex $\gamma^*\to3\pi$ current amplitude, phase conventions, and correlated variations of fit and phase inputs. The cited 2023/2025 records do not supply a complete complex partial-wave grid or correlated replicas. |
| Three-pion contribution to the doubly virtual TFF | **Not started.** | Derive single- and double-spectral pieces, including exchanged photon ordering and subtraction terms; establish the unitarity normalization for the chosen ALP source and tensor convention. A prescription must separate explicit omega/phi poles from continuum strength. |
| $K\bar K$ cut and phi region | **Partly complete:** real-photon odd-resonance $a\to\gamma K\bar K$ amplitudes and neutral-kaon projections exist; the dispersive TFF construction is not started. | The real-photon [ALP photon-pseudoscalar-pair amplitude](../photon-two-pseudoscalar/alp-photon-pp-amplitude.md) is a starting point, not the required virtual-current isoscalar amplitude. Derive the virtual-current extension and isoscalar projection; obtain complex charged/neutral kaon vector form factors, phases, and correlated uncertainties; define how overlapping phi-to-three-pion and phi-to-kaon strength is handled. |
| Higher isoscalar strength | **Not started.** | Select spectral inputs or an effective-pole parameterization only after the low-energy cuts are defined. Constrain it with suitable data and compare parameterizations; report model dependence separately. |
| Short-distance matching | **Not started for the ALP source basis.** | Derive high-virtuality coefficients from source-specific QCD matrix elements and the relevant short-distance expansion. Match the isoscalar and isovector pieces at an explicit scale and vary that scale. Do not use meson asymptotic coefficients as ALP inputs. |
| Two-photon lepton loop | **Universal normalization complete; ALP-specific full-TFF result not started.** | The [loop normalization record](alp-lepton-pair-normalization.md) fixes the universal normalization, but full use requires the completed source-resolved TFF, its cuts and high-energy continuation, and a decision on the separately matched local ALP--lepton term; see the [lepton-pair calculation plan](alp-lepton-pair-calculation-plan.md). |
| Reproducibility and validation | **Complete for supporting tree derivations only; spectral implementation not started.** | The isoscalar spectral derivations, deterministic integration, and covariance propagation remain to be implemented and checked for the supported kinematics. |

# Staged work plan

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

The first LO step gauges the existing [one-ALP $a\Phi^3$ chiral
interaction](../three-pseudoscalar/alp-phi3-lo-chipt-derivation.md), retains
the charged-pion pole graphs and computes the contact term, then checks the
current Ward identity. The executed [SymPy derivation](../../../scripts/derive_alp_virtual_photon_three_pion_born.py)
finds equal shifted source vertices, a zero contact coefficient, and a zero
Ward contraction. The equations and exact scope are recorded in the [LO
current-amplitude record](alp-three-pion-current-born-amplitude.md).

In exact light-quark isospin, decompose the electromagnetic charge matrix as
$Q_{u,d}=\frac16\mathbf1+\frac12\tau_3$. The isoscalar light-quark current
commutes with pure-pion generators, so its minimal gauged LO contribution
vanishes. This is only a statement about that LO current coupling. It does not
exclude the resonance mechanism in which an ALP-vector-photon vertex produces
an isoscalar vector which then decays to three pions. The existing
[$a\to V\gamma$ record](../../rcht/odd-parity/alp-vgamma-amplitude.md) is
for an on-shell vector and a real photon; the literature $V\to3\pi$ record
also lacks the direct contact terms needed for a complete amplitude. These
records do not yet close the off-shell $a\gamma^*\to3\pi$ isoscalar
amplitude. Keep arbitrary diagonal ALP sources until the isospin projection
is explicitly applied.

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

## 2. Build the three-pion production amplitude

Derive $a\gamma^*(y)\to\pi^+\pi^-\pi^0$ for each independent ALP source,
including its Lorentz/tensor decomposition, isospin and charge factors,
subtractions, and left-hand-cut terms. Adapt the Khuri--Treiman treatment of
the three-pion state, but do not substitute the SM
$\gamma^*\to3\pi$ amplitude for the ALP transition amplitude. Obtain or
reproduce a numerical complex solution with correlated variations of the
three-pion fit, pion phase, and required vector-current inputs
[@Hoferichter:2023bjm; @Hoferichter:2025lcz].

**Acceptance:** the amplitude has the declared phase and normalization,
source responses are independently available, its current Ward identity
holds, the required three-pion partial-wave projection is explicit, and
correlated variations can be propagated to the spectral density.

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

# Immediate prerequisites

The next work should not be a numerical isoscalar spectral integral. First
obtain or produce:

1. the explicit source-resolved ALP $a\gamma^*\to3\pi$ amplitude and its
   subtraction/left-hand-cut information;
2. a complex three-pion partial-wave solution with correlated fit and phase
   variations, or the inputs and approval to reproduce its fit and dispersive
   equation;
3. the matched current-amplitude normalization and phase convention needed
   for the three-pion unitarity product;
4. a decision on the target mass/virtuality range, which controls whether a
   narrow phi approximation is adequate and which additional channels are
   open;
5. source-resolved virtual-photon kaon amplitudes and kaon-current form
   factors before claiming a $K\bar K$ spectral contribution;
6. ALP-source short-distance matrix elements before claiming a matched
   high-energy tail.

The first three items are prerequisites to the leading low-energy three-pion
component. Items 4--6 become prerequisites as the calculation expands to the
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
