---
title: Calculation plan for gluon- and quark-coupled ALPs decaying to lepton pairs
status: mixed
created: 2026-10-23
updated: 2026-10-23
bibtex_keys: [Messerli:2025rnv, Alda:2026nrd, Bauer:2020jbp]
provenance: assumption
---

# Scope

This is a staged calculation plan, not a derived ALP width. The target is the
leading two-photon contribution to $a\to\ell^+\ell^-$ for an ALP whose UV
couplings are only to gluons and/or quarks. Use Alda:2026nrd's no-mixing basis
and the project's positive Bauer photon-operator normalization. The input is
the doubly virtual ALP transition form factor (TFF); the final state is
$\ell^+\ell^-$.

The method reference is Messerli:2025rnv, which writes the two-photon
contribution in terms of a reduced amplitude and normalized TFF, separates
low-energy spectral contributions from a high-energy asymptotic contribution,
and evaluates the latter with finite pseudoscalar-mass effects. Its dispersive
construction also resolves imaginary parts from hadronic intermediate states
beyond the two-photon cut. Its $\eta^{(\prime)}$ spectral weights,
asymptotic coefficients, and numerical inputs cannot simply be reused for an
ALP. The method can be adapted if the library builds the corresponding
source-dependent hadronic amplitudes and short-distance response.

## Library interface and source dependence

The intended deliverable is a library function whose runtime physics inputs
are the ALP mass and its gluon and quark couplings, together with the lepton
species. These remain free parameters; no single benchmark ALP is selected.
The coupling scale and convention must be explicit at the function boundary.
Unless UV running is added as a separate library component, the supplied
couplings are inputs at the scale expected by the source matching. The
library supplies a declared hadronic and short-distance model. It converts
the supplied couplings to the project's hadronic-scale, no-mixing source
convention, constructs the ALP TFF from that source, evaluates the lepton
loop, and returns the width or branching fraction and a separated uncertainty
estimate. Any local ALP--lepton coupling is a separate input or matching
prescription and must interfere coherently with the two-photon amplitude.

Because the leading TFF is linear in the ALP source coefficients, implement
source responses for a documented independent basis of hadronic-scale quark and gluon
operators, then combine them with the supplied couplings. Include the
source-dependent contact term in the same combination. Keep the loop in terms
of the unnormalized TFF whenever possible: normalization by
$F_{a\gamma\gamma}=F_{a\gamma^*\gamma^*}(0,0)$ is singular for coupling
choices where source contributions cancel at the real-photon point. The
normalized Messerli:2025rnv formula remains useful for derivation and checks
when this value is nonzero.

## Spectral construction to adapt from Messerli:2025rnv

Messerli:2025rnv decomposes the normalized $\eta^{(\prime)}$ TFF into
isovector, isoscalar, effective-pole, and asymptotic pieces. The source
extract [records the reduced loop and TFF conventions](../../extracts/Messerli:2025rnv/dilepton-loop-formalism.md).
For an ALP, keep this decomposition as a starting strategy, while deriving
the weights and amplitudes for each source basis element:

- **Low-energy isovector part:** use a symmetric double-spectral
  representation in the two photon virtualities. In Messerli:2025rnv, its
  density is built from the pion vector form factor and the
  $\eta^{(\prime)}\to\pi\pi\gamma^*$ amplitude, with the two-pion cut
  treated dispersively. For an ALP, replace that decay amplitude with the
  source-dependent $a\to\pi\pi\gamma^*$ amplitude derived from the ALP
  sources. This is a promising first hadronic cut, not a complete spectrum.
- **Low-energy isoscalar part:** begin with narrow $\omega$ and $\phi$ poles,
  using source-dependent $a\to V\gamma$ residues and the electromagnetic
  vector residues. Add dispersive finite-width propagators if the desired
  precision or mass range requires them. Messerli:2025rnv found the narrow
  treatment adequate for its $\eta^{(\prime)}$ isoscalar contribution; that
  accuracy conclusion must be checked for the ALP model.
- **Higher resonances:** use one or more effective poles to represent
  spectral strength not captured by the low-energy states. Constrain their
  parameters with available single-virtual data, normalization conditions,
  and matching constraints. Compare alternative pole parameterizations to
  estimate model dependence, as Messerli:2025rnv did.
- **High virtualities:** derive the source-basis asymptotic coefficients from
  the appropriate QCD matrix elements and operator product/light-cone
  expansion. Join this tail to the low-energy representation at a declared
  matching scale, vary that scale, and retain finite-$m_a$ corrections in the
  loop kernel when needed. Do not import the $\eta$ or $\eta'$ coefficient.
- **Additional open channels:** include further hadronic cuts when their
  source amplitudes are available and their thresholds are open. This is
  needed for the full imaginary part and may invalidate a lower bound based
  only on the two-photon cut. State the supported ALP mass range and the
  omitted-channel uncertainty.

Thus a Messerli-like dispersive library is feasible in stages. A first
version can use the computed tree-level ALP TFF and explicit pole terms, with
the high-energy coefficient exposed as a model parameter and uncertainty.
The improved version should replace those assumptions with source-resolved
spectral densities and independently derived QCD matching. The two versions
must be labelled separately; the pole model alone does not establish a
controlled dispersive prediction.

## Existing inputs and conventions

- [Alda:2026nrd's no-mixing prescription](../../extracts/Alda:2026nrd/appendix-a-expansion.md)
  fixes the transformed quark sources for a one-ALP amplitude. In particular,
  $\widehat c_q^A=-\mathbb M_a^{-1}(C)$,
  $\widehat c_q^P=-c_q^{AP}-\mathbb M_a^{-1}(C)$,
  $\widehat c_q^V=c_q^{VS}$, and $\widehat c_q^S=0$. For a purely gluonic
  UV source, use the source-dependent $C$ from that prescription; do not
  replace it with a meson mixing angle without deriving the matching.
- The [ALP photon invariant](../photon/alp-photon-invariant-from-mixing.md)
  already includes the leading WZW/contact contribution in its specified
  approximation. The [vector-conversion note](../photon/alp-photon-vector-mixing.md)
  supplies a real-photon resonance result only. Neither by itself determines
  the doubly virtual TFF required in the loop.
- Use the photon operator normalization
  $\mathcal L\supset (\alpha_{\rm em}/4\pi f_a)
  C_{a\gamma\gamma}aF_{\mu\nu}\widetilde F^{\mu\nu}$ and the canonical
  $\epsilon^{0123}=+1$ convention recorded in the photon branch. Keep explicit
  whether the $e^2$, $f_a^{-1}$, and $C_{a\gamma\gamma}$ factors are inside
  or outside the scalar TFF and reduced loop amplitude.

## Ordered calculation

## Progress recorded on 2026-10-23

The tree-level, source-resolved neutral-vector TFF in the no-mixing basis is
derived in [the doubly virtual TFF record](alp-doubly-virtual-tff.md). The
leading-QED width normalization and universal two-real-photon discontinuity
are derived in [the loop normalization record](alp-lepton-pair-normalization.md)
and its executable script. This establishes the normalization for the next
loop calculation, but does not provide the dispersive real part or full width.
The ALP mass and couplings are free runtime inputs by design. The library's
source-response basis, vector residues/pole prescription, ALP-specific
short-distance coefficient and matching scheme, hadronic cuts, and local
ALP--lepton matching prescription remain to be specified or derived. The next
calculation should implement a clearly labelled pole-model prototype, then
build the source-dependent two-pion spectral component and its QCD matching;
the two-photon cut alone must not be used as the full result.

1. **Construct the ALP-to-two-photon TFF.** Define
   $\mathcal F_{a\gamma^*\gamma^*}(q_1^2,q_2^2)$ by the same signed tensor
   structure as the project's photon amplitude. First fix
   $\mathcal F_{a\gamma\gamma}(0,0)$ from the complete invariant contact,
   WZW and resonance contributions in one basis. Then derive the virtuality
   dependence from the source-matched hadronic amplitudes. In the no-mixing
   basis, do not add meson-pole WZW terms already absorbed into the transformed
   contact coefficient. Extend the existing vector-conversion calculation
   away from $q_1^2=q_2^2=0$ and include the source-specific pseudoscalar and
   multiparticle cuts required by the chosen precision.
2. **Match the low- and high-virtuality regions.** Build a symmetric
   double-spectral representation or a validated alternative for the low
   region. Constrain it with relevant singly- and doubly-virtual data,
   source-specific low-energy amplitudes, and hadronic discontinuities. Derive
   the ALP's own short-distance coefficient and matching scale from the
   relevant QCD matrix elements / operator product expansion; do not transfer
   the $\eta$ or $\eta'$ asymptotic coefficient. Include finite-$m_a$ terms
   when the matching scale is not parametrically above $m_a$.
3. **Insert the TFF in the lepton-loop functional.** Adapt the symmetric
   two-photon kernel and reduced-amplitude normalization in Messerli:2025rnv
   to the project's TFF convention. Treat the $a$ mass as the external
   invariant and retain the full lepton-mass dependence. Compute the real
   dispersive part and the imaginary parts from every open cut represented in
   the TFF. The two-real-photon cut gives a useful normalization and unitarity
   check, but is not a substitute for the full imaginary part when other
   channels are open.
4. **Add independent local contributions consistently.** Match any induced
   ALP-lepton pseudoscalar/axial operator at the same scale and combine its
   amplitude coherently with the two-photon term. Keep the local term's
   renormalization convention explicit; avoid double counting contributions
   already incorporated in the TFF or its high-virtuality matching.
5. **Compute the width and uncertainties.** Derive the spin-summed width and
   branching fraction with explicit identical-particle and normalization
   conventions. Propagate uncertainties from source matching, hadronic
   spectral input, high-energy matching, and numerical integration separately.
   State clearly whether the result is leading-order QED or includes radiative
   corrections and the photon-cut convention used for any comparison.

## Reproducibility and checks

When the calculation is undertaken, derive normalization and tensor algebra
in an executable script under `scripts/`, then implement deterministic loop
and spectral integrations with units and a minimal invocation example. Pair
each script with a derivation record that cites every input. Verify the local
pointlike-TFF limit, the two-photon absorptive part, Bose symmetry under
$q_1\leftrightarrow q_2$, electromagnetic Ward identities, the $m_\ell$
helicity suppression of the two-photon contribution, and invariance under
the field redefinition to the no-mixing basis. Check the zero-virtuality limit
against the existing $C_{a\gamma\gamma}$ result and verify that the high-energy
tail matches the independently derived short-distance behavior.

The main unresolved inputs are a source-resolved doubly virtual ALP TFF and
its high-energy coefficient, the hadronic discontinuities for the chosen ALP
mass, and the status or boundary assumption for a local ALP-lepton
coefficient. The hadronic-scale ALP couplings are supplied by ALP-aca; their
UV evolution is outside this calculation. The current real-photon ALP-photon
notes do not settle the doubly virtual quantities.

## Literature source

Messerli:2025rnv, Secs. 2--4, provides the reduced-amplitude/TFF formalism,
low-energy dispersive construction, and asymptotic matching with finite-mass
terms. Its Sec. 6 also discusses local pseudoscalar, axial, and gluonic
operators and their meson matrix elements. Alda:2026nrd, Sec. III E and
Appendix A, supplies the no-mixing source conventions used here. The
photon-operator normalization follows Bauer:2020jbp as recorded in the ALP
photon branch.
