---
title: Complete the physical ALP two-vector amplitude in the retained RChT action
status: mixed
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg]
provenance: assumption
---

# Physical aVV completion plan

## Starting evidence and scope

The [two-source record](two-source-alp-rcht-odd-parity-extension.md) calculates
physical-ALP VVP insertions, but neither it nor the numerical operator kernel
assembles a complete normalized physical $a\to VV$ predictor.
The fresh [two-boson audit](alp-two-boson-blind-direction-audit.md) projects out
$d_4$ and the VJP blind directions for wedge legs; it does not replace this
completion. Work within the retained tensor nonet/single-trace truncation and
state any added singlet/OZI term separately. HLS is not part of this plan.

## Inputs

- Alda:2026nrd's arXiv-v1 $\mathcal C$, $X=\mathbb M_a^{-1}(\mathcal C)$,
  $Y=c_q^{AP}+X$, independent derivative/mass sources and physical generators.
- The [VVP basis](../../extracts/Ruiz-Femenia:2003jdx/vvp-three-pseudoscalar.md),
  canonical external tensor wavefunctions and operator normalizations.
- The [stage-4 sector](vector-stage4-amplitudes.md), with both mass branches
  and explicit source hypotheses. Near-ideal is preferred only if other
  criteria are equivalent; there is no averaging of mass scenarios.
- Explicit ALP UV coefficients and $m_a,f_a$ as inputs, not SM-fit outputs.

## Ordered calculation

1. Set the order: leading strong $e^0/f_a$ VVP contribution first; list any
   retained electromagnetic corrections separately. Insert physical-ALP
   $u_\mu$ and $\chi_-$ with independent $X,Y$ and contract both labeled vector
   assignments in `scripts/derive_physical_alp_vv.py` using SymPy.
2. Use each physical external generator and its own canonical
   $i(q_\mu\varepsilon_\nu-q_\nu\varepsilon_\mu)/M_i$ wavefunction.
   Derive the scalar coefficient multiplying
   $\epsilon_{\alpha\beta\rho\sigma}\varepsilon_1^\alpha
   \varepsilon_2^\beta q_1^\rho q_2^\sigma$, retaining unlike and identical
   vector multiplicities explicitly.
3. For electromagnetic terms, assemble photon/vector mixing and every induced
   local term in one field basis. A two-background-source anomaly seed is not
   an independently fixed hadronic resonance coupling: do not introduce a free
   hadronic $c_{AB}$ by treating a resonance as a UV gauge source. State any
   additional matching assumption and prevent double counting with VVP terms.
4. Verify equivalence of direct-plus-meson-pole LSZ assembly and Alda:2026nrd's no-mixing
   variables at $O(1/f_a)$. Check field reparametrization and generator-phase
   invariance, vector assignment symmetry and the existing projected nulls.
   Add any otherwise required invariant local terms as explicit new inputs.
5. Derive the polarization sum, two-body phase space, identical-particle factor
   and width in the executable. Verify dimensions and threshold closure.
   Narrow on-shell daughter vectors are an explicit approximation; broad or
   near-threshold channels require their physical daughter amplitudes/spectral
   treatment and a fresh off-shell sensitivity check.
6. Differentiate the completed coefficient along the stored SM blind columns
   before any OPE substitution. Test every mass/source hypothesis and independent
   $X,Y$ inputs. Record row-span coverage separately from statistical constraints
   and from the source/UV matching uncertainty.

## Deliverables and validity

Pair the executable with a stepwise memory derivation and a numerical predictor
under `scripts/implementation/`. Add independent limiting-normalization,
assignment/phase, null-projection and threshold tests. Enforce Alda:2026nrd's exact
source-overlap near-degeneracy condition point by point; a finite supplied
$X$ alone is insufficient. Do not reuse the WZW meson-pole term a second time
in the no-mixing basis. A physical prediction also needs calibrated couplings
and declared uncertainties; those are not supplied by this plan.
