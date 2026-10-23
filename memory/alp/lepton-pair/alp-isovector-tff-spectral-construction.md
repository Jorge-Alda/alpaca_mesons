---
title: Low-energy isovector double-spectral construction for the ALP TFF
status: mixed
created: 2026-10-23
updated: 2026-10-23
bibtex_keys: [Messerli:2025rnv, Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and inputs

This record specifies the first hadronic spectral component for the
doubly-virtual ALP TFF in the lepton-pair calculation. The two-pion cut
produces the low-energy isovector channel. The representation is adapted from
Messerli:2025rnv, Eqs. (5.12)--(5.14), but its transition amplitude and weights
must be calculated for the ALP sources; the eta and eta-prime amplitudes are
not inputs [@Messerli:2025rnv]. Source definitions and the scalar TFF tensor
convention are those of [the source-resolved ALP TFF record](alp-doubly-virtual-tff.md)
[@Alda:2026nrd].

Use $s_i=q_i^2$ and the unnormalized scalar $F_{a\gamma^*\gamma^*}$ defined
there by $i\mathcal M=i e^2 R F$. Let $\mathcal F_{a\pi\pi\gamma^*}(x,y)$
denote the isovector $P$-wave amplitude for $a\gamma^*(y)\to\pi^+\pi^-$,
with the partial-wave and tensor normalization chosen to match
Messerli:2025rnv's $\mathcal F_{P\pi\pi\gamma^*}$. Here $x$ is the dipion
invariant mass squared and $y$ is the second photon's virtuality. This
normalization match is required before using the numerical prefactor below.

# Spectral density and TFF

For $x\geq4M_\pi^2$, define
$$
\sigma_\pi(x)=\sqrt{1-\frac{4M_\pi^2}{x}},\qquad
\rho_a^{I=1}(x,y)=\frac{x\sigma_\pi^3(x)}{192\pi}
\operatorname{Im}\!\left\{[F_\pi^V(x)]^*\,
\mathcal F_{a\pi\pi\gamma^*}(x,y)\right\}.
$$
The $x\sigma_\pi^3/(192\pi)$ factor is the two-pion $P$-wave phase-space
factor in the normalization of Messerli:2025rnv, Eq. (5.14). $F_\pi^V$ is
the charged-pion vector form factor, normalized by $F_\pi^V(0)=1$;
the phase convention of the $P$-wave amplitude must be fixed consistently
with it [@Messerli:2025rnv].

For spacelike $s_i=-Q_i^2$ and a declared low-energy cutoff $\Lambda$, the
isovector TFF is
$$
F_{a\gamma^*\gamma^*}^{I=1}(-Q_1^2,-Q_2^2;\Lambda)=
\frac{1}{\pi^2}\int_{4M_\pi^2}^{\Lambda^2}\!dx
\int_{4M_\pi^2}^{\Lambda^2}\!dy\,
\frac{\rho_a^{I=1}(x,y)}{(x+Q_1^2)(y+Q_2^2)}
+(Q_1\leftrightarrow Q_2).
$$
This is a low-energy component, not a complete TFF or a prescription for its
high-virtuality continuation. The cutoff and any smooth replacement of the
sharp cutoff are model choices. Their variation belongs in the spectral-model
uncertainty; the omitted strength must be supplied by the higher-energy
components and matched without double counting [@Messerli:2025rnv]. For
timelike virtualities, use the Feynman boundary values of both denominators
and retain the resulting imaginary parts; the spacelike equation alone does
not prescribe those cuts.

# ALP source response

At the retained leading order the amplitude is linear in the no-mixing-basis
ALP source. Choose an explicitly documented independent source basis
$\{\mathcal O_A\}$, and calculate
$$
\mathcal F_{a\pi\pi\gamma^*}(x,y)
=\sum_A g_A(\mu)\,\mathcal F_A(x,y;\mu),\qquad
\rho_a^{I=1}(x,y)
=\sum_A g_A(\mu)\,\rho_A^{I=1}(x,y;\mu),
$$
where $g_A(\mu)$ are the source coefficients in the convention of
[@Alda:2026nrd], including the common scale factors exactly once. The
source-resolved functions $\mathcal F_A$ are the required hadronic matrix
elements, not fitted eta-like weights. Include the WZW/contact and resonance
pieces in the same amplitude only when their photon and source placements
contribute to this $P$ wave; maintain the no-mixing bookkeeping so meson-pole
terms already absorbed in $\widehat c_\gamma$ are not counted again. The
tree-level neutral-vector TFF in the linked record does not itself determine
the two-pion discontinuity.

For a first hadronic model, the needed partial wave can be organized by
two-pion unitarity and a dispersive $P$-wave representation: the right-hand
cut is fixed by the isovector pion phase, while subtraction data and any
left-hand-cut terms must be fixed from the ALP source amplitudes. A
source-dependent $a\rho\gamma^*$ production vertex followed by
$\rho\to\pi\pi$ is one resonance contribution to that partial wave; it is
not by itself the full partial wave once pion rescattering or nonresonant
terms are retained. No such decomposition or subtraction constants are
numerically specified here.

The project's existing [odd $a\to\gamma PP$ amplitude](../photon-two-pseudoscalar/alp-photon-pp-amplitude.md)
provides a tree-level starting kernel for the charged-pion channel, including
the $a\to V^*\gamma$ and crossed $a\to\pi V^*$ placements. Its formula is
for a real photon, retains only the stated odd resonance-exchange terms, and
explicitly omits local anomalous contributions [@Ecker:1988te; @Ecker:1989yg;
@Ruiz-Femenia:2003jdx]. For the required virtual current one must rederive
the amplitude with $q_\gamma^2=y$: replace the real-photon response
$\mathcal W_\gamma(0)$ by the virtual vector response where applicable, keep
the $y$-dependent VJP/VVP numerators and all crossed placements, then project
the result onto the dipion $P$ wave. In particular, a real-photon width or
its integrated Dalitz spectrum cannot determine this virtuality dependence.
The projector normalization must be matched to Messerli:2025rnv's
$\mathcal F_{P\pi\pi\gamma^*}$ convention.

Completeness of the double-spectral term must be checked after this
projection. A local $a\pi\pi\gamma^*$ term, or a subtraction contribution
with no discontinuity in $y$, can have an $x$-channel two-pion cut without
appearing in $\rho_a^{I=1}(x,y)$. Such terms require separate
single-spectral/subtraction contributions in the TFF if they survive the
chosen source basis and truncation. The displayed double-spectral integral
therefore defines the part supported by both spectral variables; it is not
declared the complete low-energy isovector TFF until this audit is done.

# Assumptions, missing inputs, and validity

- Isospin-symmetric charged-pion threshold and the two-pion $P$-wave are used
  for this first low-energy isovector component. Isospin breaking and other
  intermediate states are outside this record.
- The two-pion unitarity relation and the quoted prefactor are borrowed in
  Messerli:2025rnv's partial-wave normalization. A direct ALP amplitude
  derivation must verify the normalization, sign, and phase before numerical
  use. A real tree-level ALP kernel alone is not a spectral density: the
  density is the imaginary part of the product with the conjugate pion vector
  form factor, with the required rescattering phases treated consistently.
- Required inputs are a source-resolved $a\gamma^*\to\pi^+\pi^-$ $P$-wave
  amplitude, the pion vector form factor and its phase, a cutoff prescription,
  and their correlated uncertainties. None is numerically selected here.
- The existing $a\to\gamma\pi^+\pi^-$ result is for a real photon and excludes
  local anomalous terms; its crossing gives a starting amplitude, not the
  required doubly virtual partial wave.
- This channel covers the low-energy isovector cut only. It does not fix the
  isoscalar poles, higher resonances, the short-distance tail, local
  ALP-lepton matching, or additional open cuts. Do not call its contribution
  the full imaginary part or a lower bound on the width.

# Reproducibility

The executable [`derive_alp_isovector_tff_spectral.py`](../../../scripts/derive_alp_isovector_tff_spectral.py)
records the source-linear construction and checks Bose symmetry and the
double-dispersion denominators. It is a formal scaffold: until the physical
$\mathcal F_A$ and $F_\pi^V$ inputs are supplied, it does not evaluate a
spectral density or a numerical TFF. Run
`.venv/bin/python scripts/derive_alp_isovector_tff_spectral.py`.

## Literature source

Messerli:2025rnv, Eqs. (5.12)--(5.14), gives the isovector decomposition,
double-spectral representation, cutoff integral, and two-pion spectral
density used here [@Messerli:2025rnv]. Alda:2026nrd supplies the ALP
no-mixing source convention [@Alda:2026nrd].
