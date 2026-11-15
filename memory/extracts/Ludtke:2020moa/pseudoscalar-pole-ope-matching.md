---
title: Pseudoscalar-pole matching to the OPE constraint, mass shift of the pole, and eta-eta-prime relation
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Ludtke:2020moa, Agaev:2014wna]
provenance: literature
---

# Source and scope

Read from `memory/sources/Ludtke:2020moa/src/SDC.tex` [@Ludtke:2020moa]. Labels are
LaTeX labels. Leading order in $\alpha_s$ unless stated. This record is
about how the pseudoscalar pole is matched to the OPE, not about the Wilson
coefficient (see [the NLO record](ope-nlo-wilson-coefficient.md)).

# Equations

- Flavor weights (label `Eq:Ca`): $C_a={\rm tr}(\hat Q^2\lambda_a)/2$, with
  $C_3=1/6$, $C_8=1/(6\sqrt3)$, $C_0=2/(3\sqrt6)$.
- Longitudinal HLbL function in the OPE regime (label `Eq:OPEConstraint`):
  $\bar\Pi_1^{(a),{\rm OPE}}(Q^2,Q^2,Q_3^2)=-2N_cC_a^2/(\pi^2Q^2Q_3^2)$ for
  $a=3,8,0$, valid for $Q^2\gg Q_3^2\gg\Lambda_{\rm QCD}^2$.
- The source states that for the non-singlet components perturbative and
  nonperturbative corrections to the anomaly are absent in the chiral limit,
  so the $Q_3^2\gg\Lambda_{\rm QCD}^2$ restriction can be dropped. For the
  singlet, the gluon anomaly does not change the equation for
  $Q_3^2\gg\Lambda_{\rm QCD}^2$, but extrapolation to small $Q_3^2$ holds only
  at large $N_c$, where the anomaly vanishes [@Ludtke:2020moa].
- Leading-order pion limit used (text before `Eq:PiPoleVanQ3`):
  $\lim_{Q^2\to\infty}Q^2F_{\pi\gamma^*\gamma^*}(-Q^2,-Q^2)=4C_3F_\pi$, attributed
  to Nucl. Phys. B 237 (1984) 525 and Manohar:1990hu.
  With $C_3=1/6$ this is $2F_\pi/3$.
- Real-photon value quoted: $F_{\pi\gamma^*\gamma^*}(0,0)=3C_3/(2\pi^2F_\pi)$,
  i.e. $1/(4\pi^2F_\pi)$.
- Mass shift (label `Eq:OPEPionMass`, extended in `Eq:OPEConstraintMasses`):
  $\lim_{Q^2\to\infty}Q^2\bar\Pi_1^{\rm PS}(Q^2,Q^2,Q_3^2)=-6C_{\rm PS}^2/(\pi^2(Q_3^2+m_{\rm PS}^2))$.
  For $\eta,\eta'$ the source assumes ground-state pole dominance at $Q_3^2=0$ to
  define $C_{\eta,\eta'}$ (label `Eq:defCetaetap`); this is an assumption of the
  source, flagged as such there.
- Chiral-limit sum rule without gluon anomaly (label `Eq:mixingOPE`, then
  `Eq:CPhysFlavor`): $C_\eta^2+C_{\eta'}^2=C_8^2+C_0^2$, violated by the
  anomaly-induced running of the singlet decay constant.
- The source states that the gluon anomaly induces a running of the singlet
  decay constant, citing Leutwyler:1997yr, Kaiser:1998ds and Kaiser:2000gs,
  and that this gives an incomplete cancellation between the symmetric
  asymptotic and the real-photon limits with a sizable effect due to the
  large scale separation, citing [@Agaev:2014wna] (text after `Eq:mixingOPE`).

# Flags

- The source is an HLbL paper: the ALP-specific pole structure is not
  treated. The longitudinal OPE saturation statement is a statement about
  massless non-singlet Goldstone poles; it is not transferable to an ALP of
  mass $m_a$ without re-derivation.
