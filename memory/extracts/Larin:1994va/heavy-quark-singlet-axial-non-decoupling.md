---
title: Non-decoupling of a heavy quark from the singlet axial correlator and the alpha_s decoupling relation
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Larin:1994va, Chetyrkin:1997un]
provenance: literature
---

# Source and scope

Read from `memory/sources/Larin:1994va/src/main.tex` (arXiv hep-ph/9411260),
introduction and the paragraphs on the singlet axial part and the decoupling
relation [@Larin:1994va]. The paper treats $Z^0\to$ hadrons with a heavy top quark and five
massless flavors, order $\alpha_s^3$. It is a different observable from a
photon-photon transition form factor; only two general statements are used.

# Statements

- The top quark "does not decouple" from the axial-vector part because of
  diagrams of axial-anomaly type (introduction), citing Collins, Wilczek and Zee.
  The vector and non-singlet axial parts decouple.
- The singlet axial contribution vanishes for light doublets in the massless limit
  because the up-type and down-type quarks of a weak doublet enter with opposite
  axial coupling; only the heavy mass splitting (top-bottom) gives a nonzero
  singlet contribution, starting at $(\alpha_s/\pi)^2$ with the leading coefficient
  $-17/6$ (times the squared axial coupling) in the effective five-flavor theory.
- Decoupling relation of $\alpha_s$ (label `decouplingrelation`), with $T_F=1/2$
  for QCD and $\mu_t=m_t(\mu)$:
  $$
  \frac{\alpha_s^{(6)}}{\pi}=\frac{\alpha_s^{(5)}}{\pi}+\Bigl(\frac{\alpha_s^{(5)}}{\pi}\Bigr)^2\frac{T_F}3\ln\frac{\mu^2}{m_t^2}
  +\Bigl(\frac{\alpha_s^{(5)}}{\pi}\Bigr)^3\Bigl(\frac{T_F^2}9\ln^2\frac{\mu^2}{m_t^2}+\frac{5C_AT_F-3C_FT_F}{12}\ln\frac{\mu^2}{m_t^2}+\frac{13}{48}T_FC_F-\frac29T_FC_A\Bigr).
  $$
  The source remarks that its $\frac{13}{48}T_FC_F$ term differs slightly from
  the Z. Phys. C 11 (1981) 79 paper's value. At $\mu=m_t$ and $T_F=\frac12$, $C_F=\frac43$, $C_A=3$ the constant is
  $-\frac{11}{72}$ (hand arithmetic), the inverse of the coefficient $c_2=\frac{11}{72}$ of
  [Chetyrkin:1997un](../Chetyrkin:1997un/alphas-decoupling-constant.md) at this order,
  consistent with the cross-reference in that source.

# Flags

- The heavy-quark non-decoupling here is for a weak-current correlator with large mass splitting.
  It is not shown to apply to a photon-photon OPE with a single flavor of heavy
  quark. Applicability to the ALP/TFF problem is an open question; whether the heavy-quark
  axial current gives a power-suppressed or a constant contribution to the
  OPE coefficient is not established by this source.
