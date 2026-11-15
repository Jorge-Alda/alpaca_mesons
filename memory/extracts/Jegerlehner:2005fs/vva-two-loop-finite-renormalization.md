---
title: Two-loop VVA correlator and the finite renormalization constant of the singlet axial current
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Jegerlehner:2005fs, Larin:1993tq]
provenance: literature
---

# Source and scope

Read from `memory/sources/Jegerlehner:2005fs/src/anomaly.tex` (arXiv hep-ph/0510308),
abstract and "Results and Discussion" [@Jegerlehner:2005fs]. Order: two loops, chiral
limit, arbitrary momenta, $\overline{\rm MS}$ scheme.

# Statements

- Form factors (label `formfactor`):
  $w_{T,L}=n_fN_c\,w_{1,\cdot}+a\,n_fN_c\,C_2(R)\,w_{2,\cdot}$ with
  $a=\alpha_s/(4\pi)$ and $C_2(R)=4/3$ for QCD.
- The two-loop correction equals the one-loop result times a constant
  $C_2(R)\alpha_s/\pi$ in $\overline{\rm MS}$, for the longitudinal and
  transversal amplitudes alike, so the correlator is proportional to the
  one-loop term. This confirms the non-renormalization theorems of
  [Vainshtein:2002nv](../Vainshtein:2002nv/vva-soft-photon-nonrenormalization.md)
  [@Jegerlehner:2005fs].
- The constant is removed by the finite renormalization of the singlet axial
  current, $(J^5_\rho)_r=Z_5Z_{\overline{\rm MS}}(J^5_\rho)_0$ with
  $Z_{\overline{\rm MS}}=1+\mathcal O(a^2)$ and
  $$Z_5=1-4C_2(R)\,a,$$
  where the source states $Z_5$ is known to three loops from [@Larin:1993tq] and
  uses the one-loop term. The abstract states that the constant is removed by
  renormalizing the axial current according to Witten's algebraic/geometrical
  constraint on the anomalous Ward identity.
- Numerically $-4C_2(R)a=-\tfrac43\,\alpha_s/\pi$.

# Flags

- This is the coefficient of the full three-point function. It is not the
  Wilson coefficient $1-\alpha_s/\pi$ of the two-hard-photon OPE
  ([Ludtke:2020moa](../Ludtke:2020moa/ope-nlo-wilson-coefficient.md));
  the relation between the two (the Wilson coefficient multiplies the
  renormalized $j_5$) is not stated in this source.
- Scheme dependence: statements hold with the 't Hooft-Veltman/Larin
  $\gamma_5$ treatment for which $Z_5$ is defined.
