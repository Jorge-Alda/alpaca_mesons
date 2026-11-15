---
title: Non-renormalization of the VVA triangle with one soft vector current
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Vainshtein:2002nv]
provenance: literature
---

# Source and scope

Read from `memory/sources/Vainshtein:2002nv/src/nonren3.tex` (arXiv hep-ph/0212231)
[@Vainshtein:2002nv]. Labels are LaTeX labels. This concerns the amplitude with one
hard axial current and one hard vector current while the second vector
current is a soft photon. It is not the Wilson coefficient of the
two-hard-photon OPE; the two statements are different objects and are
both needed (see [Melnikov:2003xd](../Melnikov:2003xd/hard-ope-axial-current.md)
and [Ludtke:2020moa](../Ludtke:2020moa/ope-nlo-wilson-coefficient.md)).

# Definitions (label `invfun`)

$T_{\mu\nu}=i\int d^4x\,e^{iqx}\langle0|T\{j_\mu(x)j^5_\nu(0)\}|\gamma(k)\rangle$,
$$
T_{\mu\nu}=-\frac{i}{4\pi^2}\Bigl[w_T(q^2)\bigl(-q^2\tilde f_{\mu\nu}+q_\mu q^\sigma\tilde f_{\sigma\nu}-q_\nu q^\sigma\tilde f_{\sigma\mu}\bigr)+w_L(q^2)\,q_\nu q^\sigma\tilde f_{\sigma\mu}\Bigr],
$$
$\tilde f_{\mu\nu}=\tfrac12\epsilon_{\mu\nu\gamma\delta}f^{\gamma\delta}$,
$f_{\mu\nu}=k_\mu e_\nu-k_\nu e_\mu$ [@Vainshtein:2002nv].

# Result

- Label `wtwl`: $w_L[m=0]=2\,w_T[m=0]$, true at one loop, receives no
  perturbative corrections from gluon exchange. The longitudinal part
  is not renormalized (Adler-Bardeen) in the chiral limit.
- The argument uses the symmetry $\mu\leftrightarrow\nu$, $q\to-q$ of the
  singular-in-$q$ part, obtained by moving $\gamma_5$ through an even
  number of gamma matrices in the chiral limit.
- Nonperturbatively the transversal part is corrected, the longitudinal
  one is not (conclusions section). The source uses this to determine the
  magnetic susceptibility $\chi=-N_c/(4\pi^2F_\pi^2)$ [@Vainshtein:2002nv].
- Order in $\alpha_s$: all orders in perturbation theory, chiral limit.

# Flags

- Conventions for $\epsilon$ and $\gamma_5$ are not used to fix any sign
  here; the record is not used to fix the project sign.
- Applies to the non-anomalous-dimension (non-singlet) normalization.
  [Jegerlehner:2005fs](../Jegerlehner:2005fs/vva-two-loop-finite-renormalization.md)
  gives the two-loop check and the finite renormalization needed for the singlet current.
