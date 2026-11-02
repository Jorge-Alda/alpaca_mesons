---
title: Singlet-source derivative terms and the V4 field convention
status: literature
created: 2026-11-02
updated: 2026-11-02
bibtex_keys: [Kaiser:1998ds]
provenance: literature
---

# Singlet-source derivative terms and the V4 field convention

In the general effective Lagrangian, Kaiser:1998ds includes

$$
V_3D_\mu\psi D^\mu\psi+V_4D_\mu\psi D^\mu\theta+V_5D_\mu\theta D^\mu\theta
$$

alongside the usual meson derivative and mass terms (Eq. (10)). Its covariant derivatives are

$$
D_\mu\psi=\partial_\mu\psi-2\langle a_\mu\rangle,\qquad
D_\mu\theta=\partial_\mu\theta+2\langle a_\mu\rangle.
$$

The paper states that the invariant variable is $\widetilde\psi=\psi+\theta$, and the coefficient functions $V_n$ may depend on it. It then uses the field-variable freedom $U\to U\exp[i f(\widetilde\psi)]$ to choose $V_4=0$ (Eq. (13)). This is why the cross term is absent in the chosen reduced basis: it is removable by a field redefinition, not forbidden by chiral symmetry. The paper also discards $V_5$, stating that this term does not contribute to the meson masses, decay constants, or photonic transition rates under study (source lines 210–236). [@Kaiser:1998ds]

At order $p^4$, the paper says that the general Lagrangian contains 57 potentials, then retains the subset that contributes to masses or decay constants and survives at leading order in $N_c$ (source lines 309–329). The full source-dependent list includes further terms with derivatives of $\hat\theta$; in particular, its operators $O_{46}$ and $O_{52}$ are said to contribute to decay constants but are removed by a field redefinition $U\to U\exp[f_1D_\mu\theta\,U^\dagger D^\mu U+i f_2D_\mu D^\mu\theta]$ (Eq. (17), footnote 4; source lines 332–335). The full list also contains source-derivative structures with additional meson derivatives or external field strengths, whose relevance depends on the observables and sources retained. Accordingly, [the repository's retained NLO basis](../../nlo/lagrangian.md) is not a claim to reproduce every source-dependent operator in the unrestricted effective action. Its derivative structures use the $V_4=0$ convention. If one keeps $V_4D_\mu\psi D^\mu\theta$ in another field convention, its quadratic ALP–meson mixing must be included in the kinetic blocks; the alignment formula derived for the retained vielbein terms alone does not automatically cover it. [@Kaiser:1998ds]
