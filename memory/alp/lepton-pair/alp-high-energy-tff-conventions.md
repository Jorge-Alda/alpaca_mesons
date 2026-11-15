---
title: Shared conventions for the ALP high-energy TFF completion
status: mixed
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Bauer:2020jbp]
provenance: literature
---

# Scope

This record fixes, in one place, the conventions used by the short-distance
completion of the doubly virtual ALP transition form factor (TFF). It collects
conventions already stated in the [tree TFF](alp-doubly-virtual-tff.md) and the
[leading hard-current expansion](alp-hard-current-ope.md); it does not derive
anything new, except for the elementary momentum relation in the last section.
The source conventions are those of Alda:2026nrd; the TFF tensor and its
normalization are those of Messerli:2025rnv [@Alda:2026nrd; @Messerli:2025rnv].

# Conventions

| Item | Convention | Where fixed |
|---|---|---|
| Metric | $(+---)$ | [hard-current OPE](alp-hard-current-ope.md), [tree TFF](alp-doubly-virtual-tff.md) |
| Levi-Civita tensor | $\epsilon^{0123}=+1$ | same; see also [Levi-Civita convention](../../rcht/odd-parity/levi-civita-convention.md) |
| Fifth Dirac matrix | $\gamma_5=i\gamma^0\gamma^1\gamma^2\gamma^3$ | [hard-current OPE](alp-hard-current-ope.md) |
| Electromagnetic current | $J_{\rm em}^\mu=\sum_fQ_f\bar q_f\gamma^\mu q_f$, without the charge $e$ | [hard-current OPE](alp-hard-current-ope.md) |
| Photon tensor | $R=\epsilon_{\mu\nu\rho\sigma}\varepsilon_1^\mu\varepsilon_2^\nu q_1^\rho q_2^\sigma$ | [@Messerli:2025rnv, Eq. (9)] |
| Amplitude | $i\mathcal M=ie^2R\,F_{a\gamma^*\gamma^*}(s_1,s_2)$, $s_i=q_i^2$ | [@Messerli:2025rnv, Eq. (9)] |
| Current-product tensor | $T^{\mu\nu}=\epsilon^{\mu\nu\rho\sigma}q_{1\rho}q_{2\sigma}F_{a\gamma^*\gamma^*}$ | [hard-current OPE](alp-hard-current-ope.md) |
| Normalized TFF | $\widetilde F=F/F_{a\gamma\gamma}$, $F_{a\gamma\gamma}=F_{a\gamma^*\gamma^*}(0,0)$ | [@Messerli:2025rnv, Eq. (7)] |
| Dimension of $F$ | $\mathrm{GeV}^{-1}$, relativistic one-particle normalization | [tree TFF](alp-doubly-virtual-tff.md) |
| ALP sources | no-mixing matrices $U_a$, $Y_a$, contact coefficient $\widehat c_\gamma=C_{a\gamma\gamma}(m_a)$, with $U_a=-2S/f_a$, $Y_a=2\{c_q^{AP}+S,\chi\}/f_a$ | [tree TFF](alp-doubly-virtual-tff.md), [@Alda:2026nrd] |
| Source equation of motion (LO) | $Y_a+2m_a^2U_a-2Z_a\mathbf 1=0$ | [tree TFF](alp-doubly-virtual-tff.md) |
| Photon coefficient sign | positive $\widehat c_\gamma$ in the canonical dual convention | [@Bauer:2020jbp], [photon phase record](../photon/alp-photon-odd-action-phase-matching.md) |

## Sign of the tensor in the existing scripts

The intermediate source-contraction scripts use the tensor
$E=\epsilon(q_1,\varepsilon_1,q_2,\varepsilon_2)=-R$, so their scalar is
$F=-G$. Results in this change are always quoted for the tensor $R$, as in the
[tree TFF](alp-doubly-virtual-tff.md).

## Momentum labels

The tree record uses $p_a+q_1+q_2=0$ (all momenta incoming, [@Messerli:2025rnv]).
The hard-current record uses the vacuum-to-ALP matrix element
$\langle0|\cdots|a(p)\rangle$ with $p=q_1+q_2$ (ALP incoming, photon momenta
outgoing). The two labelings differ by $q_i\to-q_i$ for both photons. Since $R$
and $T^{\mu\nu}$ are bilinear in $q_1,q_2$, and $s_i=q_i^2$ is even, the scalar
$F_{a\gamma^*\gamma^*}(s_1,s_2)$ is the same in both.

With $q=(q_1-q_2)/2$ and $p=q_1+q_2$, the antisymmetric part of
$q_{1\rho}q_{2\sigma}$ contracted with $\epsilon^{\mu\nu\rho\sigma}$ equals
$\epsilon^{\mu\nu\rho\sigma}q_\rho p_\sigma$ (elementary algebra), which is the
tensor structure appearing in the local expansion. The virtualities are
$$
s_{1,2}=q^2\pm q\cdot p+\frac{m_a^2}{4},
$$
so they are nearly equal when $|q\cdot p|\ll|q^2|$.

# Related records

- [Tree TFF](alp-doubly-virtual-tff.md): resonance form of $F$ in these conventions.
- [Leading hard-current expansion](alp-hard-current-ope.md): the local short-distance coefficient.
- [Input gates](alp-isoscalar-tff-input-gates.md): open scale, scheme, and low-energy separation.
- [Dilepton loop formalism](../../extracts/Messerli:2025rnv/dilepton-loop-formalism.md): extract for the TFF normalization.

# Checks

The links above were verified to resolve relative to this file (see the change
task 1.2). No script is attached because no derivation is performed here beyond
the stated elementary identity.
