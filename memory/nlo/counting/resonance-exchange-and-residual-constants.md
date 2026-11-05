---
title: Separation of explicit resonance exchange from the local constants L_i^r(mu)
status: mixed
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Ecker:1988te, Kaiser:1998ds, Pich:2008jm, Rosell:2004mn, Rosell:2005ai, Alda:2026nrd]
provenance: assumption
---

# Resonance exchange and residual local constants

## Assumptions and inputs

- The fields are those of [field-content.md](field-content.md). The orders are those of [power-counting.md](power-counting.md): tree-level exchange of a resonance has the low-energy limit $\mathcal O(\delta)$ for its leading term, and the full propagator is kept unexpanded.
- The loop content and the conditions under which resonance loops are excluded are in [loop-classes.md](loop-classes.md). The scheme, the scale $\mu$ and the coefficients $\Gamma_i$ are in [subtraction-scheme-and-scales.md](subtraction-scheme-and-scales.md). The per-constant table is [counterterm-basis.md](counterterm-basis.md).
- The tree-level contribution of each resonance multiplet to each constant, and the coupled Goldstone–pseudoscalar-resonance system, are in [resonance-saturation-and-dynamical-masses.md](../resonance-saturation-and-dynamical-masses.md) [@Ecker:1988te].
- **Notation (approved by the user).** $L_i^r(\mu)$ is the renormalized local constant. This record fixes what it means when resonances are explicit.

This record is a definition and a rule. It contains no derivation, so it has no script of its own. The $\mu$ dependence it relies on is checked by `--scheme` in [`scripts/derive_nlo_u3_resonance_counting.py`](../../../scripts/derive_nlo_u3_resonance_counting.py).

## Definition

Integrate out a resonance at tree level. It leaves a local Lagrangian in the Goldstone fields,

$$
\Delta\mathcal L_{\rm loc}=-\frac12J_R\mathcal D_R^{-1}J_R ,
$$

with $J_R$ the linear Goldstone current of the resonance and $\mathcal D_R$ its quadratic operator ([resonance-saturation-and-dynamical-masses.md](../resonance-saturation-and-dynamical-masses.md)). In a Goldstone-only theory, the coefficient of each operator is a single constant, which contains this exchange contribution. For example, the coefficient of $\langle u_\mu u^\mu\chi_+\rangle$ includes the term $c_dc_m/M_{S_8}^2$ from the scalar octet.

In the framework of this change the resonances stay in the Lagrangian. Then:

> $L_i^r(\mu)$ is the renormalized coefficient that remains after the exchange of all explicit resonance multiplets is removed. The coefficient seen in a Goldstone-only low-energy expansion is $L_i^r(\mu)$ plus the low-energy limit of the explicit exchange.

Every record and every script that uses $L_i^r(\mu)$ must say whether it means this residual constant or the coefficient of a Goldstone-only expansion. In the alp-gii records, $L_i^r(\mu)$ always means the residual constant.

## Rule against double counting

1. **Each contribution is counted once.** An exchange kept as an explicit propagator, or as the coupled-field block $\mathcal Q(s)$ of the explicit pseudoscalar resonance, is not also included in $L_i^r(\mu)$. In particular, the induced mass term $-B^TM_R^{-2}B$ of the pseudoscalar octet, with its contribution to $L_7$ and $L_8$, is not counted in the residual constants.
2. **Literature values are not residual constants.** Values of $L_i^r(\mu)$ obtained from Goldstone-only fits to data are coefficients of a Goldstone-only expansion. They contain the exchange of the resonances in the low-energy region. Using them as $L_i^r(\mu)$ in the explicit theory counts the exchange twice. Before such a value is used, the low-energy limit of the explicit exchange, evaluated with the couplings and masses of the same calculation, must be subtracted from it. This applies to alp-gii.5.
3. **Terms beyond the leading exchange.** The series $s^n/M_R^{2n}$ of the exchange is carried entirely by the propagator. The quadratic Lagrangian keeps no local $\mathcal O(N_cp^6)$ term ([quadratic-lagrangian.md](../quadratic-lagrangian.md)), so nothing is counted twice at this order. If such a term is added later, the same subtraction applies.
4. **Constants with no resonance contribution.** For a constant to which no explicit multiplet contributes at tree level (for example $L_7$ without a pseudoscalar singlet, or $L_{47}$), the residual constant is the whole constant.
5. **Divergences.** The tree-level exchange does not depend on $\mu$. The divergences of the included loop classes, G-LO and G-ALP, are the graphs of the Goldstone-only theory, so they are absorbed with the same $\Gamma_i$ ([counterterm-basis.md](counterterm-basis.md)). The full chiral running is therefore carried by $L_i^r(\mu)$ alone. The divergences of the excluded classes, for example R-vertex, are not renormalized in this framework ([loop-classes.md](loop-classes.md)).

## Saturation hypothesis (assumption)

Resonance saturation, as proposed by Ecker:1988te, states that the Goldstone-only constants at a scale $\mu\approx M_\rho$ are given by tree-level exchange of the lightest resonances [@Ecker:1988te]. In the present definition this means

$$
L_i^r(\mu\approx M_R)\approx 0 .
$$

This is an **assumption**, not a symmetry statement, and it is used only as a starting estimate in the following sense.

- It applies only to constants with a resonance contribution in the table of [resonance-saturation-and-dynamical-masses.md](../resonance-saturation-and-dynamical-masses.md): $L_4$–$L_8$ and, outside the quadratic sector, $L_1$–$L_3$, $L_9$, $L_{10}$.
- There is no single scale. $M_V\approx0.77\,$GeV, $M_A\approx1.0\,$GeV and $M_S\approx1.09\,$GeV differ. The default is $\mu=M_\rho$, with variation over $[0.5,1.0]\,$GeV. The shifts of $L_i^r$ over this range are tabulated in [subtraction-scheme-and-scales.md](subtraction-scheme-and-scales.md) and give the uncertainty of the hypothesis.
- The residual constants are not zero in general. In the $1/N_c$ analysis of $\Pi_{VV}-\Pi_{AA}$, the NLO corrections shift the central values of $L_{10}^r$ and $C_{87}^r$ by about 17% and 28% ([loop-classes.md](loop-classes.md)) [@Pich:2008jm]. Under the exclusion of resonance loops, residual constants in the quadratic sector are therefore fixed from data (alp-gii.5), and the hypothesis serves only as a prior.

## Singlet and topological constants

The standard saturation table does not fix $\Lambda_1$, $\Lambda_2$, $L_{18}$ and $L_{25}$ ([resonance-saturation-and-dynamical-masses.md](../resonance-saturation-and-dynamical-masses.md)). They are therefore treated as free residual constants:

- No resonance contribution is assigned to them, unless a singlet/topological resonance Lagrangian is specified and matched in alp-gii.3.
- The only statement in the sources is qualitative: resonance exchange with OZI-respecting couplings gives $L_{18}=\Lambda_1=0$ at $\mu=M_\rho$, so that $L_A\simeq2L_5^r(M_\rho)$. Kaiser:1998ds warns that the result is sensitive to the scale at which $|L_{18}^r|\ll|L_5^r|$ and $|\Lambda_1|\ll1$ are assumed ([extract](../../extracts/Kaiser:1998ds/u3-loop-renormalization-and-gamma-table.md)) [@Kaiser:1998ds]. This is not used as an input.
- The QCD-scale ($Z_A$) behaviour of these constants is in [subtraction-scheme-and-scales.md](subtraction-scheme-and-scales.md). The OZI-violating difference between singlet and octet resonance couplings is $\delta^{3/2}$ ([power-counting.md](power-counting.md)), so its exchange enters at $\delta^2$.
- The ALP sources of Alda:2026nrd enter through $X_{TC}$, $\bm{u_\mu}$ and $\bm{\chi_\pm}$ ([field-content.md](field-content.md)). The ALP–meson terms of the quadratic Lagrangian use the same constants, and they need no separate residual constant at the declared order. Whether a further source operator is needed is decided in alp-gii.8 (marked row of [counterterm-basis.md](counterterm-basis.md)) [@Alda:2026nrd].

## What alp-gii.3 has to do with this record

- Compare the full quadratic NLO Lagrangian with the explicit exchange, coefficient by coefficient, and compute the low-energy limit of the exchange for each constant in the quadratic sector.
- Report, for each constant, the residual constant that remains and the matching scale $\mu$ used.
- Recover the neutral mass and mixing matrices of [neutral-meson-matrices.md](../neutral-meson-matrices.md) in the low-energy limit, with the exchange removed from the local terms.

## Consistency check against the existing record

[resonance-saturation-and-dynamical-masses.md](../resonance-saturation-and-dynamical-masses.md) states that, once a resonance is explicit, $B$ and $M_R$ are retained and the induced pieces are omitted from the local constants, and that an independent short-distance remainder stays in the local Lagrangian. This is the definition of $L_i^r(\mu)$ above. The two records agree.
