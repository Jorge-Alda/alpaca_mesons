---
title: Field content of the U(3) NLO theory with explicit resonances
status: mixed
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Alda:2026nrd, Kaiser:2000gs, Kaiser:1998ds, Ecker:1988te, Ecker:1989yg]
provenance: assumption
---

# Field content

This record declares the fields used by every calculation of the alp-gii epic. The choice of fields, the nonet assignment of each resonance multiplet and the truncation to one multiplet per channel are **assumptions** of the framework. The transformation laws and Lagrangians that are quoted come from the cited literature.

## Pseudoscalar nonet and topological source

The light pseudoscalars form a $U(3)$ nonet that includes $\eta_0$. They are collected in $\Phi$ and enter through $\bm u$, with $\bm U=\bm u^2$, in the conventions of Alda:2026nrd: $\phi=\Phi/(2f_\pi)$ and the expansion formulas of [Appendix A](../../extracts/Alda:2026nrd/appendix-a-expansion.md) [@Alda:2026nrd]. Under a chiral transformation $(g_L,g_R)$,

$$
\bm u\to g_R\bm u h^\dagger=h\bm u g_L^\dagger ,
$$

where $h(g_L,g_R,\bm U)$ is the compensator [@Ecker:1988te].

The $U(1)_A$ anomaly enters through the topological combination

$$
X_{TC}=2i\left(\operatorname{tr}\log \bm u+i c_G\frac{a}{f_a}\right),
\qquad
\mathcal L_{TC}=-\frac{f_\pi^2m_0^2}{4N_f}X_{TC}^2 ,
$$

as in [the LO topological record](../../extracts/Alda:2026nrd/lo-topological-action.md) [@Alda:2026nrd]. In the notation of Kaiser:2000gs this is the invariant combination $\psi+\theta=-X_{TC}$ ([lagrangian.md](../lagrangian.md)). A nonet theory without this term is not $U(1)_A$ covariant at leading order in the $\delta$ expansion [@Kaiser:2000gs]. The leading meson Lagrangian is the one of Alda:2026nrd. Its NLO extension is the reduced $V_4=0$ basis in [lagrangian.md](../lagrangian.md) [@Kaiser:1998ds].

## ALP

The ALP $a$ is not a dynamical field of the chiral expansion. It enters through the quark-level couplings $c_q^A$, $c_q^V$, $c_q^S$, $c_q^P$ and $c_G$ of Alda:2026nrd, which appear inside $\bm{u_\mu}$, $\bm{\chi_\pm}$ and $X_{TC}$ ([Appendix A](../../extracts/Alda:2026nrd/appendix-a-expansion.md)). Calculations are carried out at first order in $1/f_a$. Its kinetic and mass terms are external: $m_a^2$ is a free parameter, counted like $p^2$ in [the power-counting record](power-counting.md). The mixed ALP–meson two-point functions are defined in [the LSZ invariant-split record](../lsz-invariant-split.md).

## Resonance nonets

Each resonance multiplet $R$ is a $U(3)$ nonet that transforms homogeneously,

$$
R\to hRh^\dagger ,
$$

with the octet plus singlet normalization of [rcht-basics.md](../../rcht/foundations/rcht-basics.md) [@Ecker:1988te]. The retained multiplets are as follows.

| Nonet | $J^{PC}$ | Field | LO coupling to Goldstones and sources | Memory record of LO Lagrangian |
|---|---|---|---|---|
| Vector | $1^{--}$ | antisymmetric tensor $V_{\mu\nu}$ | $\frac{F_V}{2\sqrt2}\langle V_{\mu\nu}f_+^{\mu\nu}\rangle+\frac{iG_V}{\sqrt2}\langle V_{\mu\nu}u^\mu u^\nu\rangle$ | [rcht-basics.md](../../rcht/foundations/rcht-basics.md), [Ecker:1988te extract](../../extracts/Ecker:1988te/rcht-leading-lagrangian.md) |
| Axial vector | $1^{++}$ | antisymmetric tensor $A_{\mu\nu}$ | $\frac{F_A}{2\sqrt2}\langle A_{\mu\nu}f_-^{\mu\nu}\rangle$ | **missing**; input from Ecker:1988te |
| Scalar | $0^{++}$ | $S$ | $c_d\langle S u_\mu u^\mu\rangle+c_m\langle S\chi_+\rangle$ | **missing**; input from Ecker:1988te |
| Pseudoscalar | $0^{-+}$ | $P$ | $i d_m\langle P\chi_-\rangle$ | **missing**; input from Ecker:1988te; coupled pole system in [resonance-saturation-and-dynamical-masses.md](../resonance-saturation-and-dynamical-masses.md) |

The kinetic and mass terms of $V_{\mu\nu}$ are those of [the Ecker:1988te extract](../../extracts/Ecker:1988te/rcht-leading-lagrangian.md). $A_{\mu\nu}$ has the same form with $M_A$. $S$ and $P$ have canonical kinetic terms, $\frac12\langle\nabla_\mu R\nabla^\mu R-M_R^2R^2\rangle$ [@Ecker:1988te]. For the scalar and pseudoscalar singlets, Ecker:1988te uses independent couplings $\widetilde c_d,\widetilde c_m,\widetilde d_m$ and masses $M_{S_1},M_{P_1}$. Here they are kept as independent inputs. In the nonet limit they are related to the octet couplings (see [the power-counting record](power-counting.md)).

**Assumptions.**
1. One multiplet per channel. Heavier excitations, for example $\rho(1450)$, are excluded.
2. Spin-one fields are antisymmetric tensors. Ecker:1989yg shows that this representation gives the $\mathcal O(p^4)$ resonance contributions at tree level without additional local terms. This is a choice of representation, not a physical hypothesis [@Ecker:1989yg].
3. No ALP coupling is written directly to a resonance field. Any ALP–resonance coupling comes through the sources in $u_\mu$, $\chi_\pm$, $f_\pm^{\mu\nu}$ and $X_{TC}$.
4. The $S$, $P$ and $A$ Lagrangians above are quoted from Ecker:1988te without a local extract. An extract must be added before their normalizations are used in a numerical result.

## Inputs and provenance

- Goldstone, source and topological conventions: [Appendix A](../../extracts/Alda:2026nrd/appendix-a-expansion.md), [the LO topological record](../../extracts/Alda:2026nrd/lo-topological-action.md) [@Alda:2026nrd].
- $U(1)_A$ counting and topological source: [Kaiser:2000gs extract](../../extracts/Kaiser:2000gs/u1a-counting-and-topological-source.md) [@Kaiser:2000gs].
- NLO local basis: [lagrangian.md](../lagrangian.md) [@Kaiser:1998ds].
- Resonance transformation law and vector Lagrangian: [Ecker:1988te extract](../../extracts/Ecker:1988te/rcht-leading-lagrangian.md), [Ecker:1989yg extract](../../extracts/Ecker:1989yg/spin-one-representations.md).
