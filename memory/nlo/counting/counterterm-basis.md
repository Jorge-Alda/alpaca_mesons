---
title: Renormalized counterterm basis of the U(3) NLO theory with explicit resonances
status: mixed
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Kaiser:1998ds, Herrera-Siklody:1996tqr, Gasser:1984gg, Ecker:1988te, Alda:2026nrd]
provenance: analytical-derivation
---

# Counterterm basis

## Assumptions and inputs

- Orders come from [power-counting.md](power-counting.md). Loop content is from [loop-classes.md](loop-classes.md). The scheme and scales are from [subtraction-scheme-and-scales.md](subtraction-scheme-and-scales.md).
- The operators are those of [lagrangian.md](../lagrangian.md), in the reduced $V_4=0$ basis, together with the source operator $O_{47}$ of [the singlet-source derivative basis](../singlet-source-derivative-basis.md).
- The constants that enter the quadratic sector are listed in [one-loop-effects-and-renormalized-lecs.md](../one-loop-effects-and-renormalized-lecs.md).
- The resonance contributions are the tree-level assignments of [resonance-saturation-and-dynamical-masses.md](../resonance-saturation-and-dynamical-masses.md) [@Ecker:1988te]. When resonances are explicit, these contributions are **not** part of $L_i^r(\mu)$ (see the resonance-exchange record of this folder).
- Every entry in the $\delta$, $\Gamma$ and $Z_A$ columns is asserted by `--orders` and `--scheme` of [`scripts/derive_nlo_u3_resonance_counting.py`](../../../scripts/derive_nlo_u3_resonance_counting.py):

```bash
.venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --orders --scheme
```

## Table

| Constant | Operator | $\delta$ order | $\Gamma_i$ ($U(3)$) | QCD scale | Tree-level resonance contribution | Remarks |
|---|---|---|---|---|---|---|
| $F$ | $\frac{F^2}{4}\langle u_\mu u^\mu\rangle$ | $\delta^0$ | none | invariant | — | no chiral-loop pole; R-int shifts it at $\delta$ (excluded, fixed from data) |
| $B$ | $\frac{F^2}{4}\langle\chi_+\rangle$ | $\delta^0$ | $B=B^r(1+4\tau\lambda/F^4)$ | $Z_M$ | — | $\eta'$ tadpole, $\delta^2$ |
| $\tau$ | $-\frac{f_\pi^2m_0^2}{4N_f}X_{TC}^2$ | $\delta^0$ | none | $Z_A^2$ | — | |
| $L_5$ | $\langle u_\mu u^\mu\chi_+\rangle$ | $\delta$ | $3/8$ | in $2L_5+3L_{18}$ | $S_8$ | |
| $L_8$ | $\frac12\langle\chi_+^2+\chi_-^2\rangle$ | $\delta$ | $3/16$ | in $2L_8-3L_{25}$ | $S_8$, $P_8$ | leading value; next term $-v_{31}/2$ is $\delta^3$ |
| $\Lambda_1$ | $\frac{f_\pi^2}{12}\langle u_\mu\rangle\langle u^\mu\rangle$ | $\delta$ | none | $1+\Lambda_1\to Z_A^2(1+\Lambda_1)$ | not fixed by standard saturation | singlet kinetic term |
| $\Lambda_2$ | $i\frac{f_\pi^2}{12}X_{TC}\langle\chi_-\rangle$ | $\delta$ | none | $1+\Lambda_2\to Z_A(1+\Lambda_2)$ | not fixed by standard saturation | singlet/topological mass term |
| $L_4$ | $\langle u_\mu u^\mu\rangle\langle\chi_+\rangle$ | $\delta^2$ | $1/8$ | invariant | $S_8$, $S_1$ (cancel in the nonet limit) | |
| $L_6$ | $\langle\chi_+\rangle^2$ | $\delta^2$ | $1/16$ ($SU(3)$: $11/144$) | invariant | $S_8$, $S_1$ (cancel in the nonet limit) | |
| $L_7$ | $\langle\chi_-\rangle^2$ | $\delta^2$ | $0$ | invariant | $P_8$, $P_1$ | $\mathcal O(N_c^0)$ in $U(3)$; $\eta'$ is not integrated out |
| $L_{18}$ | $\langle u_\mu\rangle\langle u^\mu\chi_+\rangle$ | $\delta^2$ | $-1/4$ | $2L_5+3L_{18}\to Z_A(2L_5+3L_{18})$ | not fixed by standard saturation | |
| $L_{25}$ | $-iX_{TC}\langle\chi_+\chi_-\rangle$ | $\delta^2$ | $0$ | $2L_8-3L_{25}\to Z_A(2L_8-3L_{25})$ | not fixed by standard saturation | $\Gamma_{25}=0$ quoted by Kaiser:1998ds, not reproduced from Herrera-Siklody:1996tqr |
| $L_{47}$ | $D_\mu\hat\theta\langle C^\mu\rangle\langle M\rangle$ | $\delta^3$ ($N_c^{-1}p^4$) | $0$ at the retained $1/N_c$ order | — | none known | **beyond the declared order**; $B_{47}(0)=v_{12}(v_{31}-3)=\mathcal O(N_c^{-2})$ |
| *further ALP-source operators* | *to be fixed in alp-gii.8* | — | — | — | — | **marked row**: needed only if the source specialization of alp-gii.8 generates operators outside this basis |

Vector and axial-vector nonets contribute at tree level to $L_1$, $L_2$, $L_3$, $L_9$ and $L_{10}$. None of these enters the quadratic sector, but they enter ALP vertices with external vector sources [@Ecker:1988te].

## How the table is used

1. **Absorbing divergences.** The UV divergences of the included loop classes, G-LO and G-ALP, are absorbed by the constants with $\Gamma_i\neq0$ and by $B$. In a two-point function, $\mathcal Q_{IJ}(s)$ of [one-loop-effects-and-renormalized-lecs.md](../one-loop-effects-and-renormalized-lecs.md), the $\mu$ dependence of $\Sigma^{\rm loop}_{IJ}(s;\mu)$ cancels against that of the local terms built from $L_i^r(\mu)$ and $B^r(\mu)$ (alp-gii.2).
2. **Constants without a chiral pole.**
   - $L_7$ and $\Lambda_{1,2}$: zero at this order without further qualification.
   - $L_{25}$: zero as quoted by Kaiser:1998ds.
   - $L_{47}$: zero only under the large-$N_c$ truncation of Herrera-Siklody:1996tqr.

   These zeros fail at the next $1/N_c$ order, which is $\delta^3$, beyond the declared order.
3. **$L_{47}$ is beyond the declared order.** $O_{47}$ has two flavor traces and one factor of $\hat\theta$, and $D_\mu\hat\theta$ counts as $p$. So $L_{47}\sim N_c^{-1}$, in agreement with the general rule (4.1) of Herrera-Siklody:1996tqr, and the operator is $N_c^{-1}p^4$, that is, $\delta^3$. It is dropped at $\mathcal O(\delta^2)$, together with its counterterm. The quadratic records keep it only as a marker for a calculation at higher order ([singlet-source derivative basis](../singlet-source-derivative-basis.md)). The counting assumes that rule (4.1) also applies to derivatives of $\hat\theta$; the paper states the rule generally but does not apply it to $L_{47}$ explicitly.
4. **Inputs.** Under the exclusion of resonance loops, every constant in the table, and the resonance masses and couplings, must be fixed from data ([loop-classes.md](loop-classes.md)). The QCD scale of each singlet input must be stated ([subtraction-scheme-and-scales.md](subtraction-scheme-and-scales.md)).

## Verification

Each $\delta$ order is asserted by `--orders`. Each $\Gamma_i$, each flag and each QCD-scale law is asserted by `--scheme`. The resonance column is copied from [resonance-saturation-and-dynamical-masses.md](../resonance-saturation-and-dynamical-masses.md), and its coefficient-by-coefficient check is part of alp-gii.3.
