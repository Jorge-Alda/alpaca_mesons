---
title: Operator kinematics for the neutral pseudoscalar transition form factor
status: derived
created: 2026-10-05
updated: 2026-10-13
bibtex_keys: [Ruiz-Femenia:2003jdx, CLEO:1997fho, BaBar:2009rrj]
provenance: analytical-derivation
---

# Operator kinematics for $P\to\gamma^*\gamma$

## Conventions and inputs

Use the RChT VJP/VVP basis and Fourier convention in [the normalized amplitude record](radiative-tff-amplitude.md). All momenta are incoming, $p+q_1+q_2=0$, $p^2=m_P^2$, $q_1^2=s$, and $q_2^2=t$. The photon on the vector line has virtuality $s$; the other photon has virtuality $t$. The pseudoscalar and electromagnetic flavor inputs are $F_Q=\operatorname{tr}(\{Q,Q\}T_P)$ and $F_H=\operatorname{tr}(\{Q,Q\}\{T_P,2B_0m_q\})$.

## VJP numerator

After Fourier derivatives, building-block signs, and contraction of the operator indices, the numerator associated with the photon on the vector line is

$$
C(s,t)=F_Q\big[(c_1+c_2-c_5)m_P^2+(-c_1+c_2+c_5-2c_6)s+(c_1-c_2+c_5)t\big]+4c_3F_H.
$$

Equivalently, its individual coefficients are:

| Coupling | Coefficient in $C(s,t)$ |
|---|---|
| $c_1$ | $F_Q(m_P^2-s+t)$ |
| $c_2$ | $F_Q(m_P^2+s-t)$ |
| $c_3$ | $4F_H$ |
| $c_4$ | $0$ |
| $c_5$ | $F_Q(-m_P^2+s+t)$ |
| $c_6$ | $-2F_Qs$ |
| $c_7$ | $0$ |

The second photon placement contributes $C(t,s)$ with the other vector denominator. The VVP numerator is also symmetric in the photon virtualities and is given in the full amplitude record. Consequently a single operator placement does not determine the complete TFF sensitivity.

## Rank implication

The complete amplitude contains both photon placements, WZW, and VVP exchange. With $F_H=2m_P^2F_Q$, its TFF rows lie in the OPE span. For independent $F_H$, neutral TFF rows can resolve the additional $c_3,d_2$ combination left by the action-derived common-pole radiative rows. The executable rank comparison is in [`scripts/diagnose_radiative_rank.py`](../../../scripts/diagnose_radiative_rank.py); assumptions and numerical ranks are recorded in [radiative-rank-saturation.md](radiative-rank-saturation.md).

## Reproduction

```text
.venv/bin/python scripts/derive_radiative_tff_amplitude.py
.venv/bin/python scripts/diagnose_radiative_rank.py
```
