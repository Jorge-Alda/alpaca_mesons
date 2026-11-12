---
title: Minimal LO neutral block and the exact $K^0$ pole
status: derived
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Minimal LO neutral block and the exact $K^0$ pole

## Assumptions and inputs

- The neutral $6\times6$ system on $(P_3,P_8,P_0,P_6,P_7,a)$ of the [charged-decoupling record](alp-lo-charged-meson-decoupling.md), with the LO blocks of the [finite-mixing record](finite-alp-meson-mixing.md) and the conventions of Alda:2026nrd [@Alda:2026nrd].
- Sources $c_q^A$, $c_q^P$, $\chi_0$ Hermitian and commuting with $Q$. Their $d$-$s$ entries are written $(c)_{ds}=\operatorname{Re}+i\operatorname{Im}$.
- A flavour rotation $V$ that commutes with $Q$ acts on the $d$-$s$ plane as $X\to VXV^\dagger$ on $c_q^A$, $c_q^P$, $\chi_0$. It is vector-like and does not change the anomaly term.
- Assumption: $m_d\neq m_s$.
- New notation (proposed, pending approval): $m_{K^0,0}^2\equiv\tfrac12\big[(\chi_0)_{dd}+(\chi_0)_{ss}\big]$, the diagonal entry of the $(P_6,P_7)$ block of $M^2$. It is the mass squared of the neutral kaon before mixing with the ALP, and equals $B_0(m_d+m_s)$ only if $\chi_0=2B_0m_q$ is not modified by the ALP couplings. It is not the physical kaon mass at finite mixing: the $K^0$ combination that couples to $a$ is shifted, and so can $\chi_0$ be when the ALP has flavour-violating couplings (see alp-sjp.6). Only the combination of $(P_6,P_7)$ that decouples from $a$ keeps a pole exactly at $m_{K^0,0}^2$.

## Results

(a) *Diagonal $\chi_0$.* Then $K_{ij}=\delta_{ij}$ and $M^2$ has no entries between $(P_3,P_8,P_0)$ and $(P_6,P_7)$. The $(P_6,P_7)$ block of $M^2$ is $m_{K^0,0}^2\,\mathbb 1$. The flavour-diagonal mesons and the $K^0$ doublet mix only through $a$. With a $d$-$s$ entry in $\chi_0$ the entries between $(P_3,P_8,P_0)$ and $(P_6,P_7)$ are nonzero, but the $(P_6,P_7)$ block stays proportional to the unit matrix.

(b) *CP condition.* The ALP couples to the $K^0$ doublet through $K_{ia}$ and $M^2_{ia}$ with $i=6,7$. Their $2\times2$ determinant is
$$
K_{6a}M^2_{7a}-K_{7a}M^2_{6a}
=2\Big(\frac{f_\pi}{f_a}\Big)^2\big[(\chi_0)_{dd}+(\chi_0)_{ss}\big]\operatorname{Im}\!\big(c^A_{ds}\,c^{P*}_{ds}\big).
$$
A $K^0$ combination with zero $K_{ia}$ and zero $M^2_{ia}$ exists if and only if this vanishes. For real $d$-$s$ entries (CP conservation in the convention $\Phi\to-\Phi^{\mathsf T}$, $a\to-a$, under which $P_7$ is CP-even and $P_3,P_8,P_0,P_6,a$ are CP-odd) one has $K_{7a}=M^2_{7a}=0$ and $P_7$ decouples. For aligned complex phases the $s$-quark rephasing $V=\operatorname{diag}(1,1,e^{i\varphi})$ is a constant rotation of $(P_6,P_7)$, keeps $\chi_0$ diagonal and makes both entries real.

(c) *Exact pole.* For any phases, $s=m_{K^0,0}^2$ is a pole. With $D_{ja}(m_{K^0,0}^2)=m_{K^0,0}^2K_{ja}-M^2_{ja}$, the vector with components $(-D_{7a},\,D_{6a})$ on $(P_6,P_7)$ and zero elsewhere is a null vector of $D(m_{K^0,0}^2)$, because $D_{PP}(m_{K^0,0}^2)$ vanishes on $(P_6,P_7)$ and has no entries to the other mesons, while the $a$ row gives $-D_{6a}D_{7a}+D_{7a}D_{6a}=0$. Its ALP component is zero, so it does not contribute to the $a$-$a$ propagator.

(d) *Rotation condition.* A $d$-$s$ rotation separates $(P_3,P_8,P_0)$ from $(P_6,P_7)$ (coupled only through $a$) if and only if, in the rotated basis, the $d$-$s$ entries of $\chi_0$, $c_q^A$, $c_q^P$ all vanish, that is, if and only if their $d$-$s$ blocks commute pairwise. For $m_d\neq m_s$ the only rotations keeping $\chi_0$ diagonal are phases. Diagonalizing the $d$-$s$ block of $c_q^A$ when it does not commute with $\chi_0$ produces a $d$-$s$ entry in $\chi_0$ and therefore direct $M^2$ entries between $(P_3,P_8,P_0)$ and $(P_6,P_7)$.

## Minimal block size ($\chi_0$ diagonal, $m_d\neq m_s$)

| Sources | Block |
|---|---|
| generic, $\operatorname{Im}(c^A_{ds}c^{P*}_{ds})\neq0$ | $6\times6$ on $(P_3,P_8,P_0,P_6,P_7,a)$; five poles are zeros of $\det D(s)/(s-m_{K^0,0}^2)$ and one pole is at $m_{K^0,0}^2$ |
| $\operatorname{Im}(c^A_{ds}c^{P*}_{ds})=0$ (real entries, or $c^A_{ds}=0$ or $c^P_{ds}=0$) | $5\times5$ on $(P_3,P_8,P_0,P_6',a)$ after rephasing, plus a decoupled combination $P_7'$ at $m_{K^0,0}^2$ |
| $c^A_{ds}=c^P_{ds}=0$ | $4\times4$ on $(P_3,P_8,P_0,a)$, plus the degenerate pair $(P_6,P_7)$ at $m_{K^0,0}^2$ |

## Script and verification

[`scripts/derive_alp_neutral_block_reduction.py`](../../../scripts/derive_alp_neutral_block_reduction.py) checks (a)-(d) symbolically or with exact rationals; the exact pole is also found numerically by [`scripts/verify_alp_lo_spectrum.py`](../../../scripts/verify_alp_lo_spectrum.py). Units: $f_\pi,f_a$ in GeV; $\chi_0,m_0^2,m_{a,0}^2,M^2,s$ in GeV$^2$.

```bash
.venv/bin/python scripts/derive_alp_neutral_block_reduction.py
```

All checks pass. The exact pole with zero ALP component is also the reason why the single-pole formula of the finite-mixing record is singular there; see [the degeneracy record](alp-lo-near-degeneracy-two-level.md).
