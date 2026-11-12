---
title: Weak-mixing expansion of the LO poles and pole vectors in $f_\pi/f_a$
status: mixed
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Weak-mixing expansion of the LO poles and pole vectors in $f_\pi/f_a$

## Assumptions and inputs

- Exact LO blocks of the [finite-mixing record](finite-alp-meson-mixing.md), solved by the [LO pole solver](alp-lo-pole-solver.md). Conventions of Alda:2026nrd [@Alda:2026nrd].
- Expansion parameter $\epsilon=f_\pi/f_a$ (new symbol, proposed, pending approval), varied at fixed $f_\pi$, $c_q^{A,P}$, $c_G$, $\chi_0=2B_0m_q$, $m_0^2$, $m_{a,0}^2$. Then $K_{ia}$ and $M^2_{ia}$ are proportional to $\epsilon$, $K_{aa}-1$ and $M^2_{aa}-m_{a,0}^2$ to $\epsilon^2$, and $K_{ij}=\delta_{ij}$ and $M^2_{ij}$ do not depend on $\epsilon$.
- A constant orthogonal rotation of the mesons (an allowed field change, see the [congruence record](alp-inverse-propagator-congruence.md)) makes $M^2_{PP}=\operatorname{diag}(m_{j}^2)$ and keeps this structure. Below, the index $j$ labels this meson mass basis and $m_{j,0}^2$ is the $\epsilon=0$ mass.
- Assumption: the meson masses $m_{j,0}^2$ are distinct from $m_{a,0}^2$ (no degeneracy with the ALP). Degenerate meson poles that are not mixed with each other, such as the $K^0$ doublet, are handled as a group: the ALP couples to one combination only.
- Mixed inverse-propagator entry: $D_{ja}(s)=sK_{ja}-M^2_{ja}=O(\epsilon)$.
- New notation (proposed, pending approval): $D^{\rm eff}_{aa}(s)=sK_{aa}-M^2_{aa}-\sum_jD_{ja}(s)^2/(s-m_{j,0}^2)$, the ALP entry of the inverse propagator after the mesons are eliminated (Schur complement), and $\delta D^{\rm eff}_{aa}(s)\equiv D^{\rm eff}_{aa}(s)-(s-m_{a,0}^2)$, its deviation from the bare ALP inverse propagator. $\delta D^{\rm eff}_{aa}(s)=s(K_{aa}-1)-2m_0^2c_G^2/N_f-\sum_jD_{ja}(s)^2/(s-m_{j,0}^2)$ is $O(\epsilon^2)$; here $M^2_{aa}-m_{a,0}^2=2m_0^2c_G^2/N_f$ (this $c_G$ term is already $O(\epsilon^2)$ since it carries $f_\pi^2/f_a^2$).

## Results

Order $\epsilon$ (first order; literature). Alda:2026nrd's no-mixing prescription gives, for the ALP pole, $m_a^2=m_{a,0}^2+O(\epsilon^2)$, meson components of the pole vector $(w_i^{(a)},1)$ solving the first-order equation $\sum_i w_i\lambda_i=2\epsilon\,[c_q^A+Y]$ with $\{\chi_0,Y\}/2+m_0^2\operatorname{tr}(Y)/N_f-m_{a,0}^2Y=C$, and a unit LSZ factor. See the [baryon-coupling record](../baryon/finite-mixing-baryon-couplings.md) and the [Appendix-A extract](../../extracts/Alda:2026nrd/appendix-a-expansion.md). In the meson mass basis this is $(m_{j,0}^2-s)\,w_j=D_{ja}(s)$ at $s=m_{a,0}^2$, that is $w_j=D_{ja}(m_{a,0}^2)/(m_{j,0}^2-m_{a,0}^2)$ (verified in Part B below). Here $C=-\{\chi_0,c_q^A-c_q^P\}/2-m_0^2(c_G+\operatorname{tr}c_q^A)/N_f$.

Higher orders (derived here; not given in Alda:2026nrd).

- Poles are even in $\epsilon$ ($a\to-a$ is a congruence that flips the sign of $\epsilon$). ALP-like pole:
$$
m_a^2=m_{a,0}^2-\delta D^{\rm eff}_{aa}(m_{a,0}^2)+\delta D^{\rm eff}_{aa}(m_{a,0}^2)\,\delta D^{\rm eff\,\prime}_{aa}(m_{a,0}^2)+O(\epsilon^6),
$$
with $\delta D^{\rm eff}_{aa}$ as defined above and the prime denoting $d/ds$. This follows from the zero of $D^{\rm eff}_{aa}(s)=(s-m_{a,0}^2)+\delta D^{\rm eff}_{aa}(s)$ by iteration.
- Meson-like pole $n$: $m_n^2=m_{n,0}^2+D_{na}(m_n^2)^2/(m_n^2-m_{a,0}^2)+O(\epsilon^4)$, and its ALP component is $v^{(n)}_a=D_{na}(m_n^2)/(m_{a,0}^2-m_n^2)\,(1+O(\epsilon^2))$ up to the overall sign, so it starts at $O(\epsilon)$; other meson components start at $O(\epsilon^2)$.
- ALP-like pole vector: $v_a=1+O(\epsilon^2)$.

Residual exponents. Compared with the first-order result, the residual of the ALP and meson poles is $O(\epsilon^2)$, that of the meson components of the ALP vector and of the ALP components of meson vectors is $O(\epsilon^3)$, and that of $v_a-1$ is $O(\epsilon^2)$. Compared with the second-order pole formula the residual is $O(\epsilon^4)$. Fitted exponents on $\epsilon=0.02,\dots,0.00125$ (illustrative flavour-diagonal sources): $1.97$ (ALP pole minus first order), $3.93$ (minus second order), $2.91$ (meson components of the ALP vector), $1.90$ ($v_a-1$), $1.96$ (meson poles minus first order); neutral sources give $1.96,3.92,2.90,1.89,1.92$ and $2.98$ for the ALP components of meson vectors.

## Derivation steps

1. Eliminate the meson rows of $D(s)\Phi=0$: $\Phi_j=-D_{ja}\Phi_a/(sK_{jj}-m_{j,0}^2)$ with $K_{jj}=1$; the remaining equation is $D^{\rm eff}_{aa}(s)\Phi_a=0$.
2. Expand $D^{\rm eff}_{aa}$ in $\epsilon$ at fixed $s$ to get the ALP pole condition above. Meson poles: $(s-m_{n,0}^2)D^{\rm eff}_{aa}(s)-D_{na}(s)^2=0$ near $s=m_{n,0}^2$ (other meson terms enter at higher order).
3. Build the null vectors from the same elimination and normalize with $K$.

## Script and verification

[`scripts/verify_alp_weak_mixing_limit.py`](../../../scripts/verify_alp_weak_mixing_limit.py). Part A (symbolic, three mesons): series of the pole conditions with symbolic $m_{j,0}^2$, $K_{ja}$, $M^2_{ja}$, confirming the absence of odd terms and of an $O(\epsilon)$ ALP-pole shift, the second-order formulas, the $\epsilon^4$ coefficient $\delta D^{\rm eff}_{aa}\,\delta D^{\rm eff\prime}_{aa}$, and nonzero next coefficients (so the residual exponents are exact). Part B (symbolic, all nine generators, symbolic Hermitian sources): the first-order pole equation equals the Alda:2026nrd equation $M_s(Y)=C$, with $M_s(X)=\{\chi_0,X\}/2+m_0^2\operatorname{tr}(X)/N_f-sX$ and $C=-\{\chi_0,c_q^A-c_q^P\}/2-m_0^2(c_G+\operatorname{tr}c_q^A)/N_f$. Part C (numerical, solver on the LO $6\times6$ blocks): the exponents above and the vanishing of first-order components on charged generators. Units: $f_\pi,f_a$ in GeV, masses squared in GeV$^2$, $\epsilon$ dimensionless.

```bash
.venv/bin/python scripts/verify_alp_weak_mixing_limit.py
```

All checks pass (about 20 s). Inputs are illustrative rationals.
