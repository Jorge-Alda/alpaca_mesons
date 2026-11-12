---
title: Decoupling of charged mesons from the neutral ALP system at LO
status: derived
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Decoupling of charged mesons from the neutral ALP system at LO

## Assumptions and inputs

- LO blocks $K$, $M^2$ of the [finite-mixing record](finite-alp-meson-mixing.md) in the conventions of Alda:2026nrd [@Alda:2026nrd]: $N_f=3$, generators $\lambda_0,\dots,\lambda_8$ with $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, fields $\Phi_I=(P_0,\dots,P_8,a)$, $\mathcal L_2=\tfrac12\partial\Phi^{\mathsf T}K\partial\Phi-\tfrac12\Phi^{\mathsf T}M^2\Phi$.
- Assumption (charge conservation of the sources): $c_q^A$, $c_q^P$ and $\chi_0$ are Hermitian and commute with the quark charge matrix $Q=\operatorname{diag}(2/3,-1/3,-1/3)$. Equivalently their $u$-$d$ and $u$-$s$ entries vanish; the diagonal entries and the $d$-$s$ entry are arbitrary, and $\chi_0$ need not be diagonal.
- No NLO terms, no loops, canonical bare ALP term (assumptions of the finite-mixing record).

## Statement

Under these assumptions every entry of $K$ and $M^2$ between a charged meson ($P_1,P_2=\pi^\pm$; $P_4,P_5=K^\pm$) and any of $P_0,P_3,P_6,P_7,P_8,a$ vanishes, including the $m_0^2$ and $c_G$ terms. The nine-generator system therefore splits into the charged block on $(P_1,P_2,P_4,P_5)$ and the neutral block on $(P_3,P_8,P_0,P_6,P_7,a)$. The neutral block is the default basis of the [LO pole solver](alp-lo-pole-solver.md).

## Derivation

1. The generators commuting with $Q$ are $\lambda_0,\lambda_3,\lambda_6,\lambda_7,\lambda_8$ (neutral). The generators $\lambda_1,\lambda_2,\lambda_4,\lambda_5$ lie in the image of $\operatorname{ad}_Q=[Q,\cdot\,]$.
2. For $Y$ commuting with $Q$ and any $Z$, $\operatorname{tr}([Q,Z]Y)=\operatorname{tr}(Z[Y,Q])=0$. Hence $\operatorname{tr}(\lambda_cY)=0$ for charged $\lambda_c$. The entries $K_{ia}\propto\operatorname{tr}(c_q^A\lambda_i)$ and $M^2_{ia}\supset\operatorname{tr}(\lambda_i\{c_q^P,\chi_0\})$ are of this form, because $c_q^A$ and $\{c_q^P,\chi_0\}$ commute with $Q$. The $c_G$ and $m_0^2$ terms contain $\operatorname{tr}(\lambda_c)=0$.
3. Because $\operatorname{ad}_Q$ is a derivation, $\{[Q,Z],\lambda_n\}=[Q,\{Z,\lambda_n\}]$ for neutral $\lambda_n$. So $\{\lambda_c,\lambda_n\}$ is also in the image of $\operatorname{ad}_Q$ and $M^2_{cn}=\tfrac14\operatorname{tr}(\{\lambda_c,\lambda_n\}\chi_0)=0$.
4. $K_{aa}=1+\sum_iK_{ia}K_{ia}$ contains all nine generators. For charge-conserving $c_q^A$ the charged $K_{ia}$ vanish, so the $6\times6$ $K_{aa}$ equals $1+\sum$ over the five neutral $K_{ia}$.

Without $Q$-commutation the statement fails: the charged-neutral entries are linear in the $u$-$d$ and $u$-$s$ entries of $c_q^A$, $\{c_q^P,\chi_0\}$ and $\chi_0$, and vanish if and only if those entries vanish. Single-entry counterexamples are checked in the script (a $u$-$d$ entry in $c_q^A$ gives $K_{P_1a}\neq0$; a $u$-$s$ entry in $c_q^P$ gives $M^2_{P_4a}\neq0$; a $u$-$d$ entry in $\chi_0$ gives $M^2_{P_1P_8}\neq0$). Charge-violating sources are outside the scope of the neutral reduction; the solver must reject them instead of dropping the charged mixing.

## Script and verification

[`scripts/derive_alp_charged_meson_decoupling.py`](../../../scripts/derive_alp_charged_meson_decoupling.py) builds the ten-field Hessian of the explicit LO Lagrangian with fully symbolic sources and compares it with the closed forms. The blocks and their $6\times6$ restriction are built and cross-checked by [`scripts/derive_alp_lo_neutral_blocks.py`](../../../scripts/derive_alp_lo_neutral_blocks.py), which also checks that the charged-neutral entries are nonzero for a charge-violating example. Units: $f_\pi,f_a$ in GeV; $\chi_0,m_0^2,m_{a,0}^2$ in GeV$^2$; couplings dimensionless.

```bash
.venv/bin/python scripts/derive_alp_charged_meson_decoupling.py
.venv/bin/python scripts/derive_alp_lo_neutral_blocks.py
```

All checks pass. The numerical inputs of the second script are exact rationals chosen for testing; they are not fitted values.
