---
title: Second-order dressed photon coefficient $\widehat\Theta_{\gamma,2}$
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd, Bauer:2020jbp]
provenance: analytical-derivation
---

# Second-order dressed photon coefficient $\widehat\Theta_{\gamma,2}$

## Convention

Taylor coefficients are derivatives at $x=0$, $F_n=d^nF/dx^n|_0$ with $x=a/f_a$, as in the [dressed photon record](alp-dressed-photon-invariant.md) and the [second-order mass record](../shared/alp-dressed-mass-second-order-invariants.md).

## Inputs and assumptions

- Invariant $\widehat\Theta_\gamma(x)=\Theta_\gamma(x)+2N_c\int_0^x\operatorname{tr}[Q^2C_A(x')]\,dx'$, normalization $\mathcal L_{\gamma\gamma}=\frac{\alpha_{\rm em}}{4\pi}\Theta_\gamma F_{\mu\nu}\widetilde F^{\mu\nu}$ [@Bauer:2020jbp; @Alda:2026nrd]: [dressed photon record](alp-dressed-photon-invariant.md).
- Assumptions carried over from that record (not results): charge conservation, $[Q,M]=[Q,C_{R,L}]=0$; redefinitions commuting with $Q$, $[Q,g_{R,L}]=0$; the finite anomaly law $\Theta_\gamma'=\Theta_\gamma-N_c\operatorname{tr}[Q^2(G_R-G_L)]$ for $g_\eta=e^{iG_\eta}$. $Q=\operatorname{diag}(2/3,-1/3,-1/3)$, $N_c=3$.

## Result

The Taylor law $\widehat\Theta_{\gamma,n}=\Theta_{\gamma,n}+2N_c\operatorname{tr}(Q^2C_{A,n-1})$ of the dressed photon record gives at $n=2$
$$
\boxed{\widehat\Theta_{\gamma,2}=\Theta_{\gamma,2}+2N_c\operatorname{tr}\!\left(Q^2C_{A,1}\right)
=\Theta_{\gamma,2}+N_c\operatorname{tr}\!\left[Q^2\left(C_{R,1}-C_{L,1}\right)\right].}
$$
With $N_c=3$ this is $\Theta_{\gamma,2}+\tfrac83(C_{A,1})_{uu}+\tfrac23\big[(C_{A,1})_{dd}+(C_{A,1})_{ss}\big]$; only diagonal entries enter. As at first order, it is independent of the gluonic $\widehat\Theta_2=\Theta_2+\operatorname{tr}C_{A,1}$ of the [second-order topological record](../shared/alp-dressed-topological-second-order.md), because the weights $2N_cQ_b^2$ are not all equal. The first-order coefficient is $\widehat\Theta_{\gamma,1}=c_\gamma+2N_c\operatorname{tr}(Q^2c_q^A)=c_{\gamma A}$.

## Reproducible calculation and checks

[`derive_alp_second_order_invariants.py`](../../../scripts/derive_alp_second_order_invariants.py) checks the coefficient, its $C_{R,1}$, $C_{L,1}$ form and the $N_c=3$ component form; invariance under $x$-dependent Cayley-transform redefinitions commuting with $Q$ (noncommuting generators in the $d$-$s$ block, nonzero block traces so the raw $\Theta_{\gamma,2}$ changes), with the anomalous shift computed from block determinants and verified to equal $N_c\operatorname{tr}[Q^2(G_R-G_L)]$ through $x^2$; and the first-order regression. Redefinitions that do not commute with $Q$ are outside the assumptions and were not tested. Series checks with fixed seeds, not all-orders proofs.

```sh
.venv/bin/python scripts/derive_alp_second_order_invariants.py
```
