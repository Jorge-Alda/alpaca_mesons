---
title: Second-order dressed gluonic coefficient $\widehat\Theta_2$
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Second-order dressed gluonic coefficient $\widehat\Theta_2$

## Convention

Taylor coefficients are derivatives at $x=0$, $F_n=d^nF/dx^n|_0$ with $x=a/f_a$, as in the [dressed topological record](alp-dressed-topological-invariant.md) and the [second-order mass record](alp-dressed-mass-second-order-invariants.md).

## Inputs and assumptions

- Invariant $\widehat\Theta(x)=\Theta(x)+\int_0^x\operatorname{tr}C_A(x')\,dx'$, with $\Theta$ in the normalization of Alda:2026nrd [@Alda:2026nrd], and the anomaly law $\Theta'=\Theta-\frac1{2i}\log(\det g_R/\det g_L)$: [dressed topological record](alp-dressed-topological-invariant.md), [QCD topological record](alp-qcd-topological-cp.md).
- $C_{A,n}=\tfrac12(C_{R,n}-C_{L,n})$, Hermitian; $g_\eta(0)=\mathbf1$; same truncation as the [general record](alp-general-chiral-redefinition-invariants.md).

## Result

The general Taylor law $\widehat\Theta_n=\Theta_n+\operatorname{tr}C_{A,n-1}$ of the dressed topological record gives at $n=2$
$$
\boxed{\widehat\Theta_2=\Theta_2+\operatorname{tr}C_{A,1}=\Theta_2+\tfrac12\operatorname{tr}\left(C_{R,1}-C_{L,1}\right).}
$$
Only the trace of the first derivative of the axial current coefficient enters; $\widehat\Theta_2$ does not depend on the mass sources or on $C_{\eta,0}$. The first-order coefficient is $\widehat\Theta_1=c_G+\operatorname{tr}c_q^A=c_{GA}$.

## Relation to the dressed mass source

The vacuum-angle combination $\bar\theta=\Theta+\tfrac12\arg\det M$ satisfies $\widehat\Theta+\tfrac12\arg\det\widehat M=\bar\theta$ (dressed topological record). At second order this links $\widehat\Theta_2$ to $\widehat S_2$, $\widehat P_2$ of the [second-order mass record](alp-dressed-mass-second-order-invariants.md); the identity was checked at $O(x^2)$ with exact rational matrices.

## Reproducible calculation and checks

[`derive_alp_second_order_invariants.py`](../../../scripts/derive_alp_second_order_invariants.py) checks: the coefficient and its $C_{R,1}$, $C_{L,1}$ form; invariance under general $x$-dependent Cayley-transform redefinitions with $\operatorname{tr}H_1\neq0$, so that the anomalous shift is nonzero and the raw $\Theta_2$ changes; the first-order regression $\widehat\Theta_1=c_G+\operatorname{tr}c_q^A$; and the vacuum-angle identity at $O(x)$ and $O(x^2)$. These are series checks with fixed seeds; $\log\det g$ is expanded with $\log\det g(0)=0$.

```sh
.venv/bin/python scripts/derive_alp_second_order_invariants.py
```
