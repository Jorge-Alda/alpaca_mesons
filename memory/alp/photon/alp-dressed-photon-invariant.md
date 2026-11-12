---
title: Dressed photon coefficient and its expansion in $1/f_a$
status: derived
created: 2026-11-11
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd, Bauer:2020jbp, Bai:2025fvl]
provenance: analytical-derivation
---

# Dressed photon coefficient and its expansion in $1/f_a$

## Scope and inputs

This record adds the photon to the restricted partonic Lagrangian of the [general record](../shared/alp-general-chiral-redefinition-invariants.md) and generalizes the invariant $c_{\gamma A}\equiv c_\gamma^{\rm inv}=c_\gamma+2N_c\operatorname{tr}(Q^2c_q^A)$ of the [photon-anomaly extract](../../extracts/Alda:2026nrd/photon-anomaly-normalization.md) to arbitrary $g_{R,L}(x)$. It is the photon analogue of the [gluonic record](../shared/alp-dressed-topological-invariant.md). Normalization: $\mathcal L_{\gamma\gamma}=\frac{\alpha_{\rm em}}{4\pi}\Theta_\gamma(x)F_{\mu\nu}\widetilde F^{\mu\nu}$ with $\Theta_{\gamma,1}=c_\gamma$ [@Bauer:2020jbp; @Alda:2026nrd]. The photon couples to quarks minimally with fixed charge $Q$. No photon dipoles, no CP-even $F^2$ coefficient (it has no anomalous shift), no weak or $Z$ terms, no vector mesons.

## Assumptions (not results)

1. **Charge conservation.** The ALP is neutral, so $[Q,M(x)]=[Q,C_R(x)]=[Q,C_L(x)]=0$.
2. **Redefinitions preserve the photon coupling.** $[Q,g_{R,L}(x)]=0$. For $Q=\operatorname{diag}(2/3,-1/3,-1/3)$ this is $U(1)_u\times U(2)_{ds}$, with noncommuting generators inside the $ds$ block. A rotation mixing $u$ with $d$ or $s$ turns the photon coupling into an $x$-dependent chiral charge $g_{R,L}Qg_{R,L}^\dagger$; that case needs chiral charge matrices and counterterms and is not treated here.
3. **Finite anomaly law.** The extract gives only the linear law $c_\gamma\to c_\gamma-2N_c\operatorname{tr}(Q^2\kappa)$. Integrating the infinitesimal anomaly $d\Theta_\gamma=-N_c\operatorname{tr}[Q^2(d\alpha_R-d\alpha_L)]$ along $x$, with $\operatorname{tr}_b(g\,dg^\dagger)=d\log\det_b g^\dagger$ in each charge block $b$, gives
$$
\Theta_\gamma'(x)=\Theta_\gamma(x)-2N_c\,\frac1{2i}\sum_bQ_b^2\log\frac{\det_bg_R}{\det_bg_L}
=\Theta_\gamma(x)-N_c\operatorname{tr}\!\big[Q^2(G_R-G_L)\big]
$$
for $g_\eta=e^{iG_\eta}$.

## Derivation

Because $[Q,C_\eta]=0$, the transport matrices commute with $Q$, so $W_\eta^\dagger QW_\eta=Q$: the photon connection $eQA_\mu$ is unchanged by the dressing. With $[g_\eta,Q]=0$,
$$
\operatorname{tr}(Q^2C_\eta')=\operatorname{tr}(Q^2C_\eta)-i\,\partial_x\sum_bQ_b^2\log\det{}_bg_\eta ,
$$
which cancels the anomalous shift of $\Theta_\gamma$. Hence

$$
\boxed{\widehat\Theta_\gamma(x)=\Theta_\gamma(x)+2N_c\int_0^x\operatorname{tr}\!\big[Q^2C_A(x')\big]dx'}
$$

is invariant, and

$$
\widehat\Theta_{\gamma,n}=\Theta_{\gamma,n}+2N_c\operatorname{tr}\!\big(Q^2C_{A,n-1}\big).
$$

## Expansion in $1/f_a$

At order $1/f_a$,

$$
\widehat\Theta_{\gamma,1}=c_\gamma+2N_c\operatorname{tr}(Q^2c_q^A)=c_{\gamma A}
=c_\gamma+\tfrac83c^A_{uu}+\tfrac23\left(c^A_{dd}+c^A_{ss}\right)\quad(N_c=3).
$$

Only the diagonal entries of $c_q^A$ enter, since $Q^2$ is diagonal. The next coefficient is $\widehat\Theta_{\gamma,2}=\Theta_{\gamma,2}+2N_c\operatorname{tr}(Q^2C_{A,1})$, where $C_{A,1}=ic_q^{A,2}$-type terms are the order-$x$ axial coefficient of the general record; see the [second-order record](alp-dressed-photon-second-order.md).

## Relations to the other dressed sources

1. **Independence from $c_{GA}$.** On the diagonal of $c_q^A$, $\widehat\Theta_1$ has weights $(1,1,1)$ and $\widehat\Theta_{\gamma,1}$ has weights $2N_c(Q_u^2,Q_d^2,Q_s^2)=(8/3,2/3,2/3)$. The map has rank 2, so $c_{\gamma A}$ and $c_{GA}$ are independent invariants; the rank drops to 1 only if all $Q_b^2$ are equal. If $C_A$ lives only in the $ds$ block, $\widehat\Theta_\gamma-\Theta_\gamma=\tfrac23(\widehat\Theta-\Theta)$.
2. **Mixing contribution.** For any $\delta,\kappa$-invariant matrix $X$, $c_\gamma+2N_c\operatorname{tr}[Q^2(c_q^A+X)]=\widehat\Theta_{\gamma,1}+2N_c\operatorname{tr}(Q^2X)$. With Alda:2026nrd's $X=\mathbb M_a^{-1}(C)$ and $C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf1$ (gluonic record), the small-mixing result of the [photon mixing record](alp-photon-invariant-from-mixing.md) reads
$$
C_{a\gamma\gamma}(m_a)=\widehat\Theta_{\gamma,1}+2N_c\operatorname{tr}\!\big[Q^2\,\mathbb M_a^{-1}(C)\big],
$$
i.e. it depends only on the first Taylor coefficients of $\widehat M$, $\widehat\Theta$ and $\widehat\Theta_\gamma$. This uses that record's LSZ framework; it is not a finite-mixing statement.
3. **Meson-level photon term.** The WZW photon vertex of the mixing record is not re-derived under $g(x)$ here; matching to $\widehat\Theta_\gamma$ relies on the anomaly consistency argued in [@Bai:2025fvl] and is an assumption beyond this record.

## Reproducible calculation

[`derive_alp_dressed_photon_invariant.py`](../../../scripts/derive_alp_dressed_photon_invariant.py) uses exact rational $3\times3$ matrices and series through $O(x^4)$. It checks: $[Q,W_\eta]=0$ and $W_R^\dagger QW_R=Q$; invariance of $\widehat\Theta_\gamma$ (and $\widehat\Theta$) under nonlinear redefinitions in the commutant of $Q$, including generators traceless in each block (no $\Theta_\gamma$ shift); the Taylor form and $c_{\gamma A}$; the weight-map rank; and the identity in item 2. Also checked: for a $u$–$d$ rotation, $gQg^\dagger\neq Q$, the new charges differ for R and L, and the naive $\widehat\Theta_\gamma$ is not invariant, with residual starting at $O(x^3)$ for the tested seeds (so $\widehat\Theta_{\gamma,1}$, $\widehat\Theta_{\gamma,2}$ are unchanged there). These are series checks, not all-orders proofs; $\widehat\Theta_{\gamma,1}=c_{\gamma A}$ holds by definition, so it is a consistency check.

```sh
.venv/bin/python scripts/derive_alp_dressed_photon_invariant.py
```
