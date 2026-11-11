---
title: Dressed topological coefficient and recovery of $c_{GA}$
status: derived
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Dressed topological coefficient and recovery of $c_{GA}$

## Scope and inputs

The dressed mass function $\widehat M$ does not contain the gluonic coefficient. This record adds the singlet counterpart and shows how $c_{GA}=c_G+\operatorname{tr}c_q^A$ of Alda:2026nrd [@Alda:2026nrd] is recovered. Inputs: the anomaly law $\Theta'(x)=\Theta(x)-\frac1{2i}\log(\det g_R/\det g_L)$ from the [topological-sector record](alp-qcd-topological-cp.md); the connection law, transport matrices and $\widehat M$ from the [general record](alp-general-chiral-redefinition-invariants.md); the first-order $\widehat M_1$ from the [first-order record](alp-dressed-mass-first-order-invariants.md); and $X_{TC}=2i(\operatorname{tr}\log u+ic_Ga/f_a)$ with the tadpole source $C$ from the [Alda:2026nrd extracts](../../extracts/Alda:2026nrd/appendix-a-expansion.md). Same restricted truncation as the general record; $\Theta_1=c_G$.

## Derivation

Since $\partial_x\det W_\eta=i\operatorname{tr}C_\eta\det W_\eta$,
$$
\det W_\eta=\exp\!\Big(i\!\int_0^x\!\operatorname{tr}C_\eta\Big),\qquad
\frac1{2i}\log\frac{\det W_R}{\det W_L}=\int_0^x\operatorname{tr}C_A .
$$
Because $W_\eta'=g_\eta W_\eta$, the anomaly law gives $\Theta'=\Theta-\int\operatorname{tr}C_A\big|_{\rm shift}$, so the combination

$$
\boxed{\widehat\Theta(x)=\Theta(x)+\int_0^x\operatorname{tr}C_A(x')\,dx'}
$$

is invariant under every $g_{R,L}(x)$ with $g(0)=\mathbf1$, with Taylor coefficients $\widehat\Theta_n=\Theta_n+\operatorname{tr}C_{A,n-1}$. At first order

$$
\widehat\Theta_1=c_G+\operatorname{tr}c_q^A=c_{GA}.
$$

So $c_{GA}$ is the first Taylor coefficient of $\widehat\Theta$, just as $c_q^{VS},c_q^{AP}$ come from $\widehat M_1$. It cannot be read from $\widehat M$: $\widehat M$ gives only $\operatorname{tr}c_q^{AP}=-\tfrac12\operatorname{tr}(m_q^{-1}\widehat P_1)$.

## Relations to the dressed mass function

1. **Vacuum angle.** $\bar\theta(x)=\Theta(x)+\tfrac12\arg\det M(x)$ is invariant, and $\arg\det\widehat M=\arg\det M-2\int\operatorname{tr}C_A$, so $\widehat\Theta+\tfrac12\arg\det\widehat M=\bar\theta$. At first order $\bar\theta_1=c_G+\operatorname{tr}c_q^P=c_{GA}-\operatorname{tr}c_q^{AP}$. Hence $c_{GA}$ is the vacuum-angle slope plus the trace of the $\widehat M$ invariant. $c_{GA}$ and $\operatorname{tr}c_q^{AP}$ are independent invariants; $\bar\theta_1$ is their difference.
2. **Meson topological term.** With $\widehat U=W_R^\dagger UW_L$, $\operatorname{tr}\log\hat u=\operatorname{tr}\log u-i\int\operatorname{tr}C_A$, so
$$
\operatorname{tr}\log u+i\Theta=\operatorname{tr}\log\hat u+i\widehat\Theta .
$$
In the representative with vanishing connections, $X_{TC}=2i(\operatorname{tr}\log\hat u+i\widehat\Theta(x))$. Alda:2026nrd defines $X_{TC}$ only with the linear term $c_Ga/f_a$; replacing it by $\Theta(x)$ is an extension (assumption).
3. **Tadpole source.** Alda:2026nrd's $C=-B_0\{m_q,c_q^{AP}\}-\frac{m_0^2}{N_f}c_{GA}\mathbf1$ equals
$$
C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\,\mathbf1 ,
$$
i.e. it is built only from the first Taylor coefficients of the two dressed sources, $\widehat M$ and $\widehat\Theta$.

The photon analogue, $\widehat\Theta_\gamma$ with $\widehat\Theta_{\gamma,1}=c_{\gamma A}$, is in the [dressed photon record](../photon/alp-dressed-photon-invariant.md).

## Reproducible calculation

[`derive_alp_dressed_topological_invariant.py`](../../../scripts/derive_alp_dressed_topological_invariant.py) checks (exact rational $3\times3$ matrices, series through $O(x^4)$): the determinant of $W_\eta$; invariance of $\widehat\Theta$ under finite redefinitions with a nonlinear singlet part and with traceless generators; the Taylor form; the vacuum-angle relations; the $\operatorname{tr}\log$ identity; and the tadpole identity. These are series checks, not all-orders proofs. The identity $\widehat\Theta_1=c_{GA}$ holds by definition of $c_{GA}$ and is a consistency check, not an independent test. A unit-test file accompanies it. It imports helpers from `derive_alp_dressed_mass_first_order.py`.

```sh
.venv/bin/python scripts/derive_alp_dressed_topological_invariant.py
```
