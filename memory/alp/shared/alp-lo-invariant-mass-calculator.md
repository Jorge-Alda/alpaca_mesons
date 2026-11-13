---
title: Finite-mixing LO masses from invariant dressed sources
status: derived
created: 2026-11-13
updated: 2026-11-13
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Finite-mixing LO masses from invariant dressed sources

## Assumptions and inputs

- Conventions and generator normalization of Alda:2026nrd [@Alda:2026nrd], $N_f=3$, $U=exp(i\Phi/f_\pi)$, $x=a/f_a$, and source arrays are derivatives at $x=0$.
- The aligned vacuum has diagonal $\widehat S_0$, $\widehat P_0=0$, and $\widehat\Theta_0=0$. The ALP tadpole is absent or canceled. Sources commute with electric charge; in particular $\widehat P_1$ may have a complex $d$-$s$ entry but no $u$-to-$d,s$ entries.
- Inputs are $f_a$, $m_{a,0}$, $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$, plus the independent LO constants $f_\pi$, $B_0$, and $m_0^2$. The latter constants cannot be determined from the requested invariant source lists. $\widehat S_0$ and $\widehat S_2$ have mass dimension one; $\widehat P_n$ and $\widehat\Theta_n$ are dimensionless; $B_0,f_a,f_\pi,m_{a,0}$ are in GeV and $m_0^2$ in GeV$^2$.
- The executable implementation is [`scripts/calculate_alp_meson_masses.py`](../../../scripts/calculate_alp_meson_masses.py). It requires $\widehat S_0,\widehat S_2,\widehat P_0,\widehat P_1,\widehat\Theta_0,\widehat\Theta_1$; other list entries do not enter the quadratic Lagrangian.

## Construction

1. The fixed zeroth-order mass spurion is $\chi_0=2B_0\widehat S_0$. The neutral meson basis is $(P_3,P_8,P_0,P_6,P_7)$ with the normalized Hermitian generators $\lambda_i$ satisfying $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$. The meson mass block is
$$
(M^2)_{ij}=\frac14\operatorname{tr}(\{\lambda_i,\lambda_j\}\chi_0)+\frac{m_0^2}{6}\operatorname{tr}\lambda_i\operatorname{tr}\lambda_j.
$$
2. The invariant first-order sources are $\widehat P_1=-\{c_q^{AP},m_q\}$ and $\widehat\Theta_1=c_{GA}$. The mixing column in the canonical field basis is
$$
(M^2_{\rm can})_{ia}=-\frac{f_\pi}{f_a}B_0\operatorname{tr}(\lambda_i\widehat P_1)+\frac{f_\pi}{f_a}\frac{m_0^2\widehat\Theta_1}{3}\operatorname{tr}\lambda_i.
$$
The kinetic mixing uses the corresponding axial source, but it is removed by the unit-determinant linear field change in the [exact pole-equation record](alp-lo-exact-pole-equations.md). The canonical meson block remains the one in Step 1, and its ALP mixing is equivalently expressed through the invariant source $\mathcal C=-\tfrac12\{\chi_0,c_q^{AP}\}-\tfrac{m_0^2}{3}c_{GA}\mathbb 1$. Direct diagonalization of the original generalized eigenvalue problem or of this canonical mass matrix gives the same poles.
3. The ALP diagonal entry in the invariant canonical frame is
$$
M^2_{\hat a\hat a,\mathrm{can}}=m_{a,0}^2+\left(\frac{f_\pi}{f_a}\right)^2\left[\frac{2m_0^2}{3}\widehat\Theta_1^2-B_0\operatorname{tr}\widehat S_2\right].
$$
This uses the $\widehat S_2$ derivative, not the coefficient of $x^2$; the Taylor expansion contains $x^2\widehat S_2/2!$.
4. Diagonalize the real symmetric canonical $6\times6$ neutral mass-squared matrix. Its six eigenvalues are the finite-mixing neutral poles, valid without expansion in $f_\pi/f_a$. Charged mesons decouple at LO; for diagonal $\widehat S_0$, their squared masses are $B_0[(\widehat S_0)_{uu}+(\widehat S_0)_{dd}]$ and $B_0[(\widehat S_0)_{uu}+(\widehat S_0)_{ss}]$ for $\pi^\pm$ and $K^\pm$ respectively.
5. The script returns square roots of nonnegative pole values in GeV. Neutral labels are assigned by maximizing total squared overlap with an unmixed reference basis; the two neutral-kaon labels denote the coupled and orthogonal real combinations. Individual labels are conventional at degeneracies or strong avoided crossings; the unordered pole set is invariant.

## Invocation and limits

```bash
PYTHONPATH=scripts .venv/bin/python -c 'import numpy as np; from calculate_alp_meson_masses import calculate_masses; z=np.zeros((3,3)); print(calculate_masses(1., .6, [np.diag([.003,.005,.1]),z,z], [z,z], [0.,0.], f_pi=.092, B_0=2.5, m_0_squared=.75))'
```

The function rejects a misaligned vacuum, non-Hermitian required source matrices, charge-violating $\widehat P_1$, and negative squared masses. It assumes the LO truncation and does not include loops, NLO terms, complex poles, or higher-derivative ALP effects. Source list indices are derivative orders. The formulas follow the [exact pole-equation record](alp-lo-exact-pole-equations.md), the [neutral-block record](alp-lo-neutral-block-reduction.md), the [mass sum-rule record](alp-lo-mass-sum-rules.md), and the [joint-vertex record](alp-lo-joint-vertex-expansion.md).
