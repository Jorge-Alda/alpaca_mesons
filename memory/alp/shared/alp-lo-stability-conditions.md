---
title: Stability conditions of the LO ALP-meson quadratic theory
status: mixed
created: 2026-11-11
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Stability conditions of the LO ALP-meson quadratic theory

## Assumptions and inputs

- LO blocks $K$, $M^2$ on $\Phi=(P_i,a)$ of the [finite-mixing record](finite-alp-meson-mixing.md) [@Alda:2026nrd]; $K_{ij}=\delta_{ij}$, $K_{aa}=1+K_{ia}K_{ia}$, real $K_{ia}$ (Hermitian $c_q^A$). The poles come from the [pole solver](alp-lo-pole-solver.md).
- A stationary CP-conserving vacuum is assumed, so that the quadratic terms are the Hessian at the stationary point.

## Conditions for a physical spectrum

(i) $K$ is positive definite (no ghost and no non-propagating field); (ii) $M^2$ is positive semidefinite (the stationary point is a minimum); (iii) all poles are real and nonnegative.

For $K>0$ the poles are the eigenvalues of $L^{-1}M^2L^{-\mathsf T}$ with $K=LL^{\mathsf T}$, so they are real. By Sylvester's law of inertia the number of negative poles equals the number of negative eigenvalues of $M^2$. Hence (ii) and (iii) are equivalent once (i) holds.

## $K>0$ at LO

1. With $K_{PP}=\mathbb 1$, the Schur complement gives $\det K=K_{aa}-K_{ia}K_{ia}=1$.
2. Write $\kappa=K_{ia}K_{ia}$. Then $\det(K-t)=(1-t)^{N-1}\big(t^2-(2+\kappa)t+1\big)$, with $N$ the number of fields. The two nontrivial eigenvalues have product $1$ and sum $2+\kappa>0$, so both are positive. The other eigenvalues equal $1$.
3. Consequently $K>0$ for every real $K_{ia}$. This is also the statement that $\mathcal L^{(2)}_{\rm kin}=\tfrac12(\partial P_i+K_{ia}\partial a)^2+\tfrac12(\partial a)^2$ of the finite-mixing record is a sum of squares.

At LO, stability therefore reduces to $M^2\ge0$.

## Bound on the bare ALP mass

If $M^2_{PP}>0$, then $M^2\ge0$ if and only if the Schur condition $M^2_{aa}-M^2_{aP}(M^2_{PP})^{-1}M^2_{Pa}\ge0$ holds. This is a lower bound on $m_{a,0}^2$ of order $f_\pi^2/f_a^2$ times a meson mass squared. Its closed form, its scaling, its dependence on the assumption that the quark spurion is linear in $a$, and the inputs of any numerical value are in the [bare-ALP-mass bound record](alp-lo-bare-alp-mass-bound.md). The unstable branch of the solver (slightly below the bound $M^2$ has a negative eigenvalue) raises a diagnostic that names the violated conditions instead of returning a spectrum. No numerical value of the bound is quoted here, because it depends on the sources and on $f_a$.

## Script and verification

[`scripts/diagnose_alp_meson_spectrum.py`](../../../scripts/diagnose_alp_meson_spectrum.py) verifies $\det K=1$ and the characteristic polynomial symbolically (also for symbolic charge-conserving sources, $N=5$), evaluates conditions (i)-(iii) independently, and tests the unstable branch just below the bound, for illustrative inputs only (the script prints them). Units: $K$ dimensionless; $M^2$, $m_n^2$ in GeV$^2$.

```bash
.venv/bin/python scripts/diagnose_alp_meson_spectrum.py
```

All checks pass.
