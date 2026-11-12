---
title: Near-degeneracy of ALP-like and meson poles at LO, two-level formula
status: mixed
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Near-degeneracy of ALP-like and meson poles at LO, two-level formula

## Assumptions and inputs

- LO blocks and the two-level formula of the [finite-mixing record](finite-alp-meson-mixing.md) [@Alda:2026nrd], with canonical blocks defined there; poles and vectors from the [pole solver](alp-lo-pole-solver.md).
- Assumption: the working tolerance for "near-degenerate" is a stated, adjustable relative tolerance (default $10^{-2}$ of the largest $|m_n^2|$), not a physical criterion.
- New notation (proposed, pending approval): $\theta=\arctan\big(|v^{\rm lower}_a|/|v^{\rm upper}_a|\big)\in[0,\pi/2]$, with $v^{\rm lower}$, $v^{\rm upper}$ the pole vectors of the lower and upper of two ALP-carrying poles; $\theta=0$ means that the lower state has no ALP component, $\pi/4$ is maximal mixing.

## Method

1. Sort the poles and chain neighbouring poles into groups when they differ by at most the tolerance.
2. For each group report the spread, the minimal gap to the rest of the spectrum, a $K$-orthonormal basis of the group, and the ALP weight $\sum_{n\in G}(v^{(n)}_a)^2$, which does not depend on the basis. For an exactly degenerate group the basis is rotated so that at most one vector has a nonzero ALP component. Individual vectors are not unique; the sum of residues $\sum_{n\in G}v^{(n)}v^{(n)\mathsf T}$ is.
3. Flag a group as an ALP-meson near-degeneracy, where the small-mixing expansion of the [weak-mixing record](alp-lo-weak-mixing-expansion.md) fails, only if it contains the ALP-like pole (largest $|v_a|$) and at least two states with a nonzero ALP component. Groups with exactly one ALP-carrying state, such as the $K^0$ doublet for flavour-diagonal sources or the exact $K^0$ pole of the [reduction record](alp-lo-neutral-block-reduction.md), are reported without the flag.
4. Two-level estimate. For $K_{PP}=\mathbb 1$, let $u$ be the eigenvector of $M^2_{PP}$ with eigenvalue $m_P^2$ closest to the group (a degenerate meson eigenspace is rotated so that one combination couples to the ALP), and write $P=P'-K_{Pa}\,a$. The canonical blocks are $M^2_{P'a,\rm can}=u_i(M^2_{ia}-M^2_{ij}K_{ja})$ and $M^2_{aa,\rm can}=M^2_{aa}-2K_{ia}M^2_{ia}+K_{ia}M^2_{ij}K_{ja}$, and
$$
m_\pm^2=\frac{m_P^2+M^2_{aa,\rm can}}{2}\pm\frac12\sqrt{(M^2_{aa,\rm can}-m_P^2)^2+4(M^2_{P'a,\rm can})^2}.
$$
At $M^2_{aa,\rm can}=m_P^2$ the splitting is $2|M^2_{P'a,\rm can}|$ and $\theta=\pi/4$; in general $\theta=\tfrac12\operatorname{atan2}\big(2|M^2_{P'a,\rm can}|,\,M^2_{aa,\rm can}-m_P^2\big)$ for an isolated two-level system. With other mesons present the formula is an estimate; the full solution is the solver's.
5. The single-pole formula $w^{(a)}(s)=(M^2_{PP}-sK_{PP})^{-1}(sK_{Pa}-M^2_{Pa})$ is singular where $M^2_{PP}-m_n^2K_{PP}$ is singular (the exact $K^0$ pole; an ALP-like pole exactly degenerate with a decoupled meson). The solver reports this for every pole. There, use the full generalized eigenproblem, as stated in the finite-mixing record.

## Script and verification

[`scripts/diagnose_alp_meson_spectrum.py`](../../../scripts/diagnose_alp_meson_spectrum.py). For an isolated two-level system at degeneracy (exact rational test) the exact poles equal $m_P^2\mp|M^2_{P'a,\rm can}|$, the numerical splitting is $2|M^2_{P'a,\rm can}|$ ($3/500\ \mathrm{GeV}^2$ for the test input), $\theta=\pi/4$, exactly one ALP-meson near-degeneracy is found, and the $K^0$ doublet is reported without the flag. For the illustrative flavour-diagonal LO inputs the $K^0$ doublet is an exact degeneracy with zero ALP weight. Units: GeV$^2$ for gaps and splittings; angles in radians.

```bash
.venv/bin/python scripts/diagnose_alp_meson_spectrum.py
```

All checks pass.
