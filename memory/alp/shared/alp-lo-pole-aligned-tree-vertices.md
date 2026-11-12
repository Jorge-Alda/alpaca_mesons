---
title: Pole-aligned leading-order tree vertices
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Pole-aligned leading-order tree vertices

## Inputs and assumptions

Use natural units and the CCWZ Lagrangian, aligned vacuum, arbitrary source Taylor derivatives and all-incoming vertex convention of the [joint expansion](alp-lo-joint-vertex-expansion.md), following [@Alda:2026nrd]. The external states and full internal propagator are those of the [pole-projected construction](alp-lo-pole-projected-tree-amplitudes.md). The parent pole is simple, has a nonzero ALP component and has an invertible meson block at its pole; the [single-pole record](alp-single-pole-alignment.md) states the exceptional cases. The field expansion uses the local logarithm around $U=1$ and a complete orthonormal Hermitian $U(N)$ generator basis. The same approach applies to $U(3)$; charged fields can be formed as complex combinations of the real generator directions. Tadpoles vanish or are cancelled. The result is at LO and tree level, without expansion in the mixing; the limitations of the [spectrum](alp-lo-truncation-limits.md) apply.

## Derivation

1. For the pole vector $v^{(n)}$, set $w_i=v_i^{(n)}/v_a^{(n)}$ and use the shift $P_i=P'_i+w_i a'$ at linear order. The inverse of the exact axial transformation in the [generator record](alp-pole-alignment-chiral-generator.md) is
$$
 e^{i\Phi/f_\pi}=
 e^{i a'\sum_i w_i\lambda_i/(2f_\pi)}
 e^{i\Phi'/f_\pi}
 e^{i a'\sum_i w_i\lambda_i/(2f_\pi)},\qquad a=a'.
$$
The matrices need not commute. This substitution in the original CCWZ building blocks is equivalent to transforming the sources and fields together, by the identities in the [finite-redefinition record](alp-ccwz-finite-redefinitions.md). It retains the full source Taylor series rather than assuming linear sources.

2. Expand the matrix logarithm near the identity. Its linear term is $\Phi'+a'\sum_iw_i\lambda_i$. There is no quadratic coordinate term. Through cubic field order the nonlinear part is
$$
 \frac{a'^2}{24f_\pi^2}\left[\sum_iw_i\lambda_i,\left[\sum_jw_j\lambda_j,\Phi'\right]\right]
 +\frac{a'}{12f_\pi^2}\left[\Phi',\left[\sum_iw_i\lambda_i,\Phi'\right]\right].
$$
The script derives these terms by expanding the three matrix exponentials and the matrix logarithm, not by assuming commuting matrices. At a stationary vacuum, terms beyond cubic order in the coordinate map cannot contribute to vertices with at most four legs: they would multiply a one-point vertex. This is exact at that external-field order, not a small-mixing expansion.

3. The transformed quadratic and cubic vertices therefore contain the linear shift only. The transformed quartic vertex also receives the quadratic-Lagrangian vertex evaluated once with a cubic-map derivative and once with the remaining linear-map direction, summed over the four choices of that leg. Its kinetic and mass pieces are separately returned by the script. Every daughter pole vector is transformed with the inverse shift, including its ALP component; the exchange propagator is obtained from the transformed quadratic matrices. The full result includes all three exchange channels from the [pole-projected construction](alp-lo-pole-projected-tree-amplitudes.md).

4. For external pole vectors at their own on-shell momenta and conserved total momentum, the nonlinear quartic contact vanishes: each term contains the inverse quadratic operator acting on the remaining external pole vector. The kinetic and mass pieces need not vanish separately. This statement concerns this symmetric axial map; it does not imply that the full contact vertex or the exchange diagrams vanish. Off shell the nonlinear contact is generally nonzero.

## Verification and limitations

The tests use a complete $U(2)$ basis, noncommuting generators and an axial pole alignment, with the quadratic matrices and vertices generated from the same CCWZ sources. The source inputs include nonlinear mass and current coefficients; a separate linear counterterm cancels the illustrative ALP tadpole. The parent and one daughter both have nonzero ALP components. At two seeded, exactly on-shell algebraic momentum configurations, the full transformed four-leg amplitude, including nonzero exchange contributions, equals the original pole projection; the three-leg result also agrees. These configurations test the polynomial identities and are not asserted to lie in physical real-momentum phase space. This is a concrete noncommuting example, not an exhaustive numerical scan of all $U(3)$ sources.

A further three-leg $U(3)$ check uses charged generator directions and a Hermitian vector source supported only on the $d$–$s$ entries, with independent real and imaginary coefficients. It verifies the transformed vertex for an illustrative axial shift and retains a nonzero imaginary-source contribution. Its shift and kinematics are chosen algebraically rather than derived from the $U(3)$ spectrum, so this is a flavour and vertex check, not a second full pole-projected amplitude comparison. Independent checks of the cubic map cover both $a'P'^2$ and $a'^2P'$ terms.

For one exact four-leg configuration the nonlinear kinetic and mass contacts are respectively $-256/40545$ and $+256/40545$ in the stated input units, so their sum is zero. The linear-map contact is $188/901$, the combined exchanges are $-685579/54535728$, and the full reduced amplitude is $10693685/54535728$. With an independent off-shell change preserving momentum conservation, the nonlinear contact becomes $-16/40545$. These numbers illustrate the decomposition and are not phenomenological inputs.

[`derive_alp_pole_aligned_vertices.py`](../../../scripts/derive_alp_pole_aligned_vertices.py) calculates the cubic coordinate derivative, transformed vertices and complete transformed amplitude. [`test_derive_alp_pole_aligned_vertices.py`](../../../scripts/test_derive_alp_pole_aligned_vertices.py) checks pole alignment, a nonzero noncommuting map derivative, the exceptional zero-ALP parent, the separate contacts, and the comparisons stated above.

```sh
.venv/bin/python scripts/derive_alp_pole_aligned_vertices.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_pole_aligned_vertices.py
```
