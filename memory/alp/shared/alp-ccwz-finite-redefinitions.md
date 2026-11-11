---
title: Finite chiral redefinitions of CCWZ building blocks
status: derived
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Finite chiral redefinitions of CCWZ building blocks

## Inputs and assumptions

Use natural units, $x=a/f_a$, $U=u^2=\exp(i\Phi/f_\pi)$ and the CCWZ definitions in the [arbitrary-order expansion](alp-ccwz-arbitrary-order-expansion.md), following `Alda:2026nrd` [@Alda:2026nrd]. The local source transformations and restricted partonic Lagrangian are given in the [source-transformation record](alp-chiral-source-transformations.md); omitted operators and loops are not included. Flavor matrices are unitary and need not commute. Choose the square root $u'$ continuously near the original vacuum; the following expressions hold on that local branch. Derivatives may act on all spacetime-dependent matrices. $U$, $u$ and the chiral matrices are dimensionless; $r_\mu$ and $\ell_\mu$ have mass dimension one, and $\chi$ has mass dimension two.

## Derivation

1. The nonlinear field and its square root transform as
$$
U'=g_R U g_L^\dagger,\qquad
u'=g_R u h^\dagger=h u g_L^\dagger,\qquad
h=u'^\dagger g_R u=u'g_Lu^\dagger.
$$
The two forms of $h$ agree because $u'^2=g_Ru^2g_L^\dagger$. This defines the compensator to all orders in the fields without expanding the matrix exponential.

2. Apply the product rule to $u'=g_Ruh^\dagger$ and insert $r_\mu'=g_Rr_\mu g_R^\dagger+ig_R\partial_\mu g_R^\dagger$. The right one-form becomes
$$
u'^\dagger(\partial_\mu-ir'_\mu)u'
=h\bigl[u^\dagger(\partial_\mu-ir_\mu)u\bigr]h^\dagger+h\partial_\mu h^\dagger.
$$
Using $u'^\dagger=g_Lu^\dagger h^\dagger$ gives the same inhomogeneous term for the left one-form. Their difference and half-sum therefore give
$$
u_\mu'=hu_\mu h^\dagger,\qquad
\Gamma_\mu'=h\Gamma_\mu h^\dagger+h\partial_\mu h^\dagger.
$$

3. Substitute $\chi'=g_R\chi g_L^\dagger$ in each term of $\chi_\pm$. The two identities
$$
u'^\dagger\chi'u'^\dagger=h(u^\dagger\chi u^\dagger)h^\dagger,\qquad
u'\chi'^\dagger u'=h(u\chi^\dagger u)h^\dagger
$$
give $\chi_\pm'=h\chi_\pm h^\dagger$. Cyclicity and $h^\dagger h=1$ imply
$$
\langle u'_\mu u'^{\mu}\rangle=\langle u_\mu u^\mu\rangle,\qquad
\langle\chi'_+\rangle=\langle\chi_+\rangle.
$$
The covariant expressions are exact locally, without a meson-number or $1/f_a$ expansion.

For the symmetric axial rotation used in the [finite-mixing record](finite-alp-meson-mixing.md), expand both sides of $U'=g_RUg_L^\dagger$ around the vacuum. The terms linear in $a$ and linear in the meson coordinates agree with the exponential obtained from the linear shift $P'_i=P_i-w_i^{(a)}a$: there is no $aP_i$ term in the coordinate map. The first possible nonlinear term linear in $a$ is of order $aP_iP_j$. Consequently that nonlinear term first contributes to cubic and higher interactions when substituted into a stationary quadratic Lagrangian at order $aP^3$, not $aP^2$. This observation is specific to the symmetric axial rotation; a general chiral transformation with a vector part can have an $aP_i$ term. It does not by itself establish the complete physical amplitude, which also requires the other external states and exchange diagrams where present.

The absence of that axial-map contribution to $aPP$ does **not** imply that the physical vertex vanishes. As established in the [flavour-violating amplitude record](../two-pseudoscalar/alp-flavour-violating-two-pseudoscalar-amplitudes.md), a flavour-changing $c_q^{VS}$ gives, for distinct LO meson mass eigenstates, the generally nonzero result
$$
i\mathcal M(a\to P_iP_j)=-\frac{m_i^2-m_j^2}{2f_a}
\left\langle c_q^{VS}[\lambda_i^\dagger,\lambda_j^\dagger]\right\rangle.
$$
This is a first-order $1/f_a$ benchmark with the qualifications in that record, not a finite-mixing result.

For a separate tree-level $U(2)$ example, use two meson coordinates and a vector rotation with generator $\sigma_2/2$. Its coordinate map is $P'_1=P_1+xP_2$, $P'_2=P_2-xP_1$ to first order in $x=a/f_a$. Substitute its inverse into a diagonal quadratic Lagrangian with squared masses $m_1^2,m_2^2$. Momentum conservation converts the resulting contact contribution to the three-point vertex into
$$
-(p_1^2-m_1^2)+(p_2^2-m_2^2).
$$
It is generally nonzero off shell and vanishes when both meson legs are on shell. An independent vector-source interaction remains nonzero on shell when its coefficient and $m_1^2-m_2^2$ are nonzero. This example uses a local two-mode quadratic Lagrangian; it is not a claim that arbitrary unequal masses for the displayed $U(2)$ generators follow from the LO QCD mass operator. At finite mixing, all external states must be projected onto the physical pole vectors before applying the same argument.

To include a mixed incoming pole in this formal example, let $D(s)=sK-M^2$ be the inverse propagator in the two-dimensional ALP and neutral-meson subspace, as in the [finite-mixing record](finite-alp-meson-mixing.md). Take a pole vector $(1,w_n^{(a)})$ and write the old neutral coordinate as $P_n=P'_n+w_n^{(a)} a'$. Its transformed off-diagonal inverse propagator is $D_{na}(s)+w_n^{(a)}D_{nn}(s)$; it vanishes at the selected pole when $(1,w_n^{(a)})$ is a null vector. The linear part of the off-shell amputated three-point vertex is then $\mathcal V_{a12}^{(3)}+w_n^{(a)}\mathcal V_{n12}^{(3)}$. Composing this shift with the same formal vector rotation gives the additional nonlinear-map term proportional to
$$
-(p_1^2-m_1^2)+(p_2^2-m_2^2).
$$
It is generically nonzero off shell, while the coherent pole-vector contribution can remain nonzero on shell. This extends the formal algebraic check to a mixed incoming pole; it does not establish a physical charge-conserving $U(3)$ comparison or an $aPPP$ calculation.

## Verification

[`derive_alp_ccwz_redefinitions.py`](../../../scripts/derive_alp_ccwz_redefinitions.py) applies the product rule using abstract noncommuting matrices and derives both mass-block identities. It also checks the two traces with symbolic $2\times2$ matrices and an arbitrary invertible compensator, and expands the symmetric axial map through one ALP and two meson coordinates. The abstract product-rule result is not restricted to two flavors; the concrete trace check illustrates cyclicity for $U(2)$. The test exercises the calculation and its assertions.

[`derive_alp_two_meson_vector_redefinition.py`](../../../scripts/derive_alp_two_meson_vector_redefinition.py) derives the $U(2)$ vector map from explicit Pauli matrices, substitutes it in the quadratic Lagrangian, and compares the off-shell contact term with a nonzero on-shell vector-source vertex. The existing channel-specific first-order calculation is linked in the flavour-violating record.

```sh
.venv/bin/python scripts/derive_alp_ccwz_redefinitions.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_ccwz_redefinitions.py
.venv/bin/python scripts/derive_alp_two_meson_vector_redefinition.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_two_meson_vector_redefinition.py
```
