---
title: Joint leading-order meson and ALP vertex expansion
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Joint leading-order meson and ALP vertex expansion

## Assumptions and inputs

Use natural units, the $U(3)$ generators and source conventions of [@Alda:2026nrd], $U=\exp(i\Phi/f_\pi)$, and $x=a/f_a$. The [arbitrary-order CCWZ record](alp-ccwz-arbitrary-order-expansion.md) gives the source definitions and the [general-source record](alp-general-chiral-redefinition-invariants.md) defines $M$, $C_R$, and $C_L$. Use the **undressed** sources together, $\chi(x)=2B_0M(x)$ and $r_\mu=(\partial_\mu x)C_R(x)$, $\ell_\mu=(\partial_\mu x)C_L(x)$. The topological normalization and LO singlet potential follow the [finite-mixing record](finite-alp-meson-mixing.md). The bare ALP potential is taken quadratic. No expansion in $f_\pi/f_a$ is made.

Expand about the aligned vacuum $U=1$, $x=0$, with $\chi(0)=\chi(0)^\dagger$ and $\Theta(0)=0$. An ALP tadpole must be absent or cancelled by a separate linear term. General Taylor coefficients are independent inputs, with $\chi_n=d^n\chi/dx^n|_0$ and similarly for $C_{R,L}$ and $\Theta$. Charged generators may be supplied as non-Hermitian pairs; the conjugate exponential is the formal exponential with the sign of $i$ reversed in the same real-field expression.

## Method

1. For each labeled leg, insert its selected meson-generator matrix, its ALP component, and its incoming momentum into $\Phi$, $a$ and their derivatives. A derivative contributes $-ip_{i\mu}$. The executable treats products of its momentum symbols as $p_i\cdot p_j$.
2. Expand $u=\exp(i\Phi/(2f_\pi))$ and $u^\dagger$ to the number of legs, differentiating their matrix powers by the product rule. Form directly the CCWZ building blocks
   $$
   u_\mu=i(u^\dagger\partial_\mu u-u\partial_\mu u^\dagger)+u^\dagger r_\mu u-u\ell_\mu u^\dagger,
   \qquad \chi_+=u^\dagger\chi u^\dagger+u\chi^\dagger u.
   $$
   The source arrays contain Taylor **derivatives**; the expansion divides by $n!$. All products are truncated at the requested total field degree.
3. Form $f_\pi^2\operatorname{tr}(u_\mu u^\mu+\chi_+)/4$, and add the singlet and bare-ALP terms. Differentiate once with respect to each independent leg variable and set all of them to zero. This generates all permutations, including those of repeated species and derivative legs. The returned amputated expression omits the overall Feynman-rule factor $i$ and has all momenta incoming.

## Quadratic result and checks

For the stated aligned vacuum, the general kinetic blocks and meson/ALP mass-mixing blocks depend only on the zeroth current sources and the first mass-source derivative. Relative to the linear-spurion result of the [finite-mixing record](finite-alp-meson-mixing.md), the only change to the mass Hessian is
$$
\Delta M^2_{aa}=-\frac{f_\pi^2}{2f_a^2}\operatorname{Re}\operatorname{tr}\chi_2.
$$
$\Theta_2$ has no quadratic contribution at $\Theta(0)=0$; it enters higher vertices. This correction can be supplied directly to the $(K,M^2)$ interface in the [spectrum solver](alp-lo-pole-solver.md). It changes the eigenvalues and pole vectors in general, but does not change the solver itself. The same correction, written in terms of the coefficient of $x^2$ rather than the second derivative, is independently recorded in the [stability-bound record](alp-lo-bare-alp-mass-bound.md).

The unit tests compare all entries of the neutral quadratic blocks with `closed_form_blocks` for exact non-diagonal complex sources. They also compare the three types of quadratic entry with the closed forms using generic Hermitian meson generators, axial and pseudoscalar source matrices, and symbolic diagonal masses; the second-derivative correction is checked for a general complex matrix. The `quadratic_blocks` function constructs both matrices directly and supplies them to the existing numerical spectrum solver in a test: for a decoupled example, the ALP squared pole shifts by the stated correction. The general coupled spectrum requires diagonalizing the corrected matrices; individual pole shifts need not equal the $aa$ entry.

The tests also compare linear-source two- and three-meson vertices with the expressions in `derive_alp_phi2_lo_chipt.py` and `derive_alp_phi3_lo_chipt.py`. The three-meson comparison uses generic real meson matrices and general Hermitian axial and pseudoscalar source matrices, with diagonal symbolic masses; the two source sectors are independently variable in the symbolic comparison. The two-meson comparison uses generic real meson and real symmetric source matrices. A separate charged-meson check retains the independent real and imaginary parts of the Hermitian $d$–$s$ vector source. Other checks cover source-free kinetic normalization, permutations and nonlinear sources with a charged pair.

## Invocation

[`derive_alp_lo_joint_vertices.py`](../../../scripts/derive_alp_lo_joint_vertices.py) performs the SymPy expansion. Its minimal example prints the two-ALP vertex with a symbolic second mass-source derivative. [`test_derive_alp_lo_joint_vertices.py`](../../../scripts/test_derive_alp_lo_joint_vertices.py) contains the checks.

```sh
.venv/bin/python scripts/derive_alp_lo_joint_vertices.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_joint_vertices.py
```
