---
title: KK pi vector-current coefficient amplitudes and electromagnetic integration
status: derived
created: 2026-10-30
updated: 2026-10-30
bibtex_keys: [GomezDumm:2010ra, BaBar:2007ceh, Alda:2026nrd]
provenance: analytical-derivation
---

# Scope and inputs

The executable [module](../../../scripts/implementation/kkpi_vector_current.py)
provides the resonance part of the charged-current $F_3$ for
$K^+K^-\pi^-$ and $K^-K^0\pi^0$, using Eqs. (30), (35) and Appendix A
of [@GomezDumm:2010ra]. It also integrates an explicitly supplied
electromagnetic scalar amplitude into a Born cross section and energy-bin
averages. This is a reusable calculation for a future inclusive or Dalitz
likelihood. It does not predict the experimental component in Table 7 of
[@BaBar:2007ceh]. A separate direct crossed $K^*K$ current calculation now
fits that quasi-two-body component in
[babar-kkpi-rcht-coefficient-fit.md](babar-kkpi-rcht-coefficient-fit.md);
it does not use the inclusive $F_3$ integration implemented here. The full
inclusive electromagnetic coefficient fit is documented separately in
[babar-kkpi-unified-vector-coefficient-fit.md](babar-kkpi-unified-vector-coefficient-fit.md).

Inputs are the [coefficient polynomials and their invariant ordering](kkpi-vector-current-coefficient-rank.md),
the [canonical tensor and coefficient conversion](levi-civita-convention.md),
and the [measured-component definition](kkpi-isovector-spectrum-input-and-extraction-readiness.md).
All masses and decay constants are in GeV; invariants are in GeV$^2$.
The sixteen real coefficients are ordered as
$(c_1,\ldots,c_7,d_1,\ldots,d_4,g_1,\ldots,g_5)$ and are dimensionless.
No short-distance relation is imposed.

# Form factors and assumptions

Assume tree-level production, exact strong isospin and ideal neutral-vector
mixing, as in the source formula. The caller supplies the three complex
propagators, including all masses, self energies, subtractions and upper-rim
prescriptions. There is no default width model. The charged and neutral
resonance scalars have the common factors
$-4G_V/(M_VF^3)$ and $2\sqrt{2}G_V/(M_VF^3)$, respectively, multiplying
the expressions in the coefficient audit. The canonical linear scalar
gradients retain these reference signs because reversing the odd
coefficients and the epsilon tensor cancels. The WZW constant requires its
own tensor conversion and must be supplied explicitly.

`resonance_scalar_basis` returns the sixteen gradients of this weak-current
scalar. It does not identify this scalar with a neutral electromagnetic
current. The required electromagnetic charge/isospin map, the anomaly term,
excited-vector contributions and charge-conjugate or $K_S$ factors remain
explicit external inputs. In particular, no physical coefficient value is
selected by this module.

The optional `kstar` selection retains the $K^*(892)$ denominators in the
source formula and removes the local VPPP and charged-mode omega terms.
All five $g_i$ gradients then vanish exactly. This selection of terms is
an algebraic diagnostic; it is not a reconstruction of BaBar's Dalitz fit.
Its energy-dependent fitted component cannot inherit the rank of the full
$F_3$.

# Electromagnetic normalization and integration

For a supplied electromagnetic current in the convention
$$
H^\mu=i\epsilon^{\mu\nu\rho\sigma}p_{1\nu}p_{2\rho}p_{3\sigma}
\mathcal F(Q^2,s_{12},s_{13}),
$$
the scalar has units GeV$^{-3}$. Neglecting the electron mass, the spin
average and orientation integral give
$$
\sigma(Q^2)=\frac{e^4}{6(Q^2)^2}\int d\Phi_3\,|\boldsymbol H|^2
=\frac{\alpha^2}{48\pi(Q^2)^3}
\int ds_{12}\,ds_{13}\,
|\epsilon p_1p_2p_3|^2|\mathcal F|^2.
$$
Here $|\epsilon p_1p_2p_3|^2$ denotes the positive spatial norm in the
center-of-mass frame. With $z$ the pair-rest-frame angle,
$$
\begin{aligned}
|\epsilon p_1p_2p_3|^2&=
\frac{\lambda(Q^2,s_{12},m_3^2)\lambda(s_{12},m_1^2,m_2^2)}{16s_{12}}
(1-z^2),\\
\left|\frac{ds_{13}}{dz}\right|&=
\frac{\sqrt{\lambda(Q^2,s_{12},m_3^2)\lambda(s_{12},m_1^2,m_2^2)}}{2s_{12}}.
\end{aligned}
$$
Gauss-Legendre quadrature integrates these expressions. A supplied set of
complex scalar coefficients is integrated coherently into a real symmetric
matrix. Contracting it with the real coefficient vector on both sides
returns the cross section in nb, including interference. A fixed amplitude
can be represented by a column whose coefficient is fixed to one.
The conversion is $1\,\mathrm{GeV}^{-2}=389379.3656\,\mathrm{nb}$.

`energy_bin_matrices` averages uniformly in $E=\sqrt{Q^2}$ over the supplied
published bin edges. Uniform averaging is an explicit approximation: the
ISR luminosity and acceptance response within the BaBar bins have not been
provided. A bin below threshold needs an explicit scalar dimension; the
single-energy integral returns `None` below threshold. Quadrature refinement
for physical narrow propagators remains necessary before a fit.

# Verification and remaining inputs

Invocation:
```text
.venv/bin/python -m unittest scripts.test_kkpi_vector_current
```
Four focused tests passed. They reproduce the source charged-mode
polynomials and relative signs, verify that every VPPP gradient vanishes
under the $K^*$ selection, check the coherent interference and energy-bin
average, and check the absolute normalization against the analytic massless
constant-scalar integral
$$
\int ds_{12}ds_{13}\,|\epsilon p_1p_2p_3|^2=(Q^2)^5/480,
\qquad \sigma=\alpha^2(Q^2)^2|\mathcal F|^2/(23040\pi).
$$
These mathematical fixtures are not physical predictions for the data.

The `babar_table7_prediction` entry explicitly refuses an inclusive $F_3$
integral. The separate full-current fit uses the source-specific
electromagnetic form factors and the inclusive BaBar Tables 2/3 data,
including the local VPPP terms, and reports the coefficient map and rate
rank in
[babar-kkpi-unified-vector-coefficient-fit.md](babar-kkpi-unified-vector-coefficient-fit.md).
The weak-current module here still has no combined tau likelihood.
Published bin-covariance limitations remain relevant.
