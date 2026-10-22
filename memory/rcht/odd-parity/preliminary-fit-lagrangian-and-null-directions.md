---
title: Preliminary vector-sector fit Lagrangian and odd-coefficient null directions
status: mixed
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: [Ecker:1988te, Ruiz-Femenia:2003jdx, Alda:2026nrd, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# Preliminary fit Lagrangian and coefficient directions

## Scope and conventions

This record spells out the Lagrangian underlying the selected
`central_single_trace_negative` preliminary point, including the conditional
leading $G_V$ value. The vector fields are antisymmetric tensors in the
canonical convention of [the RChT basics](../foundations/rcht-basics.md)
[@Ecker:1988te]. The odd tensors use metric $(+---)$ and
$\epsilon^{0123}=+1$; their coefficients are the convention-preserving canonical
ones defined in [the convention record](levi-civita-convention.md). The
original odd basis is that of [@Ruiz-Femenia:2003jdx]. The fit point and its
conditional status are documented in [the focused fit record](odd-preliminary-fit.md).

## Even-intrinsic-parity sector

Let $\langle\cdots\rangle$ denote a flavor trace and
$\chi_+=\operatorname{diag}(x_l,x_l,x_s)$ in the strong-isospin vacuum.
With $\nabla_\mu$ the chiral covariant derivative, the retained even sector is

$$
\begin{aligned}
\mathcal L_{\mathrm{even}}={}&
\frac{F^2}{4}\langle u_\mu u^\mu+\chi_+\rangle
-\frac12\left\langle\nabla^\lambda V_{\lambda\mu}\nabla_\nu V^{\nu\mu}
-\frac{M_V^2}{2}V_{\mu\nu}V^{\mu\nu}\right\rangle\\
&-e_m^V\langle V_{\mu\nu}V^{\mu\nu}\chi_+\rangle
+\frac{\Delta_V}{12}\langle V_{\mu\nu}\rangle\langle V^{\mu\nu}\rangle
+\frac{\kappa_V}{4}\langle V_{\mu\nu}\rangle\langle V^{\mu\nu}\chi_+\rangle\\
&+\frac{F_V}{2\sqrt2}\langle V_{\mu\nu}f_+^{\mu\nu}\rangle
+\frac{z_V}{2\sqrt2}\langle V_{\mu\nu}\{\chi_+,f_+^{\mu\nu}\}\rangle
+\frac{iG_V}{\sqrt2}\langle V_{\mu\nu}u^\mu u^\nu\rangle.
\end{aligned}
$$

The $e_m^V$ term is the single-trace leading quark-mass correction to the
vector mass; $\Delta_V$ shifts the normalized singlet mass. The retained
$\kappa_V$ term is the additional double-trace mass insertion with one
$\chi_+$, using exactly the normalization
$(\kappa_V/4)\langle V_{\mu\nu}\rangle\langle V^{\mu\nu}\chi_+\rangle$
of [the mass-extension derivation](../foundations/vector-mass-extension.md).
The $z_V$ term is the single-trace $\chi_+$ correction to the vector-current
source coupling; it is an explicit fit extension, not part of the minimal
leading RChT Lagrangian. The fit chooses the single-trace source hypothesis,
so its additional double-trace source coefficient is fixed to $y_V=0$.
For electromagnetism $f_+^{\mu\nu}=-2eQF^{\mu\nu}$, giving the fitted source
matrix $F_VQ+z_V\{\chi_+,Q\}$.

The displayed mass parameters may equivalently use $m_\rho$ as a fitted
coordinate, with $M_V^2=m_\rho^2+4e_m^Vx_l$ after applying the charged mass
relation. The neutral mass matrix and $\omega$--$\phi$ angle are reconstructed
from the mass sector and the selected root; the angle is not an independent
Lagrangian coefficient. The coefficient $G_V$ is not determined by the odd
fit: its appended value is the separate conditional four-row fit described in
[the $G_V$ record](../foundations/vector-gv-fit.md). Widths fix $|G_V|$; the
positive representative is a convention.

## Exported central parameter representative

The production file
[`odd_preliminary_fit.json`](../../../data/vector_sector/odd_preliminary_fit.json)
contains the selected parameter representative below. The even and odd
entries are the stored preliminary-fit point; $G_V$ is appended from its
separate conditional fit. Masses and widths are in GeV, $\Delta_V$ is in
GeV$^2$, and all $c_i,d_i$ are dimensionless.

| Parameter | Value | Parameter | Value |
|---|---:|---|---:|
| $m_\rho$ | 0.7750812133 | $m_{K^*}$ | 0.8916532038 |
| $\Delta_V$ | 0.01819699987 | $\kappa_V$ | 0.04284647010 |
| $e_m^V$ (reconstructed) | -0.1071003214 | $M_V$ (reconstructed) | 0.7648951208 |
| $F_V$ | 0.1689527776 | $z_V$ | -0.1702699573 |
| $G_V$ (conditional) | 0.05925137853 | $\theta_V$ | 0.08241351807 rad |

| Odd coefficient | Value | Odd coefficient | Value |
|---|---:|---|---:|
| $c_1$ | -0.0035140102158 | $c_2$ | -0.0240705857790 |
| $c_3$ | 0.00563895764777 | $c_4$ | 0.0195211429354 |
| $c_5$ | 0.00351401021582 | $c_6$ | 0.0205565755632 |
| $c_7$ | 0 (representative) | $d_1$ | -0.129639753068 |
| $d_2$ | 0.00335172479195 | $d_3$ | 0.0653589455103 |
| $d_4$ | 0 (representative) | | |

The shared fitted total-width nuisances are
$\Gamma_{K^{*0}}=0.04729821441$,
$\Gamma_{K^{*+}}=0.05139889416$,
$\Gamma_{\eta'}=0.0002042322937$,
$\Gamma_\omega=0.008658598423$,
$\Gamma_\phi=0.004249086249$,
$\Gamma_{\rho^0}=0.1475047700$, and
$\Gamma_{\rho^+}=0.1490875366$ GeV. The zeros shown for $c_7,d_4$ are
representative coordinates along the null space; they do not determine the
corresponding coefficient values. Do not interpret the exported coefficients
as a covariance-bearing determination of all eleven couplings.

## Original eleven odd operators

The odd vector Lagrangian retained in the fit is

$$
\mathcal L_{\mathrm{odd}}=
\sum_{a=1}^{7}\frac{c_a}{M_V}\mathcal O^a_{VJP}
+\sum_{a=1}^{4}d_a\mathcal O^a_{VVP}.
$$

The seven VJP operators are

$$
\begin{aligned}
\mathcal O^1_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\nu},f_+^{\rho\alpha}\}\nabla_\alpha u^\sigma\rangle,\\
\mathcal O^2_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\alpha},f_+^{\rho\sigma}\}\nabla_\alpha u^\nu\rangle,\\
\mathcal O^3_{VJP}&=i\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\nu},f_+^{\rho\sigma}\}\chi_-\rangle,\\
\mathcal O^4_{VJP}&=i\epsilon_{\mu\nu\rho\sigma}\langle V^{\mu\nu}[f_-^{\rho\sigma},\chi_+]\rangle,\\
\mathcal O^5_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla_\alpha V^{\mu\nu},f_+^{\rho\alpha}\}u^\sigma\rangle,\\
\mathcal O^6_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla_\alpha V^{\mu\alpha},f_+^{\rho\sigma}\}u^\nu\rangle,\\
\mathcal O^7_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla^\sigma V^{\mu\nu},f_+^{\rho\alpha}\}u_\alpha\rangle.
\end{aligned}
$$

The four VVP operators are

$$
\begin{aligned}
\mathcal O^1_{VVP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\nu},V^{\rho\alpha}\}\nabla_\alpha u^\sigma\rangle,\\
\mathcal O^2_{VVP}&=i\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\nu},V^{\rho\sigma}\}\chi_-\rangle,\\
\mathcal O^3_{VVP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla_\alpha V^{\mu\nu},V^{\rho\alpha}\}u^\sigma\rangle,\\
\mathcal O^4_{VVP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla^\sigma V^{\mu\nu},V^{\rho\alpha}\}u_\alpha\rangle.
\end{aligned}
$$

These are the original seven VJP and four VVP operators. The VVP basis is
complete for vertices with one pseudoscalar, as stated in
[@Ruiz-Femenia:2003jdx]; this statement does not establish completeness for
three-pseudoscalar vertices.

## Measured image and null directions

At the selected fixed even/source point, the seven fitted odd amplitudes have
rank seven in the eleven-dimensional coefficient space. Hence the fit leaves
four independent coefficient directions unmeasured. The following columns
give one numerically orthonormal basis for that null space, in coefficient
order $(c_1,\ldots,c_7,d_1,\ldots,d_4)$; entries are rounded, and the basis
depends on the channel map and Euclidean coefficient normalization used by
the SVD:

$$
N\simeq\begin{pmatrix}
 0.623446&-0.459682& 0.002738&0\\
 0.134811& 0.617880&-0.007078&0\\
 0&0&0&0\\
 0&0&0&0\\
 0.758257&0.158199&-0.004340&0\\
 0.134811&0.617880&-0.007078&0\\
-0.003493&-0.010692&-0.999937&0\\
 0&0&0&0\\
 0&0&0&0\\
 0&0&0&0\\
 0&0&0&1
\end{pmatrix},\qquad \delta\mathbf g=N\boldsymbol\alpha.
$$

Here $\boldsymbol\alpha$ is any four-vector and $\mathbf g$ denotes the
ordered coefficient vector. The first three directions mix VJP coefficients;
the fourth is $d_4$. This is a local fixed-nuisance null basis for the
selected measured channels, not a universal operator identity or a statement
that these directions vanish in other amplitudes. The JSON stores a particular
zero-null representative, not measurements or uncertainties for the null
coordinates. The fitted values of the individual $c_i,d_i$ therefore depend
on this representative choice; channel predictions in the measured image do
not. Null-sensitive held-out or ALP amplitudes require a separate projection
and cannot be assigned a value by setting $\boldsymbol\alpha=0$ without an
additional assumption.

The null matrix was obtained with the existing executable projection in
`scripts/summarize_odd_preliminary_fit.py` by calling its `main(write_output=False)`
path at the selected stored fit. The stored calculation reports the original
eleven-operator measured rank as seven. The fit's diagonal quoted-summary
covariance, fixed neutral inputs, conditional source hypothesis and tree-level
truncation remain in force; the result is not a unique coefficient
determination.
