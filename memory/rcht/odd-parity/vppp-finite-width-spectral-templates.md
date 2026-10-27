---
title: Finite-width resonance and conformal spectral templates for the VPPP profile
status: derived
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Hoferichter:2019mqg, Hoferichter:2025lcz, Niecknig:2012sj]
provenance: analytical-derivation
---

# Source convention and inputs

The dispersive normalization uses the resonance function and conformal
polynomial in [@Hoferichter:2019mqg, Eqs. (15)--(17)], repeated in
[@Hoferichter:2025lcz]. This calculation produces their absorptive templates
for an explicitly supplied grid $s$ and energy-dependent width functions
$\Gamma_V(s)$. The inputs must also state the quadrature weights, high-energy
tail prescription, resonance masses, conformal powers, and width provenance.
No physical width grid or tail is bundled.

The [two-body partial-width helper](vppp-two-body-partial-width-inputs.md)
provides conditional channel shapes from explicitly supplied on-pole partial
widths. The [three-pion width helper](vppp-three-pion-running-width.md)
converts an externally supplied full KT amplitude grid into a normalized
running $3\pi$ width. Neither helper supplies the required physical
rescattered amplitude grid or pole partial-width inputs.

# Executed spectral formulas

For real residue $c_V$, the source denominator is
$$
\mathcal A_V(s)=\frac{c_V}{M_V^2-s-i\sqrt{s}\,\Gamma_V(s)}.
$$
SymPy derives the per-unit-residue absorptive template
$$
\operatorname{Im}\frac{1}{M_V^2-s-i\sqrt{s}\,\Gamma_V(s)}
=\frac{\sqrt{s}\,\Gamma_V(s)}
{(M_V^2-s)^2+s\Gamma_V(s)^2}.
$$
It has units $\mathrm{GeV}^{-2}$; multiplying by $c_V$ in
$\mathrm{GeV}^{-1}$ gives $\operatorname{Im}\mathcal A_V$ in
$\mathrm{GeV}^{-3}$. The existing profile calculator inserts these templates
into the once-subtracted dispersion kernel.

For the conformal map, define
$$
A=\sqrt{s_{\mathrm{inel}}-s_1},\qquad
B=\sqrt{s-s_{\mathrm{inel}}}.
$$
On the upper edge of the cut, $s>s_{\mathrm{inel}}$,
$\sqrt{s_{\mathrm{inel}}-(s+i0)}=-iB$, hence
$$
z(s+i0)=\frac{A+iB}{A-iB}=e^{2i\theta},\qquad
\theta=\arctan\frac{B}{A}.
$$
Since $z(0)$ is real, the spectral template for conformal power $n\geq1$ is
$$
\operatorname{Im}\big[z(s+i0)^n-z(0)^n\big]
=\sin(2n\theta),\qquad s>s_{\mathrm{inel}},
$$
and it vanishes below threshold. These templates are dimensionless; the
conformal coefficients carry units $\mathrm{GeV}^{-3}$.

# Reproducibility and limitations

The executable [template builder](../../../scripts/build_vppp_spectral_templates.py)
prints the SymPy derivation with:

    .venv/bin/python scripts/build_vppp_spectral_templates.py --derive-formulas

For actual templates it requires an NPZ with explicit widths, resonance
masses, spectral grid, integration weights, conformal powers, and tail/source
metadata. It refuses to extrapolate a tail or invent a width. Its real-residue
templates do not include $\rho$--$\omega$ mixing or complex residues.

[@Hoferichter:2019mqg] states that the $\omega$ and $\phi$ widths include the
main decay modes; the three-pion component must include final-state
rescattering [@Niecknig:2012sj]. The fitted paper does not publish numerical
width arrays. Therefore the builder presently supplies the functional form,
not physical spectral values. The required width grids, in particular the
rescattered $3\pi$ widths, remain to be calculated before the finite-width
profile can be evaluated.
