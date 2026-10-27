---
title: Running three-pion width from a supplied unit-normalized KT amplitude
status: derived
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Hoferichter:2019mqg, Hoferichter:2025lcz]
provenance: analytical-derivation
---

# Inputs and assumptions

The calculator converts a supplied complex full scalar
$\mathcal F(s,t,u;q^2)$ into a three-pion width shape for the spectral
templates. It assumes exact isospin, a fixed vector-current coupling, and a
unit-subtraction KT solution whose full scalar excludes the fitted production
normalization $a(q^2)$. The caller supplies the amplitude and quadrature grid;
this record does not assert that a physical KT solution is available.

# Cross-section integral and width ratio

The $e^+e^-\to3\pi$ cross-section formula is
$$
\sigma_{3\pi}(q^2)=\alpha^2\int ds\,dt\,
\frac{s\kappa(s,q^2)^2(1-z_s^2)}{768\pi(q^2)^3}
|\mathcal F(s,t,u;q^2)|^2,
$$
with $s$ the pair invariant, $z_s=(t-u)/\kappa$, and
$\kappa=\sqrt{1-4m_\pi^2/s}\,\lambda^{1/2}(q^2,m_\pi^2,s)$
[@Hoferichter:2019mqg, Eq. (5); @Hoferichter:2025lcz, Eq. (9)]. Since
$dt=\kappa\,dz_s/2$, the executable integrates this expression on
caller-supplied $(s,z_s)$ nodes.

Contracting the transverse current spectral tensor with the lepton tensor
and with a spin-averaged vector source gives, for fixed vector-current
coupling, $\Gamma_{3\pi}(q^2)/\sigma_{3\pi}(q^2)\propto(q^2)^{3/2}$.
Normalizing to an explicit on-pole three-pion partial width therefore gives
$$
\Gamma_{3\pi}(q^2)=\Gamma_{3\pi}(M_V^2)
\left(\frac{q^2}{M_V^2}\right)^{3/2}
\frac{\sigma_{\mathrm{unit}}(q^2)}
{\sigma_{\mathrm{unit}}(M_V^2)}.
$$
The common $\alpha^2$ factor cancels in the ratio. Here
$\sigma_{\mathrm{unit}}$ denotes the phase-space integral divided by
$\alpha^2$ for the supplied unit-subtraction full KT scalar; it is not a
measured cross section.

# Reproducibility and limits

Run:

    .venv/bin/python scripts/build_vppp_three_pion_width.py --derive-formulas

The default mode requires a complex full-amplitude grid, its Dalitz
quadrature, resonance pole masses, and on-pole three-pion partial widths. It
returns a normalized width grid for combination with the two-body widths in
the [partial-width helper](vppp-two-body-partial-width-inputs.md). It does
not solve the Khuri--Treiman equations, choose a pion phase, or invent
on-pole widths. The published width prescription includes rescattering
[@Hoferichter:2019mqg; @Niecknig:2012sj], so replacing the missing KT input
with a tree-level amplitude would define a different conditional model.
