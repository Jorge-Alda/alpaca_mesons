---
title: Conditional two-body partial-width shapes for VPPP spectral templates
status: derived
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Ecker:1989yg, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk, Hoferichter:2019mqg]
provenance: analytical-derivation
---

# Assumptions and inputs

The width helper uses the two-body Lorentz structures for $V\to P\gamma$ in
the [radiative-width record](radiative-width-kinematics.md) and
$V\to P_1P_2$ from the even-parity RChT vector vertex
[@Ecker:1989yg; @Ruiz-Femenia:2003jdx]. It assumes a constant reduced vertex
for each channel. Pole masses, daughter masses, and on-pole partial widths
are explicit inputs; PDG values are available in
[@ParticleDataGroup:2024cfk], but this helper does not silently combine their
conventions with VP-subtracted fit masses.

# Derived running-width shapes

For invariant mass squared $s$ and two-body momentum
$$
k(s)=\frac{\sqrt{\lambda(s,m_1^2,m_2^2)}}{2\sqrt{s}},\qquad
\lambda(s,m_1^2,m_2^2)=[s-(m_1+m_2)^2][s-(m_1-m_2)^2],
$$
the width is exactly zero at and below threshold. With a constant reduced
$V P\gamma$ vertex, the polarization-summed amplitude and two-body phase
space give $\Gamma_{V\to P\gamma}(s)\propto k(s)^3$. Normalized to an
explicit on-pole partial width,
$$
\Gamma_{V\to P\gamma}(s)=\Gamma_{V\to P\gamma}(M_V^2)
\left[\frac{k(s)}{k(M_V^2)}\right]^3.
$$
For the vector current vertex $g\,\varepsilon_V\cdot(p_1-p_2)$,
$\Gamma_{V\to P_1P_2}(s)\propto k(s)^3/s$, hence
$$
\Gamma_{V\to P_1P_2}(s)=\Gamma_{V\to P_1P_2}(M_V^2)
\frac{M_V^2}{s}
\left[\frac{k(s)}{k(M_V^2)}\right]^3.
$$
The SymPy derivation and threshold implementation are in the
[width builder](../../../scripts/build_vppp_partial_widths.py).

# Role in the dispersive spectra

The builder accepts distinct channels and thresholds, including separate
$K^+K^-$ and $K^0\bar K^0$ inputs for $\phi$. It returns the sum of only the
channels supplied by the caller and requires their sources and omitted
channels to be declared. It also requires an external, rescattered
$V\to3\pi$ width grid for every resonance; its output must not be described
as a total width unless that grid and all other relevant channels are
included.

This is a phase-space model with constant reduced vertices. It adds no
barrier factors, off-shell vertex corrections, or photon radiative
corrections. These choices are explicit approximations, not PDG line-shape
prescriptions. The $3\pi$ width remains the main missing contribution in
[@Hoferichter:2019mqg], which requires the normalized Khuri--Treiman solution
[@Niecknig:2012sj].

# Reproducibility

Run:

    .venv/bin/python scripts/build_vppp_partial_widths.py --derive-formulas

The default executable mode requires a caller NPZ containing resonance and
daughter masses, on-pole partial widths, and the external three-pion width
grid. It writes an NPZ with the summed supplied partial widths and a
provenance sidecar; no physical width values are bundled.
