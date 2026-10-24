---
title: Numerical input contract for the pion TFF spectral representation
status: derived
created: 2026-10-24
updated: 2026-10-24
bibtex_keys: [Hoferichter:2018kwz]
provenance: analytical-derivation
---

# Scope

This record defines a reproducible interface for the low-energy double-
spectral contribution to the physical pion transition form factor. It does
not provide physical spectral amplitudes or an odd-coefficient fit. The
dispersive representation and its source amplitudes follow
Hoferichter:2018kwz, Eqs. (3.11)--(3.14) and the discussion of the
$\gamma^*\to3\pi$ normalization [@Hoferichter:2018kwz].

The executable [pion TFF spectral input script](../../../scripts/pion_tff_spectral_input.py)
accepts explicitly supplied complex arrays for the pion vector form factor
$F_\pi^V(x)$ and the $P$-wave $f_1(x,y)$, with $x$ the isovector and $y$ the
isoscalar invariant mass squared. It prints the missing-input contract when
no metadata file is supplied. With data, it evaluates the central-grid
moments using composite Simpson quadrature. It does not append an explicit
tree $\rho$ contribution, effective pole, or asymptotic term.

# Convention and spectral moments

The source paper uses a positive real-photon pion TFF, whereas the repository
uses $F_{\pi\gamma\gamma}(0,0)<0$. After converting the complete scalar by
the external-pion phase, the double-spectral density is
$$
\rho^{\mathrm{disp}}(x,y)=-\frac{q_\pi^3(x)}{12\pi\sqrt{x}}
\operatorname{Im}\!\left([F_\pi^V(x)]^*f_1(x,y)\right),\qquad
q_\pi(x)=\sqrt{\frac{x}{4}-M_\pi^2}.
$$
Here $F_\pi^V$ is dimensionless, $f_1$ has units $\mathrm{GeV}^{-3}$,
and $\rho^{\mathrm{disp}}$ has units $\mathrm{GeV}^{-1}$. The source's
$F_\pi=92.28\,\mathrm{MeV}$ normalization is recorded as input metadata;
the code does not silently replace it with the repository's chiral-limit
$f_\pi=93\,\mathrm{MeV}$ [@Hoferichter:2018kwz].

For a finite spectral rectangle with the Bose-exchanged term included, the
script calculates the correlated moments
$$
\begin{aligned}
F_{\mathrm{disp}}(0,0)&=\frac{1}{\pi^2}\int dx\,dy\,
\rho^{\mathrm{disp}}(x,y)\frac{2}{xy},\\
\left.\frac{\partial F_{\mathrm{disp}}}{\partial s_1}\right|_{0,0}
&=\frac{1}{\pi^2}\int dx\,dy\,\rho^{\mathrm{disp}}(x,y)
\left(\frac{1}{x^2y}+\frac{1}{xy^2}\right),\\
\left.\frac{\partial^2 F_{\mathrm{disp}}}{\partial s_1\partial s_2}\right|_{0,0}
&=\frac{1}{\pi^2}\int dx\,dy\,\rho^{\mathrm{disp}}(x,y)\frac{2}{x^2y^2},\\
C_{\mathrm{disp}}^{\mathrm{single}}
&=\lim_{Q^2\to\infty}Q^2F_{\mathrm{disp}}(-Q^2,0)
=\frac{1}{\pi^2}\int dx\,dy\,\rho^{\mathrm{disp}}(x,y)
\left(\frac1x+\frac1y\right).
\end{aligned}
$$
The units are respectively $\mathrm{GeV}^{-1}$, $\mathrm{GeV}^{-3}$,
$\mathrm{GeV}^{-5}$, and $\mathrm{GeV}$. These are finite-rectangle
moments, not the full pion TFF unless the declared cutoffs and complementary
effective and asymptotic terms are also supplied.

# Required numerical input and limitations

The metadata must declare source keys, the source convention, decay-constant
normalization, charged-pion mass, isoscalar threshold, source provenance,
correlated uncertainty prescription, both amplitude models, and the
cutoff/matching prescription. The NPZ input must contain increasing $x,y$
grids, complex $F_\pi^V(x)$, and complex $f_1(x,y)$ with the shapes and units
specified in the script. The isoscalar threshold is explicit because the
included channel set determines it; the hadronic $3\pi$ threshold and the
$\pi^0\gamma$ threshold are not interchangeable.

The repository currently has no selected numerical $P$-wave phase/input
replicas for $F_\pi^V$, nor a numerical complex $f_1(x,y)$ grid with
correlated fit variations. A total $e^+e^-\to3\pi$ cross section constrains
the amplitude modulus but does not by itself supply this complex partial
wave. Consequently, the script reports no physical density, moments,
uncertainty, or coefficient constraint until these inputs are sourced and
matched. A standalone tree $\rho$ term must not be added to a spectral input
that already contains the same resonance strength [@Hoferichter:2018kwz].

The next calculation is to provide phase-dependent correlated grids or
replicas for both amplitudes, validate their threshold and normalization
conventions, then evaluate grid and cutoff convergence. Only after that may
the moments be matched to the anomaly/effective-pole/asymptotic pieces and
used to constrain the original odd-parity coefficient combinations.
