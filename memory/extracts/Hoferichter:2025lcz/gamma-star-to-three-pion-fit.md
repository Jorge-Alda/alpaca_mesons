---
title: Dispersive gamma-star to three-pion fit for the pion-TFF input
status: literature
created: 2026-10-25
updated: 2026-10-27
bibtex_keys: [Hoferichter:2025lcz]
provenance: literature
---

# Amplitude and partial wave

Hoferichter:2025lcz writes the current matrix element as
$$
\langle0|j_\mu(0)|\pi^+(p_+)\pi^-(p_-)\pi^0(p_0)\rangle
=-\epsilon_{\mu\nu\alpha\beta}p_+^\nu p_-^\alpha p_0^\beta
\mathcal F(s,t,u;q^2),
$$
with $q=p_++p_-+p_0$ and $s+t+u=3M_\pi^2+q^2$. Below the region where
higher partial waves become important, the reconstruction is
$$
\mathcal F(s,t,u;q^2)=\mathcal F(s,q^2)+\mathcal F(t,q^2)
+\mathcal F(u,q^2),
$$
and the Khuri--Treiman solution gives
$$
f_1(s,q^2)=\mathcal F(s,q^2)+\widehat{\mathcal F}(s,q^2).
$$
The function $\mathcal F(s,q^2)$ obeys the dispersive integral equation with
the elastic $\pi\pi$ P-wave phase, while the crossed-channel angular average
$\widehat{\mathcal F}$ supplies the left-hand-cut contribution
[@Hoferichter:2025lcz, Eqs. (2.1)--(2.8)].

# Fitted normalization and data choices

The partial-wave solution is multiplied by a fitted normalization $a(q^2)$.
Its dispersive parameterization contains explicit vector-resonance terms and
a conformal polynomial for inelastic effects, with a sum rule tying the
subtraction constant to the resonance and inelastic spectral terms
[@Hoferichter:2025lcz, Eq. (2.9) and following discussion]. The global fit to
SND, CMD-2, and BaBar data gives $F_{3\pi}=33.1(1.7)\,\mathrm{GeV}^{-3}$;
the paper also reports separate variants in conformal-polynomial order and
the phase of $\rho$--$\omega$ mixing [@Hoferichter:2025lcz, Sec. 3 and
Tables 1--3]. Belle II fits show poorer fit quality and tensions in the
resonance widths and the inferred anomaly, so those fits should not be mixed
into the global-fit central input without a separate model for the tension
[@Hoferichter:2025lcz, Sec. 4].

The global fit uses the cross-section phase-space integral of the full
complex three-pion amplitude [@Hoferichter:2025lcz, Eq. (9)]. Its
normalization is written $a(q^2)=\alpha_A+$ a dispersive vector-resonance
term $+C_p(q^2)$, with $\alpha_A=C_A/3$ and a sum rule relating $\alpha_A$
to the resonance and inelastic spectral functions [@Hoferichter:2025lcz,
Eqs. (11), (12), and (14)]. The paper reports
$C_A/C_A^{\mathrm{WZW}}=1.028(15)(48)(17)$ for the global fit, where the
quoted components are statistical, conformal-truncation, and
$\rho$--$\omega$-phase uncertainties; it also tabulates the vector residues
such as $c_\omega$ and conformal coefficients for
$p_{\mathrm{conf}}=2,3,4$ [@Hoferichter:2025lcz, Tables 1--3]. These tables
give marginal uncertainties, not the joint covariance of the fit parameters.
Reproducing their likelihood requires the underlying experiment data and
covariance treatment.

The reconstruction writes the full scalar amplitude as a sum of three
single-variable functions, and its P-wave solution is the single-variable
function plus the crossed-channel average [@Hoferichter:2025lcz,
Eqs. (6)--(8); Hoferichter:2023bjm, Eqs. (4), (16), and (17)]. Before
rescattering, a Dalitz-independent full scalar shift is divided equally
among the three single-variable terms. After Khuri--Treiman rescattering,
the corresponding subtraction shift multiplies the normalized dispersive
solution and its crossed-channel average. Therefore the bare RChT scalar
cannot simply be added to the physical P-wave amplitude; its matching to
$a(q^2)$ requires a declared subtraction and contact prescription. The
cited papers do not provide an RChT $g_i$-to-$a(q^2)$ map.

The detailed amplitude and fitted-normalization framework is recorded in
[Hoferichter:2023bjm](../Hoferichter:2023bjm/three-pion-amplitude-and-normalization.md)
[@Hoferichter:2023bjm].

# Limit for the pion-TFF spectral grid

This is a useful updated source for the normalization $a(q^2)$ and its
systematic variants. The paper does not publish a complex grid of
$f_1(s,q^2)$ or complete correlated fit replicas/covariance for all
parameters. Its tables of central fit parameters and marginal errors alone
are therefore insufficient to manufacture correlated spectral replicas.
Constructing the required $f_1(x,y)$ still needs the numerical Khuri--Treiman
solution, the chosen $\pi\pi$ P-wave phase and pion vector form factor, and
correlated variations of the fitted normalization and model inputs.

# Provenance

The complete arXiv v2 source is archived in
[the source index](../../sources/Hoferichter:2025lcz/INDEX.md). The original
pion-TFF dispersive construction is summarized in
[Hoferichter:2018kwz](../Hoferichter:2018kwz/pion-tff-dispersive-asymptotic-representation.md)
[@Hoferichter:2018kwz].
