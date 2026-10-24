---
title: Dispersive three-pion amplitude and fitted normalization
status: literature
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Hoferichter:2023bjm]
provenance: literature
---

# Amplitude and partial wave

Hoferichter:2023bjm decomposes the current matrix element as
$$
\langle0|j_\mu(0)|\pi^+(p_+)\pi^-(p_-)\pi^0(p_0)\rangle
=-\epsilon_{\mu\nu\alpha\beta}p_+^\nu p_-^\alpha p_0^\beta
\mathcal F(s,t,u;q^2),
$$
with $\mathcal F(s,t,u;q^2)=\mathcal F(s,q^2)+\mathcal F(t,q^2)+\mathcal F(u,q^2)$.
The dominant partial wave is
$$
f_1(s,q^2)=\mathcal F(s,q^2)+\widehat{\mathcal F}(s,q^2),
$$
where $\widehat{\mathcal F}$ is the crossed-channel angular average.
The single-variable function obeys the Khuri--Treiman integral equation with
the elastic $\pi\pi$ P-wave phase as input [@Hoferichter:2023bjm,
Sec. 2, Eqs. (2.1)--(2.8)].

# Overall normalization and fit

The Khuri--Treiman solution fixes the momentum dependence up to the function
$a(q^2)$. Its dispersive parameterization combines a Wess--Zumino--Witten
subtraction, explicit vector-resonance contributions, and a conformal
polynomial for non-resonant strength. The global fit includes the anomaly
normalization and fitted resonance and conformal parameters for
$e^+e^-\to3\pi$ data [@Hoferichter:2023bjm, Sec. 2 and Tables 1--3]. This is
the amplitude model subsequently used by Hoferichter:2025lcz.

# Scope and limits for the TFF input

The archived `anc/bareta3pi.dat` file tabulates the radiative-correction
function $\bar\eta_{3\pi}$, not the complex partial wave. The paper gives
the dispersive definitions and fit results, but its local source archive
does not supply a numerical complex grid of $f_1(s,q^2)$ or fit covariance
replicas. A spectral calculation still requires a numerical Khuri--Treiman
solution, specified $\pi\pi$ phase and continuation, and correlated fit
variations or a reproduced likelihood fit. Do not treat independent
marginal parameter errors as correlated replicas.

# Provenance

The source equations and ancillary file are indexed in
[the source archive](../../sources/Hoferichter:2023bjm/INDEX.md). The updated
global-fit choices and limitations are in
[Hoferichter:2025lcz](../Hoferichter:2025lcz/gamma-star-to-three-pion-fit.md)
[@Hoferichter:2025lcz].
