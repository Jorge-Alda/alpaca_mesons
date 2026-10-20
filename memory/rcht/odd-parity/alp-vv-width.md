---
title: Analytic widths for ALP decays to neutral vector pairs
status: derived
created: 2026-10-20
updated: 2026-10-20
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Kinematics

Use the amplitude and conventions in [the no-mixing ALP vector-pair
derivation](alp-vv-amplitude.md). For $a(p)\to V_i(q_1)V_j(q_2)$ define
$$
\lambda(m_a^2,M_i^2,M_j^2)=
[m_a^2-(M_i+M_j)^2][m_a^2-(M_i-M_j)^2].
$$
For the tensor amplitude
$$
\mathcal M_{ij}=G_{ij}\epsilon_{\mu\nu\rho\sigma}
\varepsilon_1^\mu\varepsilon_2^\nu q_1^\rho q_2^\sigma,
$$
the polarization sum is $\sum|\mathcal M_{ij}|^2=|G_{ij}|^2\lambda/2$.
The integrated two-body phase space then gives
$$
\Gamma(a\to V_iV_j)=\frac{|G_{ij}|^2}{32\pi m_a^3}
\lambda(m_a^2,M_i^2,M_j^2)^{3/2}\frac{1}{1+\delta_{ij}},
\quad m_a>M_i+M_j,
$$
and zero for a closed channel. The identical-particle factor is $1/2!$ for
$\omega\omega$ and $\phi\phi$, and one for $\omega\phi$.

Thus
$$
\Gamma_{\omega\omega}=\frac{|G_{\omega\omega}|^2
\lambda(m_a^2,M_\omega^2,M_\omega^2)^{3/2}}{64\pi m_a^3},
$$
$$
\Gamma_{\omega\phi}=\frac{|G_{\omega\phi}|^2
\lambda(m_a^2,M_\omega^2,M_\phi^2)^{3/2}}{32\pi m_a^3},
$$
$$
\Gamma_{\phi\phi}=\frac{|G_{\phi\phi}|^2
\lambda(m_a^2,M_\phi^2,M_\phi^2)^{3/2}}{64\pi m_a^3}.
$$
The dimensions are $[G_{ij}]=-1$ and $[\Gamma]=1$ in mass units.

# Domain

These are narrow-vector tree-level widths. Near threshold, broad-vector
spectral integrals and the daughter decay amplitudes must replace stable
two-body phase space. No such spectral treatment is included here.
