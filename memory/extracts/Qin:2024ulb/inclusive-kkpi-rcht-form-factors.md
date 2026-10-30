---
title: Inclusive KK pi form factors and VPPP sensitivity in RChT
status: literature
created: 2026-10-30
updated: 2026-10-30
bibtex_keys: [Qin:2024ulb, BaBar:2007ceh]
provenance: literature
---

# Amplitudes and cross-section normalization

For $e^+e^-\to K\bar K\pi$, [@Qin:2024ulb] writes
$$
\mathcal M=-\frac{4\pi\alpha}{q^2}\,iF_R(q^2,s,t)
\varepsilon_{\mu\nu\alpha\beta}p_1^\nu p_2^\alpha p_3^\beta
\bar v(q_1)\gamma^\mu u(q_2),
\qquad F_R=F_a+F_b+F_c+F_d.
$$
The inclusive Born cross section is
$$
\sigma(q^2)=\frac{\alpha^2}{192\pi q^6}
\int ds\,dt\,\phi(q^2,s,t)|F_R(q^2,s,t)|^2.
$$
The source gives $\phi$ and the physical Dalitz boundaries in its Eq. (12)
and the following definitions. The masses in form-factor numerators use
averaged pseudoscalar masses, while phase space uses the channel's physical
masses [@Qin:2024ulb].

# VPPP combinations

The local VPPP contribution is part of $F_b$. For $K^+K^-\pi^0$,
$$
F_b^1=\frac{2\sqrt2F_V}{F^3M_V}
\left[
\left(D_\rho(q^2)+\frac13D_\omega(q^2)\right)G_{R1}(q^2,s)
-\frac23D_\phi(q^2)G_{R2}(q^2,s)\right],
$$
where $D_R(x)=1/(M_R^2-x)$ at the source's zero-width level. For
$K_S^0K^\pm\pi^\mp$,
$$
F_b^3=\frac{2\sqrt2F_V}{F^3M_V}
\left[
(g_1+2g_2-g_3)(t-u)D_\rho(q^2)
-\frac13D_\omega(q^2)G_{R1}(q^2,s)
+\frac23D_\phi(q^2)G_{R2}(q^2,s)\right].
$$
The source definitions (Appendix Eq. CRDR) are
$$
\begin{aligned}
G_{R1}(q^2,s)={}&(g_1+2g_2-g_3)(q^2+s-4m_K^2-m_\pi^2)
+4(2g_4+g_5)m_K^2\\
&-4g_2(q^2-2m_K^2-m_\pi^2)-4g_4(m_K^2-m_\pi^2),\\
G_{R2}(q^2,s)={}&(g_1+2g_2-g_3)(q^2-s-m_\pi^2)
+2(2g_4+g_5)m_\pi^2\\
&-2g_2(q^2-2m_K^2-m_\pi^2)+4g_4(m_K^2-m_\pi^2).
\end{aligned}
$$
Thus the inclusive $K_S K^\pm\pi^\mp$ and $K^+K^-\pi^0$ rates retain VPPP
terms. The isobar-only Table 7 component of [@BaBar:2007ceh] does not
include this full amplitude.

# Other terms and higher resonances

The same source supplies $F_c$ from VJP and $K^*$ exchange and $F_d$ from
VVP exchange, with their complete channel formulae in Eqs.
(kpkmpi0), (ksklpi0), and (kskpi). The WZW term $F_a$ is nonzero
for the neutral-pion charged-kaon channels and vanishes for
$K_S K^\pm\pi^\mp$ [@Qin:2024ulb].

To model higher states, [@Qin:2024ulb] extends a ground-state propagator by
$$
D_R(x)\longrightarrow D_R(x)
+\frac{\beta^{R'}_j}{M_{R'}^2-x}
+\frac{\beta^{R''}_j}{M_{R''}^2-x},
$$
with channel-dependent residues and energy-dependent widths for the rho and
kaon vector families. The analysis also applies
$\exp[-q^6/\Lambda^6]$ to regulate high-energy growth and compares
$\Lambda=2.3$--$2.7\,\mathrm{GeV}$ [@Qin:2024ulb]. Its reported VPPP values
come from a joint fit to several experiments, invariant-mass spectra and
resonance widths, with short-distance constraints and tensor-resonance
terms. They are not a BaBar-only determination in the repository's
preliminary coefficient scheme. The paper explicitly notes that fitted
$g_i$ may compensate for incomplete final-state interactions
[@Qin:2024ulb].
