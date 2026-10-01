---
title: On-shell vector-to-pseudoscalar-photon amplitudes and widths
status: literature
created: 2026-09-30
updated: 2026-09-30
bibtex_keys: [Ruiz-Femenia:2003jdx]
provenance: literature
---

# Exact published formulas

The local auditable source is `memory/sources/Ruiz-Femenia:2003jdx/vvp.tex`; source
line numbers below refer to that file. The antisymmetric tensor convention is
`epsilon_{0123}=+1` (lines 352--354), and the octet matrices are normalized as
`Phi=lambda^a phi^a/sqrt(2)` and `V=lambda^a V^a/sqrt(2)` (Eqs. (4),(8), lines
206--219 and 291--304). Ideal mixing is
$$
|\omega\rangle=\sqrt{2/3}|\omega_1\rangle+\sqrt{1/3}|\omega_8\rangle,
\qquad
|\phi\rangle=-\sqrt{1/3}|\omega_1\rangle+\sqrt{2/3}|\omega_8\rangle
$$
(lines 784--797). The singlet is added as `omega_1/sqrt(3)` to the octet
matrix. These mixing statements are needed for any singlet-containing channel.

The even-parity tensor normalization is
$$
\mathcal L_2(V)=\frac{F_V}{2\sqrt2}\langle V_{\mu\nu}f_+^{\mu\nu}\rangle
+\frac{iG_V}{\sqrt2}\langle V_{\mu\nu}u^\mu u^\nu\rangle,
$$
with `f_+` defined from the left/right external field strengths (Eq. (7),
lines 278--290). The odd terms use
$$
\mathcal L_{VJP}=\sum_{a=1}^7\frac{c_a}{M_V}\mathcal O^a_{VJP},
\qquad
\mathcal L_{VVP}=\sum_{a=1}^4d_a\mathcal O^a_{VVP}
$$
(Eq. (12), lines 415--430). Thus $c_a$ are dimensionless while $d_a$ carry
the normalization of the displayed operators. The source does not print a
separate $V$--photon mixing equation or an explicit `v_mu=... A_mu`
substitution; the factors of $e$ in the following published amplitudes are
the operative photon normalization.

For $\omega(q)\to\pi(k)\gamma$, the uneliminated direct and $\rho$-exchange
amplitudes are (Eq. (27), source lines 810--831):
$$
i\mathcal M^{\rm direct}=i\epsilon_{\alpha\beta\rho\sigma}
\epsilon_\omega^\alpha\epsilon_\gamma^\beta q^\rho k^\sigma
\frac{2\sqrt2e}{M_\omega M_VF}
\left[(c_2-c_1+c_5-2c_6)M_\omega^2
+(c_1+c_2+8c_3-c_5)m_\pi^2\right],
$$
$$
i\mathcal M^\rho=-i\epsilon_{\alpha\beta\rho\sigma}
\epsilon_\omega^\alpha\epsilon_\gamma^\beta q^\rho k^\sigma
\frac{4e}{M_V^2M_\omega}\frac{F_V}{F}
\left[d_3M_\omega^2+(d_1+8d_2-d_3)m_\pi^2\right].
$$
After the OPE constraints, the published on-shell result is (Eq. (28),
source lines 832--843):
$$
i\mathcal M_{\omega\to\pi\gamma}=i\epsilon_{\alpha\beta\rho\sigma}
\epsilon_\omega^\alpha\epsilon_\gamma^\beta q^\rho k^\sigma
\frac{e}{F_V}\left[\frac{N_C}{8\pi^2}\frac{M_\omega}{F}
-\frac F2\frac{M_\omega}{M_V^2}\left(1+\frac{m_\pi^2}{M_\omega^2}\right)\right].
$$
Using $F_V=\sqrt2F$, its width is (Eq. (29), source lines 856--868):
$$
\Gamma(\omega\to\pi\gamma)=\frac{\alpha}{192}M_\omega
\left(1-\frac{m_\pi^2}{M_\omega^2}\right)^3
\left[\frac{N_C}{4\pi^2}\frac{M_\omega^2}{F^2}
-\frac{M_\omega^2}{M_V^2}\left(1+\frac{m_\pi^2}{M_\omega^2}\right)\right]^2.
$$

The tensor/Lorentz form and the two displayed coupling combinations generalize
channel-by-channel after replacing the flavor trace, external masses, vector
pole, and mixing factors. The numerical OPE prediction does **not** generalize
automatically: the authors explicitly say that $K^*\to K\gamma$ and
$\phi\to\eta\gamma$ require flavor-resolved short-distance matching (source
lines 880--890).
