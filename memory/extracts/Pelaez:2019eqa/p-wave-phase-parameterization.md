---
title: Isovector P-wave phase from the global scattering parameterization
status: literature
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Pelaez:2019eqa]
provenance: literature
---

# Scope and provenance

The arXiv source `main-v3.tex`, Sec. III.B, Eqs. (20)–(24) and
the P-wave coefficient tables (source labels `tab:para11`, `tab:para112`, `tab:para113`) of [@Pelaez:2019eqa] provide the following physical-axis
parameterization below $\sqrt{s}=1.4$ GeV. The fit reproduces the Madrid
CFD/GKPY low-energy amplitudes. Solutions I, II and III represent distinct
higher-energy data choices; they are alternatives, not correlated samples.
The [source index](../../sources/Pelaez:2019eqa/INDEX.md) records the archive.

# Equations

Use $q_i=\sqrt{s/4-m_i^2}$ and $\sigma_i=2q_i/\sqrt{s}$, with the physical
upper-lip kaon loop

$$
\bar J(s,m_K)=\frac2\pi+\frac{\sigma_K(s)}\pi
\log\frac{\sigma_K(s)-1}{\sigma_K(s)+1}.
$$

The logarithm above kaon threshold has imaginary part $+\pi$.
Below threshold the loop is real. The conformal amplitude is

$$
t^1_{1,\mathrm{conf}}=\frac1{\sigma_\pi(\Phi^1_1-i)},\qquad
\Phi^1_1=\frac{\sqrt{s}(m_\rho^2-s)}{2q_\pi^3}
\left[\frac{2m_\pi^3}{m_\rho^2\sqrt{s}}+
\sum_{n=0}^4 B_n\omega(s)^n\right],
$$

$$
\omega(s)=\frac{\sqrt{s}-\alpha\sqrt{s_0-s}}
{\sqrt{s}+\alpha\sqrt{s_0-s}},\qquad
\alpha=0.3,\qquad s_0=(1.43\,\mathrm{GeV})^2.
$$

Include the inelastic factor even below the kaon threshold:

$$
t^1_{1,\mathrm{in}}=\frac{e^{i\delta^1_{1,\mathrm{in}}}-1}{2i\sigma_\pi},
\quad
\delta^1_{1,\mathrm{in}}=\bar J(s,m_K)
\left(K_0+K_1\frac{s}{m_K^2}\right)
\frac{q_\pi^3}{\sqrt{s}m_\pi^2}\frac{q_K^2}{m_K^2}.
$$

The full amplitude is

$$
t^1_1=t^1_{1,\mathrm{conf}}+t^1_{1,\mathrm{in}}
+2i\sigma_\pi t^1_{1,\mathrm{conf}}t^1_{1,\mathrm{in}},
\quad S^1_1=1+2i\sigma_\pi t^1_1=\eta^1_1e^{2i\delta^1_1}.
$$

Thus the physical phase is half the continuously unwrapped argument of
$S^1_1$. The quantity named $\delta^1_{1,\mathrm{in}}$ in the paper
contributes **half** its real part to the physical phase, since its
$S$ factor is $e^{i\delta^1_{1,\mathrm{in}}}$. Taking the conformal phase
alone would omit a small published contribution.

# Central coefficients with 2024 erratum

| Solution | $B_0$ | $B_1$ | $B_2$ | $B_3$ | $B_4$ | $K_0$ | $K_1$ |
|---|---|---|---|---|---|---|---|
| I | 0.96 | 0.09 | -0.07 | 0.58 | 1.39 | 0.049 | -0.0051 |
| II | 0.96 | 0.09 | -0.03 | 0.64 | 1.1 | 0.054 | -0.0062 |
| III | 0.97 | 0.07 | -0.003 | 0.67 | 1.0 | 0.045 | -0.0037 |

All solutions use $m_\rho=0.7749\pm0.0012$ GeV.
The erratum corrects $K_1$ for I/II and $B_0,B_2$ for III relative to the
arXiv v3 tables. Table errors are marginal errors only; no covariance is
published in these tables. The executable input manifest preserves these
errors separately and does not construct a multivariate probability model.

# Applicability to KT input

This phase is a documented Madrid alternative to the exact Bern input
used by the retained three-pion baseline. It is not a numerical reproduction
of that Bern input. Above kaon threshold the scattering $S$ matrix is
inelastic, whereas the single-channel KT solver uses only its phase.
The separate [phase preparation record](../../rcht/odd-parity/vppp-pipi-phase-input.md)
states the mass approximation, finite table domain and solver continuation.

Erratum provenance: [official 2024 erratum](https://link.springer.com/article/10.1140/epjc/s10052-024-13589-6),
DOI `10.1140/epjc/s10052-024-13589-6`, cited through the erratum note in
`Pelaez:2019eqa`. This erratum, rather than the unchanged arXiv v3 table,
is the source of the four corrected central entries above.
