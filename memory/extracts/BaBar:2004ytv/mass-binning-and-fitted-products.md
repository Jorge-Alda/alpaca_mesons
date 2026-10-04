---
title: BaBar ISR mass binning, threshold, and fitted resonance products
status: literature
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [BaBar:2004ytv]
provenance: literature
---

# BaBar ISR Mass Binning, Threshold, and Fitted Resonance Products

## Mass Binning and Range

The published Born cross-section spectrum (Table~\ref{sumtab}) covers the $e^+e^-$ c.m. energy range **1.05 to 3.00 GeV** [@BaBar:2004ytv, Sec. 7, line 1934; Abstract line 745]. The table contains **78 rows** with uniform **25 MeV bin width** [@BaBar:2004ytv, Table~\ref{sumtab} caption]. Bin centers run from $\sqrt{s'}=1.0625$ GeV to 2.9875 GeV in steps of 0.025 GeV.

The **fit to the $3\pi$ mass spectrum** (Sec. 6) uses a different binning: **2.5 MeV/$c^2$ for the $\omega-\phi$ mass range** (approximately 0.70--1.05 GeV) and **25 MeV/$c^2$ for masses above the $\phi$** [@BaBar:2004ytv, Sec. 6, lines 1754--1755]. The fitted mass region is restricted to **masses below 1.8 GeV** [@BaBar:2004ytv, Sec. 6, line 1777]; above this, "data in the higher mass range would require a more complicated function... There are no reliable models available in the literature" [@BaBar:2004ytv, Sec. 6, lines 1778--1783].

## Relation to Dalitz/Detected Mass

The invariant mass $m \equiv \sqrt{s'}$ is the **$3\pi$ invariant mass** ($M_{3\pi}$) reconstructed from the detected $\pi^+\pi^-\pi^0$ final state [@BaBar:2004ytv, Eq. (16) text: "$m\equiv \sqrt{s^\prime}$ is the $3\pi$ invariant mass"]. The Dalitz plot variables ($s_{12}, s_{13}, s_{23}$) are not used in the cross-section measurement; the analysis integrates over the full Dalitz plot. The event generator for simulation uses a $\rho\pi$ intermediate state model [@BaBar:2004ytv, Sec. 2, line 915--916; Sec. 5, line 1562--1566], which determines the Dalitz distribution used for efficiency calculation.

## $\omega$ and $\phi$ Bin Edges

**No bin edge in Table~\ref{sumtab} coincides with the $\omega$ or $\phi$ masses.** The table starts at 1.05 GeV, which is **above the $\phi$ mass (1.020 GeV)**. The $\omega$ (0.782 GeV) and $\phi$ (1.020 GeV) lie **below the published spectrum range**. The fit in Sec. 6 covers the $\omega-\phi$ region (0.70--1.05 GeV) with fine 2.5 MeV binning, but this region is **not included in the published Table~\ref{sumtab}**.

## Fitted Branching-Fraction Products

The fit to the $3\pi$ mass spectrum (below 1.8 GeV) extracts the following products of branching fractions [@BaBar:2004ytv, Sec. 6, Eqs. (17)--(18); Summary Eqs. (19)--(20)]:

| Resonance | $\mathcal B(V\to e^+e^-)\mathcal B(V\to 3\pi)$ | Statistical | Systematic |
|-----------|-----------------------------------------------|-------------|------------|
| $\omega$  | $(6.70 \pm 0.06 \pm 0.27) \times 10^{-5}$     | $\pm 0.06$  | $\pm 0.27$ |
| $\phi$    | $(4.30 \pm 0.08 \pm 0.21) \times 10^{-5}$     | $\pm 0.08$  | $\pm 0.21$ |
| $\omega'$ | $(0.82 \pm 0.05 \pm 0.06) \times 10^{-6}$     | $\pm 0.05$  | $\pm 0.06$ |
| $\omega''$| $(1.3 \pm 0.1 \pm 0.1) \times 10^{-6}$        | $\pm 0.1$   | $\pm 0.1$  |

The systematic error for $\mathcal B(V\to e^+e^-)\mathcal B(V\to 3\pi)$ **includes**:
- Statistical error from simulation
- Error on the efficiency correction (Table~\ref{effcorr})
- **1.2% uncertainty in the luminosity**
- 1% theoretical uncertainty on the radiative correction
- Background subtraction uncertainty (0.4% at $\omega$, 0.6% at $\phi$)
- Uncertainty from errors on $\Gamma_\omega$, $\Gamma_\phi$, and $\phi_\phi$ (1% at $\omega$, 2.8% at $\phi$) [@BaBar:2004ytv, Sec. 6, lines 1815--1822]

The **1.2% luminosity uncertainty** is explicitly listed as a separate component of the systematic error on the fitted products. It is **not folded into the per-bin cross-section systematic errors** in Table~\ref{sumtab} (which only include event-count, efficiency, and radiative-correction systematics) [@BaBar:2004ytv, Sec. 7, lines 1992--1994; metadata.json confirms this separation].

## Fit Model for Resonance Amplitudes

The Born cross section in the fit is modeled as a coherent sum of four resonances ($\omega, \phi, \omega', \omega''$) [@BaBar:2004ytv, Eq. (14)]:

$$
\sigma_{3\pi}(m)=\frac{12\pi}{m^3}F_{\rho\pi}(m)
\left|\sum_{V=\omega,\phi,\omega',\omega''}
\frac{\Gamma_V m_V^{3/2}\sqrt{{\cal B}(V\to e^+e^-){\cal B}(V\to 3\pi)}}{D_V(m)}
\frac{e^{i\phi_V}}{\sqrt{F_{\rho\pi}(m_V)}}\right|^2,
$$

where $D_V(m)=m_V^2 - m^2 - i m\Gamma_V(m)$, $\Gamma_V(m)=\sum_f \Gamma_f(m)$, and $F_{\rho\pi}(m)$ is the $3\pi$ phase space volume assuming $V\to 3\pi$ proceeds via $\rho\pi$ [@BaBar:2004ytv, Sec. 6, lines 1687--1693]. The $\omega$-$\phi$ relative phase is fixed to $\phi_\phi = (163\pm7)^\circ$ from SND [@BaBar:2004ytv, Sec. 6, line 1767--1768]. Phases of $\omega, \omega', \omega''$ are fixed at $0^\circ, 180^\circ, 0^\circ$ [@BaBar:2004ytv, Sec. 6, line 1769--1770, citing @Clegg].

The fit **does not include** the $e^+e^-\to\omega\pi\gamma\to 3\pi\gamma$ process via excited $\rho$ states [@BaBar:2004ytv, Sec. 6, lines 1771--1776].