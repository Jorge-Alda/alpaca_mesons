---
title: BaBar ISR model ingredients for e+e- -> pi+ pi- pi0 line shape
status: literature
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [BaBar:2004ytv]
provenance: literature
---

# Model Ingredients in the BaBar ISR $e^+e^-\to\pi^+\pi^-\pi^0$ Analysis

This extract separates what the paper **explicitly includes** in its line-shape model from what it **does not include** or **approximates**.

## $\rho$ Exchange Contribution

**Included as the sole hadronic amplitude model.** The event generator for $e^+e^-\to 3\pi\gamma$ uses a model with an **intermediate $\rho\pi$ state** [@BaBar:2004ytv, Sec. 2, lines 915--916; Sec. 5, lines 1562--1566]. The $\rho\pi$ model is checked against high-statistics $e^+e^-$ experiments (SND, DAFNE) in the $\omega-\phi$ region; small deviations found in $\phi$ decays are negligible at BaBar's statistical level [@BaBar:2004ytv, Sec. 5, lines 1564--1567].

In the **fit to the mass spectrum (Sec. 6)**, the $V\to 3\pi$ decay is **assumed to proceed via the $\rho\pi$ intermediate state**, and $F_{\rho\pi}(m)$ is the $3\pi$ phase space volume calculated under this hypothesis [@BaBar:2004ytv, Sec. 6, lines 1689--1693]. The formulae for $F_{\rho\pi}$ are taken from Ref. [SND2003].

**Above the $\phi$**, the $\omega\pi$ intermediate state becomes noticeable. SND established that its contribution to the total $e^+e^-\to 3\pi$ cross section does not exceed 10% up to 1.4 GeV [@BaBar:2004ytv, Sec. 5, lines 1573--1578]. A narrow $\omega$ peak is seen in the $\pi^+\pi^-$ mass spectra for $1.1<M_{3\pi}<1.8$ GeV, with fraction ~10% (1.1--1.4 GeV) and ~6% (1.4--1.8 GeV), consistent with SND [@BaBar:2004ytv, Sec. 5, lines 1585--1589]. The detection efficiency difference for events with $0.77<M_{\pi^+\pi^-}<0.80$ GeV vs. full range is $(1.6\pm3.9)\%$ [@BaBar:2004ytv, Sec. 5, lines 1590--1592].

**Conclusion:** The measured line shape is **generated and fitted assuming $\rho\pi$ dominance** for $V\to 3\pi$. The $\omega\pi$ contribution is acknowledged but not included in the fit model. No other intermediate states (e.g., $f_2\pi$, $\rho'\pi$, direct contact) are included in the fit.

## Timelike $\omega/\phi$ Propagator with Mixing

**Included as relativistic Breit-Wigner propagators with mass-dependent widths and a fixed relative phase.** The fit uses [@BaBar:2004ytv, Eq. (14), lines 1684--1686]:

$$
D_V(m) = m_V^2 - m^2 - i m \Gamma_V(m), \quad \Gamma_V(m) = \sum_f \Gamma_f(m),
$$

with mass-dependent partial widths $\Gamma_f(m)$ for $\omega$ and $\phi$ taken from SND2003 [@BaBar:2004ytv, Sec. 6, lines 1687--1690]. The **$\omega$-$\phi$ relative phase is fixed to $\phi_\phi = (163\pm7)^\circ$ from SND2003** [@BaBar:2004ytv, Sec. 6, lines 1767--1768]. The phases of $\omega, \omega', \omega''$ are fixed at $0^\circ, 180^\circ, 0^\circ$ [@BaBar:2004ytv, Sec. 6, line 1769--1770, citing Clegg].

**No dynamical $\omega$-$\phi$ mixing** is implemented in the propagator; the mixing effect is absorbed into the fixed relative phase and the fitted branching-fraction products. The propagators are standard Breit-Wigners with energy-dependent widths, not a coupled-channel mixing matrix.

## Direct VVP Contact Terms

**NOT included in the fit model.** The fit assumes $V\to 3\pi$ proceeds **exclusively via $\rho\pi$** [@BaBar:2004ytv, Sec. 6, line 1690: "We assume that $V\to 3\pi$ decay proceeds via the $\rho\pi$ intermediate state"]. There is no direct $V\to 3\pi$ contact term (VVP operator) in the amplitude. The paper does not discuss VVP operators or chiral Lagrangians; it is a phenomenological resonance fit.

The paper **does note** that vector-meson dominance (VMD) alone gives $\Gamma(\omega\to 3\pi) \approx 1.4$ MeV, while the measured value is $7.52\pm0.06$ MeV, requiring a "sizeable direct contribution" [@BaBar:2004ytv, Sec. 1, lines 856--860, citing Ruiz-Femenia:2003jdx]. However, this statement is in the Introduction as motivation/context; **the actual fit does not include a direct term**—the $\rho\pi$ model with fitted resonance parameters absorbs the total width.

## Final-State and Radiative Corrections

### Radiative Corrections (ISR)
- **Included:** The radiative correction factor $R=0.9994$ accounts for multiple photon emission by initial $e^\pm$ (structure function method) [@BaBar:2004ytv, Sec. 6, lines 1695--1707].
- **Excluded from $R$:** Leptonic and hadronic vacuum polarization corrections [@BaBar:2004ytv, Sec. 6, lines 1710--1713]. These are conventionally absorbed into the "dressed" resonance electronic widths.
- **Theoretical uncertainty on $R$:** $\le 1\%$ [@BaBar:2004ytv, Sec. 6, line 1709].

### Final-State Radiation (FSR)
- **Not explicitly discussed** as a separate correction. The Phokhara generator (used for $e^+e^-\to\pi^+\pi^-\gamma$ background) includes NLO QED corrections [@BaBar:2004ytv, Sec. 2, lines 932--934], but the signal generator (based on Czyz-Kuhn) uses the structure function method for extra soft-photon radiation [@BaBar:2004ytv, Sec. 2, lines 920--922]. FSR from final-state pions is not separately quantified.

### Vacuum Polarization
- **Not applied to the cross section.** The paper measures the "dressed" cross section (including VP) by convention [@BaBar:2004ytv, Sec. 7, lines 1998--2002]. The radiative correction factor $R$ explicitly excludes VP.

## Summary: What the Measured Line Shape Constrains

| Ingredient | In Fit Model? | In Simulation? | Notes |
|------------|---------------|----------------|-------|
| $\rho\pi$ intermediate state ($V\to\rho\pi\to 3\pi$) | Yes (sole hadronic amplitude) | Yes (event generator) | Phase space $F_{\rho\pi}(m)$ from SND2003 |
| $\omega\pi$ intermediate state | No | Acknowledged $\le 10\%$ up to 1.4 GeV | Seen in $\pi^+\pi^-$ mass spectra above $\phi$ |
| Direct $V\to 3\pi$ contact term (VVP) | No | No | VMD alone fails for $\omega\to 3\pi$ width (Intro only) |
| $\omega$-$\phi$ mixing (dynamical) | No | No | Fixed relative phase $\phi_\phi=163^\circ$ from SND |
| Relativistic BW with mass-dependent width | Yes | N/A (fit only) | $\Gamma_V(m)$ from SND2003 for $\omega,\phi$ |
| ISR radiative corrections (multiple photons) | Yes (factor $R=0.9994$) | Yes (structure function) | Generator-level, no detector sim |
| Vacuum polarization | No (absorbed in dressed widths) | No | Convention: "dressed" cross section |
| FSR from final pions | Not quantified | Part of Phokhara for background | Not separately corrected for signal |

**Key implication for RChT fitting:** The BaBar ISR Born cross section **does not directly constrain VVP contact terms** ($d_i$). It constrains the **product** of the $\rho\pi$-mediated amplitude (which depends on $G_V$ and the VVP combination $d_1+8d_2-d_3$ per Ruiz-Femenia:2003jdx Eq. (30)) and the timelike vector propagators. A direct VVP term would require a dedicated amplitude analysis (e.g., Dalitz plot) or a line-shape model that includes it explicitly.