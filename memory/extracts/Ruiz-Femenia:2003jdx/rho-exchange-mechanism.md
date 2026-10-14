---
title: Rho exchange mechanism for pi+ pi- pair in V -> 3pi and e+e- -> 3pi
status: literature
created: 2026-10-04
updated: 2026-10-15
bibtex_keys: [Ruiz-Femenia:2003jdx, Niecknig:2012sj, BaBar:2004ytv]
provenance: literature
---

Displayed tensors and odd couplings use $\epsilon^{0123}=+1$ and
$\epsilon_{0123}=-1$ for $(+---)$. The canonical coefficients are
$c_i=-c_i^R$, $d_i=-d_i^R$ relative to the opposite source convention.
The [conversion audit](../../rcht/odd-parity/levi-civita-convention.md)
derives the action-preserving map; archived source text remains literal.


# $\rho$ Exchange Mechanism for the $\pi^+\pi^-$ Pair

## Ruiz-Femenia:2003jdx — $\rho$-Mediated $\omega\to 3\pi$ (Eq. 30)

The paper gives the $\rho$-exchange amplitude for $\omega\to\pi^+\pi^-\pi^0$ (including cyclic permutations) [@Ruiz-Femenia:2003jdx, Eq. (30), lines 913--924]:

$$
i\mathcal M_{\omega\to3\pi}=i\epsilon_{\alpha\beta\rho\sigma}
k_1^\alpha k_2^\beta k_3^\rho\epsilon_\omega^\sigma
\frac{8G_V}{M_\omega F^3}\left[
\frac{m_\pi^2(d_1+8d_2-d_3)+(M_\omega^2+s_{12})d_3}{M_V^2-s_{12}}
+\{s_{12}\to s_{13}\}+\{s_{12}\to s_{23}\}\right].
$$

**Key features:**
- The $\rho$ propagator is $1/(M_V^2-s_{ij})$ with the **octet mass $M_V$** (not the physical $\rho$ mass) [@Ruiz-Femenia:2003jdx, lines 827--830].
- The amplitude depends on $d_1+8d_2-d_3$ (multiplying $m_\pi^2$) and $d_3$ (multiplying $M_\omega^2+s_{ij}$).
- $G_V = F/\sqrt{2}$ from short-distance constraints [@Ruiz-Femenia:2003jdx, line 949].
- The VMD-only width is $\sim 1.4$ MeV, far below the experimental $7.52\pm0.06$ MeV, requiring a **direct $\omega 3\pi$ term** [@Ruiz-Femenia:2003jdx, lines 950--956].

## Niecknig:2012sj — Dispersive View: $\rho\pi$ is Not a Distinct State

In dispersion theory, **there is no meaningful $\rho\pi$ two-body intermediate state** distinct from $3\pi$ [@Niecknig:2012sj, lines 267--272]:

> "Note finally that in the sense of a dispersion-theory analysis, there is no $\rho\pi$ two-body intermediate state contributing to the $\phi$ decay in a similar fashion as $K\bar K$... $\rho\pi$ is no distinct state from $3\pi$, the effect conventionally encoded that way shows up in our analysis as the resonant two-body $\pi\pi$ P-wave interaction. From the point of view of dispersion theory, there is no meaningful way to differentiate between $\rho\pi$ and $3\pi$ final states in $\phi$ decays."

The $\rho$ resonance appears as the **$\pi\pi$ P-wave phase shift** $\delta_1^1(s)$ in the Omnès function. The dispersive amplitude automatically includes crossed-channel rescattering (three-particle cuts) that go beyond a simple $\rho$-exchange isobar model.

## BaBar:2004ytv — Event Generator Uses $\rho\pi$ Model

The BaBar ISR analysis uses an event generator with an **intermediate $\rho\pi$ state** [@BaBar:2004ytv, Sec. 2, lines 915--916; Sec. 5, lines 1562--1566]. The $\rho\pi$ model is checked against SND/DAFNE data in the $\omega-\phi$ region; small deviations in $\phi$ decays are negligible at BaBar's statistics [@BaBar:2004ytv, Sec. 5, lines 1564--1567].

In the fit to the mass spectrum (Sec. 6), the $V\to 3\pi$ decay is **assumed to proceed via $\rho\pi$**, and $F_{\rho\pi}(m)$ is the $3\pi$ phase space volume under this hypothesis [@BaBar:2004ytv, Sec. 6, lines 1689--1693].

**Above the $\phi$**, an $\omega\pi$ intermediate state contributes $\le 10\%$ up to 1.4 GeV [@BaBar:2004ytv, Sec. 5, lines 1573--1578]. A narrow $\omega$ peak is seen in $\pi^+\pi^-$ mass spectra for $1.1<M_{3\pi}<1.8$ GeV (fraction $\sim 10\%$ and $\sim 6\%$), consistent with SND [@BaBar:2004ytv, Sec. 5, lines 1585--1589].

## Role in the Measured Line Shape

| Framework | $\rho$ Exchange Role |
|-----------|---------------------|
| **Ruiz-Femenia (RChT)** | Explicit $\rho$ propagator $1/(M_V^2-s_{ij})$; amplitude $\propto d_3$ and $d_1+8d_2-d_3$; VMD-only width too small by factor $\sim 4$ |
| **Niecknig (Dispersive)** | $\rho$ is the $\pi\pi$ P-wave phase shift; no separate $\rho\pi$ state; crossed-channel effects modify the $\rho$ peak (shift + attenuation) |
| **BaBar (Experimental)** | $\rho\pi$ model used for simulation efficiency; fit assumes $V\to 3\pi$ via $\rho\pi$ phase space $F_{\rho\pi}(m)$; $\omega\pi$ contribution acknowledged but not in fit |

## $\omega$ Exchange in the $\pi^+\pi^-$ Pair

- **Ruiz-Femenia:** Not discussed for $e^+e^-\to 3\pi$; $\omega$ appears as external decaying particle in $\omega\to 3\pi$.
- **Niecknig:** $\omega\pi$ intermediate states contribute to $\pi\pi$ P-wave inelasticity above $\omega\pi$ threshold ($\sqrt{s}\simeq 0.922$ GeV) [@Niecknig:2012sj, Sec. 6]. The analytic continuation of $\omega\to 3\pi$ gives $\pi\pi\to\omega\pi$ amplitude.
- **BaBar:** $\omega\pi$ intermediate state visible in $\pi^+\pi^-$ mass spectra above $\phi$; efficiency difference $(1.6\pm3.9)\%$ [@BaBar:2004ytv, Sec. 5, lines 1590--1592].

## Summary for RChT Fitting

The BaBar ISR Born cross section is **generated and fitted assuming $\rho\pi$ dominance** for $V\to 3\pi$. In RChT language, this corresponds to the $\rho$-exchange diagram with the $VVP$ vertex proportional to $d_3$ and $d_1+8d_2-d_3$. The direct $VVP$ contact term (which Ruiz-Femenia argues is necessary) is **not included** in the BaBar model. A dispersive analysis (Niecknig) shows that crossed-channel rescattering modifies the $\rho$ peak shape significantly, which a simple $\rho$-exchange model with fixed width does not capture.
