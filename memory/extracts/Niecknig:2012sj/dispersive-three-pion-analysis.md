---
title: Dispersive analysis of phi -> 3pi and omega -> 3pi from Niecknig:2012sj and GomezDumm:2013sib
status: literature
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [Niecknig:2012sj, GomezDumm:2013sib]
provenance: literature
---

# Dispersive Analysis of $\phi\to 3\pi$ and $\omega\to 3\pi$

## Niecknig:2012sj — Dispersive Framework for $V\to 3\pi$

### Amplitude Decomposition

The decay amplitude for $V=\omega,\phi$ is decomposed as [@Niecknig:2012sj, Eq. (1)]:

$$
\mathcal M(s,t,u)=i\epsilon_{\mu\nu\alpha\beta}n^\mu p_+^\nu p_-^\alpha p_0^\beta \mathcal F(s,t,u),
$$

where $n^\mu$ is the polarization vector of the decaying vector meson. The scalar function $\mathcal F(s,t,u)$ is decomposed into single-variable functions [@Niecknig:2012sj, Eq. (7)]:

$$
\mathcal F(s,t,u) = \mathcal F(s) + \mathcal F(t) + \mathcal F(u),
$$

where $\mathcal F(s)$ has only a right-hand cut. This corresponds to a symmetrized P-wave expansion.

### Unitarity and Integral Equation

The discontinuity of $\mathcal F(s)$ is related to the $\pi\pi$ P-wave phase shift $\delta(s)$ via a Watson-like unitarity relation [@Niecknig:2012sj, Eq. (11)]:

$$
\mathrm{disc}\,\mathcal F(s) = 2i\bigl(\mathcal F(s)+\hat{\mathcal F}(s)\bigr)\theta(s-4m_\pi^2)\sin\delta(s)e^{-i\delta(s)},
$$

where $\hat{\mathcal F}(s)$ is the inhomogeneity from crossed-channel rescattering (angular averages over $\mathcal F$). The solution is [@Niecknig:2012sj, Eq. (14)]:

$$
\mathcal F(s) = \Omega(s)\left\{ a + \frac{s}{\pi}\int_{4m_\pi^2}^\infty \frac{ds'}{s'} \frac{\sin\delta(s')\hat{\mathcal F}(s')}{|\Omega(s')|(s'-s)} \right\},
$$

with the Omnès function $\Omega(s)$ [@Niecknig:2012sj, Eq. (12)]:

$$
\Omega(s) = \exp\left\{ \frac{s}{\pi}\int_{4m_\pi^2}^\infty ds' \frac{\delta(s')}{s'(s'-s)} \right\}.
$$

The subtraction constant $a$ fixes the overall normalization (fitted to $\Gamma_{V\to 3\pi}$). The Dalitz plot shape is a **pure prediction** from the $\pi\pi$ phase shift.

### Oversubtraction (Two Subtractions)

To suppress high-energy uncertainties, a twice-subtracted representation is used [@Niecknig:2012sj, Eq. (18)]:

$$
\mathcal F(s) = \Omega(s)\left\{ a + b's + \frac{s^2}{\pi}\int_{4m_\pi^2}^\infty \frac{ds'}{s'^2} \frac{\sin\delta(s')\hat{\mathcal F}(s')}{|\Omega(s')|(s'-s)} \right\},
$$

where $b'$ is a complex subtraction constant. If $b'$ satisfies the sum rule [@Niecknig:2012sj, Eq. (19)], this is equivalent to the once-subtracted form. In fits to data, $b'$ is allowed to vary, introducing two real parameters (modulus and phase).

### Fit to KLOE and CMD-2 $\phi\to 3\pi$ Dalitz Data

- **KLOE:** $\sim 2\times 10^6$ events, bin size 8.75 MeV, efficiency-corrected [@Niecknig:2012sj, Sec. 4.1].
- **CMD-2:** $\sim 8\times 10^4$ events, convoluted with efficiency matrices [@Niecknig:2012sj, Sec. 4.1].
- **Isospin breaking:** Charged pion mass used in phase shift; physical phase-space factor with $m_{\pi^0}$ applied to amplitude squared [@Niecknig:2012sj, Eq. (23)].
- **$\omega\pi^0$ background:** Added as a resonance term $a_\omega e^{i\phi_\omega} M_\omega^2/(M_\omega^2 - i\sqrt{s}\Gamma_\omega - s)$ [@Niecknig:2012sj, Eq. (24)].

**Key results:**
- Once-subtracted: $\chi^2/\mathrm{ndof} = 1.17$ (best, upper boundary) to $2.06$ (Omnès only) [@Niecknig:2012sj, Table 1].
- Twice-subtracted: $\chi^2/\mathrm{ndof} \approx 1.02$ (both Bern and Madrid phases) [@Niecknig:2012sj, Table 2].
- Crossed-channel rescattering **significantly improves** the fit ($\chi^2$ drops from $\sim 2.0$ to $\sim 1.2$ for once-subtracted) [@Niecknig:2012sj, Sec. 4.1].
- The fitted $b$ parameter deviates from its sum-rule value, suggesting high-energy/inelastic effects [@Niecknig:2012sj, Sec. 4.2].

### $\omega\to 3\pi$ Predictions

Dalitz plot parameters predicted for future measurements (WASA, KLOE-2) [@Niecknig:2012sj, Sec. 4.3, Table 3]. At least two parameters ($\alpha,\beta$) needed for 1% precision. Leading slope $\alpha > 0$ fixed by $\rho$ peak kinematics.

### Warning: Fixed-Width Models vs. Dispersive Data

> "descriptions in terms of vector-meson dominance and improved tree-level models (mostly incorporating a finite width of the intermediate $\rho$ resonances)... which do not fully respect analyticity and unitarity constraints, still seem to be state of the art." [@Niecknig:2012sj, Sec. 1, lines 117--121]

The paper explicitly contrasts the dispersive approach (which respects analyticity and unitarity exactly) with isobar models (Breit-Wigner $\rho$ exchange). The dispersive amplitude includes **crossed-channel three-particle cuts** that modify the $\rho$ peak shape (shift + attenuation) in a way that fixed-width models cannot capture. The $\rho$ peak in $\phi\to 3\pi$ is **attenuated** by crossed-channel effects, while in $\omega\to 3\pi$ it is **enhanced and shifted** [@Niecknig:2012sj, Sec. 3.2, lines 582--588].

**Explicit warning:** "it is obvious that these [isobar models] are doomed to fail for the decays of $\omega$ or $\phi$ into three pions, where the influence of the $\rho$ resonance is already strongly felt ($\omega$) or even shows up in the form of resonance bands inside the Dalitz plot ($\phi$)." [@Niecknig:2012sj, Sec. 1, lines 106--109]

## GomezDumm:2013sib — Dispersive Pion Vector Form Factor

While focused on $F_V^\pi(s)$ for $\tau\to\pi\pi\nu_\tau$, this work uses a similar dispersive framework [@GomezDumm:2013sib, Sec. 2]. The form factor is written as an $n$-subtracted Omnès representation [@GomezDumm:2013sib, Eq. (15)]:

$$
F_V^\pi(s) = \exp\left\{ \alpha_1 s + \frac{\alpha_2}{2}s^2 + \frac{s^3}{\pi}\int_{s_{\rm thr}}^\infty ds' \frac{\delta_1^1(s')}{(s')^3(s'-s-i\epsilon)} \right\}.
$$

The phase shift $\delta_1^1(s)$ is taken from the $\mathcal O(p^4)$ resonance chiral theory expression [@GomezDumm:2013sib, Eq. (10)]:

$$
\tan\delta_1^1(s) = \frac{\Im m F_V^{\pi(0)}(s)}{\Re e F_V^{\pi(0)}(s)},
$$

with $F_V^{\pi(0)}(s)$ having the $\rho$ propagator with energy-dependent width from chiral loops [@GomezDumm:2013sib, Eqs. (8)--(9)]. **Three subtractions** ($n=3$) are needed for a good fit to Belle data up to $s\simeq 1.5$ GeV$^2$ [@GomezDumm:2013sib, Sec. 2].

**Key statement:** "In order to account for the resonance width, here one should not simply replace $M^2-s$ by $M^2-s-iM\Gamma(s)$ in the effective propagator, since this would double count $\Im m[A_P(s)]$ and analyticity would be violated at $\mathcal O(p^6)$ in the chiral expansion." [@GomezDumm:2013sib, Sec. 2, lines 247--250]

This is a direct warning against using a fixed-width Breit-Wigner propagator in a chiral/dispersive framework.

## Comparison: Dispersive vs. RChT for $e^+e^-\to 3\pi$

| Aspect | Niecknig (Dispersive) | GomezDumm (Dispersive $F_V^\pi$) | RChT (Ruiz-Femenia) |
|--------|----------------------|----------------------------------|---------------------|
| **Input** | $\pi\pi$ P-wave phase shift (Roy eqs.) | $\pi\pi$ phase shift from RChT $\mathcal O(p^4)$ | VVP/VVP operators + short-distance constraints |
| **$\rho$ propagator** | Omnès function $\Omega(s)$ | Omnès with RChT phase shift | $1/(M_V^2-s_{ij})$ (tree-level) |
| **Crossed-channel** | Included via integral equation | Not in form factor (2-body) | Not in tree-level amplitude |
| **$\omega$-$\phi$ mixing** | Not in formalism | $\rho$-$\omega$ mixing in $F_V^\pi$ | Ideal mixing fixed |
| **Direct $V\to 3\pi$** | Absorbed in subtraction constants | N/A (2-body form factor) | Required but not in published basis |
| **Warning on BW** | Isobar models "doomed to fail" | Fixed-width BW violates analyticity | Tree-level only, no width |

## Implications for BaBar ISR Fit

The BaBar ISR Born cross section is fitted with a **coherent sum of fixed-width Breit-Wigner resonances** [@BaBar:2004ytv, Eq. (14)]. Niecknig:2012sj demonstrates that:
1. This model does not respect analyticity/unitarity.
2. Crossed-channel rescattering modifies the $\rho$ peak shape significantly.
3. The dispersive Dalitz plot shape is a **parameter-free prediction** (up to normalization) from the $\pi\pi$ phase shift.

**For an RChT predictor:** Using the BaBar ISR spectrum to constrain $d_i$ couplings requires a line-shape model that includes at least the dispersive $\rho$ propagator (Omnès function) and crossed-channel effects. A simple Breit-Wigner sum with fixed widths (as in BaBar's fit) will not correctly map the RChT amplitude parameters to the measured cross section, especially in the $\phi$ region where the $\rho$ bands are visible.