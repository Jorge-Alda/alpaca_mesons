---
title: Dispersive TFF and asymptotic matching for dilepton decays
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Messerli:2025rnv, Hoferichter:2020lap]
provenance: literature
---

# Low-energy representation and normalization

Messerli:2025rnv decomposes each normalized $η^{(\prime)}$ TFF into isovector, isoscalar, effective-pole, and asymptotic pieces [@Messerli:2025rnv, Eq. (22)]. The isovector contribution is a finite-cutoff double-spectral integral, while narrow isoscalar poles describe the isoscalar part [@Messerli:2025rnv, Eqs. (16), (28)]. Effective poles restore the measured TFF normalization not saturated by these low-energy pieces and represent higher resonances; their parameters are constrained by normalization and singly-virtual spacelike data [@Messerli:2025rnv, Sec. 3.3, Eqs. (37)--(38)]. Thus the pole terms are part of the low-energy interpolation and are not, by themselves, a subtraction prescription for the asymptotic continuum.

# Massless asymptotic spectral term

Starting from the leading-twist light-cone convolution, the paper gives equal-virtuality and singly-virtual limits $Q^2F(-Q^2,-Q^2)\to\bar F^P_{\rm asym}/3$ and $Q^2F(-Q^2,0)\to\bar F^P_{\rm asym}$ [@Messerli:2025rnv, Eqs. (39)--(41)]. It notes that the singly-virtual limit is more sensitive to the distribution amplitude and is not a strict local OPE result. The asymptotic contribution is rewritten as a double-spectral term with a lower threshold $s_{\rm m}$:
$$
F^{\rm asym}_{P\gamma^*\gamma^*}(q_1^2,q_2^2)=\bar F^P_{\rm asym}\int_{s_{\rm m}}^\infty dx\,\frac{q_1^2q_2^2}{(x-q_1^2)^2(x-q_2^2)^2}.
$$
The corresponding density is proportional to $-\pi^2\bar F^P_{\rm asym}q_1^2q_2^2\delta''(q_1^2-q_2^2)$, with boundary terms chosen to produce this representation [@Messerli:2025rnv, Eqs. (42)--(43)]. The explicit virtuality factors make this term vanish if either photon is real, so it leaves $F(0,0)$ unchanged; it supplies the doubly-virtual asymptotic limit. The singly-virtual behavior in their construction is instead fixed by low-energy dispersive/effective-pole contributions constrained by spacelike data. The coefficient is inferred through a superconvergence relation of that representation, not borrowed independently from a meson decay constant [@Messerli:2025rnv, Sec. 4.1].

# Mass corrections and subtraction

For nonzero pseudoscalar mass the hard convolution becomes
$$
F_{P\gamma^*\gamma^*}(q_1^2,q_2^2)=-\frac{\bar F^P_{\rm asym}}{3}\int_0^1du\,\frac{\phi_P(u)}{u q_1^2+(1-u)q_2^2-u(1-u)M_P^2}.
$$
Messerli:2025rnv states that this kinematic mass correction is not a complete higher-order light-cone result. It cites the generalized double-spectral construction and subtractions that enforce $F_P^{\rm asym}(q_1^2,q_2^2)=\mathcal O(q_1^2q_2^2)$ near the origin [@Messerli:2025rnv, Eq. (44)--(45), citing Hoferichter:2020lap, Zanke:2021wiq, and Hoferichter:2024bae]. This small-virtuality behavior protects the low-energy normalization. For the $η$ and $η'$, the matching thresholds are varied around $1.4$ and $1.5\,\mathrm{GeV}^2$, respectively, and their variation contributes to the asymptotic uncertainty [@Messerli:2025rnv, Secs. 4.1 and 5.1]. These numerical values and fitted coefficients are meson-specific and do not determine an ALP threshold or coefficient.

# Scope and overlap

The added asymptotic spectral term represents strength above $s_{\rm m}$; it does not assert that all low-energy cuts are unchanged if an existing low-energy density already extends into that region. A construction must define the low/high spectral partition and avoid double counting. The effective poles are separately refit when other contributions or normalization constraints change [@Messerli:2025rnv, Sec. 5.1]. For an ALP application, the source-specific coefficient, threshold, and shared spectral support remain inputs to be derived; the numerical $η^{(\prime)}$ coefficients and thresholds cannot be transplanted.
