---
title: GomezDumm:2013sib dispersive pion vector form factor and Breit-Wigner warning
status: literature
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [GomezDumm:2013sib]
provenance: literature
---

# GomezDumm:2013sib — Dispersive Pion Vector Form Factor and Breit-Wigner Warning

## Dispersive Representation of $F_V^\pi(s)$

The paper constructs a dispersive representation of the charged pion vector form factor $F_V^\pi(s)$ for $\tau\to\pi\pi\nu_\tau$ decays [@GomezDumm:2013sib, Sec. 2]. The form factor is written as an $n$-subtracted Omnès representation [@GomezDumm:2013sib, Eq. (15)]:

$$
F_V^\pi(s) = \exp\left\{ \alpha_1 s + \frac{\alpha_2}{2}s^2 + \frac{s^3}{\pi}\int_{s_{\rm thr}}^\infty ds' \frac{\delta_1^1(s')}{(s')^3(s'-s-i\epsilon)} \right\}.
$$

The $\pi\pi$ P-wave phase shift $\delta_1^1(s)$ is obtained from the $\mathcal O(p^4)$ resonance chiral theory expression [@GomezDumm:2013sib, Eq. (10)]:

$$
\tan\delta_1^1(s) = \frac{\Im m F_V^{\pi(0)}(s)}{\Re e F_V^{\pi(0)}(s)},
$$

where $F_V^{\pi(0)}(s)$ includes the $\rho$ propagator with an energy-dependent width from chiral loops [@GomezDumm:2013sib, Eqs. (8)--(9)]. **Three subtractions** ($n=3$) yield a good fit to Belle data up to $s\simeq 1.5$ GeV$^2$ [@GomezDumm:2013sib, Sec. 2].

## Explicit Warning Against Fixed-Width Breit-Wigner

> "In order to account for the resonance width, here one should not simply replace $M^2-s$ by $M^2-s-iM\Gamma(s)$ in the effective propagator, since this would double count $\Im m[A_P(s)]$ and analyticity would be violated at $\mathcal O(p^6)$ in the chiral expansion." [@GomezDumm:2013sib, Sec. 2, lines 247--250]

The paper explains that the energy-dependent width from chiral loops already contains the imaginary part of the loop functions; adding it by hand in the propagator double-counts the absorptive part and violates analyticity. The correct procedure is to resum the loop functions into the denominator (or use the Omnès exponentiation).

## Relevance to $V\to 3\pi$ and $e^+e^-\to 3\pi$

While this paper focuses on the **two-body** form factor $F_V^\pi(s)$, the warning applies directly to any RChT/dispersive calculation of $V\to 3\pi$ or $e^+e^-\to 3\pi$ that uses a $\rho$ propagator. The $\rho$ exchange in $V\to 3\pi$ (Ruiz-Femenia Eq. 30) uses a tree-level propagator $1/(M_V^2-s_{ij})$. If one dresses this with a width, the same double-counting/analyticity issue arises unless the dispersive/Omnès resummation is used.

The dispersive approach of Niecknig:2012sj for $V\to 3\pi$ avoids this by using the Omnès function built from the empirical $\pi\pi$ phase shift, which automatically includes the correct analytic structure.