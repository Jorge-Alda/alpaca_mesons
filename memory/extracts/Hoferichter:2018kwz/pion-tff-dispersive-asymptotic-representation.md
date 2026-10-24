---
title: Pion TFF dispersive and asymptotic representation
status: literature
created: 2026-10-24
updated: 2026-10-24
bibtex_keys: [Hoferichter:2018kwz]
provenance: literature
---

# Pion TFF dispersive and asymptotic representation

## Convention and decomposition

The pion transition form factor (TFF) is defined from the QCD matrix element
in Eq. (2.1), with $\epsilon^{0123}=+1$. In this source convention its real-
photon normalization is positive,
$$
F_{pi^0\gamma^*\gamma^*}(0,0)=F_{\pi\gamma\gamma}
=\frac{1}{4\pi^2 F_\pi},
$$
where the paper uses $F_\pi=92.28(9)\,\mathrm{MeV}$ including its stated
quark-mass renormalization [@Hoferichter:2018kwz, Eqs. (2.1), (2.2)]. A sign
conversion is required before combining this scalar with conventions in which
$F_\pi(0,0)<0$; the magnitude of the decay constant also follows the source's
convention and should not silently replace $f_\pi=93\,\mathrm{MeV}$.

Assuming exact isospin, one photon has isovector quantum numbers and the
other isoscalar quantum numbers. The full TFF is the Bose-symmetrized sum
$$
F_{\pi^0\gamma^*\gamma^*}(q_1^2,q_2^2)
=F_{vs}(q_1^2,q_2^2)+F_{vs}(q_2^2,q_1^2)
$$
[@Hoferichter:2018kwz, Eq. (2.3)]. The paper's final representation is
$$
F_{\pi^0\gamma^*\gamma^*}
=F_{\pi^0\gamma^*\gamma^*}^{\mathrm{disp}}
+F_{\pi^0\gamma^*\gamma^*}^{\mathrm{eff}}
+F_{\pi^0\gamma^*\gamma^*}^{\mathrm{asym}},
$$
where the first term contains low-energy intermediate states, the second
represents higher intermediate states and normalization/data not supplied by
the low-energy part, and the third imposes the large-virtuality behavior
[@Hoferichter:2018kwz, Eq. (1.1) and surrounding text].

## Low-energy spectral contribution

The source starts from an unsubtracted dispersion relation for one isovector
virtuality,
$$
F_{vs}(q_1^2,q_2^2)=\frac{1}{12\pi^2}
\int_{4M_\pi^2}^{\infty}dx\,
\frac{q_\pi^3(x)\,[F_\pi^V(x)]^*f_1(x,q_2^2)}
{\sqrt{x}\,(x-q_1^2)},\qquad
q_\pi(x)=\sqrt{x/4-M_\pi^2},
$$
with the pion vector form factor $F_\pi^V$ and the $P$-wave partial wave
$f_1$ for $\gamma_s^*\to3\pi$ as inputs [@Hoferichter:2018kwz,
Eq. (3.11)]. To make both virtualities explicit, the low-energy double-
spectral piece is
$$
F_{\pi^0\gamma^*\gamma^*}^{\mathrm{disp}}(-Q_1^2,-Q_2^2)
=\frac{1}{\pi^2}\int_{4M_\pi^2}^{s_v}dx
\int_{s_{\mathrm{thr}}}^{s_s}dy\,
\frac{\rho^{\mathrm{disp}}(x,y)}{(x+Q_1^2)(y+Q_2^2)}
+(Q_1\leftrightarrow Q_2),
$$
$$
\rho^{\mathrm{disp}}(x,y)=\frac{q_\pi^3(x)}{12\pi\sqrt{x}}
\operatorname{Im}\!\left([F_\pi^V(x)]^*f_1(x,y)\right),
$$
as in Eq. (3.14) [@Hoferichter:2018kwz]. The cutoffs $s_v$ and $s_s$
bound the explicit isovector and isoscalar spectral regions. The construction
uses elastic $2\pi$ rescattering in the isovector channel and a Khuri--Treiman
description of the $3\pi$ amplitude; exact full three-body unitarity is not
implemented, and the crossed-channel left-hand cuts are approximated by
pion-pole terms [@Hoferichter:2018kwz, discussion below Eq. (2.3),
Eqs. (2.8), (3.14)].

## Effective and asymptotic pieces

The effective contribution is a symmetric pole,
$$
F_{\pi^0\gamma^*\gamma^*}^{\mathrm{eff}}(q_1^2,q_2^2)
=\frac{g_{\rho'}}{4\pi^2F_\pi}
\frac{M_{\rho'}^4}{(M_{\rho'}^2-q_1^2)(M_{\rho'}^2-q_2^2)},
$$
whose coupling restores the anomaly sum rule and whose mass is fit to
singly-virtual spacelike data [@Hoferichter:2018kwz, Eq. (4.10)]. It
collects higher intermediate states and high-energy contributions omitted
from the explicit low-energy spectral integral; the paper reports typical
$g_{\rho'}\sim10\%$ and $M_{\rho'}\sim1.5$--$2\,\mathrm{GeV}$, and notes
that the explicit low-energy representation already extends partly into the
1--2 GeV region [@Hoferichter:2018kwz, discussion around Eq. (4.10)].

At leading order with the asymptotic pion distribution amplitude
$\phi_\pi(u)=6u(1-u)$, the hard expression is
$$
F_{\pi^0\gamma^*\gamma^*}(q_1^2,q_2^2)
=-\frac{2F_\pi}{3}\int_0^1du\,
\frac{\phi_\pi(u)}{u q_1^2+(1-u)q_2^2}
+O(q_i^{-4})
$$
[@Hoferichter:2018kwz, Eq. (4.1)]. Writing
$\omega=(q_1^2-q_2^2)/(q_1^2+q_2^2)$ gives
$$
F_{\pi^0\gamma^*\gamma^*}
=-\frac{4F_\pi}{3}\frac{f(\omega)}{q_1^2+q_2^2}
+O(q_i^{-4}),\qquad
f(\omega)=\int_0^1du\,
\frac{\phi_\pi(u)}{u(1-\omega)+(1-u)(1+\omega)},
$$
as in Eqs. (4.2), (4.3) [@Hoferichter:2018kwz]. The source obtains the
diagonal spacelike coefficient $Q^2F(-Q^2,-Q^2)=2F_\pi/3$ and formally
the single-virtual coefficient $Q^2F(-Q^2,0)=2F_\pi$. It explicitly cautions
that this OPE expression is justified only for $|\omega|<1/2$; the
$|\omega|=1$ Brodsky--Lepage result is a formal endpoint evaluation of that
formula, not a rigorous consequence of that OPE domain [@Hoferichter:2018kwz,
Eqs. (4.4), (4.5) and text following]. In particular, continuing the
fixed-ratio leading-twist formula all the way to a singly-virtual endpoint is
not a strict OPE result.

For the asymptotic term, the source uses the double-spectral density
$\rho^{\mathrm{asym}}(x,y)=-2\pi^2F_\pi xy\,\delta''(x-y)$ and keeps the
doubly-high-energy region above a continuum threshold $s_m$. This gives
$$
F_{\pi^0\gamma^*\gamma^*}^{\mathrm{asym}}(q_1^2,q_2^2)
=2F_\pi\int_{s_m}^{\infty}dx\,
\frac{q_1^2q_2^2}{(x-q_1^2)^2(x-q_2^2)^2},
$$
which yields the leading-twist result for two nonzero hard virtualities
[@Hoferichter:2018kwz, Eqs. (4.7)--(4.9)]. Mixed low/high spectral regions
are omitted because the asymptotic density used there vanishes; the effective
pole is assumed to represent their net effect. Thus the source provides a
constructive dispersive-plus-continuum model, but it does not uniquely
determine mixed-region spectral strength [@Hoferichter:2018kwz, discussion
before Eq. (4.9)]. Since this term is proportional to $q_1^2q_2^2$, it
vanishes when either photon is real. The singly-virtual coefficient is
therefore supplied by the low-energy dispersive and effective terms: the
low-energy part alone accounts for about $55\%$ of the Brodsky--Lepage value
in their representation, and the effective-pole parameters are fit to
singly-virtual spacelike data [@Hoferichter:2018kwz, discussion around
Eq. (4.10)].

## Limits for using this representation in the RChT fit

This is a pion-specific TFF reconstruction for the pion-pole contribution
to hadronic light-by-light scattering. Its spectral functions are tied to
the pion vector form factor, the $\gamma^*\to3\pi$ partial wave, and fits to
$e^+e^-\to3\pi$ and singly-virtual pion-TFF data [@Hoferichter:2018kwz,
Eqs. (3.11), (3.14), (4.10)]. It supplies neither the source-resolved ALP
short-distance coefficient nor a direct mapping of its spectral parameters
to the RChT $c_i,d_i$ combinations. Adding its explicit $2\pi/3\pi$ strength
to an existing resonance TFF would require a matching prescription to avoid
counting the same rho-region strength twice. These limits are applications
to the present fit, not claims made by the paper.
