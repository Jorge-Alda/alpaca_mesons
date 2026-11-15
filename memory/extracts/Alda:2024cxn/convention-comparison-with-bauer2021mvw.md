---
title: Sign and normalization comparison of Alda:2024cxn and Bauer:2021mvw for the ALP-photon loop functions
status: mixed
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2024cxn, Bauer:2021mvw]
provenance: literature
---

# Comparison of Alda:2024cxn and Bauer:2021mvw

Literature facts (equation numbers in the two sibling records): the loop functions $B_1$, $B_3$ and $f$ are identical in both papers (Bauer:2021mvw Eqs. (2.71), (2.82); Alda:2024cxn Eqs. (A1), (A3), (A5)).

Differences in coupling sign (literature facts):

- Photon operator: Bauer:2021mvw Eq. (2.10) has $+c_{\gamma\gamma}\frac{\alpha}{4\pi}\frac af F\tilde F$; Alda:2024cxn Eq. (1) has $-c^0_{\gamma\gamma}\frac{\alpha_{\rm em}}{4\pi}\frac{a}{f_a}F\tilde F$.
- Derivative coupling: Bauer:2021mvw Eq. (2.19) has $+\frac{\partial^\mu a}{2f}c_{ff}\bar f\gamma_\mu\gamma_5f$; Alda:2024cxn Eq. (1) has $-\frac{c_\tau}{2f_a}\partial_\mu a\,\bar\tau\gamma^\mu\gamma_5\tau$.
- Hence the dictionary $c_{\gamma\gamma}=-c^0_{\gamma\gamma}$, $c_{\tau\tau}=-c_\tau$, $f=f_a$ maps $C^{\rm eff}_{\gamma\gamma}$ of Bauer:2021mvw Eq. (2.70) to $-c^{{\rm eff},0}_{\gamma\gamma}$ of Alda:2024cxn Eq. (4), consistent with $\Gamma\propto|C^{\rm eff}|^2$, where $\frac{\alpha^2m_a^3}{64\pi^3}=\frac{\alpha_{\rm em}^2}{(4\pi)^3}m_a^3$ (same normalization). This dictionary is an inference from the displayed Lagrangians, not a statement of either paper.
- The common overall sign flip is a matter of the sign convention of $\epsilon^{0123}$ or of $a$. Bauer:2021mvw states $\epsilon^{0123}=1$; Alda:2024cxn does not state it. Resolution with the project convention (Alda:2026nrd) is open.

Derivative vs. pseudoscalar coupling (Alda:2024cxn Eqs. (1), (2)): integrating $-\frac{c_\tau}{2f_a}\partial_\mu a\,\bar\tau\gamma^\mu\gamma_5\tau$ by parts and using $\partial_\mu(\bar\tau\gamma^\mu\gamma_5\tau)=2im_\tau\bar\tau\gamma_5\tau+$ anomaly gives $+i\,m_\tau\frac{c_\tau}{f_a}a\bar\tau\gamma_5\tau$ plus a shift of the photon coefficient $c^0_{\gamma\gamma}\to c^0_{\gamma\gamma}+c_\tau$. This matches Eq. (2) (the derivation step is my check of the stated relation; the anomaly coefficient is taken from the source). Consequently in the pseudoscalar basis the loop gives $c_\tau[B_1-1]$, as in Alda:2024cxn Eq. (21), and $B_1\to1$ for $m_\tau\to0$ cancels the shift, as expected for a derivative coupling to a massless fermion.

Argument ordering: both papers use (ALP, other boson): Bauer:2021mvw $B_3(4m_f^2/m_a^2,4m_f^2/m_Z^2)$; Alda:2024cxn $B_3(4m_\tau^2/m_a^2,4m_\tau^2/s)$. Bauer:2021mvw Eq. (2.74) writes the one-argument form $B_1(4m_f^2/k^2)$ with $k$ the ALP momentum. Neither paper defines a function with both photons virtual.

Limits (from the closed form, my algebra): $B_3\to1$ for $m_f\to0$ (both arguments $\to0$); for $m_f\to\infty$, $B_3\to-\frac13\big(\frac1x+\frac1y\big)=-\frac{m_a^2+s}{12m_f^2}$, reducing to $B_1(x)\to-\frac{1}{3x}$ for $y\to\infty$.
