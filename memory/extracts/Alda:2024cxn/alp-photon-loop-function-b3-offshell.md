---
title: Alda:2024cxn loop functions B_1, B_3, A_3 for a tau-philic ALP, scalar and pseudoscalar
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2024cxn]
provenance: literature
---

# Alda:2024cxn: loop functions with one real and one virtual photon

Source: arXiv:2407.18296 (`memory/sources/Alda:2024cxn/src/draft.tex`). Equation numbers are those of the compiled arXiv version (appendix equations are (A1)-(A17)).

## Conventions

Derivative basis, Eq. (1):
$$\mathcal{L}^a_{d\le5}=\tfrac12(\partial_\mu a)^2-\tfrac12m_a^2a^2-\frac{\alpha_{\rm em}}{4\pi}c^0_{\gamma\gamma}\frac{a}{f_a}F_{\mu\nu}\tilde F^{\mu\nu}-\frac{c_\tau}{2f_a}(\partial_\mu a)\bar\tau\gamma^\mu\gamma_5\tau .$$
Chirality-flipping (pseudoscalar) basis, Eq. (2):
$$\mathcal{L}^a_{\rm int}\supset i\,m_\tau\frac{c_\tau}{f_a}\,a\,\bar\tau\gamma_5\tau-\frac{\alpha_{\rm em}}{4\pi}\frac{c^0_{\gamma\gamma}+c_\tau}{f_a}\,a\,F_{\mu\nu}\tilde F^{\mu\nu}.$$
The source gives both bases but no derivation, and does not state the Levi-Civita sign or the definition of $\tilde F^{\mu\nu}$.

General spin-0 Lagrangian, Eqs. (19), (20): $\mathcal{L}=i g^\varphi_\tau\varphi\bar\tau\gamma_5\tau+g^\phi_\tau\phi\bar\tau\tau$ and $\mathcal{L}=-g^\varphi_{\gamma\gamma}\frac{\alpha_{\rm em}}{4\pi}\varphi F_{\mu\nu}\tilde F^{\mu\nu}-g^\phi_{\gamma\gamma}\frac{\alpha_{\rm em}}{4\pi}\phi F_{\mu\nu}F^{\mu\nu}$.

## Effective couplings

Two real photons, Eqs. (3), (4): $\Gamma(a\to\gamma\gamma)=\frac{\alpha_{\rm em}^2}{(4\pi)^3}m_a^3\frac{|c^{{\rm eff},0}_{\gamma\gamma}|^2}{f_a^2}$, $\;c^{{\rm eff},0}_{\gamma\gamma}=c^0_{\gamma\gamma}+c_\tau B_1\!\big(\tfrac{4m_\tau^2}{m_a^2}\big)$.

One virtual photon of squared momentum $s$, Eqs. (6), (7):
$$\sigma_{\rm NR}(e^+e^-\to\gamma a)=\frac{\alpha_{\rm em}^3}{24\pi^2}\frac{|c^{{\rm eff},s}_{\gamma\gamma}|^2}{f_a^2}\Big(1-\frac{m_a^2}{s}\Big)^3,\qquad
c^{{\rm eff},s}_{\gamma\gamma}=c^0_{\gamma\gamma}+c_\tau\,B_3\!\Big(\frac{4m_\tau^2}{m_a^2},\frac{4m_\tau^2}{s}\Big).$$
Pseudoscalar and scalar in the Yukawa basis, Eqs. (21), (22):
$$g^{\varphi,s}_{\gamma\gamma}=\frac{g^\varphi_\tau}{m_\tau}\Big[B_3\Big(\frac{4m_\tau^2}{m_a^2},\frac{4m_\tau^2}{s}\Big)-1\Big],\qquad
g^{\phi,s}_{\gamma\gamma}=\frac{g^\phi_\tau}{m_\tau}A_3\Big(\frac{4m_\tau^2}{m_\phi^2},\frac{4m_\tau^2}{s}\Big).$$
Here $s$ is the momentum squared of the off-shell photon.

## Loop functions (App. A), Eqs. (A1)-(A12)

$$B_3(x,y)=1+\frac{xy}{x-y}\big[f^2(x)-f^2(y)\big],\quad
f(x)=\begin{cases}\arcsin\frac1{\sqrt x},&x\ge1\\ \frac\pi2+\frac i2\log\frac{1+\sqrt{1-x}}{1-\sqrt{1-x}},&x<1\end{cases}$$
$$A_3(x,y)=\frac{xy}{(x-y)^2}\Big[x-y+(x-y+xy)\big(f^2(x)-f^2(y)\big)-x\big(g(x)-g(y)\big)\Big],$$
$$g(x)=\begin{cases}-\sqrt{x-1}\arccos\!\big(1-\frac2x\big),&x\ge1\\ \sqrt{1-x}\big(\log\frac{2-x-2\sqrt{1-x}}{x}+i\pi\big),&x<1.\end{cases}$$
Limits: $B_1(x)=B_3(x,y\gg1)\simeq1-xf^2(x)$ (A5), $B_1(x\gg1)\simeq-\frac1{3x}$ (A6), valid for $s\ll4m_\tau^2$ (two on-shell photons). For $s\gg4m_\tau^2$:
$B_3(x,y\ll1)\simeq1+y\arcsin^2\frac1{\sqrt x}+\frac y4\log^2\!\big(\frac{-y}4\big)$ (A9), $B_3(x\gg1,y\ll1)\simeq1+\frac y4\log^2\!\big(\frac{-y}4\big)+\frac yx$ (A10).
Text after Eq. (7): $B_3$ is almost independent of $m_a$ for $m_a\ll2m_\tau$; $|B_3|\approx0.5$ at BESIII and $\approx1.2$ at Belle II energies.

Arguments: first argument is $4m_f^2/m_a^2$ (ALP), second is $4m_f^2/s$ (virtual photon). The function is symmetric under $x\leftrightarrow y$ by its closed form. Both photons off shell is not treated.

Note on the $y\ll1$ forms (A9)-(A10): they use $\log(-y/4)$, i.e. the timelike continuation with $s+i\epsilon$.
