---
title: Bauer:2021mvw ALP-photon loop functions B_1, B_3 and the off-shell vertex
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Bauer:2021mvw]
provenance: literature
---

# Bauer:2021mvw: loop functions $B_1$, $B_3$ and the off-shell ALP-photon vertex

Source: arXiv:2110.10698 (`memory/sources/Bauer:2021mvw/src/paper.tex`, Sec. 2). Equation numbers are those of the compiled arXiv version.

## Conventions

- Lagrangian, Eq. (2.1), $D\le5$, above the weak scale: $\mathcal{L}\ni \frac{\partial^\mu a}{f}\sum_F\bar\psi_F\bm{c}_F\gamma_\mu\psi_F + c_{BB}\frac{\alpha_1}{4\pi}\frac{a}{f}B_{\mu\nu}\tilde B^{\mu\nu}+\dots$, with $\tilde B^{\mu\nu}=\tfrac12\epsilon^{\mu\nu\alpha\beta}B_{\alpha\beta}$ and $\epsilon^{0123}=1$ (stated below Eq. (2.1)).
- Weak-scale Lagrangian, Eq. (2.10): $+\,c_{\gamma\gamma}\frac{\alpha}{4\pi}\frac{a}{f}F_{\mu\nu}\tilde F^{\mu\nu}$ (positive sign), $c_{\gamma\gamma}=c_{WW}+c_{BB}$, Eq. (2.11); $c_{\gamma Z}\frac{\alpha}{2\pi s_wc_w}\frac{a}{f}F_{\mu\nu}\tilde Z^{\mu\nu}$.
- Fermion couplings, Eq. (2.14): $c_{f_if_i}(\mu)\equiv[k_f(\mu)]_{ii}-[k_F(\mu)]_{ii}$ ($f$ right-handed, $F$ left-handed). Below the weak scale, Eq. (2.19): $\mathcal{L}^{\rm diag}_{\rm fermion}=\sum_{f\ne t}\frac{\partial^\mu a}{2f}\,c_{ff}\,\bar f\gamma_\mu\gamma_5 f$.
- Pseudoscalar (Yukawa-like) form: the source states only that $(\tilde{\bm Y}_f)_{ii}=i\,y_{f_i}c_{f_if_i}$ (text after Eq. (2.17)), i.e. it does not write $a\bar f\gamma_5 f$ explicitly in this paper. The explicit derivative-to-pseudoscalar rewriting is not given here (see the comparison record under `Alda:2024cxn`).

## Loop functions

On-shell ALP, one photon pair on shell, with $\tau_f\equiv 4m_f^2/m_a^2$, Eqs. (2.70), (2.71):
$$\Gamma(a\to\gamma\gamma)=\frac{\alpha^2 m_a^3}{64\pi^3 f^2}\left|C_{\gamma\gamma}^{\rm eff}\right|^2,\qquad
C_{\gamma\gamma}^{\rm eff}(m_a)=c_{\gamma\gamma}+\sum_{f\ne t}N_c^fQ_f^2\,c_{ff}(m_a)\,B_1(\tau_f)\quad(m_a\gg\Lambda_{\rm QCD}),$$
$$B_1(\tau)=1-\tau f^2(\tau),\qquad
f(\tau)=\begin{cases}\arcsin\frac{1}{\sqrt\tau}, & \tau\ge1,\\[1mm] \frac\pi2+\frac i2\ln\frac{1+\sqrt{1-\tau}}{1-\sqrt{1-\tau}}, & \tau<1.\end{cases}$$
Limits (text after Eq. (2.71)): $B_1\approx1$ for $m_f\ll m_a$; $B_1\approx -m_a^2/(12m_f^2)$ for $m_f\gg m_a$. Below the QCD scale, Eq. (2.72), leptons enter as $\sum_\ell c_{\ell\ell}B_1(\tau_\ell)$ and $c,b$ quarks as $3Q_q^2c_{qq}(\mu_0)B_1(\tau_q)$.

Two different virtualities (real photon with momentum $q$ and a $Z$ of mass $m_Z$), Eqs. (2.81), (2.82):
$$C_{\gamma Z}^{\rm eff}=c_{\gamma Z}+\sum_fN_c^fQ_f\Big(\tfrac12T_3^f-Q_fs_w^2\Big)c_{ff}(m_Z)\,B_3\!\Big(\frac{4m_f^2}{m_a^2},\frac{4m_f^2}{m_Z^2}\Big),$$
$$B_3(\tau_1,\tau_2)=1+\frac{\tau_1\tau_2}{\tau_1-\tau_2}\big[f^2(\tau_1)-f^2(\tau_2)\big],\qquad f\text{ as in Eq. (2.71)}.$$
Text: $B_3\approx1$ for all light fermions other than the top; $B_3(4m_t^2/m_a^2,4m_t^2/m_Z^2)\approx-0.024$. This enters $\Gamma(Z\to\gamma a)=\frac{m_Z^3}{96\pi^3f^2}\frac{\alpha\,\alpha(m_Z)}{s_w^2c_w^2}|C_{\gamma Z}^{\rm eff}|^2(1-m_a^2/m_Z^2)^3$, Eq. (2.80).

## Off-shell ALP-photon vertex (one real photon), Eqs. (2.73), (2.74)

Kinematics: $q$ inflowing on-shell photon ($q^2=0$, index $\mu$); $k$ outflowing off-shell ALP ($k^2\ne m_a^2$); $p=k-q$ inflowing off-shell photon (index $\alpha$). One-loop result (finite):
$$\Gamma^{\mu\alpha}_{\gamma\gamma^*a^*}(q,k)=\frac{i\alpha}{\pi f}\epsilon^{\mu\alpha\beta\gamma}q_\beta k_\gamma\Big\{c_{\gamma\gamma}+\sum_fN_c^fQ_f^2c_{ff}\Big[1-\int_0^1\!dx\!\int_0^1\!dy\,\frac{m_f^2}{m_f^2-x(1-x)(k^2-2y\,k\cdot q)-i\epsilon}\Big]\Big\}.$$
Dropping the $k\cdot q$ term (allowed for $g-2$, not for flavor-changing dipoles), the $y$ integral gives 1 and
$$\Gamma^{\mu\alpha}_{\gamma\gamma^*a^*}\to\frac{i\alpha}{\pi f}\epsilon^{\mu\alpha\beta\gamma}q_\beta k_\gamma\Big[c_{\gamma\gamma}+\sum_fN_c^fQ_f^2c_{ff}\,B_1\!\Big(\frac{4m_f^2}{k^2}\Big)\Big].$$
For $k^2\to m_a^2$ this reduces to Eq. (2.70). Heavy fermions ($m_f^2\gg|k^2|$) decouple; light fermions give 1 in the bracket.

## Inference (not stated in the source)

Using $p^2=k^2-2k\cdot q$ for $q^2=0$, the exponent combination is $k^2-2yk\cdot q=(1-y)k^2+y\,p^2$. In the same-kinematics integral (2.73) the bracket is then
$$1-\int_0^1\!dx\!\int_0^1\!dy\,\frac{m_f^2}{m_f^2-x(1-x)[(1-y)k^2+y\,p^2]}=B_3\!\Big(\frac{4m_f^2}{k^2},\frac{4m_f^2}{p^2}\Big).$$
This identification was checked numerically (scratch computation, four kinematic points including spacelike $p^2<0$ and timelike $p^2$ below threshold; agreement to about $10^{-13}$) and is not stated in Bauer:2021mvw. The check is not stored under `scripts/`.

Caveat: the printed $f(\tau)$ for $\tau<0$ (spacelike photon) must be continued with $f(\tau)=\tfrac i2\ln\frac{\sqrt{1-\tau}+1}{\sqrt{1-\tau}-1}$ so that $f^2<0$ is real; the $\tau<1$ branch as printed is ambiguous by $\pm i\pi$ in the logarithm for $\tau<0$.

## Not covered

Both photons off shell is not treated in Bauer:2021mvw (only $q^2=0$).
