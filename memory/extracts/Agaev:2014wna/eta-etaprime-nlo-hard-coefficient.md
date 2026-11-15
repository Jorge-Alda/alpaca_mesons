---
title: NLO quark and LO gluon hard coefficients for the eta and eta-prime transition form factors, and heavy-quark treatment
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Agaev:2014wna]
provenance: literature
---

# Source and scope

Read from `memory/sources/Agaev:2014wna/src/140914.tex` (arXiv 1409.4311), section on
the leading-twist contribution [@Agaev:2014wna]. Labels are LaTeX labels. Order in
$\alpha_s$: NLO for the quark coefficient function, leading order for the gluon
coefficient function; one virtuality is large ($q_1^2=-Q^2$) and the other is
$q_2^2=-q^2\to0$ (singly virtual limit).

# Conventions

- $\gamma_5=i\gamma^0\gamma^1\gamma^2\gamma^3$ and $\epsilon_{0123}=+1$ (after the gluon DA
  definition), so $\epsilon^{0123}=-1$ in the metric $(+---)$. This is the opposite
  $\epsilon$ orientation from the project ($\epsilon^{0123}=+1$).
- Label `eq:Fgamma`: $\int d^4x\,e^{iq_1x}\langle M(p)|T\{j^{\rm em}_\mu(x)j^{\rm em}_\nu(0)\}|0\rangle=ie^2\varepsilon_{\mu\nu\alpha\beta}q_1^\alpha q_2^\beta F_{\gamma^*\gamma^*\to M}(q_1^2,q_2^2)$
  with $q_2=q_1+p$: production kinematics, overall factor $ie^2$ in the definition.
  The sign and phase of $F$ therefore differ from the project's decay convention
  $T^{\mu\nu}=\epsilon^{\mu\nu\rho\sigma}q_{1\rho}q_{2\sigma}F$. Only magnitudes
  below should be compared with project quantities.
- Decay constants: $\langle0|J^{(i)}_{\mu5}|M(p)\rangle=if^{(i)}_Mp_\mu$ with
  $J^{(1)}_{\mu5}=\frac1{\sqrt3}(\bar u\gamma_\mu\gamma_5u+\bar d\gamma_\mu\gamma_5d+\bar s\gamma_\mu\gamma_5s)$,
  $J^{(8)}_{\mu5}=\frac1{\sqrt6}(\bar u\ldots+\bar d\ldots-2\bar s\ldots)$,
  $J^{(q)}_{\mu5}=\frac1{\sqrt2}(\bar u\ldots+\bar d\ldots)$, $J^{(s)}_{\mu5}=\bar s\gamma_\mu\gamma_5s$.
  Octet and singlet couplings are related to the quark-flavor ones by the rotation
  with $\varphi_0=\arctan\sqrt2$ (label `eq:SOQF`).

# Factorization formula (label `eq:lt`)

$$
F_{\gamma^*\gamma\to M}(Q^2)=\frac{f_M^{(8)}}{3\sqrt6}\int_0^1du\,T_H^{(8)}\phi_M^{(8)}
+\frac{2f_M^{(1)}}{3\sqrt3}\int_0^1du\,\bigl(T_H^{(1)}\phi_M^{(1)}+T_H^{(g)}\phi_M^{(g)}\bigr),
$$
with a two-gluon distribution amplitude $\phi_M^{(g)}$ (antisymmetric in $u\leftrightarrow\bar u$,
expanded in $C^{5/2}_{n-1}$, label `phig`) [@Agaev:2014wna].

# Coefficient functions

- Quark, NLO $\overline{\rm MS}$, same for octet and singlet (label `eq:NLOcf`),
  attributed to delAguila:1981nk, Phys. Rev. D 28 (1983) 524 and Yad. Fiz. 44 (1986) 507:
  $$
  T_H^{\rm NLO}=\frac{2}{uQ^2}\Bigl\{1+C_F\frac{\alpha_s(\mu)}{2\pi}\Bigl[\tfrac12\ln^2u-\tfrac12\frac u{\bar u}\ln u-\tfrac92+\bigl(\tfrac32+\ln u\bigr)\ln\frac{Q^2}{\mu^2}\Bigr]\Bigr\}.
  $$
- Gluon, light quarks, LO in $\alpha_s$ (label `lightbox`), from Kroll:2002nt and Kroll:2013iwa:
  $T_H^{g}|_{\rm light}=-C_F\frac{\alpha_s(\mu)}{2\pi}\frac{2\ln u}{\bar u^2Q^2}\{\frac1u-3+\frac12\ln u+\ln\frac{Q^2}{\mu^2}\}$.
- Gluon, charm (label `eq:charm`), given as a new result in the source, with
  $\beta(Q^2)=\sqrt{1+4m_c^2/Q^2}$ and $m_c=1.42$ GeV (pole mass); the bottom
  contribution has $m_c\to m_b$ and an extra factor $1/4$ from $e_b^2/e_c^2$, and
  is said to be negligible for $Q^2\lesssim100$ GeV$^2$.

# Asymptotic value and singlet finite renormalization (label `eq:finrem`)

$$
\lim_{Q^2\to\infty}Q^2F_{\gamma^*\gamma\to M}(Q^2)=\sqrt{\tfrac23}\Bigl[f_M^{(8)}+2\sqrt2\,f_M^{(1)}(\mu_0)\Bigl(1-\frac{2n_f}{\pi\beta_0}\alpha_s(\mu_0)\Bigr)\Bigr].
$$
For $n_f=4$, $\mu_0=1$ GeV, $\alpha_s(1\,{\rm GeV})=0.5$ the source states the
finite renormalization factor is about $0.85$, a $\lesssim5\%$ effect for $\eta$ and
a $20\%$ reduction for $\eta'$ because of the octet-singlet cancellation
(label `eq:asy`: $0.173$ and $0.247$ GeV for $\eta,\eta'$ with the first FKS
parameter set) [@Agaev:2014wna]. The evolution that produces it is in the
[singlet-evolution record](singlet-decay-constant-evolution.md).

# Heavy-flavor statements

- The source adopts the decoupling (fixed-flavor-number) scheme, keeping the full
  heavy-quark mass dependence in the coefficient functions; it states that large
  logarithms $\ln Q^2/m_h^2$ are not resummed (this needs a variable-flavor-number
  scheme with ERBL evolution, left for future work), and terms
  $\alpha_s^2\ln Q/m_h$ in light-quark coefficient functions are omitted.
- Footnote to `eq:finrem`: heavy quarks must be added at their thresholds, so that
  $n_f\alpha_s(\mu_0)\mapsto3\alpha_s(\mu_0)+\alpha_s(\mu_c)+\alpha_s(\mu_b)+\ldots$;
  numerically the difference is stated to be small.

# Flags

- Singly-virtual only; no equal-virtuality or fixed-ratio coefficient is given.
- The gluon coefficient function enters at $\mathcal O(\alpha_s)$ but multiplies the
  singlet decay constant and a gluon DA that is a nonperturbative input.
