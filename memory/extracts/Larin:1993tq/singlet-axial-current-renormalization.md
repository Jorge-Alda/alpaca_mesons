---
title: Renormalization, finite constant Z5, and anomalous dimension of the singlet axial current
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Larin:1993tq]
provenance: literature
---

# Source and scope

Read from `memory/sources/Larin:1993tq/src/main.tex` (arXiv hep-ph/9302240),
Sections 1 and 3 [@Larin:1993tq]. Labels are LaTeX labels. Massless perturbative QCD,
$\overline{\rm MS}$ scheme, $D=4-2\epsilon$, coupling $a=\alpha_s/(4\pi)=g^2/(16\pi^2)$.
Order: three loops for the anomalous dimension. The two-loop term is stated
to agree with Nucl. Phys. B 165 (1980) 129; the operator mixing literature is also
associated with Nucl. Phys. B 213 (1983) 301 in the source's bibliography.

# Definitions

- Singlet current (label `sa`): $J^5_\mu=\bar\psi\gamma_\mu\gamma_5\psi=\frac{i}{3!}\varepsilon_{\mu\nu_1\nu_2\nu_3}\bar\psi\gamma_{\nu_1}\gamma_{\nu_2}\gamma_{\nu_3}\psi$
  ('t Hooft-Veltman $\gamma_5$ with four-dimensional $\varepsilon$).
- Renormalized current (label `rens`): $(J^5_\mu)_R=Z_5^sZ^s_{\overline{MS}}(J^5_\mu)_B$.
- Anomaly relation kept one-loop in dimensional regularization (label `anom`):
  $(\partial_\mu J^5_\mu)_R=a\frac{n_f}{2}(G\tilde G)_R$ with
  $G\tilde G=\varepsilon_{\mu\nu\lambda\rho}G^a_{\mu\nu}G^a_{\lambda\rho}$, i.e.
  $n_f\frac{\alpha_s}{4\pi}G^a_{\mu\nu}\tilde G^{a}_{\mu\nu}$ with $\tilde G=\frac12\epsilon G$
  (hand conversion, consistent with the per-flavor relation in
  [Agaev:2014wna](../Agaev:2014wna/singlet-decay-constant-evolution.md))
  [@Larin:1993tq].

# Results

Finite constants (labels `z5`, `z5s`):
$$
Z_5^{ns}=1-4C_Fa+a^2\bigl(22C_F^2-\tfrac{107}{9}C_FC_A+\tfrac29C_Fn_f\bigr),
$$
$$
Z_5^{s}=1-4C_Fa+a^2\bigl(22C_F^2-\tfrac{107}{9}C_AC_F+\tfrac{31}{18}C_Fn_f\bigr).
$$
They differ only in the $C_Fn_f$ term, from light-by-light-type diagrams.
Pole part (label `zas`): $Z^s_{\overline{MS}}=1+a^2\bigl[\tfrac{22}{3\epsilon}C_FC_A+\tfrac{5}{3\epsilon}C_Fn_f\bigr]+\mathcal O(a^3)$.
Anomalous dimension (label `gammaA`), $\gamma^s_J=\mu^2\,d\ln(Z_5^sZ^s_{\overline{MS}})/d\mu^2$:
$$
\gamma_J^s(a)=-6C_Fn_f\,a^2+\Bigl(-\tfrac{142}{3}C_FC_An_f+\tfrac43C_Fn_f^2+18C_F^2n_f\Bigr)a^3 .
$$
For QCD ($C_F=4/3$): leading term $-8n_fa^2=-\frac{n_f}{2\pi^2}\alpha_s^2$ in
$\mu^2d/d\mu^2$, i.e. $-\frac{n_f}{\pi^2}\alpha_s^2$ in $\mu\,d/d\mu$ (hand
conversion), matching the equation quoted in
[Agaev:2014wna](../Agaev:2014wna/singlet-decay-constant-evolution.md). The
non-singlet current has $\gamma_J^{ns}=0$ after $Z_5^{ns}$ is applied (label `zero`
region) [@Larin:1993tq].
The source states that the finite constants hold also for massive quarks
since the $\overline{\rm MS}$ constants do not depend on masses
(last paragraph of the paper).

# Flags

- Here $n_f$ counts the active (massless) flavors. The paper gives no
  threshold matching for the singlet current; see
  [Larin:1994va](../Larin:1994va/heavy-quark-singlet-axial-non-decoupling.md).
- The paper concerns current renormalization only; no sign or coefficient of an
  electromagnetic-current OPE follows from it.
