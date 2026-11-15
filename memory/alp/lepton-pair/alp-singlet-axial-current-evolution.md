---
title: Scheme, scale, and evolution of the singlet axial current in the hard ALP form factor
status: mixed
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Larin:1993tq, Larin:1994va, Agaev:2014wna, Ludtke:2020moa]
provenance: analytical-derivation
---

# Scope, inputs, and assumptions

The [axial-current matrix element](alp-axial-current-matrix-element.md) is the leading-order (LO)
chiral perturbation theory value
$\langle0|\bar q\gamma^\mu\gamma_5Xq|a(p)\rangle=\frac{2if_\pi^2}{f_a}p^\mu\langle X\,\mathbb M_a^{-1}(\bm{\mathcal C})\rangle$.
This record states in which scheme and at which scale that matrix element is used in the short-distance
expansion, how its singlet part depends on the scale, and what role the anomaly and the gluonic operator play.
Conventions are those of the [shared-conventions record](alp-high-energy-tff-conventions.md);
$q^2$ is the photon virtuality and $Q=\operatorname{diag}(\tfrac23,-\tfrac13,-\tfrac13)$ the charge matrix.

Inputs and provenance:

- Singlet current $\bar\psi\gamma_\mu\gamma_5\psi$ summed over the $n_f$ active flavors, renormalized in the
  $\overline{\rm MS}$ scheme with the finite constant $Z_5^s$, in massless QCD; three-loop anomalous dimension
  $\mu^2\,d\ln J/d\mu^2=-6C_Fn_f\,(\alpha_s/4\pi)^2+\bigl(-\tfrac{142}{3}C_FC_An_f+\tfrac43C_Fn_f^2+18C_F^2n_f\bigr)(\alpha_s/4\pi)^3$
  [@Larin:1993tq, see the [extract](../../extracts/Larin:1993tq/singlet-axial-current-renormalization.md)].
  The non-singlet current has vanishing anomalous dimension in the same scheme.
- Same result at leading order, and its integrated form at order $\alpha_s$:
  $\mu\,d\ln f^{(1)}/d\mu=-n_f\alpha_s^2/\pi^2$, $f^{(1)}(\mu)=f^{(1)}(\mu_0)\{1+\frac{2n_f}{\pi\beta_0}[\alpha_s(\mu)-\alpha_s(\mu_0)]\}$,
  $\beta_0=11-\tfrac23n_f$, evaluated from a low scale $\mu_0\sim1$ GeV
  [@Agaev:2014wna, see the [extract](../../extracts/Agaev:2014wna/singlet-decay-constant-evolution.md)].
- Next-to-leading-order (NLO) correction to the Wilson coefficient of the leading local operator,
  $1-\alpha_s(\mu)/\pi$ with $\mu^2=-q^2$, for the symmetric limit [@Ludtke:2020moa, see the
  [extract](../../extracts/Ludtke:2020moa/ope-nlo-wilson-coefficient.md)].
- The two-loop coefficient $\beta_1=102-\tfrac{38}{3}n_f$ of the beta function in $\mu^2\,d(\alpha_s/4\pi)/d\mu^2=-\beta_0(\alpha_s/4\pi)^2-\beta_1(\alpha_s/4\pi)^3$
  is the standard value; **it is not taken from a source consulted in this repository** and is used only in the next-order series check.

Assumptions (not results): massless quarks with a fixed number $n_f$ of active flavors (thresholds are task 2.3);
the LO ChPT matrix element is identified with the $\overline{\rm MS}$ current at the low scale $\mu_0$; the NLO
coefficient $1-\alpha_s/\pi$ applies to the singlet part as it does to the non-singlet part.

# Singlet and non-singlet parts

For a diagonal weight, with $\langle\cdot\rangle$ the flavor trace and $N_f=3$,
$$
\langle Q^2M\rangle=\Bigl\langle\bigl(Q^2-\tfrac29\mathbf 1\bigr)M\Bigr\rangle+\tfrac29\langle M\rangle ,
\qquad \tfrac29=\frac{\langle Q^2\rangle}{N_f},
$$
where $Q^2-\tfrac29\mathbf 1=\operatorname{diag}(\tfrac29,-\tfrac19,-\tfrac19)$ is traceless. The traceless part
contains the isovector and octet currents, which do not depend on the scale. Only the second term
involves the singlet current $\sum_f\bar q_f\gamma_\mu\gamma_5q_f$ and therefore depends on the scale.
For $M=\mathbb M_a^{-1}(\bm{\mathcal C})$ this is the splitting already shown in the matrix-element record.

# Scale dependence of the singlet matrix element

The paired [script](../../../scripts/derive_alp_singlet_axial_rg.py) integrates the anomalous dimension over the
beta function. At the leading order in each, the one-loop coupling
$$
\frac{\alpha_s(\mu)}{4\pi}=\frac{\alpha_s(\mu_0)/4\pi}{1+\beta_0\frac{\alpha_s(\mu_0)}{4\pi}\ln(\mu^2/\mu_0^2)}
$$
gives for the ratio of the singlet matrix element at $\mu$ to that at $\mu_0$
$$
\boxed{\frac{\langle0|\sum_f\bar q_f\gamma^\mu\gamma_5q_f|a(p)\rangle_\mu}{\langle0|\sum_f\bar q_f\gamma^\mu\gamma_5q_f|a(p)\rangle_{\mu_0}}
=\exp\Bigl[\frac{6C_Fn_f}{\beta_0}\Bigl(\frac{\alpha_s(\mu)}{4\pi}-\frac{\alpha_s(\mu_0)}{4\pi}\Bigr)\Bigr]
\;\simeq\;1+\frac{2n_f}{\pi\beta_0}\bigl[\alpha_s(\mu)-\alpha_s(\mu_0)\bigr],}
$$
with $C_F=\tfrac43$; the factor is unity for the non-singlet part. The right-hand form is the first-order
truncation and agrees with the cited literature expression. The exponential is the solution of
$\mu^2\,d\ln(\text{ratio})/d\mu^2=-6C_Fn_f(\alpha_s/4\pi)^2$ with the one-loop coupling; the script
verifies this identity exactly. At the next order in $\alpha_s$,
$$
\ln(\text{ratio})=\frac{6C_Fn_f}{\beta_0}\Delta-\frac12\Bigl[\frac{\gamma_3}{\beta_0}+\frac{6C_Fn_f\beta_1}{\beta_0^2}\Bigr]\Delta_2+\mathcal O(\alpha_s^3),
$$
where $\gamma_3=-\tfrac{142}{3}C_FC_An_f+\tfrac43C_Fn_f^2+18C_F^2n_f$,
$\Delta=\frac{\alpha_s(\mu)-\alpha_s(\mu_0)}{4\pi}$ and
$\Delta_2=\bigl(\frac{\alpha_s(\mu)}{4\pi}\bigr)^2-\bigl(\frac{\alpha_s(\mu_0)}{4\pi}\bigr)^2$;
the script checks this as a series, with residual starting at $\alpha_s^4$ in the evolution equation. It is not an exact
two-loop solution.

**Size.** For $n_f=3$, $\mu_0=1$ GeV and $\alpha_s(\mu_0)=0.5$ (an illustration of the effect, not a benchmark and not an
ALP input), the leading-order ratio at $\mu=2,10,100$ GeV is $0.965,\,0.936,\,0.922$ with fixed $n_f=3$ and no
thresholds. The next-order term is not small at this $\mu_0$: its relative size in the exponent is
$107\,[\alpha_s(\mu)+\alpha_s(\mu_0)]/(72\pi)$ at $n_f=3$, about $0.4$ at $\mu=2$ GeV. The perturbative expansion of the
evolution therefore converges slowly from $\mu_0\sim1$ GeV, and the choice of $\mu_0$ is an uncertainty to be reported.

# Result for the singlet-resolved operator contribution

Combining the previous record's result, the NLO coefficient with $\mu^2=-q^2$, and this evolution, the operator
contribution to the nearly symmetric hard TFF at order $\alpha_s$ is
$$
F_{a\gamma^*\gamma^*}\supset-\frac{4f_\pi^2}{f_a\,q^2}\Bigl[1-\frac{\alpha_s(\sqrt{-q^2})}{\pi}\Bigr]
\Bigl\{\Bigl\langle\bigl(Q^2-\tfrac29\mathbf 1\bigr)\mathbb M_a^{-1}(\bm{\mathcal C})\Bigr\rangle
+\tfrac29\,\bigl(\text{ratio at }\mu=\sqrt{-q^2}\bigr)\,\bigl\langle\mathbb M_a^{-1}(\bm{\mathcal C})\bigr\rangle\Bigr\}.
$$
This is only the operator contribution of the previous record; the term independent of $q^2$ from the explicit
couplings is task 3.1. The product of the coefficient and the matrix element is independent of $\mu$ at the
order kept only if the singlet coefficient has the corresponding dependence on $\ln(-q^2/\mu^2)$; that dependence is not
given in the sources read, and this is the reason for fixing $\mu^2=-q^2$.

# Anomaly and gluonic operator

- In massless $\overline{\rm MS}$ QCD with the finite constant $Z_5^s$, the renormalized divergence of the singlet current obeys
  $(\partial_\mu J^\mu_5)_R=\frac{\alpha_s}{4\pi}\frac{n_f}2(\epsilon_{\mu\nu\lambda\rho}G^{a\,\mu\nu}G^{a\,\lambda\rho})_R$ [@Larin:1993tq]; with quark masses,
  $\partial^\mu J^{(s)}_{\mu5}=2m_s\bar si\gamma_5s+\frac{\alpha_s}{4\pi}G\widetilde G$ for the strange component [@Agaev:2014wna].
  The renormalized gluonic operator mixes with the divergence of the singlet current, so
  $\langle0|\frac{\alpha_s}{4\pi}G\widetilde G|a\rangle$ is not an independent input: it equals the divergence of the singlet matrix element minus the mass terms.
- The mass-plus-topological form of the [matrix-element record](alp-axial-current-matrix-element.md) splits the divergence into
  $B_0\{m_q,c_q^{AP}+\mathbb M_a^{-1}(\bm{\mathcal C})\}$ and $\frac{m_0^2}{N_f}(c_{GA}+\langle\mathbb M_a^{-1}(\bm{\mathcal C})\rangle)\mathbf 1$.
  Matching the second term to the anomaly operator is consistent with the statement in [@Alda:2026nrd] that for $m_a^2\gg m_q^2$ the mass term is suppressed
  and the anomalous combination $c_{GA}$ survives, but **it is a hypothesis here, not a derivation**: the sign depends on the orientation
  of $\epsilon$ in $G\widetilde G$ (the cited sources use different orientations from this project), and the pseudoscalar-density
  matrix element $\langle0|\bar qi\gamma_5Xq|a\rangle$ that defines the mass part has not been derived.
- The direct coefficient of $\frac{\alpha_s}{4\pi}G\widetilde G$ in the product of two electromagnetic currents is of relative order $\alpha_s$.
  Only a singly virtual (one photon real) coefficient at leading order is available in the literature read [@Agaev:2014wna];
  none for equal or fixed-ratio virtualities. The gluonic contribution is therefore not included, and is
  carried as a named input and an uncertainty estimate of relative order $\alpha_s$.

# Remaining named inputs

1. The low scale $\mu_0$ at which the LO ChPT matrix element represents the $\overline{\rm MS}$ current, and the size of the corrections to that identification.
2. The scheme of the NLO coefficient $1-\alpha_s/\pi$ for the singlet part, and its dependence on $\ln(-q^2/\mu^2)$.
3. The gluonic-operator coefficient for equal and fixed-ratio virtualities, and the pseudoscalar-density matrix element.
4. The number of active flavors and threshold matching of the singlet current (task 2.3); [@Larin:1994va] states that a heavy quark need not decouple from the singlet axial part, which is not resolved here.
5. $\beta_1$ from a consulted source, if the next-order series is used.

# Verification

Run
```text
.venv/bin/python scripts/derive_alp_singlet_axial_rg.py
.venv/bin/python scripts/derive_alp_singlet_axial_rg.py --test
```
Six unit tests pass: exact solution of the differential equation at leading order with the explicit running coupling,
equivalence to the cited linear form and its derivative, the order of the next-order residual,
the trivial non-singlet case, the flavor projection above, and the consistency of the illustration table.
The script prints $\alpha_s/(4\pi)$ as $a(\mu)$ following the cited notation of the source; here it is always written out.
