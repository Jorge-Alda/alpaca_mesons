---
title: Matching of the strong coupling at flavor thresholds and four-loop running
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Chetyrkin:1997sg, Larin:1994va]
provenance: literature
---

# Source and scope

Read from `memory/sources/Chetyrkin:1997sg/src/a1.tex` (arXiv hep-ph/9706430)
[@Chetyrkin:1997sg]. Labels are LaTeX labels. $\overline{\rm MS}$ scheme, couplant
$a^{(n_f)}(\mu)=\alpha_s^{(n_f)}(\mu)/\pi$. Order: three-loop matching, four-loop running.

# Equations

- Matching at scale $\mu^{(n_f)}$ in units of $\mu_h=m_h(\mu_h)$, with
  $\ell_h=\ln[(\mu^{(n_f)})^2/\mu_h^2]$ and $a'=a^{(n_l)}$, $a=a^{(n_f)}$ (label `msb`, constants in `cms`):
  $$
  \frac{a'}{a}=1-a\frac{\ell_h}{6}+a^2\Bigl(\frac{\ell_h^2}{36}-\frac{19}{24}\ell_h+c_2\Bigr)
  +a^3\Bigl[-\frac{\ell_h^3}{216}-\frac{131}{576}\ell_h^2+\frac{\ell_h}{1728}(-6793+281n_l)+c_3\Bigr],
  $$
  $c_2=\frac{11}{72}$, $c_3=-\frac{82043}{27648}\zeta(3)+\frac{564731}{124416}-\frac{2633}{31104}n_l$.
- The choice $\mu^{(n_f)}=\mu_h$ eliminates $\ell_h$ (the source calls it convenient and
  consistent with the Particle Data Group practice) and gives the relation of the
  $\Lambda$ parameters in label `lam`. A pole-mass form is label `oms`.
- Four-loop running solution (label `alp`), $L=\ln(\mu^2/\Lambda^2)$, $b_N=\beta_N/\beta_0$:
  $a=\frac1{\beta_0L}-\frac{b_1\ln L}{(\beta_0L)^2}+\frac1{(\beta_0L)^3}[b_1^2(\ln^2L-\ln L-1)+b_2]+\frac1{(\beta_0L)^4}[\ldots]$.
- Statement about arbitrariness: "the actual value of $\mu^{(n_f)}$ does not matter as
  long as it is comparable to the heavy-quark mass" in higher orders (text before label `lam`).

# Cross-reference

The source states that it settles a disagreement between two earlier two-loop
results, with its $c_2$ agreeing with Larin's ([@Larin:1994va]) and not with
Z. Phys. C 11 (1981) 79 (that paper was not read; no arXiv version).

# Flags

- Matching is for $\alpha_s$ only. It is given for matching scales near
  $m_h$; no smoothing prescription across thresholds for observables appears in the
  passages read.
- With $c_2\neq0$ the coupling is discontinuous at $\mu=\mu_h$ from two loops
  ($a'/a=1+0.1528a^2$), a fact needed for a continuity check of any coefficient built
  from $\alpha_s(-q^2)$ across $c$ and $b$ thresholds (interpretation).
