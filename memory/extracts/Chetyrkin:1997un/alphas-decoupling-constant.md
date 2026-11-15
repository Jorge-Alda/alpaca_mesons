---
title: Decoupling constant of the strong coupling through three loops
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Chetyrkin:1997un, Larin:1994va]
provenance: literature
---

# Source and scope

Read from `memory/sources/Chetyrkin:1997un/src/dec.tex` (arXiv hep-ph/9708255),
Sections "Formalism" and "Decoupling relations" [@Chetyrkin:1997un]. Labels are LaTeX labels.
$\overline{\rm MS}$ scheme, $N_c=3$, one heavy quark $h$ with $n_l$ light flavors,
$n_f=n_l+1$. Order: three loops (relative $(\alpha_s/\pi)^3$ in $\zeta_g^2$).

# Definitions

- Beta function (label `eqrgea`): $\mu^2\frac{d}{d\mu^2}\frac{\alpha_s^{(n_f)}}\pi=-\sum_{N\ge1}\beta^{(n_f)}_{N-1}\bigl(\frac{\alpha_s^{(n_f)}}\pi\bigr)^{N+1}$,
  with $\beta_0^{(n_f)}=\frac14\bigl[11-\frac23n_f\bigr]$, $\beta_1^{(n_f)}=\frac1{16}\bigl[102-\frac{38}3n_f\bigr]$.
- Decoupling (label `eqdecg`): $\alpha_s^{(n_l)}(\mu)=\zeta_g^2\,\alpha_s^{(n_f)}(\mu)$, with the right-hand
  side expressed in the full-theory ($n_f$) parameters. The primed theory is the effective
  $n_l$-flavor theory.

# Result (label `eqzetag`)

With $\ell=\ln(\mu^2/m_h^2)$, $m_h=m_h(\mu)$ the $\overline{\rm MS}$ mass, and
$\alpha_s\equiv\alpha_s^{(n_f)}(\mu)$:
$$
\zeta_g^2=1-\frac{\alpha_s}{\pi}\frac{\ell}{6}
+\Bigl(\frac{\alpha_s}{\pi}\Bigr)^2\Bigl(\frac{11}{72}-\frac{11}{24}\ell+\frac1{36}\ell^2\Bigr)
+\Bigl(\frac{\alpha_s}{\pi}\Bigr)^3\Bigl[\frac{564731}{124416}-\frac{82043}{27648}\zeta(3)
-\frac{955}{576}\ell+\frac{53}{576}\ell^2-\frac1{216}\ell^3
+n_l\Bigl(-\frac{2633}{31104}+\frac{67}{576}\ell-\frac1{36}\ell^2\Bigr)\Bigr].
$$
At $\mu=m_h(m_h)$ the source gives numerically
$\zeta_g^2\approx1+0.1528(\alpha_s^{(n_f)}/\pi)^2+(0.9721-0.0847\,n_l)(\alpha_s^{(n_f)}/\pi)^3$.
There is no $\mathcal O(\alpha_s)$ term at this scale: $\alpha_s$ is continuous
across the threshold at one loop and has a finite jump from two loops [@Chetyrkin:1997un].
A pole-mass version (label `eqzetagos`) is also given in the source.

# Cross-reference and discrepancy

The source states that, apart from the $\mathcal O(\alpha_s^3)$ term, its result agrees with
Larin:1994va (registered as
[@Larin:1994va], see [its record](../Larin:1994va/heavy-quark-singlet-axial-non-decoupling.md)),
while the constant in the two-loop term differs slightly from Z. Phys. C 11 (1981) 79, whose
authors later revised their result to agree (statement in the source, citing its reference Ber97). The Z. Phys. C 11 (1981) 79 text itself was not read (no arXiv version; journal
not accessed); its erratum is listed in the INSPIRE record.

# Flags

- This is the decoupling of $\alpha_s$ only. The source gives no decoupling
  relation for the axial current or $G\tilde G$.
