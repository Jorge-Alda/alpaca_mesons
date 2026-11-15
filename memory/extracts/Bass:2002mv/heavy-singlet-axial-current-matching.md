---
title: Heavy singlet axial-current matching onto light currents
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Bass:2002mv]
provenance: literature
---

# Source and conventions

Read Eq. (f4) and the following discussion in `memory/sources/Bass:2002mv/src/main.tex` (arXiv hep-ph/0207071) [@Bass:2002mv]. The paper uses $F=f+1$- and $f$-flavor $overline{\mathrm{MS}}$ theories, scale $\bar\mu$, heavy-quark mass $m_h$, and $\alpha_F$ for the full-theory coupling. The displayed matching is at leading power in $1/m_h$; corrections indicated by the source are $O(\alpha_F^3,m_h^{-1})$.

# Operator relation

Matching the heavy-flavor axial current amplitude $\Gamma_{\mu5}$ onto an external light-quark axial current gives
$$
\Gamma_{\mu5}=\left(\frac{\alpha_F}{\pi}\right)^2\left(\ln\frac{m_h}{\bar\mu}+\frac18\right)\gamma_\mu\gamma_5+O(\alpha_F^3,m_h^{-1}).
$$
The paper states the corresponding operator contribution as
$$
h_F=\left(\frac{\alpha_F}{\pi}\right)^2\left(\ln\frac{m_h}{\bar\mu}+\frac18\right)\sum_{i=1}^{f}(\bar q_i\gamma_\mu\gamma_5q_i)_f+O(\alpha_F^3,m_h^{-1}),
$$
with the sum over the light flavors. At $\bar\mu=m_h$, the coefficient is $(\alpha_s/\pi)^2/8$. This is the leading nonzero perturbative matching contribution, at $O(\alpha_s^2)$.

# Scope and limitations

The paper uses this relation inside Witten's renormalization-group treatment of heavy-quark contributions to the nucleon's neutral weak axial charge. The operator matching itself is state-independent, but its application here requires a consistent map between the paper's $\overline{\mathrm{MS}}$ currents, the charge-weighted ALP matrix element, and the chosen ALP source basis. It is not an $O(\alpha_s)$ correction to the symmetric electromagnetic-current Wilson coefficient. It also does not provide the symmetric doubly virtual gluonic Wilson coefficient or the source-dependent light-quark mass corrections to the ALP TFF.
