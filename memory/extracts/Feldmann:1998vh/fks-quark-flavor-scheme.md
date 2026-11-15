---
title: Quark-flavor basis mixing scheme for eta and eta-prime decay constants
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Feldmann:1998vh]
provenance: literature
---

# Source and scope

Read from `memory/sources/Feldmann:1998vh/src/hep_rev.tex` (arXiv hep-ph/9802409),
introduction and theory sections [@Feldmann:1998vh]. Labels are LaTeX labels. This is a model
assumption about the decay constants, not a QCD result, and the source ignores the
scale dependence of the singlet current.

# Statements

- General octet-singlet parametrization (label `newmix`): $f^8_\eta=f_8\cos\theta_8$,
  $f^1_\eta=-f_1\sin\theta_1$, $f^8_{\eta'}=f_8\sin\theta_8$, $f^1_{\eta'}=f_1\cos\theta_1$,
  with $\theta_1\ne\theta_8$.
- Quark-flavor basis (label `phimix`), central assumption: the decay constants follow
  the state mixing with one angle $\phi$:
  $f^q_\eta=f_q\cos\phi$, $f^s_\eta=-f_s\sin\phi$, $f^q_{\eta'}=f_q\sin\phi$,
  $f^s_{\eta'}=f_s\cos\phi$, with currents $\bar qq=(\bar uu+\bar dd)/\sqrt2$ and $\bar ss$.
  The source says it has "no rigorous justification" for this assumption
  and tests it [@Feldmann:1998vh].
- First-order estimate (label `phen6`): $f_q=f_\pi$, $f_s=\sqrt{2f_K^2-f_\pi^2}$;
  Leutwyler's relation $f_8=\sqrt{\tfrac43f_K^2-\tfrac13f_\pi^2}$ (label `f8rel`).
- Anomaly relation for the $s$-current divergence, in the source normalization:
  $\partial^\mu J^s_{\mu5}=2m_s\bar si\gamma_5s+\frac{\alpha_s}{4\pi}G\tilde G$.
- Phenomenological determination quoted in the source: $\phi=39.9^\circ\pm2.9^\circ$
  from $J/\psi$ decays, and the consistency values $f_q=1.09f_\pi$, $f_s=1.38f_\pi$,
  $\phi_q=39.4^\circ$, $\phi_s=38.5^\circ$ from the earlier four-parameter analysis.

# Flags

- Scheme/scale dependence: the source does not give a renormalization scale for
  $f_q,f_s$. [Agaev:2014wna](../Agaev:2014wna/singlet-decay-constant-evolution.md)
  assumes the scheme at $\mu_0\approx1$ GeV and evolves the singlet part.
- Relevant to meson TFFs only. An ALP has no such mixing parameters;
  its octet and singlet components are fixed by the ALP source structure.
