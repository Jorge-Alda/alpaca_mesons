---
title: Scale dependence of the singlet decay constant from the U(1) anomaly and the FKS scale assumption
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Agaev:2014wna, Feldmann:1998vh]
provenance: literature
---

# Source and scope

Read from `memory/sources/Agaev:2014wna/src/140914.tex` [@Agaev:2014wna]. Labels are LaTeX
labels. The leading-order anomalous dimension is attributed by the source to
Nucl. Phys. B 165 (1980) 129. Order: $\mathcal O(\alpha_s^2)$ in the anomalous dimension.

# Equations

- Octet decay constant $f^{(8)}_M$ is scale independent; the singlet one evolves
  (before label `eq:anomaly`):
  $$
  \mu\frac{d}{d\mu}f^{(1)}_M(\mu)=-4n_f\Bigl(\frac{\alpha_s}{2\pi}\Bigr)^2f^{(1)}_M+\mathcal O(\alpha_s^3),
  $$
  equivalently $\mu\,d\ln f^{(1)}/d\mu=-n_f\alpha_s^2/\pi^2$.
- Solution at this order (label `eq:anomaly`), with $\beta_0=11-\tfrac23n_f$:
  $$
  f^{(1)}_M(\mu)=f^{(1)}_M(\mu_0)\Bigl\{1+\frac{2n_f}{\pi\beta_0}\bigl[\alpha_s(\mu)-\alpha_s(\mu_0)\bigr]\Bigr\},
  $$
  where $n_f$ is the number of light flavors [@Agaev:2014wna].
- Anomaly relation (label `anomaly`) with $a_M=\langle0|\frac{\alpha_s}{4\pi}G^A_{\mu\nu}\tilde G^{A,\mu\nu}|M\rangle$ (label `eq:aM`):
  $\partial^\mu J^{(s)}_{\mu5}=2m_s\bar si\gamma_5s+\frac{\alpha_s}{4\pi}G\tilde G$,
  and $H^{(s)}_M=m_M^2F^{(s)}_M-a_M$.
- FKS scale assumption: because octet and singlet couplings have different scale
  dependence, the relation of the FKS quark-flavor scheme ([Feldmann:1998vh](../Feldmann:1998vh/fks-quark-flavor-scheme.md)) cannot hold at all
  scales. The source assumes it at a low scale $\mu_0\sim1$ GeV and evolves
  upward, which generates OZI-violating contributions; it uses
  $\alpha_s(\mu_0)=0.494$ (text before the numerics) [@Agaev:2014wna].

# Consistency (hand check of the stated solution)

Using $\mu\,d\alpha_s/d\mu=-\beta_0\alpha_s^2/(2\pi)$, the differential equation
above integrates to the displayed solution at this order; this is
interpretation, not a result stated by the source.

# Flags

- Scheme: the source does not state the finite-renormalization scheme of the
  singlet current; see
  [Larin:1993tq](../Larin:1993tq/singlet-axial-current-renormalization.md) for the
  $\overline{\rm MS}$ definition and the three-loop anomalous dimension.
- Applies to meson decay constants, not to an ALP; for an ALP the same
  current renormalization applies but the matrix elements come from the ALP
  sources.
