---
title: NLO Wilson coefficient of the electromagnetic-current OPE into the axial current
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Ludtke:2020moa]
provenance: literature
---

# Source and scope

Read from `memory/sources/Ludtke:2020moa/src/SDC.tex` (arXiv 2006.00007),
subsection "$\alpha_s$ corrections to the OPE" and the preceding OPE
subsection [@Ludtke:2020moa]. Equation labels are the LaTeX labels. The source
quotes, but does not rederive, the NLO coefficient; it attributes it to the
Kodaira papers (Kodaira:1978sh, Kodaira:1979ib, Nucl. Phys. B 165 (1980) 129). None of
these papers was read (no arXiv version; journal not accessed).

# Conventions

The source states $\epsilon^{0123}=+1$. The tensor is
$\Pi^{\mu\nu}(q_1,q_2)=i\int d^4x\,d^4y\,e^{-i(q_1\cdot x+q_2\cdot y)}T\{j^\mu(x)j^\nu(y)\}$,
$\hat q=(q_1-q_2)/2$, with $j^\mu=\bar\psi\hat Q\gamma^\mu\psi$ for $\psi=(u,d,s)$
and $j_5^\mu=\bar\psi\hat Q^2\gamma^\mu\gamma_5\psi$ (label `Eq:EM`, text after
`Eq:OPE`). The TFF is defined by (label `Eq:TFF`)
$i\int d^4x\,e^{iq_1\cdot x}\langle0|T\{j_\mu(x)j_\nu(0)\}|{\rm PS}(q_1+q_2)\rangle=\epsilon_{\mu\nu\alpha\beta}q_1^\alpha q_2^\beta F_{{\rm PS}\gamma^*\gamma^*}(q_1^2,q_2^2)$
[@Ludtke:2020moa]. This tensor/scalar relation has the same form as the project
convention $T^{\mu\nu}=\epsilon^{\mu\nu\rho\sigma}q_{1\rho}q_{2\sigma}F$.

# Equations

Leading order (label `Eq:OPE`):
$$
\Pi^{\mu\nu}(q_1,q_2)=\int d^4z\,e^{-i(q_1+q_2)\cdot z}
\Bigl(-\frac{2i}{\hat q^2}\epsilon^{\mu\nu\alpha\beta}\hat q_\alpha j_{5\beta}(z)\Bigr)+\cdots
$$
The source says no other dimension-3 operator appears, so $\alpha_s$ corrections
change only the coefficient of $j_5$. NLO (label `Eq:OPENLO`):
$$
\Pi^{\mu\nu}(q_1,q_2)=\int d^4z\,e^{-i(q_1+q_2)\cdot z}
\Bigl(-\frac{2i}{\hat q^2}\Bigl(1-\frac{\alpha_s}{\pi}\Bigr)
\epsilon^{\mu\nu\alpha\beta}\hat q_\alpha j_{5\beta}(z)+\mathcal O(\hat q^{-2})\Bigr)
$$
[@Ludtke:2020moa]. Symmetric-limit pion TFF (equation following `Eq:OPENLO`):
$$
F_{\pi\gamma^*\gamma^*}(-Q^2,-Q^2)=\Bigl(1-\frac{\alpha_s}{\pi}+\mathcal O(\alpha_s^2)\Bigr)\frac{2F_\pi}{3Q^2}+\mathcal O(Q^{-4}).
$$
The source states this agrees with the NLO result of Phys. Rev. D 28 (1983) 524
and with the dispersive treatment it cites, and that for $Q_1^2=Q_2^2$ the
NLO hard coefficient is independent of the quark momentum fraction, hence
independent of the pion distribution amplitude (footnote to that equation).

# Order, scale, and scheme statements

- Order given: $\mathcal O(\alpha_s)$ (NLO). No $\alpha_s^2$ term is given.
- Renormalization scale: the dimensionless part is scale independent when
  the coupling is evaluated at $-\hat q^2$:
  $d(-\hat q^2/\mu^2,\alpha_s(\mu^2))=d(1,\alpha_s(-\hat q^2))$. The source
  uses $\mu^2=Q_1^2+Q_2^2$ in its numerics and notes this equals $-2\hat q^2$ up
  to the neglected $Q_3^2$ term [@Ludtke:2020moa].
- Anomalous dimension: non-singlet axial currents are conserved up to quark
  masses, so the anomalous dimension vanishes. The singlet current has an
  anomalous dimension from the gluon anomaly (attributed to
  Nucl. Phys. B 213 (1983) 301 and Bos:1992nd), which the source neglects because it starts
  at $\mathcal O(\alpha_s^2)$. See
  [Larin:1993tq](../Larin:1993tq/singlet-axial-current-renormalization.md) for
  the value.
- The source does not state the scheme or finite renormalization used for the
  singlet axial current at $\mathcal O(\alpha_s)$; for $a=0$ the factor
  $1-\alpha_s/\pi$ is taken to be the same as for $a=3,8$. This assumption is
  the source's and is unverified here.
- Heavy quarks: not discussed; the numerics use a three-flavor one-loop beta
  function matched at $\alpha_s(M_\tau^2)=0.35$ (text after `Eq:Pi1NLOMVOPE`).

# Sign and normalization flags against the project convention

- Project record `memory/alp/lepton-pair/alp-hard-current-ope.md` has
  $+2i/(q^2+i0)\,\epsilon^{\mu\nu\rho\sigma}q_\rho\sum_fQ_f^2\bar q_f\gamma_\sigma\gamma_5q_f$
  for $i\int d^4x\,e^{+iq\cdot x}T\{J(x/2)J(-x/2)\}$. The source has the
  opposite Fourier exponent, $e^{-i(q_1x+q_2y)}$, so $\hat q$ of the source
  corresponds to $-q$ of the project and $-2i/\hat q^2\,\hat q_\alpha\to+2i/q^2\,q_\alpha$.
  Interpretation by hand, not script-verified: the signs agree and the NLO
  statement is a multiplicative factor $(1-\alpha_s/\pi)$ on the project's
  coefficient.
- The source's pion limit $+2F_\pi/(3Q^2)$ is positive; the project's canonical
  fixed-phase pion scalar is negative (see
  [Knecht:2001qf](../Knecht:2001qf/pion-tff-short-distance.md)). The difference
  is the external-meson phase, not the Wilson coefficient. It does not apply to
  an ALP, whose matrix element is defined by its own source.
- Statements on $-Q^2$, Euclidean $\hat q^2<0$: the coefficient $1-\alpha_s/\pi$
  is for spacelike (Euclidean) $\hat q$. No timelike continuation is given.
