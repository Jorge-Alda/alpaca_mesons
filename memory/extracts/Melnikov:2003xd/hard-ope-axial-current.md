---
title: Leading OPE of two electromagnetic currents into the axial current, flavor weights, and VVA coefficients
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Melnikov:2003xd]
provenance: literature
---

# Source and scope

All statements are read from `memory/sources/Melnikov:2003xd/src/lbl5.tex`
(arXiv hep-ph/0312226), subsection "OPE and triangle amplitude"
[@Melnikov:2003xd]. Equation labels are the LaTeX labels in that file; the
numbers printed in the journal version were not checked. Order in
$\alpha_s$: leading order only for the Wilson coefficient. The source states no
$\alpha_s$ correction to it.

# Leading OPE (label `jj`)

For $q_1^2\approx q_2^2\gg q_3^2$, with $\hat q=(q_1-q_2)/2\approx q_1\approx-q_2$
and $q_1+q_2=-q_3$ small, the source writes
$$
i\int d^4x\,d^4y\,e^{-iq_1x-iq_2y}\,T\{j_{\mu_1}(x)j_{\mu_2}(y)\}
=\int d^4z\,e^{-i(q_1+q_2)z}\,\frac{2i}{\hat q^2}\,
\epsilon_{\mu_1\mu_2\delta\rho}\hat q^{\delta}j_5^{\rho}(z)+\cdots
$$
[@Melnikov:2003xd]. The axial current is $j_5^\rho=\bar q\,\hat Q^2\gamma^\rho\gamma_5 q$,
so flavors enter with weights $Q_f^2$. Only the leading term in $\Lambda_{\rm QCD}/\hat q$
is kept. The source notes that the matrix element of this OPE between the
pion and the vacuum gives the asymptotic $\pi^0\gamma^*\gamma^*$ amplitude
and attributes this to Nucl. Phys. B 237 (1984) 525 ("misuse" reference in the source)
[@Melnikov:2003xd].

# Flavor decomposition (label `expj`)

$$
j_{5\rho}=\sum_{a=3,8,0}\frac{{\rm Tr}[\lambda_a\hat Q^2]}{{\rm Tr}[\lambda_a^2]}\,j_{5\rho}^{(a)},
\qquad j_{5\rho}^{(a)}=\bar q\lambda_a\gamma_\rho\gamma_5 q,\ \lambda_0=\mathbb 1,
$$
[@Melnikov:2003xd]. The source normalizes $\lambda_0$ to the unit matrix here, not
to $\sqrt{2/3}\,\mathbb 1$.

# Perturbative triangle coefficients (labels `eqw`, `eqLT`)

For the amplitude with one axial current and two vector currents, one of which
has vanishing momentum, the source writes (label `eqw`)
$T^{(a)}_{\mu_3\rho}=-\frac{ieN_c{\rm Tr}[\lambda_a\hat Q^2]}{4\pi^2}\{w_L^{(a)}q_{3\rho}q_3^\sigma\tilde f_{\sigma\mu_3}+w_T^{(a)}(\dots)\}$
and, for massless quarks (label `eqLT`),
$$
w_L^{(a)}(q^2)=2w_T^{(a)}(q^2)=-\frac{2}{q^2}
$$
[@Melnikov:2003xd]. For the non-singlet currents $a=3,8$ the source states that
$w_L$ is exact in the chiral limit (no perturbative or nonperturbative
corrections), and that perturbative corrections to $w_T^{(3,8)}$ are also
absent, citing the non-renormalization theorem (see
[Vainshtein:2002nv](../Vainshtein:2002nv/vva-soft-photon-nonrenormalization.md)).
For the transversal function and for the singlet longitudinal function
$w_L^{(0)}$ higher OPE terms exist; the source says their resummation shifts
the pole $1/q^2\to1/(q^2-m^2)$, with the singlet pole moved to $m_{\eta'}^2$
[@Melnikov:2003xd].

# Pion coupling (labels `gpgg`)

With $\langle0|j^{(3)}_{5\rho}|\pi^0\rangle=2iF_\pi q_\rho$ the source obtains
$g_{\pi\gamma\gamma}=N_c{\rm Tr}[\lambda_3\hat Q^2]/(16\pi^2F_\pi)$
[@Melnikov:2003xd]. Here $F_\pi\simeq92$ MeV because $j_5^{(3)}$ carries $\lambda_3$,
not $\lambda_3/2$.

# Convention flags

- The source does not state the sign of $\epsilon^{0123}$ or the
  $\gamma_5$ orientation in the TeX file (searched for "0123" and
  "gamma_5"; no definition found). The overall sign of `jj` therefore cannot
  be tied to the project convention ($\epsilon^{0123}=+1$) from this source.
- The sign of `jj` ($+2i/\hat q^2$ with $e^{-iq_1x-iq_2y}$) is opposite to the
  one printed in [Ludtke:2020moa](../Ludtke:2020moa/ope-nlo-wilson-coefficient.md)
  ($-2i/\hat q^2$, same Fourier phases, stated $\epsilon^{0123}=+1$). This is
  an unresolved sign mismatch between sources. The project record
  `memory/alp/lepton-pair/alp-hard-current-ope.md` agrees with the Ludtke
  sign once the Fourier phase is reversed ($q\to-q$), see that extract.
  The coefficient magnitude $2/\hat q^2$ agrees in both.
