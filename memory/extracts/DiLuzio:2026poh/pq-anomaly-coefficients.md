---
title: Color and electromagnetic anomalies of the PQ symmetry
status: mixed
created: 2026-11-19
updated: 2026-11-19
bibtex_keys: [DiLuzio:2026poh, Alda:2026nrd]
provenance: literature
---

In DiLuzio:2026poh, Eqs. (3), (5), and (9), the quark mass of flavor $q$ carries the phase $e^{in_q a}$, where $a$ is the dimensionless phase of $\varphi$ with period $2\pi$. For negative $n_q$, the mass operator uses $\varphi^*$ so that the phase has the opposite sign. The color anomaly is proportional to the signed sum $\sum_q n_q$. For a nonzero sum, the electromagnetic-to-color anomaly ratio and the QCD domain-wall number are
$$
\frac{E}{N}=\frac{2}{3}\frac{4n_u+n_d+n_s}{n_u+n_d+n_s},
\qquad
N_{\mathrm{DW}}=|n_u+n_d+n_s|.
$$
The original assignment $(1,0,0)$ gives $E/N=8/3$ and $N_{\mathrm{DW}}=1$.

These coefficients describe the PQ transformations of the quarks in Eq. (3). In the mass-phase basis used in Eqs. (3) and (6), matching retains the dynamical mass matrix and the QCD singlet anomaly. Removing the mass phases by a local chiral rotation instead produces derivative quark couplings and anomalous gauge terms fixed by the same charges. The ratio in Eq. (5) gives the relative anomaly coefficients in that representation. The physical two-photon amplitude also depends on the neutral-meson components of the mass eigenstate.

## Mapping to the Alda:2026nrd invariants

Identify $x=a/f_a$ in Alda:2026nrd with the $2\pi$-periodic phase of DiLuzio:2026poh. At the radial vacuum, the source multiplying $-\bar q_R M q_L$ in Alda:2026nrd is the conjugate of the mass matrix in DiLuzio:2026poh, so $\widehat M_{qq}(x)=m_q e^{-in_qx}$ in the mass-phase basis. In that basis the axial-current and explicit ALP gauge coefficients vanish, giving $c_{GA}=c_{\gamma A}=0$. From $\widehat P_1=-\{c_q^{AP},m_q\}$ one obtains $c_q^{AP}=\operatorname{diag}(n_u,n_d,n_s)/2$. See the [mass invariant](../../alp/shared/alp-dressed-mass-first-order-invariants.md).

With $c_q^{AP}=c_q^A-c_q^P$, the invariant color combination is
$$
c_{GP}\equiv c_G+\operatorname{tr}(c_q^P)
=c_{GA}-\operatorname{tr}(c_q^{AP}).
$$
The QCD angle in the normalization of Alda:2026nrd has linear slope $c_{GP}$. In the basis with real masses, $c_q^P=0$ and $c_q^A=c_q^{AP}$, so the explicit color and photon coefficients are $c_G=c_{GP}$ and $c_\gamma=c_{\gamma A}-2N_c\operatorname{tr}(Q^2c_q^{AP})$. Thus
$$
N_{\mathrm{DW}}=2|c_{GP}|,\qquad
\frac{E}{N}=\frac{c_{\gamma A}-2N_c\operatorname{tr}(Q^2c_q^{AP})}{c_{GP}}
\quad (c_{GP}\ne 0).
$$
For $(n_u,n_d,n_s)=(1,0,0)$, $c_{GP}=-1/2$ and the numerator is $-4/3$. These values reproduce $N_{\mathrm{DW}}=1$ and $E/N=8/3$. The factor in the domain-wall formula uses the full $2\pi$ period of $x$. A common rescaling of the linear coefficients cancels in the ratio. See the [color invariant](../../alp/shared/alp-dressed-topological-invariant.md) and the [photon invariant](../../alp/photon/alp-dressed-photon-invariant.md).
