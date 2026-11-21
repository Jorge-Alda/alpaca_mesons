---
title: Murayama GeV-scale QCD axion model definition and results
status: literature
created: 2026-11-20
updated: 2026-11-20
bibtex_keys: [Murayama:2026ioh]
provenance: literature
---

All statements are taken from the source of Murayama:2026ioh (arXiv:2601.02465, file `superheavyaxion.tex`). The notation of the paper is kept; conversion to the conventions of Alda:2026nrd is not done here.

## Field content and charges

A complex scalar $\varphi$ with Peccei–Quinn charge $+1$ and the right-handed up quark $u_R$ with charge $-1$. No other field is charged at the QCD scale. In the notation $(n_u,n_d,n_s)$ of the quark charges this is $(1,0,0)$ up to the sign convention for the phase of $\varphi$; the paper does not use this notation. The $d$ and $s$ quarks keep ordinary masses.

## Lagrangian at the QCD scale

$$
\mathcal{L}_{QCD} = -\tfrac12 \mathrm{Tr}\, G_{\mu\nu}G^{\mu\nu} + \frac{\theta}{16\pi^2}\mathrm{Tr}\, G_{\mu\nu}\tilde G^{\mu\nu} - \tfrac14 F_{\mu\nu}F^{\mu\nu}
+ \bar u\, i\slashed{D}\, u + \bar d(i\slashed{D}-m_d)d + \bar s(i\slashed{D}-m_s)s
+ \partial_\mu\varphi^*\partial^\mu\varphi - m_\varphi^2\varphi^*\varphi - (\kappa\,\varphi\,\overline{u_L}u_R + \mathrm{c.c.}).
$$

There is no up-quark mass term; it is forbidden by the PQ symmetry. The symmetry is anomalous, and $\theta$ is removed by $\varphi\to e^{-i\theta}\varphi$, $u_R\to e^{i\theta}u_R$. A quartic term $(\varphi^*\varphi)^2$ is allowed but is mostly ignored.

## Chiral Lagrangian and quark-mass matrix

With $U=e^{2i\pi^aT^a/f_\pi}$,
$$
\mathcal{L}_\chi = \frac{f_\pi^2}{4}\mathrm{Tr}\,\partial_\mu U^\dagger\partial^\mu U + \partial_\mu\varphi^*\partial^\mu\varphi - m_\varphi^2\varphi^*\varphi + \mu^3\,\mathrm{Tr}\,MU + \mathrm{c.c.},
\qquad M=\mathrm{diag}(\kappa\varphi,\, m_d,\, m_s),
$$
with $\mu^3 = f_\pi^2 m_\pi^2/[2(m_u+m_d)]\approx 0.022\ \mathrm{GeV}^3$ (so $\mu^3=\langle\bar u_L u_R\rangle$). The up-quark mass is dynamical and complex through $\varphi$.

## Vacuum, up-quark mass and radial/phase structure

$$
\langle\varphi\rangle=\frac{\kappa\mu^3}{m_\varphi^2},\qquad m_u=\kappa\langle\varphi\rangle=\frac{\kappa^2\mu^3}{m_\varphi^2}.
$$
The paper writes $\varphi=(\sigma+ia)/\sqrt2$ (Cartesian components, not a polar decomposition), giving a degenerate pair $a$, $\sigma$ with mass $\approx m_\varphi$, split only if the quartic term is included. The axion mass does not vanish for $m_u\to0$ or $\kappa\to0$.

Allowing for the Kaplan–Manohar ambiguity, $m_u^{\rm eff}=R\,m_u$ with enhancement factor $R$ of order a few, the required $m_u^{\rm eff}=2.16\pm0.04$ MeV gives $m_\varphi = 3.2\ \mathrm{GeV}\,\kappa R^{1/2}$. Perturbativity of the UV completion needs $\kappa\lesssim0.5$, so $m_\varphi\lesssim 2$ GeV.

## Pseudoscalar mass matrix

In the basis $(a,\pi_3,\eta_8,\eta_1)$ the paper gives a $4\times4$ mass matrix with the mixing entries $a$–$\pi_3$: $\sqrt2\,\kappa\mu^3/f_\pi$, $a$–$\eta_8$: $\sqrt{2/3}\,\kappa\mu^3/f_\pi$, $a$–$\eta_1$: $(2/\sqrt3)\,\kappa\mu^3/f_\pi$, and $a$–$a$: $m_\varphi^2$ (full matrix in Eq. (fullmass2) of the source). The determinant is $\frac{16}{3f_\pi^4}\mu^6 m_d m_s m_1^2 m_\varphi^2$, with $m_1$ the $\eta_1$ mass from the QCD anomaly.

## Stated results

- Pion mass splitting: $\pi^0$–$a$ mixing lowers $m_{\pi^0}^2$ by $2\mu^3 m_u/f_\pi^2$. Requiring the shift not to exceed the observed 4.6 MeV splitting gives $m_\varphi\gtrsim 6.7\,\kappa$ GeV, in tension with the mass relation unless $R>1$. This is called the most significant constraint.
- Decays (rough, untrusted estimates): $\Gamma(a\to u\bar u)=\frac{3}{16\pi}\kappa^2 m_a$; $\Gamma(a\to3\pi^0)=\frac{3}{128\pi^3}\frac{\kappa^2\mu^6}{f_\pi^6}m_a$; $\Gamma(a\to\pi^0\pi^+\pi^-)=\frac{1}{64\pi^3}\frac{\kappa^2\mu^6}{f_\pi^6}m_a$. Other multimeson channels are listed without widths. $\sigma$ decays to an even number of pseudoscalars, $a$ to an odd number.
- $a\to\gamma\gamma$ is suppressed by $(\alpha\kappa m_u/4\pi m_\varphi^2)^2$.
- Nucleon coupling from the Gell-Mann–Okubo fit: $\mathcal{L}_{aNN}=1.47\,\kappa\varphi\,\overline{p_L}p_R+\mathrm{c.c.}$, none to the neutron.
- FCNC are negligible; collider and astrophysical constraints are argued to be evaded. Quantum-gravity corrections ($\varphi^5/M_{\rm Pl}$) shift $\theta$ negligibly.
- UV completion: vector-like quark $\mathcal{Q}=(\mathcal{U},\mathcal{D})$ with $\mathcal{L}\supset\kappa_i\bar Q_i\varphi\mathcal{Q}+y\bar{\mathcal Q}Hu_R$, giving the dimension-five operator $\frac{y\kappa_u}{M_{\mathcal Q}}\bar Q_u\varphi H u_R$ and $\kappa=y v\kappa_u/M_{\mathcal Q}$. The paper suggests $\sigma$, $a$ may be among the observed $f_0$ and $\eta$ resonances.
