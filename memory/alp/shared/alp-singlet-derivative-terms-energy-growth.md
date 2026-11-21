---
title: Singlet-source terms beyond LO in the dressed frame and energy-growing ALP amplitudes
status: mixed
created: 2026-11-21
updated: 2026-11-21
bibtex_keys: [Kaiser:1998ds, Alda:2026nrd]
provenance: assumption
---

# Singlet-source terms beyond LO in the dressed frame and energy-growing ALP amplitudes

## Conventions, units and assumptions

Units GeV; notation of [@Alda:2026nrd] and [the one-loop power-counting record](alp-loop-power-counting.md) (dressed frame, $x=a/f_a$, $\widehat\Theta(x)=c_{GA}x$). All coefficients below are **unknown low-energy coefficients** (assumption: they are kept as free symbols; no value or large-$N_c$ order is assigned, since neither is recorded in the repository extract).

## Inputs

- The general Lagrangian of [@Kaiser:1998ds] contains $V_3D_\mu\psi D^\mu\psi+V_4D_\mu\psi D^\mu\theta+V_5D_\mu\theta D^\mu\theta$, with $D_\mu\psi=\partial_\mu\psi-2\langle a_\mu\rangle$, $D_\mu\theta=\partial_\mu\theta+2\langle a_\mu\rangle$; the coefficient functions depend on $\widetilde\psi=\psi+\theta$; $V_4$ can be set to zero by a meson-field redefinition and $V_5$ is discarded there ([extract](../../extracts/Kaiser:1998ds/singlet-source-derivative-terms.md)).
- $\psi+\theta=-X_{TC}$ and $D_\mu\psi=-\langle u_\mu\rangle$ ([NLO Lagrangian record](../../nlo/lagrangian.md)); in the dressed frame $X_{TC}=-\operatorname{tr}\Phi/f_\pi-2\widehat\Theta(x)$ ([dressed topological record](alp-dressed-topological-invariant.md)).
- **Name not available.** The extract does not record the name of the coefficient function multiplying $\langle D_\mu U^\dagger D^\mu U\rangle=\langle u_\mu u^\mu\rangle$ in [@Kaiser:1998ds]. It is described here in words: its value at $\widetilde\psi=0$ is $f_\pi^2/4$, and the quantity used below is
$$
\big[\text{ratio}\big]_{uu}\equiv\frac{\text{second }\widetilde\psi\text{-derivative at }0\text{ of the coefficient of }\langle u_\mu u^\mu\rangle}{f_\pi^2/4}.
$$
This bracket is only a written-out abbreviation of that sentence (a name is to be agreed).
- Assumption: the coefficient functions are even in $\widetilde\psi$ (parity), so the expansion starts at $\widetilde\psi^2$.

## Method

1. Dressed frame: the currents are removed, so $\langle a_\mu\rangle\to0$ and $\theta\to2\widehat\Theta(x)$. Hence
$$
\widetilde\psi=\frac{\operatorname{tr}\Phi}{f_\pi}+\frac{2c_{GA}a}{f_a}=\frac{\sqrt6\,\eta_0}{f_\pi}+\frac{2c_{GA}a}{f_a},\qquad
D_\mu\theta\to\frac{2c_{GA}}{f_a}\partial_\mu a,\qquad D_\mu\psi\to\frac{\partial_\mu\operatorname{tr}\Phi}{f_\pi}.
$$
At $c_{GA}=0$ none of these terms couples the ALP.
2. Expand the coefficient of $\langle u_\mu u^\mu\rangle$ to $O(\widetilde\psi^2)$, use $\langle u_\mu u^\mu\rangle=\frac{2}{f_\pi^2}\sum_i(\partial P_i)^2+O(\Phi^4)$, and expand $V_5(\widetilde\psi)(D\theta)^2=[V_5(0)+\frac12V_5''(0)\widetilde\psi^2](D\theta)^2$.
3. Contact amplitudes (all momenta incoming, $-ip$ per derivative) and S-wave projection at $s\gg$ masses. For identical pairs in initial and final states $a_0=\mathcal M/(32\pi)$; the condition is $|{\rm Re}\,a_0|=\frac12$.

## Results

(i) From the $\langle u_\mu u^\mu\rangle$ coefficient, $\mathcal L\supset[\text{ratio}]_{uu}\,c_{GA}^2\,\frac{a^2}{f_a^2}\sum_i(\partial P_i)^2$ (flavour universal over the nine real fields), giving
$$
\mathcal M(P_iP_i\to aa)=-2[\text{ratio}]_{uu}\,c_{GA}^2\,\frac{s-2m_P^2}{f_a^2}
=-2[\text{ratio}]_{uu}\,c_{GA}^2\Big(\frac{f_\pi}{f_a}\Big)^2\frac{s-2m_P^2}{f_\pi^2},
$$
and from the cross term $\mathcal M(P_iP_i\to a\eta_0)=-\sqrt6[\text{ratio}]_{uu}c_{GA}\,(s-2m_P^2)/(f_\pi f_a)$. This term does not renormalize the ALP kinetic term.

(ii) From $V_5$:
$$
K_{aa}=1+\frac{8V_5(0)c_{GA}^2}{f_a^2}=1+\Big(\frac{f_\pi}{f_a}\Big)^2c_{GA}^2\frac{8V_5(0)}{f_\pi^2},
\qquad
\mathcal M(\eta_0\eta_0\to aa)=-\frac{24V_5''(0)c_{GA}^2}{f_\pi^2f_a^2}(s-2m_a^2),
$$
while the $a^2(\partial a)^2$ term gives $\mathcal M(aa\to aa)=64V_5''(0)c_{GA}^4m_a^2/f_a^4$, which does not grow with energy on shell. If $V_4\ne0$ (other field convention) there is an $\eta_0$-ALP kinetic mixing $K_{\eta_0a}=2\sqrt6V_4(0)c_{GA}/(f_\pi f_a)$, to be included in the kinetic blocks.

**Consequence of $K_{aa}$.** After canonical normalization every vertex with $n_a$ ALPs is multiplied by $K_{aa}^{-n_a/2}$, i.e. $f_a^2\to f_a^2+8V_5(0)c_{GA}^2$ in all ALP couplings. If $V_5(0)\neq0$, then for $f_\pi/f_a\to\infty$ at fixed $c_{GA}$ the effective ALP coupling scale tends to $c_{GA}\sqrt{8V_5(0)}$ rather than to zero. The LO power counting of the companion record holds only while $8V_5(0)c_{GA}^2\ll f_a^2$.

**Energy where $|{\rm Re}\,a_0|=1/2$** (high-energy limit; replace $f_a^2\to f_a^2K_{aa}$ if $K_{aa}\neq1$):
$$
s_*(P_iP_i\to aa)=\frac{8\pi f_a^2}{c_{GA}^2\,|[\text{ratio}]_{uu}|}\ \text{(one channel)},\qquad
s_*=\frac{8\pi f_a^2}{3c_{GA}^2\,|[\text{ratio}]_{uu}|}\ \text{(aa coupled to the nine } P_iP_i\text{)},
$$
$$
s_*(\eta_0\eta_0\to aa)=\frac{2\pi f_\pi^2f_a^2}{3c_{GA}^2\,|V_5''(0)|}.
$$
So $\sqrt{s_*}\propto f_a/c_{GA}=f_\pi\,(f_\pi/f_a)^{-1}/c_{GA}$: for $f_\pi/f_a\gg1$ the energy where these terms break perturbative unitarity falls below $4\pi f_\pi$ unless the unknown coefficients are small.

Not computed (stated only): higher terms of the $\widetilde\psi$-dependent potential (e.g. the $\widetilde\psi^4$ term) give an ALP self-coupling proportional to $c_{GA}^4(f_\pi/f_a)^4$, which beyond LO would correlate the topological ALP mass with a self-coupling; the $\widetilde\psi$ dependence of the coefficient of $\langle U^\dagger\chi\rangle$ gives further $c_{GA}$-dependent mass-type vertices.

## Reproducible calculation

`singlet_derivative_terms` in [`derive_alp_loop_power_counting.py`](../../../scripts/derive_alp_loop_power_counting.py) derives all amplitudes, $K_{aa}$, $K_{\eta_0a}$ and $s_*$; the unit test `SingletDerivativeTermTests` in [`test_derive_alp_loop_power_counting.py`](../../../scripts/test_derive_alp_loop_power_counting.py) checks the $P_iP_i\to aa$ and $P_iP_i\to a\eta_0$ amplitudes, the absence of energy growth in $aa\to aa$, the growth of $\eta_0\eta_0\to aa$, $K_{aa}$, and the factor $3$ between the single-channel and nonet-coupled scales.

```sh
.venv/bin/python scripts/derive_alp_loop_power_counting.py
```
