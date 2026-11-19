---
title: Vanishing amplitude without radial mode for three distinct neutral daughters with one PQ-charged quark
status: derived
created: 2026-11-19
updated: 2026-11-19
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Why $\mathcal M_{\text{no }\sigma}(a\to\pi^0\eta\eta')=0$ in the only-up model

**Result.** In polar coordinates and at tree level, take a source in which a single quark flavour carries PQ charge, $\widehat M=\operatorname{diag}[m_q(1+\sigma/f_a)^{n_q}e^{-in_qa/f_a}]$ with exactly one $n_q\neq0$, $\theta=0$ and no bare ALP mass. Then the amplitude without the $\sigma$ line vanishes when the parent and the three daughters are the four distinct poles of the diagonal neutral sector ($a$, $\pi^0$, $\eta$, $\eta'$). The cubic vertices vanish by intrinsic parity ([parity record](murayama-cubic-vertex-parity-rule.md)). The quartic contact is proportional to
$$
\sum_{q=u,d,s}\frac{1}{\prod_{p\neq q}\big[2B_0m_q\big(1+\tfrac{n_q^2f_\pi^2}{2f_a^2}\big)-2B_0m_p\big(1+\tfrac{n_p^2f_\pi^2}{2f_a^2}\big)\big]}=0,
$$
a three-term identity in the squared masses $2B_0m_u(1+n_u^2f_\pi^2/(2f_a^2))$, $2B_0m_d$, $2B_0m_s$ (written for charged $u$; for another charged flavour the factor $1+n_q^2f_\pi^2/(2f_a^2)$ moves to that flavour). These are the squared masses of the three cosine directions of the potential without the anomalous singlet term. The $a\to P_iP_jP_k$ amplitude of these final states is then pure $\sigma$ exchange. This settles the "Numerical observation" of the [radial-exchange record](alp-three-pseudoscalar-radial-exchange.md).

The statement concerns the polar coordinates. In Cartesian coordinates the split into "no $\sigma$" and "$\sigma$" parts is different (the Cartesian contact contains $-\tfrac\lambda4p^4$); the total amplitude is the same.

## Conventions and units

Conventions of Alda:2026nrd [@Alda:2026nrd] and of the [radial-exchange record](alp-three-pseudoscalar-radial-exchange.md); $\mathbb M_s$ (`eq:murayama-operator`), $N(m^2)$ (`eq:murayama-pole-norm`), pole legs `eq:murayama-pole-direction`. Model of DiLuzio:2026poh [@DiLuzio:2026poh]. GeV, GeV$^2$. Vertices omit the factor $i$.

## Assumptions

1. Tree level, LO $U(3)$, polar connection-free coordinates, $\sigma=0$ part of the amplitude.
2. Diagonal sources with PQ charges $n_q$, exactly one nonzero; the anomalous singlet term $-\tfrac{m_1^2}{12}\langle\Phi\rangle^2$ is the only $\eta_0$ mass term and does not depend on $a$.
3. Parent and daughters are four distinct, nondegenerate poles of the diagonal neutral sector; the three squared masses above are distinct and $\langle\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\rangle\neq0$.

## Inputs

- Absence of $\mathcal V^{(3)}$: [parity record](murayama-cubic-vertex-parity-rule.md).
- Quartic vertex on diagonal legs, form (4a): [quartic-vertex record](murayama-unmixed-basis-quartic-vertex.md).
- Pole vectors in the original fields: meson matrix $\tfrac{2B_0f_\pi}{f_a}\mathbb M^{-1}_{m_x^2}(\widehat P_{01})/\sqrt{N(m_x^2)}$ and ALP component $1/\sqrt{N(m_x^2)}$, `eq:murayama-physical-fields`.

## Step 1: the diagonal neutral sector

For diagonal sources the diagonal mesons $\Phi_{qq}$ and $a$ do not mix with the off-diagonal mesons at quadratic order, and for diagonal $U$ the meson kinetic term has no quartic part. At $\sigma=0$ the Lagrangian of the sector is
$$
\mathcal L=\frac14\sum_q(\partial\Phi_{qq})^2+\frac12(\partial a)^2+B_0f_\pi^2\sum_qm_q\cos\Big(\frac{\Phi_{qq}}{f_\pi}+\frac{n_qa}{f_a}\Big)-\frac{m_1^2}{12}\Big(\sum_q\Phi_{qq}\Big)^2 .
$$
The $\sigma=0$ amplitude for four poles of this sector is the quartic contact (form (4a) of the quartic-vertex record),
$$
\mathcal V^{(4)}_{aijk}=B_0f_\pi^2\sum_qm_q\prod_{x=a,i,j,k}\frac{1}{\sqrt{N(m_x^2)}}\Big[\frac{2B_0}{f_a}\big(\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\big)_{qq}+\frac{n_q}{f_a}\Big],
$$
the bracket being the cosine argument of flavour $q$ on pole $x$ (up to $1/\sqrt{N(m_x^2)}$).

## Step 2: the cosine argument on a pole

Take the linear equations of motion $\partial\mathcal L_2/\partial\Phi_{qq}$ and $\partial\mathcal L_2/\partial a$ at momentum squared $p^2$. SymPy shows, for symbolic $n_u$ and $n_d=n_s=0$, the identity
$$
\frac{2}{f_\pi}\frac{\partial\mathcal L_2}{\partial\Phi_{qq}}+\frac{n_q}{f_a}\frac{\partial\mathcal L_2}{\partial a}=\Big[p^2-2B_0m_q\Big(1+\frac{n_q^2f_\pi^2}{2f_a^2}\Big)\Big]\Big(\frac{\Phi_{qq}}{f_\pi}+\frac{n_qa}{f_a}\Big)-\frac{m_1^2}{3f_\pi}\sum_p\Phi_{pp}\qquad(q=u,d,s).
$$
It uses that only one flavour enters the ALP equation. On a pole ($p^2=m_x^2$, both equations vanish), with $\Phi_{qq}=\tfrac{2B_0f_\pi}{f_a}\big(\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\big)_{qq}\,a$ (the pole solution of the meson equations),
$$
\Big[2B_0\big(\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\big)_{qq}+n_q\Big]\Big[m_x^2-2B_0m_q\Big(1+\frac{n_q^2f_\pi^2}{2f_a^2}\Big)\Big]=\frac{2B_0m_1^2}{3}\big\langle\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\big\rangle .
$$
Inserted in Step 1 this gives form (4b) of the quartic-vertex record,
$$
\mathcal V^{(4)}_{aijk}=\frac{B_0f_\pi^2}{f_a^4}\Big(\frac{2B_0m_1^2}{3}\Big)^4\prod_{x}\frac{\big\langle\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\big\rangle}{\sqrt{N(m_x^2)}}\;\sum_{q=u,d,s}\frac{m_q}{\prod_{x}\big[m_x^2-2B_0m_q\big(1+\tfrac{n_q^2f_\pi^2}{2f_a^2}\big)\big]} ,
$$
valid also for repeated daughters.

## Step 3: secular product and the three-term identity

Let $\{a,i,j,k\}$ be the four distinct poles of the sector, so that $\prod_x(p^2-m_x^2)$ is its characteristic polynomial in $p^2$. SymPy (determinant of the $4\times4$ inverse propagator with kinetic matrix $\operatorname{diag}(\tfrac12,\tfrac12,\tfrac12,1)$, symbolic $n_u$) gives, for $q=u,d,s$,
$$
\prod_{x=a,i,j,k}\Big[2B_0m_q\Big(1+\frac{n_q^2f_\pi^2}{2f_a^2}\Big)-m_x^2\Big]=-\frac{2B_0m_1^2}{3}\,m_q\prod_{p\neq q}\Big[2B_0m_q\Big(1+\frac{n_q^2f_\pi^2}{2f_a^2}\Big)-2B_0m_p\Big(1+\frac{n_p^2f_\pi^2}{2f_a^2}\Big)\Big].
$$
The massless combination of $\Phi_{uu}$ and $a$ orthogonal to the up-quark cosine direction contributes a factor equal to the squared mass itself, which cancels against the factor $2B_0m_q(1+\dots)$ of the general rank-one formula. Hence
$$
\sum_q\frac{m_q}{\prod_x\big[m_x^2-2B_0m_q(1+\tfrac{n_q^2f_\pi^2}{2f_a^2})\big]}=-\frac{3}{2B_0m_1^2}\sum_{q=u,d,s}\frac{1}{\prod_{p\neq q}\big[2B_0m_q(1+\tfrac{n_q^2f_\pi^2}{2f_a^2})-2B_0m_p(1+\tfrac{n_p^2f_\pi^2}{2f_a^2})\big]}=0,
$$
by the three-term identity (SymPy; it is the second divided difference of the constant $1$ over the three squared masses). In the earlier form of this argument the sum ran over four squared masses including $0$ for the massless combination, with numerator equal to the squared mass; the term with $0$ vanishes and the factor cancels, giving the three-term identity.

## Predictions and their tests

| case | mechanism | prediction | numerical result (relative to the $\sigma$ term) |
|---|---|---|---|
| single charged quark: $n=(1,0,0)$, $(0,1,0)$, $(0,0,1)$, $(2,0,0)$, $(-1,0,0)$ | as above | $0$ | $1.1\times10^{-13}$, $1.5\times10^{-13}$, $4\times10^{-16}$, $9\times10^{-14}$, $1.6\times10^{-13}$ |
| two charged quarks: $n=(1,1,0)$, $(1,0,1)$ | the ALP equation contains two flavours; the identity of Step 2 fails | $\neq0$ | $3.5\times10^{-2}$, $1.6\times10^{-2}$ |
| bare ALP mass $\mu^2$ (cosine model) | the ALP equation acquires $-\mu^2a$, which does not combine with a cosine argument | $\neq0$ | $\mu^2=10^{-3}$, $0.05$ GeV$^2$: $4.5\times10^{-3}$, $-0.116$ (absolute, against $-0.120$ for identical daughters) |
| two identical daughters ($\pi^0\pi^0\eta$, $3\pi^0$) | the product over $x$ is not the characteristic polynomial | $\neq0$ | $\mathcal V^{(4)}=0.164$ and $0.0155$ |

The cosine model of Step 1 reproduces the direct polar calculation of the full eleven-field Lagrangian to $\le1.4\times10^{-15}$ (absolute) for all seven charge assignments; its pole masses agree to $10^{-10}$. Only-up values: the three squared masses are $0.01175$, $0.02444$, $0.494$ GeV$^2$; the four poles are $1.97\times10^{-4}$, $0.0179$, $0.260$, $1.002$ GeV$^2$.

## Verification

[Script](../../../scripts/derive_murayama_only_up_three_pseudoscalar_amplitude.py) (`derive_neutral_pole_relation`, `derive_neutral_secular_product`, `derive_three_term_identity`, `neutral_sector_model`, `single_flavour_neutral_quartic_vertex`) and [unit tests](../../../scripts/test_derive_murayama_only_up_three_pseudoscalar_amplitude.py) (`SymbolicTest`, `ContactVanishingTest`, `DraftDividedDifferenceTest`). Parameters as in the [parity record](murayama-cubic-vertex-parity-rule.md), Dalitz point $s_{12}=0.3$, $s_{13}=0.5$ GeV$^2$.

| check | residual |
|---|---|
| SymPy: identity of Step 2 for $q=u,d,s$; secular product of Step 3; three-term identity and the contact sum | exact |
| Step 2 relation on the four neutral poles of the full model, $n=(1,0,0)$, $(2,0,0)$, $(0,0,1)$, $(-1,0,0)$ (relative) | $\le2.9\times10^{-13}$ |
| Step 3 product on the same poles (relative) | $\le2.2\times10^{-14}$ |
| three-term sum relative to its largest term | $\le6.4\times10^{-17}$ |
| form (4b) on four distinct poles relative to a repeated-daughter value | $\le1.2\times10^{-14}$ |

```sh
.venv/bin/python scripts/derive_murayama_only_up_three_pseudoscalar_amplitude.py
.venv/bin/python -m unittest scripts.test_derive_murayama_only_up_three_pseudoscalar_amplitude
```
