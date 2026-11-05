---
title: Subtraction scheme, chiral scale and QCD scale of the U(3) NLO theory
status: mixed
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Gasser:1983yg, Gasser:1984gg, Kaiser:1998ds, Herrera-Siklody:1996tqr]
provenance: analytical-derivation
---

# Subtraction scheme and scales

## Assumptions and inputs

- The loop content is that of [loop-classes.md](loop-classes.md). Only the G-LO and G-ALP graphs are renormalized here.
- The scheme, the coefficients $\Gamma_i$ and the QCD-scale laws are taken from the sources and checked against the arXiv sources: [Kaiser:1998ds Γ table](../../extracts/Kaiser:1998ds/u3-loop-renormalization-and-gamma-table.md), [Kaiser:1998ds $Z_A$ laws](../../extracts/Kaiser:1998ds/qcd-scale-dependence-za.md), [Herrera-Siklody:1996tqr Γ table](../../extracts/Herrera-Siklody:1996tqr/u3-one-loop-counterterms-and-gamma-table.md) and [notation dictionary](../../extracts/Herrera-Siklody:1996tqr/u3-notation-dictionary-with-kaiser-1998ds.md) [@Kaiser:1998ds; @Herrera-Siklody:1996tqr].
- **Notation (approved by the user):** $L_i^r(\mu)$ is the renormalized local constant, with $\mu$ the chiral subtraction scale. Kaiser:1998ds writes it $\mu_\chi$ to distinguish it from the QCD scale. In the framework with explicit resonances, $L_i^r(\mu)$ is the constant that remains after explicit resonance exchange is removed; see the resonance-exchange record of this folder.

The executable check is [`scripts/derive_nlo_u3_resonance_counting.py`](../../../scripts/derive_nlo_u3_resonance_counting.py):

```bash
.venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --scheme
```

## Step 1: subtraction

Dimensional regularization with the subtraction of Gasser:1983yg [@Gasser:1983yg]:

$$
\lambda(\mu)=\frac{\mu^{d-4}}{16\pi^2}\left[\frac{1}{d-4}-\frac12\left(\log4\pi+\Gamma'(1)+1\right)\right],
\qquad
L_i=L_i^r(\mu)+\Gamma_i\lambda(\mu).
$$

Kaiser:1998ds and Herrera-Siklody:1996tqr (Eq. (3.20), written with $\Gamma'(1)=-\gamma$) use the same $\lambda$. Since $\mu^{d-4}/(d-4)=1/(d-4)+\log\mu+\mathcal O(d-4)$, the bare $L_i$ is $\mu$ independent if

$$
\mu\frac{dL_i^r}{d\mu}=-\frac{\Gamma_i}{16\pi^2},\qquad
L_i^r(\mu_2)=L_i^r(\mu_1)-\frac{\Gamma_i}{16\pi^2}\log\frac{\mu_2}{\mu_1}.
$$

The script checks this symbolically, and numerically for $\mu\in\{0.5,0.77,1.0\}\,$GeV.

## Step 2: the default chiral scale

The default is $\mu=M_\rho=0.77\,$GeV, the usual scale at which resonance-saturation estimates are quoted. The variation used to probe the scale dependence is $\mu\in[0.5,1.0]\,$GeV. Over this range the constants shift as follows (script output):

| Constant | $L_i^r(0.5)-L_i^r(0.77)$ | $L_i^r(1.0)-L_i^r(0.77)$ |
|---|---|---|
| $L_4$ | $+0.34\times10^{-3}$ | $-0.21\times10^{-3}$ |
| $L_5$ | $+1.03\times10^{-3}$ | $-0.62\times10^{-3}$ |
| $L_6$ | $+0.17\times10^{-3}$ | $-0.10\times10^{-3}$ |
| $L_7$ | $0$ | $0$ |
| $L_8$ | $+0.51\times10^{-3}$ | $-0.31\times10^{-3}$ |
| $L_{18}$ | $-0.68\times10^{-3}$ | $+0.41\times10^{-3}$ |
| $L_{25}$ | $0$ | $0$ |

## Step 3: what is renormalized by chiral loops

- The local constants $L_4$–$L_8$, $L_{18}$ and $L_{25}$, with the coefficients in [counterterm-basis.md](counterterm-basis.md).
- $B$: $B=B^r(\mu)\{1+4\tau\lambda(\mu)/F^4\}$, so $\mu\,dB^r/d\mu=-\tau B^r/(4\pi^2F^4)$ (Kaiser:1998ds). Since $4\tau/F^4=m_0^2/(N_fF^2)$, this is the $\eta'$ tadpole. Its relative size is $m_0^2/(16\pi^2F^2)$, which is $\delta^2$, the order of G-LO. So $B^r$ depends on $\mu$ at the declared order. The script asserts the order by three equivalent countings.
- $F$, $\tau$, $\Lambda_1$ and $\Lambda_2$ have no chiral-loop renormalization at this order (Kaiser:1998ds).

## Step 4: combinations independent of the chiral scale

A combination $\sum_ic_iL_i^r$ does not depend on $\mu$ if $\sum_ic_i\Gamma_i=0$. Among $L_4,\dots,L_8,L_{18},L_{25}$, a basis of such combinations is

$$
3L_4-L_5,\quad L_4-2L_6,\quad L_7,\quad 3L_4-2L_8,\quad 2L_4+L_{18},\quad L_{25}.
$$

In particular, $2L_5+3L_{18}$ is independent of $\mu$ ($2\cdot\tfrac38-3\cdot\tfrac14=0$), while $2L_8-3L_{25}$ is not ($\Gamma=3/8$).

## Step 5: the QCD scale

Independently of $\mu$, the singlet axial current has an anomalous dimension, $\mu_{\rm QCD}\,dZ_A/d\mu_{\rm QCD}=\gamma_AZ_A$ with $Z_A=1+\mathcal O(1/N_c)$. Under a change of the QCD scale (Kaiser:1998ds):

| Quantity | Transformation |
|---|---|
| $\tau$ | $Z_A^2\,\tau$ |
| $1+\Lambda_1$ | $Z_A^2(1+\Lambda_1)$ |
| $1+\Lambda_2$ | $Z_A(1+\Lambda_2)$ |
| $2L_5+3L_{18}$ | $Z_A(2L_5+3L_{18})$ |
| $2L_8-3L_{25}$ | $Z_A(2L_8-3L_{25})$ |
| $L_4,\dots,L_8$, $F$ | invariant |
| $B$ | $Z_MB$ (quark-mass factor) |

Hence $L_{18}$ and $L_{25}$ carry the QCD-scale dependence of these combinations. Two combinations are invariant under both scales: $L_A=(2L_5^r+3L_{18}^r)/\sqrt{1+\Lambda_1}$ and $F^4(1+\Lambda_2)^2/\tau$ (Kaiser:1998ds). The script checks both.

**Consistency with the chiral subtraction.** The $Z_A$ laws are stated for the bare couplings. For a combination with a chiral pole ($2L_8-3L_{25}$, and $B$ through $\tau$), the renormalized combination therefore picks up an extra $(Z_A-1)\times$(pole term). This is $\mathcal O(1/N_c)$ times a one-loop pole, that is, $\delta^3$, beyond the declared order. At $\mathcal O(\delta^2)$ the laws hold for the renormalized couplings. The script asserts this. A consistent treatment at $\delta^3$ would also need the next $1/N_c$ terms of $\Gamma_8$ and $\Gamma_{25}$ (see [counterterm-basis.md](counterterm-basis.md)).

**Convention for inputs.** Neither Kaiser:1998ds nor Herrera-Siklody:1996tqr quotes numerical values of $\Lambda_1$, $\Lambda_2$, $L_{18}$ or $L_{25}$ at any QCD scale. Any input for these constants (alp-gii.5) must be quoted together with the QCD scale of its source. Observables must be built only from QCD-scale-invariant combinations, such as $L_A$, or from combinations in which the $Z_A$ factors cancel against those of the ALP couplings to the singlet current. The second point is an assumption to be checked in alp-gii.2 and alp-gii.8.

## Translation between Kaiser:1998ds and Herrera-Siklody:1996tqr

| Kaiser:1998ds | Herrera-Siklody:1996tqr | Status |
|---|---|---|
| $L_4,\dots,L_8$ | $\beta_i$ of $O_4,\dots,O_8$ (Eq. (3.27)) | identical |
| $L_{18}$ | $L_{18}(0)$ | identical |
| $L_{25}$ | $L_{25}'(0)$ | identical up to a factor $i$, which depends on whether the derivative is with respect to $X$ or $\tilde\psi$; open |
| $\Lambda_1$ | $-12W_4(0)/f^2$ | derived by hand in the extract; not script-checked |
| $\Lambda_2$ | $-3v_{31}$ | derived by hand in the extract; not script-checked |
| $V_1,V_2$ | $W_1,W_2$ (leading values $f^2/4$) | matched at leading order only |
| ($L_{47}$ not kept) | $O_{47}=D_\mu\hat\theta\langle C^\mu\rangle\langle M\rangle$ | present only in Herrera-Siklody:1996tqr |

At $n_l=3$ the $\Gamma_i$ of the two papers agree for $i=4,\dots,8$ and $18$. $\Gamma_{25}$ is given only by Kaiser:1998ds.

## Verification

`--scheme` asserts every statement of Steps 1–5 and exits with code 0. The tables are copied from its output.
