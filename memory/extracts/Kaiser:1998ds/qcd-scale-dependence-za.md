---
title: QCD-scale dependence of the singlet couplings (Z_A scaling laws)
status: literature
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Kaiser:1998ds]
provenance: literature
---

# QCD-scale dependence of the singlet couplings

Source: `memory/sources/Kaiser:1998ds/src/adehep.tex` (hep-ph/9806336), Sections "Anomalous dimensions" (line 481) and "Scaling laws for the effective coupling constants" (line 596). [@Kaiser:1998ds]

## Origin

The effective action depends on the QCD running scale $\mu$ (lines 490--493). The singlet axial current has an anomalous dimension $\gamma_A$ (lines 583--588),

$$
\mu\frac{dZ_A}{d\mu}=\gamma_AZ_A,\qquad
\gamma_A=-\frac{3N_f(N_c^2-1)}{8N_c}\left(\frac g{2\pi}\right)^4+O(g^6),
$$

and $F^0_P\to Z_AF^0_P$ (line 581). The effective action is invariant under a change of $\mu$ provided the external fields transform as (lines 574--577)

$$
m_\theta\to Z_M^{-1}m_\theta,\qquad D_\mu\theta\to Z_A^{-1}D_\mu\theta,
$$

where $Z_M$ is the quark-mass factor. Because the triangle graph responsible for $\gamma_A$ is $1/N_c$ suppressed, $Z_A=1+O(1/N_c)$ (lines 590--594).

## Scaling laws

With $\mathcal L_{\rm eff}=\frac1{12}F^2(1+\Lambda_1)D_\mu\psi D^\mu\psi-\frac12\tau\psi^2+\dots$ (Eq. (Lpsi), lines 609--611), $F^0_{\eta'}=\sqrt{1+\Lambda_1}F$ and $M_{\eta'}^2=6\tau/[(1+\Lambda_1)F^2]$ at tree level (lines 613--616), the paper finds (lines 600, 621--630):

- $F$ is scale independent; $B\to Z_MB$.
- $\tau\to Z_A^2\tau$ and $1+\Lambda_1\to Z_A^2(1+\Lambda_1)$.
- $L_4,\dots,L_8$ are scale independent (they are Taylor coefficients of $F_\pi$, $F_K$, $M_\pi^2$, $M_K^2$ in the quark masses; lines 624--626).
- $$1+\Lambda_2\to Z_A(1+\Lambda_2),\quad 2L_5+3L_{18}\to Z_A\,(2L_5+3L_{18}),\quad 2L_8-3L_{25}\to Z_A\,(2L_8-3L_{25}).$$

So the scaling is multiplicative on the whole combinations $1+\Lambda_1$, $1+\Lambda_2$, $2L_5+3L_{18}$ and $2L_8-3L_{25}$, not on $\Lambda_i$ alone. Since $L_5$ and $L_8$ are scale independent, the scale dependence of $L_{18}$ and $L_{25}$ is $3L_{18}\to Z_A(2L_5+3L_{18})-2L_5$ and $3L_{25}\to2L_8-Z_A(2L_8-3L_{25})$ (consequence of the quoted laws).

Check stated by the paper: the laws for $\tau$ and $\Lambda_2$ are consistent with the scale independence of $L_7^{SU(3)}=-F^4(1+\Lambda_2)^2/(288\tau)+L_7$ (Eq. (L7), lines 631--634). The vacuum angle enters only via $\tilde\psi=\psi+\theta$, and $\mathcal L_{\rm eff}$ is invariant under $\tilde\psi\to Z_A^{-1}\tilde\psi$ together with the above (lines 652--656). The invariant $L_A=(2L_5^r+3L_{18}^r)/\sqrt{1+\Lambda_1}$ is independent of $\mu$ (line 777--779); $\bar F_0$ in Eq. (F0pi) is as well.

## Relation to the chiral renormalization

The QCD scale $\mu$ and the chiral scale $\mu_\chi$ are independent (lines 676--683). Having no pole in $\lambda$ (e.g. $\Gamma_{25}=0$, $\Lambda_{1,2}$) does not remove the $Z_A$ scaling.
