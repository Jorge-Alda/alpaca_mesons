---
title: One-loop renormalization of the U(3) chiral couplings (lambda convention and Gamma_i table)
status: literature
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Kaiser:1998ds, Herrera-Siklody:1996tqr, Gasser:1984gg, Cirigliano:2001mk]
provenance: literature
---

# One-loop renormalization of the $U(3)$ chiral couplings

Source: TeX file `memory/sources/Kaiser:1998ds/src/adehep.tex` (arXiv hep-ph/9806336); line numbers below refer to that file. The paper's Section "Loops and chiral logarithms" starts at line 658. [@Kaiser:1998ds]

## Power counting

Up to and including $O(\delta^2)$ (with $\partial_\mu=O(\sqrt\delta)$, $m=O(\delta)$, $1/N_c=O(\delta)$) only one-loop graphs built from $\mathcal L^{(0)}$ (Eq. (L0), lines 452--456) contribute (lines 660--663). The divergences are absorbed in the effective couplings, with the coefficients taken from [@Herrera-Siklody:1996tqr] (lines 665--667).

## Definition of $\lambda$

The renormalization scale of the effective theory is denoted $\mu_\chi$, to distinguish it from the QCD running scale $\mu$ (lines 676--678). The paper defines (lines 674--675)

$$
\lambda=\frac{\mu_\chi^{d-4}}{16\pi^2}\left\{\frac1{d-4}-\frac12\left(\ln4\pi+\Gamma'(1)+1\right)\right\}.
$$

Since $\Gamma'(1)=-\gamma_E$, this equals $\lambda_\infty$ of [@Herrera-Siklody:1996tqr] (Eq. (3.20), `un.tex` lines 1351--1353), which is written with $\log4\pi-\gamma+1$. The same expression, called $\Lambda(\mu)$, is quoted for $SU(3)$ with $L_i=L_i^r(\mu)+\Gamma_i\Lambda(\mu)$ in [@Cirigliano:2001mk] (`source.tex` lines 486--496), who refer to the $SU(3)$ Gasser--Leutwyler paper [@Gasser:1984gg] for the $\Gamma_i$. [@Kaiser:1998ds] cites that paper (Nucl. Phys. B250, 465) as its standard reference "GL SU(3)", not Annals Phys. 158 ([@Gasser:1983yg], not accessible here; it is the $SU(2)$ paper). The statement that $\lambda$ is the Gasser--Leutwyler convention is therefore verified only through [@Cirigliano:2001mk] and the references of [@Kaiser:1998ds], not from the Gasser--Leutwyler papers themselves.

## Couplings without chiral-loop renormalization

The constant $B$ is renormalized (lines 679--686):

$$
B=B^r\left\{1+\frac{4\tau}{F^4}\lambda\right\},\qquad
\mu_\chi\frac{\partial B^r}{\partial\mu_\chi}=-\frac{\tau B^r}{4\pi^2F^4},
$$

with the QCD scale $\mu$ kept fixed. The paper states (lines 687--689): "The remaining couplings occurring at order $p^0$ and $p^2$ ($F$, $\tau$, $\Lambda_1$ and $\Lambda_2$) do not pick up renormalization". This means $F,\tau,\Lambda_1,\Lambda_2$ receive no pole in $\lambda$; $B$ does.

## $O(p^4)$ couplings

At order $p^4$ the renormalization is "of the standard form", but $\Gamma_6$ and $\Gamma_8$ differ from the $SU(3)$ values (lines 689--696):

$$
L_n=L_n^r+\Gamma_n\lambda,
$$

$$
\Gamma_4=\frac18,\quad\Gamma_5=\frac38,\quad\Gamma_6=\frac1{16},\quad\Gamma_7=0,\quad\Gamma_8=\frac3{16},\quad\Gamma_{18}=-\frac14,\quad\Gamma_{25}=0.
$$

The paper does not quote the $SU(3)$ values $11/144$ and $5/48$; these appear in [@Herrera-Siklody:1996tqr] (see [the Herrera-Siklody counterterm extract](../Herrera-Siklody:1996tqr/u3-one-loop-counterterms-and-gamma-table.md)).

## Consistency relations stated

Because $2\Gamma_5+3\Gamma_{18}=0$ and $3\Gamma_4-\Gamma_5=0$ (lines 778--781), the combination $L_A=(2L_5^r+3L_{18}^r)/\sqrt{1+\Lambda_1}$ is independent of $\mu_\chi$, and so is $\bar F_0$ of Eq. (F0pi) (lines 770--781). $F_0$ receives no loop contribution (line 768--769), and neither does $F^8_\eta F^0_\eta+F^8_{\eta'}F^0_{\eta'}$ (lines 822--823).

## Loop expressions retained

Only $F_\pi$ and $F_K$ are given with loops (lines 697--708) with $\mu_P=\frac{M_P^2}{32\pi^2F^2}\ln\frac{M_P^2}{\mu_\chi^2}$ (line 712).

## Renormalization points

No numerical value of $\Lambda_1$, $\Lambda_2$, $L_{18}$ or $L_{25}$ is given at any scale. The only scale statement is an estimate "at the scale $\mu=M_\rho$": resonance exchange (with OZI-respecting couplings) gives $L_{18}=\Lambda_1=0$, so $L_A\simeq2L_5^r(M_\rho)$ and $\bar F_0\simeq F_\pi$ (lines 793--798). The paper warns that the result is sensitive to the scale at which $|L_{18}^r(\mu_\chi)|\ll|L_5^r(\mu_\chi)|$ and $|\Lambda_1|\ll1$ are assumed (lines 798--801). [@Kaiser:1998ds]
