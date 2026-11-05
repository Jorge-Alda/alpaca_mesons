---
title: U(3) one-loop counterterm functions B_i(X) and the Gamma_i table at n_l = 3
status: literature
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Herrera-Siklody:1996tqr, Kaiser:1998ds]
provenance: literature
---

# $U(3)$ one-loop counterterms $B_i(X)$ and $\Gamma_i$

Source: `memory/sources/Herrera-Siklody:1996tqr/src/un.tex` (hep-ph/9610549, extracted from `source.tar.gz`); line numbers refer to that file. [@Herrera-Siklody:1996tqr]

## Convention

Heat-kernel, dimensional renormalization for $U_L(n_l)\otimes U_R(n_l)$ with generic $n_l$ (abstract, lines 77--80). The one-loop Lagrangian is (Eq. (3.22), lines 1325--1328)

$$
\mathcal L^{\rm one\text{-}loop}=-W_0^r(X,\mu)+\sum_{i=1}^6W_i^r(X,\mu)E_i+\sum_{i=0}^{57}\beta_i^r(X,\mu)O_i+\text{finite nonlocal},
$$

with $\beta_i(X)=\beta_i^r(X,\mu)+\hbar B_i(X,\mu)\lambda_\infty$ and $\delta W_i=\hbar W_i^{(1)}+\hbar\Omega_i\lambda_\infty$ (Eq. (3.19), lines 1344--1346), and (Eq. (3.20), lines 1351--1353)

$$
\lambda_\infty=\frac{\mu^{D-4}}{(4\pi)^2}\left(\frac1{D-4}-\frac12\left(\log4\pi-\gamma+1\right)\right).
$$

This is the expression used for $\lambda$ in [@Kaiser:1998ds], since $\Gamma'(1)=-\gamma$. The functions are to be read as power series in $X$; renormalizing a function means renormalizing the coefficients of its expansion (lines 1371--1374). Here $X=i\frac{\sqrt{2n_l}}{f}\eta^0+\hat\theta$ with $\hat\theta=i\theta$ (lines 1874, 1190).

The paper notation: $\omega_i=W_i/W_1$, $\omega_i'=W_i'/W_1$, $\omega_i''=W_i''/W_1$; the prime differentiates $W_i$, not $\omega_i$ (Appendix D, lines 2593--2600). $W_1(0)=W_2(0)=f^2/4$.

## Relation between $L_i$ and $\beta_i$

Eq. (3.27) (lines 1426--1437): $L_1=\beta_1+\tfrac12\beta_0$, $L_2=\beta_2+\beta_0$, $L_3=\beta_3-2\beta_0$, $L_{13}=\beta_{13}-4\beta_0$, $L_{14}=\beta_{14}+\beta_0$, $L_{15}=\beta_{15}+2\beta_0$, $L_{16}=\beta_{16}+\frac12\beta_0$, and $L_i=\beta_i$ for all other $i$. So the divergences of $L_4,\dots,L_8,L_{18},L_{25},L_{47}$ are $B_4,\dots,B_8,B_{18},B_{25},B_{47}$.

## $B_i$ with the $U_A(1)$ functions frozen

Disregarding all coefficients associated with $U_A(1)$ (Eq. (llista), lines 1454--1462), for generic $n_l$:

$$
B_4=\tfrac18,\ B_5=\tfrac{n_l}8,\ B_6=\tfrac1{16},\ B_7=0,\ B_8=\tfrac{n_l}{16},\ B_{12}=\tfrac{n_l}8,\ B_{18}=-\tfrac14 ,
$$

(also $B_0=n_l/48$, $B_1=1/16$, $B_2=1/8$, $B_3=n_l/24$, $B_9=n_l/12$, $B_{10}=-n_l/12$, $B_{11}=-n_l/24$, $B_{13}=1/4$, $B_{14}=\dots=B_{17}=0$, $B_{19}=1/24$, $B_{20}=1/12$).

## $\Gamma_i$ at $n_l=3$

Renaming the constants as $\Gamma_i$ (lines 1467--1491):

$$
\begin{aligned}
&\Gamma_1=\tfrac3{32},\ \Gamma_2=\tfrac3{16},\ \Gamma_3=0,\ \Gamma_4=\tfrac18,\ \Gamma_5=\tfrac38,\ \Gamma_6=\tfrac1{16},\ \Gamma_7=0,\\
&\Gamma_8=\tfrac3{16},\ \Gamma_9=\tfrac14,\ \Gamma_{10}=-\tfrac14,\ \Gamma_{11}=-\tfrac18,\ \Gamma_{12}=\tfrac38,\ \Gamma_{13}=0,\ \Gamma_{14}=\tfrac1{16},\\
&\Gamma_{15}=\tfrac18,\ \Gamma_{16}=\tfrac1{32},\ \Gamma_{17}=0,\ \Gamma_{18}=-\tfrac14,\ \Gamma_{19}=\tfrac1{24},\ \Gamma_{20}=\tfrac1{12}.
\end{aligned}
$$

The table stops at $i=20$: no $\Gamma_{25}$ and no $\Gamma_{47}$ are listed. They coincide with $SU(3)$ except for $\Gamma_6$, $\Gamma_8$, $\Gamma_{12}$ (the operators with $\langle M^2\rangle$, $\langle M\rangle^2$, $\langle N^2\rangle$, $\langle N\rangle^2$), because $C_\mu$ and $F_{L,R}^{\mu\nu}$ are traceless for $SU(3)$ but $M$, $N$ are not. The $SU(3)$ values quoted for comparison (lines 1505--1509):

$$
\Gamma_6^{[SU]}=\frac{11}{144},\qquad\Gamma_8^{[SU]}=\frac5{48},\qquad\Gamma_{12}^{[SU]}=\frac5{24}.
$$

Checks stated (lines 1520--1525): $\frac98(12\Gamma_6^{[SU]}+2\Gamma_8^{[SU]}+\Gamma_{12}^{[SU]})=12\Gamma_6+2\Gamma_8+\Gamma_{12}$ and $\frac98(3\Gamma_6^{[SU]}+\Gamma_8^{[SU]})=3\Gamma_6+\Gamma_8$; general $n_l$: $\Gamma_6^{[SU]}=\frac{2+n_l^2}{16n_l^2}$, $\Gamma_8^{[SU]}=\frac{n_l^2-4}{16n_l}$ (attributed to a Bijnens--Bramon--Cornet reference, lines 1536--1538).

These values agree with Eq. (scaling) of [@Kaiser:1998ds] for $i=4,\dots,8,18$.

## $1/N_c$ expansion of the $B_i$

Expanding Appendix D in $X$ and then in $1/N_c$ (keeping $1/N_c$ for the $O(p^4)$ terms, lines 2022--2030), "the $B_i$'s for $i=0$ to 20 are the same as in (llista) except for two new contributions" (lines 2044--2050):

$$
B_8=\frac{n_l}{16}-\frac12v_{31}+O(N_c^{-2}),\qquad B_{12}=\frac{n_l}8-v_{31}+O(N_c^{-2}).
$$

"The rest of coefficients either vanish exactly or do not contribute to this order" (lines 2052--2053). Large-$N_c$ orders of the couplings (lines 1993--2008): $L_{5}(0),L_8(0)=O(N_c)$; $L_4(0),L_6(0),L_7(0),L_{18}(0)=O(1)$; $L'_{21}(0),L'_{25}(0),L'_{28}(0)=O(1)$.

So the number $\Gamma_8=3/16$ is the leading term at $n_l=3$; the paper's own next term $-\frac12v_{31}$ is dropped by [@Kaiser:1998ds] (see [the basis dictionary](u3-notation-dictionary-with-kaiser-1998ds.md) for $v_{31}\leftrightarrow\Lambda_2$).

## Full functions for $\beta_{18}$, $\beta_{25}$, $\beta_{47}$ (Appendix D)

$B_{18}$ (lines 2724--2729), $B_{25}$ (lines 2756--2760) and $B_{47}$ (lines 2891--2897) are functions of $\omega_i^{(k)}$ and $n_l$; $B_{47}$ is as transcribed in [the NLO record](../../nlo/one-loop-effects-and-renormalized-lecs.md). $O_{47}=D_\mu\hat\theta\langle C^\mu\rangle\langle M\rangle$ (line 2516), $O_{18}=\langle U^\dagger D_\mu U\rangle\langle D^\mu U^\dagger\chi-D^\mu U\chi^\dagger\rangle$ (line 1263), $O_{25}=i\langle U^\dagger\chi U^\dagger\chi-\chi^\dagger U\chi^\dagger U\rangle=i\langle NM\rangle$ (lines 1286--1287).

Evaluating $B_{25}'(0)$ (needed for the $L'_{25}(0)$ constant) from the Appendix D expression has not been done in a script; $\Gamma_{25}=0$ is therefore stated here only as quoted by [@Kaiser:1998ds], not as reproduced from the $B_{25}$ expression.
