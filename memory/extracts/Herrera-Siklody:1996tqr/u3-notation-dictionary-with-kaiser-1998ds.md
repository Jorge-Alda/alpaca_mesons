---
title: Dictionary between Herrera-Siklody:1996tqr and Kaiser:1998ds notation
status: mixed
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Herrera-Siklody:1996tqr, Kaiser:1998ds]
provenance: literature
---

# Dictionary of constants: $O_i,W_i,\beta_i$ of [@Herrera-Siklody:1996tqr] versus $L_i,\Lambda_i,V_i$ of [@Kaiser:1998ds]

Items marked (stated) are written in the sources; items marked (matched) were obtained by comparing the two Lagrangians term by term, by hand (no script), and are not stated in either paper.

## Fields

- $U$: $U=\exp(i\Phi/f)$ with $\det U=\exp(i\sqrt{2n_l}\,\eta^0/f)$ (`un.tex` lines 577--587); in [@Kaiser:1998ds], $U=e^{i\psi/3}\hat U$, so $\langle\ln U\rangle=i\psi$ (`adehep.tex` line 644). Hence $\psi=\sqrt{2n_l}\,\eta^0/f$ (matched).
- Source combination: $X=i\sqrt{2n_l}\eta^0/f+\hat\theta=i(\psi+\theta)=i\tilde\psi$ (matched; $X$ in `un.tex` line 1874, $\tilde\psi$ in `adehep.tex` line 314), with $\hat\theta=i\theta$.
- Building blocks: $C_\mu=U^\dagger D_\mu U$, $M=U^\dagger\chi+\chi^\dagger U$, $N=U^\dagger\chi-\chi^\dagger U$ (`un.tex` lines 1208--1210).

## Leading Lagrangian

- $\mathcal L_{(0+2)}$ of Eq. (3.7), `un.tex` lines 755--762: $-W_0+W_1\langle D_\mu U^\dagger D^\mu U\rangle+W_2\langle M\rangle+iW_3\langle N\rangle+W_4\langle U^\dagger D_\mu U\rangle^2+W_5\langle U^\dagger D_\mu U\rangle D^\mu\hat\theta+W_6D_\mu\hat\theta D^\mu\hat\theta$.
- [@Kaiser:1998ds] Eq. (10), lines 301--312: $V_0,\dots,V_5$ with $V_4D_\mu\psi D^\mu\theta$, $V_5D_\mu\theta D^\mu\theta$ and $V_4=0$ chosen (line 335).
- $W_0\leftrightarrow V_0$ (matched, up to the sign convention for derivatives in $X=i\tilde\psi$ that is not fixed here): $V_0=\text{const}+\frac12\tau\tilde\psi^2$ in [@Kaiser:1998ds] (lines 349--354); the singlet mass $m^2_{\eta^0}|_{U_A(1)}=-\frac{2n_l}{f^2}W_0''(0)$ is stated in [@Herrera-Siklody:1996tqr] (lines 776--778). The relation between $\tau$ and $W_0''(0)$ is not worked out here.
- $F^2/4\leftrightarrow W_1(0)=W_2(0)=f^2/4$ (stated, line 768; line 768 also sets $W_4(0)=0$, while lines 805--811 keep $W_4(0)$ as a normalization freedom for $f_0$), $F=f$ (matched).
- $\Lambda_1$: the singlet kinetic coefficient is $\frac{f^2-4n_lW_4(0)}{f_0^2}=1$ (stated, lines 805--811); with $n_l=3$ and $f_0^2=F_0^2=(1+\Lambda_1)F^2$ (stated in [@Kaiser:1998ds], line 614), this gives $\Lambda_1=-12W_4(0)/f^2$ (matched).
- $\Lambda_2$: [@Kaiser:1998ds] has $\frac1{12}F^2\Lambda_2\,i\tilde\psi\langle\chi^\dagger U-U^\dagger\chi\rangle$ (line 466), and Eq. (13), lines 349--354, $V_2=\frac12F^2B(1-\frac i3\Lambda_2\tilde\psi)$. In [@Herrera-Siklody:1996tqr], $W_3(X)=-i\frac{f^2}4(v_{31}X+\dots)$ (lines 1969--1970). Matching the coefficient of $\tilde\psi\langle U^\dagger\chi-\chi^\dagger U\rangle$ gives $v_{31}=-\Lambda_2/3$ (matched).
- Other functions: $v_{12}$ is the $X^2$ coefficient of $W_1$; $[@Kaiser:1998ds]$ keeps $V_1=F^2/4+O(1/N_c)$, so $v_{12}$ has no counterpart in its list.

## $O(p^4)$ constants

| [@Kaiser:1998ds] | [@Herrera-Siklody:1996tqr] operator | Constant of [@Herrera-Siklody:1996tqr] |
|---|---|---|
| $L_4$ | $O_4=\langle D_\mu U^\dagger D^\mu U\rangle\langle U^\dagger\chi+\chi^\dagger U\rangle$ | $L_4(0)=\beta_4(0)$ |
| $L_5$ | $O_5$ | $L_5(0)$ |
| $L_6$ | $O_6=\langle M\rangle^2$ | $L_6(0)$ |
| $L_7$ | $O_7=\langle N\rangle^2$ | $L_7(0)$ |
| $L_8$ | $O_8$ | $L_8(0)$ |
| $L_{18}$ | $O_{18}=\langle U^\dagger D_\mu U\rangle\langle D^\mu U^\dagger\chi-D^\mu U\chi^\dagger\rangle$ | $L_{18}(0)$; stated: $L_{18}\equiv L_{18}(0)$ (`adehep.tex` line 430) |
| $L_{25}$ | $O_{25}=i\langle U^\dagger\chi U^\dagger\chi-\chi^\dagger U\chi^\dagger U\rangle$ | $L_{25}'(0)$; stated: $L_{25}\equiv L_{25}'(0)$ (line 431) |
| not retained | $O_{47}=D_\mu\hat\theta\langle C^\mu\rangle\langle M\rangle$ | $L_{47}(X)=\beta_{47}$ |
| $L_{46}=L_{52}=0$ | $O_{46}$, $O_{52}$ | removed by field redefinition (`adehep.tex` footnote, lines 332--336) |

Notes (matched): $\langle U^\dagger D_\mu U\rangle=iD_\mu\psi$ for vanishing sources, so $O_{18}=iD_\mu\psi\langle D^\mu U^\dagger\chi-D^\mu U\chi^\dagger\rangle$ is exactly the $L_{18}$ operator of Eq. (L4). Since $O_{25}=i\langle NM\rangle$ and $X=i\tilde\psi$, the term $L_{25}(X)O_{25}\approx L_{25}'(0)\,i\tilde\psi\,O_{25}$ equals $-L_{25}'(0)\tilde\psi\langle U^\dagger\chi U^\dagger\chi-\chi^\dagger U\chi^\dagger U\rangle$, whereas the $L_{25}$ term of [@Kaiser:1998ds] is $i\tilde\psi L_{25}\langle\cdots\rangle$; a factor of $i$ remains between $L_{25}$ and $L_{25}'(0)$ unless $L_{25}'(0)$ denotes the derivative with respect to $\tilde\psi$ (open: not settled in either paper).

## Number of flavours

[@Herrera-Siklody:1996tqr] keeps $n_l$ generic and sets $n_l=3$ in the $\Gamma_i$ table (lines 1467--1491). [@Kaiser:1998ds] works with three flavours throughout.

## Agreement of $\Gamma_i$

At $n_l=3$: $\Gamma_4=1/8$, $\Gamma_5=3/8$, $\Gamma_6=1/16$, $\Gamma_7=0$, $\Gamma_8=3/16$, $\Gamma_{18}=-1/4$ agree between the two papers. $\Gamma_{25}=0$ appears only in [@Kaiser:1998ds] (`adehep.tex` line 696). $L_{47}$ does not appear in [@Kaiser:1998ds].
