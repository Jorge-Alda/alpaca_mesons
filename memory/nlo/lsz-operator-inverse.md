---
title: NLO extension of the ALP–meson linear inverse operator
status: derived
created: 2026-11-02
updated: 2026-11-02
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds]
provenance: analytical-derivation
---

# NLO extension of the ALP–meson linear inverse operator

## Inputs and conventions

Use the full $U(3)$ pseudoscalar matrix space, the general noncommuting source in the [LSZ invariant-source record](lsz-invariant-split.md), and the local quadratic operators in the [NLO quadratic Lagrangian](quadratic-lagrangian.md) [@Alda:2026nrd; @Kaiser:1998ds]. Let $w=m_a^2$, $\chi_0=2B_0m_q$ be diagonal, and $g=m_0^2/3$. The LEC corrections are retained linearly for a strict NLO result. Loops and the omitted higher-order operators are not included.

For any Hermitian flavor matrix $X$, define the LO operator

$$
\mathbb M_{a,0}[X]=\frac12\{\chi_0,X\}+g\operatorname{tr}(X)\mathbf1-wX.
$$

This is the full-matrix form of Alda:2026nrd's neutral equation $(\chi_i-w)S_i+gT_S=C_i$. The source is $C$ from [the invariant-source result](lsz-invariant-split.md), and $S_0=\mathbb M_{a,0}^{-1}[C]$. Its matrix inverse can be written without component sums:

$$
\mathbb M_{a,0}^{-1}[Y]
=\mathbb D_{a,0}^{-1}[Y]
-\frac{g\,\operatorname{tr}(\mathbb D_{a,0}^{-1}[Y])}
{1+g\,\operatorname{tr}(\mathbb D_{a,0}^{-1}[\mathbf1])}
\mathbb D_{a,0}^{-1}[\mathbf1],
\qquad
\mathbb D_{a,0}[X]=\frac12\{\chi_0,X\}-wX.
$$

## Explicit LO inverse for a general source

In the quark-mass basis write $\chi_0=\operatorname{diag}(\chi_u,\chi_d,\chi_s)$. For an arbitrary source matrix $Y$, the off-diagonal equations have no anomaly trace term:

$$
(\mathbb D_{a,0}^{-1}[Y])_{ij}
=\frac{Y_{ij}}{(\chi_i+\chi_j)/2-w},\qquad i\ne j,
$$

where $\mathbb D_{a,0}[X]=\frac12\{\chi_0,X\}-wX$. For the diagonal components set

$$
d_i=\chi_i-w,\qquad
T_Y=\frac{\displaystyle\sum_iY_{ii}/d_i}
{\displaystyle 1+g\sum_i1/d_i}.
$$

Then

$$
(\mathbb M_{a,0}^{-1}[Y])_{ii}=\frac{Y_{ii}-gT_Y}{d_i}.
$$

These expressions follow by first applying $\mathbb D_{a,0}^{-1}$ and then solving the scalar trace equation $T_Y=\operatorname{tr}\mathbb D_{a,0}^{-1}[Y]-gT_Y\operatorname{tr}\mathbb D_{a,0}^{-1}[\mathbf1]$. Apparent individual denominators can be removable; the full operator must be nonsingular, and the ALP–meson overlaps must remain perturbative.

## NLO operator

The operator obtained from the meson two-point function at $p^2=w$ is the following matrix expression:

$$
\begin{aligned}
\mathbb M_a^{\rm loc,NLO}[X]={}&\mathbb D_a[X]
+\frac{m_0^2-w\Lambda_1}{3}\operatorname{tr}(X)\mathbf1\\
&+\left(\frac{\Lambda_2}{3}-\frac{4wL_{18}}{f_\pi^2}\right)
\left[\operatorname{tr}(\chi_0X)\mathbf1+\operatorname{tr}(X)\chi_0\right]
\\
&+\frac{16L_7}{f_\pi^2}\operatorname{tr}(\chi_0X)\chi_0
-\frac{8L_{25}}{f_\pi^2}
\left[\operatorname{tr}(\chi_0^2X)\mathbf1+\operatorname{tr}(X)\chi_0^2\right].
\end{aligned}
$$

where

$$
\begin{aligned}
\mathbb D_a[X]={}&\frac12\{\chi_0,X\}-wX
+\frac{4L_8}{f_\pi^2}\{\chi_0,\{\chi_0,X\}\}
+\frac{8L_6}{f_\pi^2}\operatorname{tr}(\chi_0)\{\chi_0,X\}\\
&-w\left[
\frac{4L_5}{f_\pi^2}\{\chi_0,X\}
+\frac{8L_4}{f_\pi^2}\operatorname{tr}(\chi_0)X
\right],
\end{aligned}
$$

and the trace coefficient matrix, in the ordered trace basis $(\mathbf1,\chi_0,\chi_0^2)$, is

$$
H=\begin{pmatrix}
(m_0^2-w\Lambda_1)/3 & \Lambda_2/3-4wL_{18}/f_\pi^2 & -8L_{25}/f_\pi^2\\
\Lambda_2/3-4wL_{18}/f_\pi^2 & 16L_7/f_\pi^2 & 0\\
-8L_{25}/f_\pi^2 & 0 & 0
\end{pmatrix}.
$$

Thus all NLO terms beyond the anticommutator polynomial enter through at most three trace moments. Since $\chi_0$ is diagonal, $\mathbb D_a$ is diagonal on matrix units:

$$
(\mathbb D_a^{-1}[Y])_{ij}=\frac{Y_{ij}}{d_{ij}},
$$

with

$$
\begin{aligned}
d_{ij}={}&\frac{\chi_i+\chi_j}{2}-w
+\frac{4L_8}{f_\pi^2}(\chi_i+\chi_j)^2
+\frac{8L_6\operatorname{tr}(\chi_0)}{f_\pi^2}(\chi_i+\chi_j)\\
&-\frac{4wL_5}{f_\pi^2}(\chi_i+\chi_j)
-\frac{8wL_4\operatorname{tr}(\chi_0)}{f_\pi^2}.
\end{aligned}
$$

For the component-level trace system, define the three trace moments and its $3\times3$ matrices by

$$
b_r=\operatorname{tr}\!\left(\chi_0^r\mathbb D_a^{-1}[Y]\right),\qquad
G_{rt}=\operatorname{tr}\!\left(\chi_0^r\mathbb D_a^{-1}[\chi_0^t]\right),
\qquad r,t=0,1,2.
$$

The trace moments of the solution are $T=(I_3+GH)^{-1}b$, and set $q=HT$. The exact inverse of this **truncated local operator**, written as a matrix expression, is

$$
\boxed{(\mathbb M_a^{\rm loc,NLO})^{-1}[Y]
=\mathbb D_a^{-1}\!\left[Y-q_0\mathbf1-q_1\chi_0-q_2\chi_0^2\right].}
$$

No inverse of $H$ is needed. At LO, $H_{00}=g$ is the only nonzero entry; the formula reduces to the rank-one trace inverse above. Exact inversion of the truncated operator resums products of LECs, so it should be distinguished from a strict NLO expansion.

## Strict NLO source solution

Write $\mathbb M_a^{\rm loc,NLO}=\mathbb M_{a,0}+\Delta\mathbb M_a$ and $C_{\rm NLO}=C+\Delta C_{\rm NLO}$, with $C$ and $\Delta C_{\rm NLO}$ given in [the invariant-source record](lsz-invariant-split.md). Expanding the operator equation $\mathbb M_a^{\rm loc,NLO}[S]=C_{\rm NLO}$ to first order gives

$$
S_0=\mathbb M_{a,0}^{-1}[C],\qquad
S_{\rm NLO}=S_0+\mathbb M_{a,0}^{-1}
\left[\Delta C_{\rm NLO}-\Delta\mathbb M_a[S_0]\right].
$$

The correction $\Delta\mathbb M_a[S_0]$ includes both the NLO mass operator and $-w$ times the NLO kinetic operator. This is the strict-order generalization of Alda:2026nrd's $S=\mathbb M_a^{-1}(C)$ and is the quantity to use in the ALP overlap $r_A+(M^2-wK)^{-1}\mathcal I$.

## Reproducibility

The executable derivation is [`scripts/derive_nlo_lsz_operator_inverse.py`](../../scripts/derive_nlo_lsz_operator_inverse.py). Run

```text
.venv/bin/python scripts/derive_nlo_lsz_operator_inverse.py
.venv/bin/python scripts/derive_nlo_lsz_operator_inverse.py --test
```

It differentiates the full-matrix trace operators, derives $\mathbb D_a$ and $H$, constructs the trace-system inverse, and verifies the LO limit, exact solution for a non-diagonal Hermitian source, and strict first-order residual. The pole domain excludes singularities of the full $\mathbb M_a$; intermediate singularities of $\mathbb D_a$ may be removable in the trace-system expression.
