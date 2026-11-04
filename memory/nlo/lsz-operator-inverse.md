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

## Step-by-step construction of the NLO operator

Use the flavor-matrix inner product $\langle X,Y\rangle=\frac12\operatorname{tr}(XY)$ for Hermitian fields. The quadratic meson Lagrangian is organized as $\frac12\langle\partial\Phi,K\partial\Phi\rangle-\frac12\langle\Phi,M^2\Phi\rangle$. At the ALP pole, $w=m_a^2$, the linear operator is $M^2-wK$.

### 1. Read the kinetic and mass Hessians from the bilinears

The local NLO kinetic bilinears in [quadratic-lagrangian.md](quadratic-lagrangian.md) are

$$
\mathcal L_{\rm kin,NLO}^{(2)}=
\frac{2L_5}{f_\pi^2}\operatorname{tr}(\chi_0\partial\Phi\partial\Phi)
+\frac{\Lambda_1}{12}\operatorname{tr}(\partial\Phi)^2
+\frac{2L_4}{f_\pi^2}\operatorname{tr}(\partial\Phi\partial\Phi)\operatorname{tr}(\chi_0)
+\frac{2L_{18}}{f_\pi^2}\operatorname{tr}(\chi_0\partial\Phi)\operatorname{tr}(\partial\Phi).
$$

Differentiating these quadratic forms with respect to the meson matrix gives the kinetic correction operator

$$
\begin{aligned}
\delta K[X]={}&\frac{4L_5}{f_\pi^2}\{\chi_0,X\}
+\frac{\Lambda_1}{3}\operatorname{tr}(X)\mathbf1
+\frac{8L_4}{f_\pi^2}\operatorname{tr}(\chi_0)X\\
&+\frac{4L_{18}}{f_\pi^2}
\left[\operatorname{tr}(\chi_0X)\mathbf1+\operatorname{tr}(X)\chi_0\right].
\end{aligned}
$$

The nondifferentiated NLO bilinears give

$$
\begin{aligned}
\delta M^2[X]={}&\frac{4L_8}{f_\pi^2}\{\chi_0,\{\chi_0,X\}\}
+\frac{8L_6}{f_\pi^2}\operatorname{tr}(\chi_0)\{\chi_0,X\}
+\frac{16L_7}{f_\pi^2}\operatorname{tr}(\chi_0X)\chi_0\\
&+\frac{\Lambda_2}{3}
\left[\operatorname{tr}(\chi_0X)\mathbf1+\operatorname{tr}(X)\chi_0\right]
-\frac{8L_{25}}{f_\pi^2}
\left[\operatorname{tr}(\chi_0^2X)\mathbf1+\operatorname{tr}(X)\chi_0^2\right].
\end{aligned}
$$

These are matrix identities for arbitrary Hermitian $X$; no coupling matrix is commuted through $\chi_0$.

### 2. Form the pole operator $M^2-wK$

The LO Hessians give $\frac12\{\chi_0,X\}+g\operatorname{tr}(X)\mathbf1$ for $M^2$ and $X$ for $K$. Add the NLO mass operator and subtract $w$ times the NLO kinetic operator. This produces

$$
\begin{aligned}
\mathbb M_a^{\rm loc,NLO}[X]={}&\mathbb D_a[X]
+\frac{m_0^2-w\Lambda_1}{3}\operatorname{tr}(X)\mathbf1\\
&+\left(\frac{\Lambda_2}{3}-\frac{4wL_{18}}{f_\pi^2}\right)
\left[\operatorname{tr}(\chi_0X)\mathbf1+\operatorname{tr}(X)\chi_0\right]
\\
&+\frac{16L_7}{f_\pi^2}\operatorname{tr}(\chi_0X)\chi_0
-\frac{8L_{25}}{f_\pi^2}
\left[\operatorname{tr}(\chi_0^2X)\mathbf1+\operatorname{tr}(X)\chi_0^2\right],
\end{aligned}
$$

where the non-trace part is grouped into

$$
\begin{aligned}
\mathbb D_a[X]={}&\frac12\{\chi_0,X\}-wX
+\frac{4L_8}{f_\pi^2}\{\chi_0,\{\chi_0,X\}\}
+\frac{8L_6}{f_\pi^2}\operatorname{tr}(\chi_0)\{\chi_0,X\}\\
&-w\left[
\frac{4L_5}{f_\pi^2}\{\chi_0,X\}
+\frac{8L_4}{f_\pi^2}\operatorname{tr}(\chi_0)X
\right].
\end{aligned}
$$

The operator $O_{47}$ contributes to the ALP–meson mixed kinetic column $\Delta_i$ and hence to the invariant source $\mathcal C_i^{\rm pole}$, as specified in the [LSZ invariant-source record](lsz-invariant-split.md). It has no meson-only $\Phi^2$ quadratic term, so it adds no term to $K_{ij}$, $M^2_{ij}$, or the linear map $\mathbb M_a^{\rm loc,NLO}:X\mapsto(M^2-wK)[X]$ derived here. The corrected NLO solution is obtained by applying this same operator inverse to the source including $\delta_{47}\bm{\mathcal C}_{\rm NLO}$.

The trace terms left outside $\mathbb D_a$ only depend on $X$ through the three scalar moments $t_0=\operatorname{tr}(X)$, $t_1=\operatorname{tr}(\chi_0X)$, and $t_2=\operatorname{tr}(\chi_0^2X)$. Substituting these moments into the remaining operator gives

$$
\begin{aligned}
\mathbb R_a[X]={}&
\left[\frac{m_0^2-w\Lambda_1}{3}t_0
+\left(\frac{\Lambda_2}{3}-\frac{4wL_{18}}{f_\pi^2}\right)t_1
-\frac{8L_{25}}{f_\pi^2}t_2\right]\mathbf1\\
&+\left[\left(\frac{\Lambda_2}{3}-\frac{4wL_{18}}{f_\pi^2}\right)t_0
+\frac{16L_7}{f_\pi^2}t_1\right]\chi_0
-\frac{8L_{25}}{f_\pi^2}t_0\chi_0^2.
\end{aligned}
$$

Thus $\mathbb R_a[X]$ always lies in the three-dimensional span of $\mathbf1,\chi_0,\chi_0^2$, even though $X$ ranges over the full nine-dimensional flavor-matrix space. This is why $\mathbb D_a$ is useful: it contains the part that can be inverted element by element, leaving only three scalar moments to determine. It is a computational decomposition of $\mathbb M_a^{\rm loc,NLO}$, not an additional physical operator.

### 3. Invert the elementwise part

Since $\chi_0$ is diagonal, write $E_{ij}$ for a flavor matrix unit. Then $\{\chi_0,E_{ij}\}=(\chi_i+\chi_j)E_{ij}$ and $\{\chi_0,\{\chi_0,E_{ij}\}\}=(\chi_i+\chi_j)^2E_{ij}$. Therefore $\mathbb D_a[E_{ij}]=d_{ij}E_{ij}$ and

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

This is a component-level equation; the preceding operator equations remain in matrix form. In particular, $\mathbb D_a$ does not mix different $E_{ij}$ directions. For each $E_{ij}$ the two anticommutators multiply by $\chi_i+\chi_j$ and $(\chi_i+\chi_j)^2$, respectively; the displayed $d_{ij}$ is the resulting scalar eigenvalue. Inverting $\mathbb D_a$ therefore means dividing each matrix entry by its own $d_{ij}$. The trace operator $\mathbb R_a$ is what couples those entries back together.

### 4. Solve the remaining trace system

The trace operator just written has output coefficients $h_0,h_1,h_2$ multiplying $\mathbf1,\chi_0,\chi_0^2$. These coefficients are linear in the input moments:

$$
\begin{pmatrix}h_0\\h_1\\h_2\end{pmatrix}
=H\begin{pmatrix}t_0\\t_1\\t_2\end{pmatrix}.
$$

In components of this three-dimensional trace space, this says

$$
\begin{aligned}
h_0&=\frac{m_0^2-w\Lambda_1}{3}t_0
+\left(\frac{\Lambda_2}{3}-\frac{4wL_{18}}{f_\pi^2}\right)t_1
-\frac{8L_{25}}{f_\pi^2}t_2,\\
h_1&=\left(\frac{\Lambda_2}{3}-\frac{4wL_{18}}{f_\pi^2}\right)t_0
+\frac{16L_7}{f_\pi^2}t_1,\\
h_2&=-\frac{8L_{25}}{f_\pi^2}t_0.
\end{aligned}
$$

The coefficient matrix, in the ordered trace basis $(\mathbf1,\chi_0,\chi_0^2)$, is

$$
H=\begin{pmatrix}
(m_0^2-w\Lambda_1)/3 & \Lambda_2/3-4wL_{18}/f_\pi^2 & -8L_{25}/f_\pi^2\\
\Lambda_2/3-4wL_{18}/f_\pi^2 & 16L_7/f_\pi^2 & 0\\
-8L_{25}/f_\pi^2 & 0 & 0
\end{pmatrix}.
$$

To solve $\mathbb M_a^{\rm loc,NLO}[X]=Y$, first move the trace part to the right and apply the elementwise inverse:

$$
X=\mathbb D_a^{-1}\!\left[Y-h_0\mathbf1-h_1\chi_0-h_2\chi_0^2\right].
$$

Now take the three traces of this equation. The input source gives the vector $b$, while the response to each basis matrix gives $G$:

$$
b_r=\operatorname{tr}\!\left(\chi_0^r\mathbb D_a^{-1}[Y]\right),\qquad
G_{rt}=\operatorname{tr}\!\left(\chi_0^r\mathbb D_a^{-1}[\chi_0^t]\right),
\qquad r,t=0,1,2.
$$

In this component-level notation $T=(t_0,t_1,t_2)^T$, $h=HT$, and the three trace equations read

$$
T=b-Gh=b-GHT,
\qquad (I_3+GH)T=b.
$$

$I_3$ is the identity in this three-dimensional trace-moment space; it is distinct from the flavor identity $\mathbf1$. Solving gives $T=(I_3+GH)^{-1}b$ and then $q=h=HT$. Substituting these three scalar coefficients back into the matrix equation gives

$$
\boxed{(\mathbb M_a^{\rm loc,NLO})^{-1}[Y]
=\mathbb D_a^{-1}\!\left[Y-q_0\mathbf1-q_1\chi_0-q_2\chi_0^2\right].}
$$

At LO, $H_{00}=g$ is the only nonzero entry; the formula reduces to the rank-one trace inverse above. No inverse of $H$ is needed. This exact inverse is for the truncated local operator and resums products of retained LECs, so it should be distinguished from a strict NLO expansion.

## Strict NLO source solution

Write $\mathbb M_a^{\rm loc,NLO}=\mathbb M_{a,0}+\Delta\mathbb M_a$, $S=S_0+\Delta S$, and $C_{\rm NLO}=C+\Delta C_{\rm NLO}$, with $C$ and $\Delta C_{\rm NLO}$ given in [the invariant-source record](lsz-invariant-split.md). Substitute these into $\mathbb M_a^{\rm loc,NLO}[S]=C_{\rm NLO}$ and retain first order:

$$
\mathbb M_{a,0}[S_0]+\mathbb M_{a,0}[\Delta S]
+\Delta\mathbb M_a[S_0]=C+\Delta C_{\rm NLO}.
$$

The LO equation is $\mathbb M_{a,0}[S_0]=C$. Subtracting it leaves $\mathbb M_{a,0}[\Delta S]=\Delta C_{\rm NLO}-\Delta\mathbb M_a[S_0]$, and applying the LO inverse gives

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
