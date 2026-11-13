---
title: Exact LO trace and determinant of the finite-mixing spectrum in invariant form
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Exact LO trace and determinant of the finite-mixing spectrum in invariant form

## Assumptions, inputs and notation

- Quadratic Lagrangian $\mathcal L_2=\tfrac12\partial\Phi^{\mathsf T}K\partial\Phi-\tfrac12\Phi^{\mathsf T}M^2\Phi$ on $\Phi=(P_i,a)$, with $P_i$ the neutral mesons $i\in N=\{0,3,8,6,7\}$, in the conventions of Alda:2026nrd [@Alda:2026nrd]: $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, $N_f=3$. Blocks as in the [finite-mixing record](finite-alp-meson-mixing.md) and the [charged-decoupling record](alp-lo-charged-meson-decoupling.md).
- Sources $c_q^A$, $c_q^P$, $\chi_0=2B_0m_q$ Hermitian and commuting with $Q=\operatorname{diag}(2/3,-1/3,-1/3)$. Their $u$-$d$ and $u$-$s$ entries vanish.
- $\epsilon=f_\pi/f_a$. $m_{K^0,0}^2=\tfrac12\big[(\chi_0)_{dd}+(\chi_0)_{ss}\big]$ as in the [neutral-block record](alp-lo-neutral-block-reduction.md). $c_G$, $m_0^2$, $m_{a,0}^2$ as in Alda:2026nrd.
- $\widehat S_2$ is the second-order dressed scalar source of the [second-order mass record](alp-dressed-mass-second-order-invariants.md). $c_{GA}\equiv c_G+\operatorname{tr}c_q^A$ and $c_q^{AP}\equiv c_q^A-c_q^P$ are the first-order invariants.
- Quark spurion linear in $a$, plus the $a^2$ term $\Delta M^2_{aa}=-\epsilon^2B_0\operatorname{tr}S_2$ from the [joint-vertex record](alp-lo-joint-vertex-expansion.md). The quadratic blocks depend on no other higher-order coefficient (Step 6).

## Step 1. Entries of $K$ and $M^2$

Entries of the neutral blocks, from the finite-mixing record:
$$
K_{ij}=\delta_{ij},\qquad K_{ia}=-v_i,\qquad K_{aa}=1+v^{\mathsf T}v,\qquad v_i=\epsilon\operatorname{tr}(c_q^A\lambda_i),
$$
$$
(M^2)_{ij}=\tfrac14\operatorname{tr}\big(\{\lambda_i,\lambda_j\}\chi_0\big)+\tfrac{m_0^2}{2N_f}\operatorname{tr}\lambda_i\operatorname{tr}\lambda_j,\qquad
(M^2)_{ia}=-\tfrac{\epsilon}{2}\operatorname{tr}\big(\lambda_i\{c_q^P,\chi_0\}\big)+\tfrac{\epsilon m_0^2c_G}{N_f}\operatorname{tr}\lambda_i,
$$
$$
(M^2)_{aa}=m_{a,0}^2+\tfrac{2\epsilon^2m_0^2}{N_f}c_G^2-\epsilon^2B_0\operatorname{tr}S_2 .
$$
Here $v$ is the column with entries $v_i$, and $K_{PP}=\mathbb 1$. With $K=\begin{pmatrix}\mathbb 1&-v\\-v^{\mathsf T}&1+v^{\mathsf T}v\end{pmatrix}$, direct multiplication gives
$$
K^{-1}=\begin{pmatrix}\mathbb 1+vv^{\mathsf T}&v\\ v^{\mathsf T}&1\end{pmatrix},\qquad \det K=1 .
$$

## Step 2. Two identities for neutral flavour traces

The neutral generators $\lambda_i$, $i\in N$, form an orthonormal basis (with respect to $\tfrac12\operatorname{tr}(XY)$) of the five-dimensional space of Hermitian matrices commuting with $Q$. This space is $\{\operatorname{diag}(x_u)\}\oplus\{\text{Hermitian }2\times2\text{ on }(d,s)\}$. Therefore, for any $X,Y$ in this space,
$$
\sum_{i\in N}\operatorname{tr}(\lambda_iX)\operatorname{tr}(\lambda_iY)=2\operatorname{tr}(XY). \tag{C1}
$$
Explicitly, using the diagonal entries of the generators,
$$
\sum_{i\in N}\lambda_i^2=\operatorname{diag}(2,4,4). \tag{C2}
$$
Also $\operatorname{tr}\lambda_i=\sqrt6\,\delta_{i0}$.

## Step 3. Trace of $K^{-1}M^2$ in flavour traces

By block multiplication,
$$
\operatorname{tr}(K^{-1}M^2)=\operatorname{tr}M^2_{PP}+v^{\mathsf T}M^2_{PP}v+2v^{\mathsf T}M^2_{Pa}+M^2_{aa}. \tag{3.1}
$$

(a) Using $\sum_i\tfrac14\operatorname{tr}(\{\lambda_i,\lambda_i\}\chi_0)=\tfrac12\operatorname{tr}\big(\chi_0\sum_i\lambda_i^2\big)$ and (C2), and $\tfrac{m_0^2}{2N_f}\sum_i(\operatorname{tr}\lambda_i)^2=m_0^2$:
$$
\operatorname{tr}M^2_{PP}=\operatorname{tr}\big(\chi_0\operatorname{diag}(1,2,2)\big)+m_0^2 .
$$

(b) Let $Y\equiv\sum_{i\in N}\operatorname{tr}(\lambda_ic_q^A)\lambda_i=2c_q^A$ by (C1)'s basis expansion. Then
$$
v^{\mathsf T}M^2_{PP}v=\epsilon^2\cdot\tfrac14\operatorname{tr}\big(\{Y,Y\}\chi_0\big)+\tfrac{\epsilon^2m_0^2}{2N_f}\Big(\sum_i\operatorname{tr}(\lambda_ic_q^A)\operatorname{tr}\lambda_i\Big)^2 .
$$
Since $\tfrac14\operatorname{tr}(\{2c_q^A,2c_q^A\}\chi_0)=2\operatorname{tr}\big((c_q^A)^2\chi_0\big)$ and $\sum_i\operatorname{tr}(\lambda_ic_q^A)\operatorname{tr}\lambda_i=\sqrt6\operatorname{tr}(\lambda_0c_q^A)=2\operatorname{tr}c_q^A$,
$$
v^{\mathsf T}M^2_{PP}v=2\epsilon^2\operatorname{tr}\big((c_q^A)^2\chi_0\big)+\tfrac{2\epsilon^2m_0^2}{N_f}\big(\operatorname{tr}c_q^A\big)^2 .
$$

(c) Using (C1) with $X=c_q^A$, $Y=\{c_q^P,\chi_0\}$-part, and the same trace of $\lambda_0$:
$$
v^{\mathsf T}M^2_{Pa}=-\epsilon^2\operatorname{tr}\big(c_q^A\{c_q^P,\chi_0\}\big)+\tfrac{2\epsilon^2m_0^2}{N_f}c_G\operatorname{tr}c_q^A .
$$

Inserting (a)-(c) and $M^2_{aa}$ into (3.1), the $m_0^2$ terms combine as $\tfrac{2\epsilon^2m_0^2}{N_f}\big(\operatorname{tr}c_q^A+c_G\big)^2=\tfrac{2\epsilon^2m_0^2}{N_f}c_{GA}^2$:
$$
\operatorname{tr}(K^{-1}M^2)=\operatorname{tr}\big(\chi_0\operatorname{diag}(1,2,2)\big)+m_0^2+m_{a,0}^2+\epsilon^2\Big[2\operatorname{tr}\big((c_q^A)^2\chi_0\big)-2\operatorname{tr}\big(c_q^A\{c_q^P,\chi_0\}\big)+\tfrac{2m_0^2}{N_f}c_{GA}^2-B_0\operatorname{tr}S_2\Big]. \tag{3.2}
$$

## Step 4. Replacing $\operatorname{tr}S_2$ by $\operatorname{tr}\widehat S_2$

The specialisation of $\widehat S_2$ to $S_0=m_q$, $P_0=0$, $P_1=\{c_q^P,m_q\}$, $C_{A,0}=c_q^A$ in the [second-order mass record](alp-dressed-mass-second-order-invariants.md) reads
$$
\widehat S_2=S_2-2i[C_{V,0},S_1]+2\{c_q^A,P_1\}-i[C_{V,1},m_q]-[C_{V,0},[C_{V,0},m_q]]-\{c_q^A,\{c_q^A,m_q\}\}.
$$
Take the trace term by term:
- every commutator is traceless, so $\operatorname{tr}\widehat S_2$ receives no contribution from $[C_{V,0},S_1]$, $[C_{V,1},m_q]$, $[C_{V,0},[C_{V,0},m_q]]$;
- $\operatorname{tr}\{c_q^A,P_1\}=2\operatorname{tr}(c_q^AP_1)$;
- $\operatorname{tr}\{c_q^A,\{c_q^A,m_q\}\}=4\operatorname{tr}\big((c_q^A)^2m_q\big)$, since $\{C,\{C,m\}\}=C^2m+2CmC+mC^2$.

Hence
$$
\operatorname{tr}\widehat S_2=\operatorname{tr}S_2+4\operatorname{tr}\big(c_q^A\{c_q^P,m_q\}\big)-4\operatorname{tr}\big((c_q^A)^2m_q\big).
$$
With $\chi_0=2B_0m_q$, multiplying by $B_0$:
$$
B_0\operatorname{tr}S_2=B_0\operatorname{tr}\widehat S_2-2\operatorname{tr}\big(c_q^A\{c_q^P,\chi_0\}\big)+2\operatorname{tr}\big((c_q^A)^2\chi_0\big). \tag{4.1}
$$
Substituting (4.1) in (3.2), the $\operatorname{tr}((c_q^A)^2\chi_0)$ and $\operatorname{tr}(c_q^A\{c_q^P,\chi_0\})$ terms cancel exactly:
$$
\boxed{\operatorname{tr}(K^{-1}M^2)=\operatorname{tr}\big(\chi_0\operatorname{diag}(1,2,2)\big)+m_0^2+m_{a,0}^2+\epsilon^2\Big[\tfrac{2m_0^2}{N_f}\,c_{GA}^2-B_0\operatorname{tr}\widehat S_2\Big]}
$$
Every term is now a zeroth-order source or an invariant under constant chiral redefinitions ($c_q^{AP}$, $c_{GA}$, $\widehat S_2$).

## Step 5. Determinant of $M^2$

Take first $\chi_0$ diagonal, $\chi_0=\operatorname{diag}(x_u,x_d,x_s)$. The general case follows from a $Q$-commuting unitary $\operatorname{diag}(1,V)$ that diagonalises the $(d,s)$ block of $\chi_0$. This is a constant orthogonal rotation of the mesons, which leaves $\det M^2$ and $\operatorname{tr}K^{-1}M^2$ invariant.

**Neutral meson block.** For diagonal $\chi_0$ the generators $\lambda_6,\lambda_7$ are off-diagonal and $\lambda_0,\lambda_3,\lambda_8$ are diagonal. Therefore $M^2_{PP}$ splits into a $3\times3$ block on $(\lambda_0,\lambda_3,\lambda_8)$ and a $2\times2$ block on $(\lambda_6,\lambda_7)$, with no mixing between them. The $2\times2$ block is
$$
\tfrac14\operatorname{tr}\big(\{\lambda_6,\lambda_6\}\chi_0\big)=\tfrac12(x_d+x_s)=m_{K^0,0}^2,\qquad \tfrac14\operatorname{tr}\big(\{\lambda_6,\lambda_7\}\chi_0\big)=0,
$$
and $\lambda_7$ gives the same value. So the $2\times2$ block is $m_{K^0,0}^2\mathbb 1_2$.

The $3\times3$ block is $G+m_0^2\,e_0e_0^{\mathsf T}$, with $e_0=(1,0,0)^{\mathsf T}$ and
$$
G=\tfrac12\Lambda^{\mathsf T}\operatorname{diag}(x_u,x_d,x_s)\Lambda,\qquad
\Lambda=\begin{pmatrix}\sqrt{2/3}&1&1/\sqrt3\\ \sqrt{2/3}&-1&1/\sqrt3\\ \sqrt{2/3}&0&-2/\sqrt3\end{pmatrix}.
$$
Here the columns of $\Lambda$ are the diagonal entries of $\lambda_0,\lambda_3,\lambda_8$. One checks $\Lambda^{\mathsf T}\Lambda=2\mathbb 1_3$, hence $\det\Lambda^2=8$ and $\det G=\tfrac18\cdot8\,x_ux_dx_s=x_ux_dx_s$.

So $\det M^2_{PP}=m_{K^0,0}^4\,\det\big(G+m_0^2e_0e_0^{\mathsf T}\big)$. Using the matrix determinant lemma and $G^{-1}=2\Lambda^{-1}\operatorname{diag}(x)^{-1}\Lambda^{-\mathsf T}$ with $\Lambda^{-\mathsf T}e_0=\tfrac12\Lambda e_0=\tfrac12\sqrt{2/3}(1,1,1)^{\mathsf T}$:
$$
e_0^{\mathsf T}G^{-1}e_0=\tfrac13\operatorname{tr}\chi_0^{-1},\qquad
\det M^2_{PP}=m_{K^0,0}^4\,D_\chi,\qquad D_\chi\equiv\det\chi_0\Big(1+\tfrac{m_0^2}{N_f}\operatorname{tr}\chi_0^{-1}\Big). \tag{5.1}
$$

**Schur complement.** The mixing column is $M^2_{Pa}=w$ with $w_i=-\tfrac{\epsilon}{2}\operatorname{tr}(\lambda_i\{c_q^P,\chi_0\})+\tfrac{\epsilon m_0^2c_G}{N_f}\operatorname{tr}\lambda_i$. Split $w$ into its $(\lambda_0,\lambda_3,\lambda_8)$ part $w_D$ and its $(\lambda_6,\lambda_7)$ part $w_{ds}$.

- $(d,s)$ part. For diagonal $\chi_0$, $\{c_q^P,\chi_0\}$ restricted to $(d,s)$ equals $(x_d+x_s)c_q^P|_{ds}$. Writing $\alpha$ for the $(d,s)$ entry of $c_q^P$, $\operatorname{tr}(\lambda_6c_q^P)=2\operatorname{Re}\alpha$ and $\operatorname{tr}(\lambda_7c_q^P)=-2\operatorname{Im}\alpha$. Hence $w_{ds}=-\epsilon\,m_{K^0,0}^2\,(\operatorname{tr}\lambda_6c_q^P,\operatorname{tr}\lambda_7c_q^P)^{\mathsf T}$, and
$$
w_{ds}^{\mathsf T}\big(m_{K^0,0}^2\mathbb 1_2\big)^{-1}w_{ds}=4\epsilon^2m_{K^0,0}^2|\alpha|^2=2\epsilon^2\operatorname{tr}\big((c_q^P)^2\chi_0\big)\big|_{ds}.
$$
The last equality uses $\operatorname{tr}((c_q^P)^2\chi_0)|_{ds}=(x_d+x_s)|\alpha|^2=2m_{K^0,0}^2|\alpha|^2$.

- $(u,d,s)$-diagonal part. Write $w_D=-\epsilon\Lambda^{\mathsf T}\operatorname{diag}(x)\,\hat c+\beta e_0$ with $\hat c$ the diagonal entries of $c_q^P$ and $\beta=\sqrt6\,\epsilon m_0^2c_G/N_f$. Then $G^{-1}\Lambda^{\mathsf T}\operatorname{diag}(x)\hat c=\Lambda^{\mathsf T}\hat c$, which gives $w_D^{\mathsf T}G^{-1}w_D$ term by term:
$$
w_D^{\mathsf T}G^{-1}w_D=2\epsilon^2\operatorname{tr}\big((c_q^P)^2\chi_0\big)\big|_{\rm diag}-\tfrac{4\epsilon^2m_0^2}{N_f}c_G\operatorname{tr}c_q^P+\tfrac{2\epsilon^2m_0^4c_G^2}{9}\operatorname{tr}\chi_0^{-1}.
$$
Using the Sherman-Morrison formula for $G+m_0^2e_0e_0^{\mathsf T}$, the rank-one correction is
$$
-\frac{m_0^2\big(e_0^{\mathsf T}G^{-1}w_D\big)^2}{1+\tfrac{m_0^2}{N_f}\operatorname{tr}\chi_0^{-1}},\qquad e_0^{\mathsf T}G^{-1}w_D=\sqrt{\tfrac23}\,\epsilon\Big(-\operatorname{tr}c_q^P+\tfrac{m_0^2c_G}{N_f}\operatorname{tr}\chi_0^{-1}\Big).
$$

Together with the $(d,s)$ part, the Schur term $S=w^{\mathsf T}(M^2_{PP})^{-1}w$ is
$$
S=2\epsilon^2\operatorname{tr}\big((c_q^P)^2\chi_0\big)-\tfrac{4\epsilon^2m_0^2}{N_f}c_G\operatorname{tr}c_q^P+\tfrac{2\epsilon^2m_0^4c_G^2}{9}\operatorname{tr}\chi_0^{-1}-\tfrac{2\epsilon^2m_0^2}{3}\,\frac{\big(\operatorname{tr}c_q^P-\tfrac{m_0^2c_G}{N_f}\operatorname{tr}\chi_0^{-1}\big)^2}{1+\tfrac{m_0^2}{N_f}\operatorname{tr}\chi_0^{-1}}.
$$

**Determinant.** $\det M^2=\det M^2_{PP}\,\big(M^2_{aa}-S\big)=m_{K^0,0}^4D_\chi\big(M^2_{aa}-S\big)$. Multiply out with $D_\chi=\det\chi_0\,(1+\tfrac{m_0^2}{N_f}\operatorname{tr}\chi_0^{-1})$. Write $y\equiv\tfrac{m_0^2}{N_f}\operatorname{tr}\chi_0^{-1}$; the $c_G^2$ terms give $\tfrac{2\epsilon^2m_0^2}{N_f}\det\chi_0\,c_G^2\big[(1+y)-y\big]$, and the $c_G\operatorname{tr}c_q^P$ and $(\operatorname{tr}c_q^P)^2$ terms, in the same way, combine into $\tfrac{2\epsilon^2m_0^2}{N_f}\det\chi_0\,(c_G+\operatorname{tr}c_q^P)^2$. The $y$-dependence cancels against the Sherman-Morrison denominator, and
$$
\det M^2=m_{K^0,0}^4\Big\{D_\chi\Big[m_{a,0}^2-2\epsilon^2\operatorname{tr}\big((c_q^P)^2\chi_0\big)-\epsilon^2B_0\operatorname{tr}S_2\Big]+\tfrac{2\epsilon^2m_0^2}{N_f}\det\chi_0\,\big(c_G+\operatorname{tr}c_q^P\big)^2\Big\}. \tag{5.2}
$$

## Step 6. Invariant form of the determinant

Substitute (4.1) into (5.2). The combination in the bracket becomes
$$
2\operatorname{tr}\big((c_q^P)^2\chi_0\big)+B_0\operatorname{tr}S_2=2\operatorname{tr}\big((c_q^A-c_q^P)^2\chi_0\big)+B_0\operatorname{tr}\widehat S_2,
$$
because $\operatorname{tr}\big((c_q^A-c_q^P)^2\chi_0\big)=\operatorname{tr}\big((c_q^A)^2\chi_0\big)-\operatorname{tr}\big(c_q^A\{c_q^P,\chi_0\}\big)+\operatorname{tr}\big((c_q^P)^2\chi_0\big)$. Also $c_G+\operatorname{tr}c_q^P=c_{GA}-\operatorname{tr}c_q^{AP}$. Therefore
$$
\boxed{\det M^2=m_{K^0,0}^4\Big\{D_\chi\Big[m_{a,0}^2-\epsilon^2\Big(2\operatorname{tr}\big((c_q^{AP})^2\chi_0\big)+B_0\operatorname{tr}\widehat S_2\Big)\Big]+\tfrac{2\epsilon^2m_0^2}{N_f}\det\chi_0\,\big(c_{GA}-\operatorname{tr}c_q^{AP}\big)^2\Big\}}
$$
with $D_\chi=\det\chi_0+\tfrac{m_0^2}{N_f}\det\chi_0\operatorname{tr}\chi_0^{-1}$.

Every quantity here is invariant under constant chiral redefinitions: $\chi_0$ (zeroth order, fixed), $c_q^{AP}$, $c_{GA}$, $\widehat S_2$. The shift $c_G\to c_G-\operatorname{tr}\kappa$ that makes $c_{GA}$ invariant is taken from the task statement, not derived here.

## Dependence on higher Taylor coefficients

The exact quadratic blocks depend on the sources through $C_{A,0}$, $S_0=m_q$, $P_1$, $\operatorname{tr}S_2$ and $\Theta_1=c_G$ only. They do not depend on $C_{R,n\ge1}$, $C_{L,n\ge1}$, $M_{n\ge3}$, $\Theta_{n\ge2}$, $S_1$, $P_0$ or $P_2$. The check is [`scripts/check_alp_quadratic_taylor_dependence.py`](../../../scripts/check_alp_quadratic_taylor_dependence.py). Higher coefficients enter only cubic and higher vertices and loop or tadpole effects, which are outside this record.

## Verification

[`scripts/derive_alp_lo_mass_sum_rules.py`](../../../scripts/derive_alp_lo_mass_sum_rules.py) checks (3.2), (5.2) and Step 6 against exact rational matrices and against the roots of $\det(sK-M^2)$ from [`scripts/solve_alp_meson_spectrum.py`](../../../scripts/solve_alp_meson_spectrum.py). The script checks the derivation; it does not derive anything. The Step 5 algebra was carried out by hand after the rewrite; the script result is from before the rewrite and has not been re-run on this text.

```bash
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_mass_sum_rules.py
.venv/bin/python scripts/derive_alp_lo_mass_sum_rules.py
```

The independent exact-rational comparison in [`scripts/check_alp_arrow_mass_sum_rules.py`](../../../scripts/check_alp_arrow_mass_sum_rules.py) checks both identities against the eigenvalues of the canonical arrow matrix and against `solve_alp_meson_spectrum.py`. The illustrative input uses an SU(3)-symmetric spurion, rational GeV$^2$ mass entries and nonzero singlet-ALP mixing (including nonzero $S_2$); this restriction permits closed-form exact roots and is not a physical specialization of the sum rules. All four trace and product comparisons agree exactly. Run it with `.venv/bin/python scripts/check_alp_arrow_mass_sum_rules.py` and `.venv/bin/python -m unittest discover -s scripts -p test_check_alp_arrow_mass_sum_rules.py`.
