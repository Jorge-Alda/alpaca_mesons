---
title: Analytic LO finite-mixing overlap, secular equation and normalization from the inverse of $\mathbb M_{a,0}$
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Analytic LO finite-mixing overlap, secular equation and normalization

## Assumptions and inputs

- LO Hessian blocks $K$, $M^2$, $D(s)=sK-M^2$ with $\Phi=(P_i,a)$, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, $N_f=3$, linear quark spurion, no independent $a^2$ term: [finite-mixing record](finite-alp-meson-mixing.md).
- LO operator $\mathbb M_{a,0}$ and its inverse (elementwise part plus one scalar trace equation): [LSZ operator inverse](../../nlo/lsz-operator-inverse.md), notation of Alda:2026nrd [@Alda:2026nrd].
- Pole-aligned shift and Schur complement: [single-pole record](alp-single-pole-alignment.md).
- Notation of Alda:2026nrd [@Alda:2026nrd]: $\epsilon=f_\pi/f_a$, $\chi_0=2B_0m_q=\operatorname{diag}(\chi_u,\chi_d,\chi_s)$, $c_q^{AP}=c_q^A-c_q^P$, $c_{GA}=c_G+\operatorname{tr}c_q^A$ (and $c_{GP}=c_{GA}-\operatorname{tr}c_q^{AP}=c_G+\operatorname{tr}c_q^P$, used in the [bare-mass-bound record](alp-lo-bare-alp-mass-bound.md)), the invariant flavour matrix
$$
\mathcal C=-\tfrac12\{\chi_0,c_q^{AP}\}-\frac{m_0^2}{N_f}c_{GA}\mathbf1
$$
(see the [LSZ invariant-split record](../../nlo/lsz-invariant-split.md)), and the operator
$$
\mathbb M_a(X)=\tfrac12\{\chi_0,X\}+\frac{m_0^2}{N_f}\operatorname{tr}(X)\mathbf1-m_a^2X .
$$
Here $\mathbb M_a$ is evaluated at a general squared-mass argument $s$ instead of $m_a^2$ and written $\mathbb M_a(X)\big|_{m_a^2\to s}$; its inverse is $\mathbb M_a^{-1}$ at the same $s$. $c_q^A,c_q^P$ are Hermitian and may have off-diagonal entries. No expansion in $\epsilon$ and no truncation of the mixing.

## Method

1. **Pole equation as an operator equation.** The mixed-field equation $(M^2-sK)_{ij}w^{(a)}_j=sK_{i\hat a}-M^2_{i\hat a}$ of the finite-mixing record, multiplied by $\tfrac12\lambda_i$ and summed, is
$$
\mathbb M_a\big(w^{(a)}_i\lambda_i\big)=\epsilon\big(\{c_q^P,\chi_0\}-2s\,c_q^A-2\tfrac{m_0^2}{N_f}c_G\mathbf1\big),
$$
because $M^2_{ij}w_j=\tfrac12\operatorname{tr}\big(\lambda_i[\tfrac12\{\chi_0,w_j\lambda_j\}+\tfrac{m_0^2}{N_f}\operatorname{tr}(w_j\lambda_j)\mathbf1]\big)$. This is Alda:2026nrd's operator equation with the pole condition $s=m_a^2$ not yet imposed.

2. **Remove the kinetic mixing.** The canonical shift of the finite-mixing record is $v_i=\epsilon\operatorname{tr}(c_q^A\lambda_i)$, i.e. $v_i\lambda_i=2\epsilon c_q^A$. Since $\mathbb M_a(2\epsilon c_q^A)$ cancels every $s$-dependent term of the source and replaces $c_q^P-c_q^A$, $c_G+\operatorname{tr}c_q^A$ by $-c_q^{AP}$, $c_{GA}$,
$$
\mathbb M_a\big((w^{(a)}_i-v_i)\lambda_i\big)=2\epsilon\,\mathcal C ,
\qquad\text{i.e.}\qquad
w^{(a)}_i\lambda_i=2\epsilon\big[c_q^A+\mathbb M_a^{-1}(\mathcal C)\big].
$$
The source $\mathcal C$ does not depend on $s$; all $s$-dependence is in $\mathbb M_a^{-1}$. At $\epsilon\to0$ and $s=m_{a,0}^2$ this is Alda:2026nrd's solution, up to the first-order normalization.

3. **Explicit inverse.** In the quark-mass basis (the trace-system inverse of the [operator-inverse record](../../nlo/lsz-operator-inverse.md) at LO):
$$
\big[\mathbb M_a^{-1}(\mathcal C)\big]_{ij}=\frac{\mathcal C_{ij}}{(\chi_i+\chi_j)/2-s}\quad(i\ne j),
$$
$$
\big[\mathbb M_a^{-1}(\mathcal C)\big]_{ii}=\frac{\mathcal C_{ii}-\frac{m_0^2}{N_f}\operatorname{tr}\mathbb M_a^{-1}(\mathcal C)}{\chi_i-s},\qquad
\operatorname{tr}\mathbb M_a^{-1}(\mathcal C)=\frac{\sum_k\mathcal C_{kk}/(\chi_k-s)}{1+\frac{m_0^2}{N_f}\sum_k1/(\chi_k-s)},
$$
with $\mathcal C_{ij}=-\tfrac12(\chi_i+\chi_j)(c_q^{AP})_{ij}-\tfrac{m_0^2}{N_f}c_{GA}\delta_{ij}$.

4. **Pole normalization.** $K_{\hat a\hat a}+2K_{i\hat a}w^{(a)}_i+K_{ij}w^{(a)}_iw^{(a)}_j=1+\sum_i(w^{(a)}_i-v_i)^2$, so
$$
N(s)=1+2\epsilon^2\operatorname{tr}\big[\mathbb M_a^{-1}(\mathcal C)^2\big],
$$
where the trace of the square of the Hermitian matrix is the sum over all ordered index pairs.

5. **ALP entry in the canonical frame.** After $P_i=P'_i+v_ia$ (the canonical kinetic frame of the finite-mixing record, a fixed constant shift, not the pole-dependent one),
$$
M^2_{\hat a\hat a,\rm can}=m_{a,0}^2+2\epsilon^2\frac{m_0^2}{N_f}c_{GA}^2-2\epsilon^2\big[\operatorname{tr}(c_q^A\{c_q^P,\chi_0\})-\operatorname{tr}(\chi_0\,c_q^Ac_q^A)\big].
$$

6. **Secular equation.** The Schur complement of $D(s)=sK-M^2$ with respect to the meson block gives the scalar condition for the ALP-like pole $s=m_a^2$:
$$
D_a(s)\equiv s-M^2_{\hat a\hat a,\rm can}+2\epsilon^2\operatorname{tr}\big[\mathcal C\,\mathbb M_a^{-1}(\mathcal C)\big]=0,
$$
$$
\operatorname{tr}\big[\mathcal C\,\mathbb M_a^{-1}(\mathcal C)\big]=\sum_{i\ne j}\frac{|\mathcal C_{ij}|^2}{(\chi_i+\chi_j)/2-s}+\sum_i\frac{\mathcal C_{ii}\big[\mathcal C_{ii}-\frac{m_0^2}{N_f}\operatorname{tr}\mathbb M_a^{-1}(\mathcal C)\big]}{\chi_i-s}.
$$
Then $\det D(s)=-\det(M^2_{PP}-s)\,D_a(s)$. For $0\le s<\min\chi_i$ the operator $\mathbb M_a$ is positive, so the meson term is positive and the ALP-like root lies below $M^2_{\hat a\hat a,\rm can}$ (level repulsion). Here $D_a(s)$ is the ALP inverse propagator after the mesons are integrated out (the Schur complement of $D(s)$ with respect to the meson block); its zero is the ALP-like pole.

7. **Normalization as a derivative.** Since $d\,\mathbb M_a^{-1}/ds=\mathbb M_a^{-2}$,
$$
N(s)=\frac{dD_a}{ds},
$$
so the external-leg factor is $1/D_a'(m_a^2)$.

8. **Amplitude.** For an amputated one-leg vertex with components $\mathcal V_{\hat a},\mathcal V_i$ and $\mathcal V=\mathcal V_i\lambda_i$,
$$
\mathcal M(a_{\rm phys},X)=\frac{\mathcal V_{\hat a}+\epsilon\operatorname{tr}\!\big\{\mathcal V\big[c_q^A+\mathbb M_a^{-1}(\mathcal C)\big]\big\}}{\sqrt{N(m_a^2)}},
$$
evaluated at $s=m_a^2$.

## Result: how the problem is solved

The finite-mixing problem splits into a scalar equation and an explicit matrix.

- Solve $D_a(s)=0$ for $m_a^2$. For flavor-diagonal sources on $(P_3,P_8,P_0,a)$ it is a monic quartic, $\det D=-D_a\prod_i(\chi_i-s)\,\big(1+\tfrac{m_0^2}{N_f}\sum_i1/(\chi_i-s)\big)$; for generic sources $\det D$ is a monic polynomial of degree 10.
- Insert $s=m_a^2$ into $\mathbb M_a^{-1}(\mathcal C)$ and $N$ above. No linear system is solved numerically.
- $\mathbb M_a^{-1}(\mathcal C)$, $N$ and $D_a+M^2_{\hat a\hat a,\rm can}$ depend on $c_q^A,c_q^P,c_G$ only through $c_q^{AP}$ and $c_{GA}$, i.e. through $\mathcal C$. $M^2_{\hat a\hat a,\rm can}$ also depends on $c_q^A$ separately. Under $(c_q^A,c_q^P,c_G)\to(c_q^A+\delta,c_q^P+\delta,c_G-\operatorname{tr}\delta)$ it changes by $-2\epsilon^2[\operatorname{tr}(\delta\{c_q^P,\chi_0\})+\operatorname{tr}(\chi_0\delta^2)]$, which is compensated by the $a^2$ term that the transformed spurion acquires (see the [bare-mass-bound record](alp-lo-bare-alp-mass-bound.md)); this record's frame has no such term, so the spectrum computed here is not invariant under the shift.
- For $\epsilon\to0$: $D_a\to s-m_{a,0}^2$, $N\to1$, and $w^{(a)}_i\lambda_i/\epsilon\to2[c_q^A+\mathbb M_a^{-1}(\mathcal C)]$ at $s=m_{a,0}^2$, the first-order solution of Alda:2026nrd.
- Decoupled states: for generic sources each pair $i<j$ has one meson combination orthogonal to $\mathcal C_{ij}$ that decouples exactly at $s=m^2_{ij}$ with no ALP component, so $w$ is not defined there (see exceptional poles in the [single-pole record](alp-single-pole-alignment.md)).
- Nearly degenerate poles: $\mathbb M_a^{-1}(\mathcal C)$ has simple poles at $s=(\chi_i+\chi_j)/2$, so near a meson mass the overlap is large. The formulas stay exact, but $D_a(s)=0$ must be solved without expanding.

## Verification

[`scripts/derive_alp_finite_mixing_operator_inverse.py`](../../../scripts/derive_alp_finite_mixing_operator_inverse.py) checks steps 1–9 (92 checks), symbolically for generic Hermitian $c_q^A,c_q^P$ where tractable and with exact rational data otherwise, and compares the roots of $D_a$, $w^{(a)}_i$ and $N$ with the [pole solver](alp-lo-pole-solver.md) to $10^{-12}$ relative (excluding the decoupled states). Unit tests: `scripts/test_derive_alp_finite_mixing_operator_inverse.py`.

```sh
.venv/bin/python scripts/derive_alp_finite_mixing_operator_inverse.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_finite_mixing_operator_inverse.py
```
