---
title: Frame-invariant one-loop ALP contributions to the meson pole masses and their crossing values of f_pi/f_a
status: mixed
created: 2026-11-21
updated: 2026-11-21
bibtex_keys: [Alda:2026nrd, Gasser:1983yg]
provenance: analytical-derivation
---

# Frame-invariant one-loop ALP contributions to the meson pole masses and their crossing values of $f_\pi/f_a$

The crossing is defined with an observable: the one-loop shift of a physical meson pole mass, $\delta m_n^2$. It is split by loop function, and every coefficient is written with the invariant dressed sources only. No values or representatives are used for $c_q^A$, $c_q^P$ or any other source coefficient.

## Conventions and units

- Units GeV. Notation of [@Alda:2026nrd]: $U=u^2=\exp(i\Phi/f_\pi)$, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, $N_f=3$, $x=a/f_a$, $\epsilon=f_\pi/f_a$ (approved notation of the [weak-mixing record](../../../../meson-mixing/alp-lo-weak-mixing-expansion.md)).
- Invariant sources: $\widehat M(x)=\sum_n\frac{x^n}{n!}(\widehat S_n+i\widehat P_n)$ and $\widehat\Theta(x)=\sum_n\frac{x^n}{n!}\widehat\Theta_n$ (Taylor coefficients are derivatives at $x=0$).
- Alda:2026nrd's tadpole source $C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1$ and the operator $\mathbb M_s(X)=\tfrac12\{\chi_0,X\}+\tfrac{m_0^2}{N_f}\operatorname{tr}(X)\mathbf 1-sX$, $\chi_0=2B_0m_q$ ([mass-shift record](../../../../meson-mixing/alp-lo-meson-mass-shift-second-order.md)). $Y$ denotes the solution of Alda:2026nrd's first-order equation, $\mathbb M_{m_a^2}(Y)=C$, evaluated at the ALP mass.
- $\lambda_n$: generator of the meson $n$. For $\pi^0,\eta,\eta'$ it is the eigenvector of $\mathbb M_0$ on the diagonal neutral matrices. For $\pi^\pm$, $K^\pm$, $K^0$ it is $\lambda_1$, $\lambda_4$, $\lambda_6$ (one real component of the pair).
- Pole vectors $v^{(j)}$ are $K$-normalized, as in the [pole solver](../../../../meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-pole-solver.md). In the dressed frame $K=\mathbf 1$ and $D_{ja}=\epsilon\operatorname{tr}(\lambda_jC)$ does not depend on $s$.
- Loop functions ($\overline{\rm MS}$, dimensional regularization, $\int d^dk\,1=0$):
$$
A_0(m^2)=\frac{m^2}{16\pi^2}\Big[1-\ln\frac{m^2}{\mu^2}\Big],\qquad
\bar B(s;m_1,m_2)=-\int_0^1\!dz\,\ln\frac{z m_1^2+(1-z)m_2^2-z(1-z)s-i0}{\mu^2}.
$$
  **Notation flagged for approval:** $\bar B$ (the subtracted bubble function). The bar distinguishes it from the low-energy constant $B_0$.
- One-loop pole shift: $\delta m_n^2=v^{(n)\mathsf T}\big[\delta M^2(s_n)-s_n\,\delta K\big]v^{(n)}$ at $s_n=m_n^2$. Here $\delta M^2$ and $\delta K$ are the one-loop mass and kinetic terms, coming from tadpoles, bubbles and loop-induced one-point functions.

## Assumptions (all explicit)

1. **Flavour-diagonal dressed sources.** All $\widehat S_n$, $\widehat P_n$ commute with a diagonal $m_q$. Then the [first-order record](../../../../chiral-transformations/source-invariants/dressed-sources/alp-dressed-mass-first-order-invariants.md) gives $\widehat S_1=-i[c_q^{VS},m_q]=0$, $\widehat S_0=m_q$, $\widehat P_0=0$. This restricts the flavour structure but fixes no value.
2. LO $U(3)$ Lagrangian in the dressed frame of the [power-counting record](../loop-counting/alp-loop-power-counting.md): $\mathcal L=\frac{f_\pi^2}{4}\langle\partial U\partial U^\dagger\rangle+\frac{f_\pi^2B_0}{2}\langle\widehat MU^\dagger+U\widehat M^\dagger\rangle-\frac{f_\pi^2m_0^2}{4N_f}X_{TC}^2+\frac12(\partial a)^2-\frac12m_{a,0}^2a^2$, with $X_{TC}=-\operatorname{tr}\Phi/f_\pi-2\widehat\Theta(x)$. Using a general $\widehat\Theta(x)$ is the extension stated in the [dressed topological record](../../../../chiral-transformations/source-invariants/dressed-sources/alp-dressed-topological-invariant.md). No ALP self-interactions are included beyond those from the sources.
3. One loop, and fields up to fourth order. Counterterms of the next order are not included, so the "with logarithm" entries depend on $\mu$ and on the subtraction.
4. Each class is given at its leading order in $\epsilon$. Validity requires the meson pole to be far from the ALP pole, $|m_n^2-m_a^2|\gg\epsilon|\operatorname{tr}(\lambda_nC)|$, and $m_a^2$ to be far from the unmixed neutral masses (otherwise $Y$ is singular). Near a degeneracy the closed forms fail, and the exact finite-mixing expression computed by the script must be used.
5. The crossing is defined class by class: $|\delta m_n^2|_{\rm class}/m_n^2=1$.
6. Table only: $\lambda_{\pi^0}\simeq\lambda_3$, i.e. $\pi^0$–$\eta$–$\eta'$ mixing from $m_u\neq m_d$ is neglected. Then $m_{\pi^0}^2=B_0(m_u+m_d)$, with $B_0$ defined from this relation.

## Inputs

- Dressed sources and their invariance: [general redefinition record](../../../../chiral-transformations/source-invariants/finite-transformations/alp-general-chiral-redefinition-invariants.md), [source transformations](../../../../chiral-transformations/source-invariants/finite-transformations/alp-chiral-source-transformations.md), [first-order](../../../../chiral-transformations/source-invariants/dressed-sources/alp-dressed-mass-first-order-invariants.md) and [second-order](../../../../chiral-transformations/source-invariants/dressed-sources/alp-dressed-mass-second-order-invariants.md) mass invariants, [topological invariant](../../../../chiral-transformations/source-invariants/dressed-sources/alp-dressed-topological-invariant.md) and its [second order](../../../../chiral-transformations/source-invariants/dressed-sources/alp-dressed-topological-second-order.md).
- Tree-level invariance under constant axial redefinitions and $D_{na}(m_n^2)=\epsilon\operatorname{tr}(\lambda_nC)$: [mass-shift record](../../../../meson-mixing/alp-lo-meson-mass-shift-second-order.md). First- and second-order pole-vector components: [weak-mixing record](../../../../meson-mixing/alp-lo-weak-mixing-expansion.md).
- Tadpole combinatorics and the charged-meson self-energies: [power-counting record](../loop-counting/alp-loop-power-counting.md), section L2.
- Physical inputs (table only): $m_{\pi^0}=0.13498$ GeV, $f_\pi=0.09207$ GeV, $(m_u,m_d,m_s)=(2.2,4.7,95)$ MeV (defaults of `charge_model_parameters` in `scripts/derive_murayama_only_up_three_pseudoscalar_amplitude.py`), $B_0=m_{\pi^0}^2/(m_u+m_d)=2.6405$ GeV [@Gasser:1983yg], $\mu=0.77$ GeV, $m_a\in\{0.1,0.3,0.5,1,2\}$ GeV.

## Method

1. Build $\mathcal L$ to fourth order in the ten real fields ($P_u,P_d,P_s,a$ and the six charged components). Do this in the dressed frame and in the frame reached by the constant axial redefinition $g_R=g_L^\dagger=e^{ix\kappa}$, i.e. the change of variables $\widehat U=g^\dagger Ug^\dagger$. That frame has $C_A=\kappa$, $M'=e^{2ix\kappa}\widehat M$ and $\Theta'=\widehat\Theta-x\operatorname{tr}\kappa$. Derivative vertices are kept exactly (kinetic term with its $\partial x$ terms).
2. Extract $K$, $M^2$, the cubic vertices and the quartic vertices (with and without derivatives). Tadpoles are Wick contractions: $\langle F_KF_L\rangle=-\sum_rv^{(r)}_Kv^{(r)}_LA_0(m_r^2)$ and $\langle\partial F_K\cdot\partial F_L\rangle=-\sum_rm_r^2v^{(r)}_Kv^{(r)}_LA_0(m_r^2)$. Bubbles give $-\frac1{32\pi^2}V_{nrt}^2\bar B(s_n;m_r,m_t)$, and one-point functions $-V_{nnK}\langle F_K\rangle$.
3. Write $\delta m_n^2=\sum_rc_{n,r}A_0(m_r^2)+\sum_{r\le t}b_{n,rt}\bar B(m_n^2;m_r,m_t)$. The loop functions of different masses are independent functions, so every coefficient is separately an observable quantity. A class is defined as one such coefficient.
4. Expand each coefficient in $\epsilon$ using the first-order ALP component of the meson poles, $v^{(j)}_a=D_{ja}/(m_a^2-m_j^2)$, the second-order meson components $v^{(j)}_k=v^{(j)}_aD_{ka}/(m_k^2-m_j^2)$, and the meson components of the ALP pole, $\sum_kv^{(a)}_k\lambda_k=2\epsilon Y$.

## Result (a): surviving classes and closed formulas

Invariants used: $\widehat S_2$, $\widehat P_1$, $\widehat P_2$, $\widehat\Theta_1$, $\widehat\Theta_2$, together with $m_q$, $B_0$, $m_0^2$, $f_\pi$ and the pole masses. Meson poles not degenerate with the ALP.

**Class A: ALP line, coefficient of $A_0(m_a^2)$, order $\epsilon^2$.** It applies to every meson that is unmixed at leading order:
$$
\boxed{\;c_{n,a}=-\frac{\epsilon^2}{2f_\pi^2}\Big[B_0\big\langle\lambda_n^2\big(\widehat S_2+4Y\widehat P_1-4Y^2m_q\big)\big\rangle-m_a^2\big\langle[Y,\lambda_n]^2\big\rangle\Big]+O(\epsilon^4),\qquad \mathbb M_{m_a^2}(Y)=C .\;}
$$
For diagonal $C$ the solution is $Y_q=\big[C_q-\tfrac{m_0^2}{3}\operatorname{tr}Y\big]/(2B_0m_q-m_a^2)$ with $\operatorname{tr}Y=\dfrac{\sum_qC_q/(2B_0m_q-m_a^2)}{1+\frac{m_0^2}{3}\sum_q1/(2B_0m_q-m_a^2)}$.
- For $\pi^0,\eta,\eta'$ the commutator vanishes.
- For a meson of flavour $q\bar r$: $\langle\lambda^2X\rangle=X_{qq}+X_{rr}$ and $-\langle[Y,\lambda]^2\rangle=2(Y_{qq}-Y_{rr})^2$.
- Origin of the terms. The $\widehat S_2$ term is the contact tadpole with a pure ALP line. The terms with $Y$ come from the meson components $2\epsilon Y$ of the ALP-like pole, contracted with the meson–meson and meson–ALP quartic vertices (including kinetic vertices for the charged mesons).
- The same result is the pure-ALP tadpole in the frame where the ALP-like pole has no meson component at $O(\epsilon)$ (Alda:2026nrd's prescription $\kappa=-c_q^A-Y$). There, $\widehat S_2+4Y\widehat P_1-4Y^2m_q$ is the raw second-order scalar source, and the $m_a^2$ term comes from the axial current.
- Crossing: $\epsilon_*=f_\pi\Big[\dfrac{2m_n^2}{|A_0(m_a^2)|\,\big|B_0\langle\lambda_n^2(\widehat S_2+4Y\widehat P_1-4Y^2m_q)\rangle-m_a^2\langle[Y,\lambda_n]^2\rangle\big|}\Big]^{1/2}$. Size estimate: $|A_0|\to m_a^2/16\pi^2$.

**Class B: ALP pair, coefficient of $\bar B(m_n^2;m_a,m_a)$, order $\epsilon^4$.** Only for $\pi^0,\eta,\eta'$:
$$
\boxed{\;b_{n,aa}=-\frac{\epsilon^4}{32\pi^2f_\pi^2}\Big\langle\lambda_n\Big(B_0\widehat P_2-\frac{m_0^2}{N_f}\widehat\Theta_2\mathbf 1\Big)\Big\rangle^2+O(\epsilon^6).\;}
$$
- The vertex is $V_{naa}=\frac{\epsilon^2}{f_\pi}\langle\lambda_n(B_0\widehat P_2-\frac{m_0^2}{N_f}\widehat\Theta_2)\rangle$. The meson components of the internal ALP lines do not change it at this order, because $V_{naP}=V_{nPP'}=0$.
- It vanishes for charged mesons and $K^0$. It needs the CP-odd combinations $\widehat P_2$, $\widehat\Theta_2$ (decay $n\to aa$), so it is absent for a CP-conserving ALP.
- Crossing: $\epsilon_*=\Big[\dfrac{32\pi^2f_\pi^2m_n^2}{\langle\lambda_n(B_0\widehat P_2-\frac{m_0^2}{N_f}\widehat\Theta_2)\rangle^2|\bar B(m_n^2;m_a,m_a)|}\Big]^{1/4}$. Here $\bar B$ is evaluated on shell, so no evaluation point has to be chosen.

**Class C: meson lines changed by ALP–meson mixing, coefficient of $A_0(m_P^2)$, order $\epsilon^2$.** For a neutral external meson $n$ and an internal neutral meson $P$:
$$
c_{n,P}=\tfrac12\Big[V_{nnPP}\big(1-(v^{(n)}_a)^2-(v^{(P)}_a)^2\big)+2v^{(n)}_aV_{naPP}+2v^{(P)}_aV_{nnPa}+2\!\sum_{k\neq n}\!v^{(n)}_kV_{knPP}+2\!\sum_{k\neq P}\!v^{(P)}_kV_{nnPk}\Big]+O(\epsilon^4),
$$
- Vertices: $V_{ijkl}=\frac{B_0}{f_\pi^2}\langle\lambda_i\lambda_j\lambda_k\lambda_lm_q\rangle$ and $V_{ijka}=-\frac{\epsilon B_0}{f_\pi^2}\langle\lambda_i\lambda_j\lambda_k\widehat P_1\rangle$.
- Charged lines, pair $q\bar r$ with both real components: $c_{n,qr}=v^{(n)\mathsf T}\Sigma^{qr}(m_n^2)v^{(n)}$, with $v^{(n)}$ through $O(\epsilon^2)$ in $(P_u,P_d,P_s,a)$ and the nonzero entries
$$
\Sigma^{qr}_{qq}=-\frac{s-2B_0m_q}{3f_\pi^2},\quad
\Sigma^{qr}_{qr}=\frac{s+2B_0(m_q+m_r)}{3f_\pi^2},\quad
\Sigma^{qr}_{qa}=-\frac{2\sqrt2\,\epsilon B_0}{3f_\pi^2}\Big(\widehat P_{1,qq}+\tfrac12\widehat P_{1,rr}\Big),
$$
  plus the same with $q\leftrightarrow r$. These are the charged self-energies of the power-counting record (with $\delta M^2=-\Sigma$ there), rewritten with $\widehat P_1$. They contain the kinetic vertices.
- The $O(\epsilon^0)$ part is the pure meson loop. It does not depend on $\epsilon$ and gives no crossing. The $O(\epsilon^2)$ part scales as $\epsilon^2[\operatorname{tr}(\lambda C)/(m_a^2-m_P^2)]^2m_P^2/(16\pi^2f_\pi^2)$. It does not factorize into a loop factor times one invariant, so it has no entry in the table.

**Absent at these orders:** the ALP–meson bubble (it would need $V_{naP}\propto$ a diagonal $\widehat S_1=0$), the meson–meson bubble ($V_{nPP'}=0$), and loop-induced one-point functions (they start beyond $O(\epsilon^2)$, numerically $\ge\epsilon^4$).

## Result (b): invariance check

- **Constant axial redefinition, explicit one loop, all orders in $\epsilon$.** For every pole and every loop function, the coefficients agree between the dressed frame and the $\kappa$ frame. This was tested at $\epsilon=1/3$, which is finite mixing, with largest difference $10^{-14}$ relative. The $\kappa$-frame tadpoles use the full momentum-dependent vertices. Exact rational check: the pulled-back self-energy differs from the dressed one by $D(s)X+X^{\mathsf T}D(s)$, where $X$ is linear in the charged-meson tadpoles only. The neutral-line part is exactly congruent, and the cubic vertices are exactly congruent. Since $D(m_n^2)v^{(n)}=0$, the difference vanishes on shell. The off-shell entries differ (7 nonzero entries, nonzero at $s=0$).
- **Second-order generator** $G=\kappa_2x^2/2$ (neutral sector, tree vertex only). The momentum-dependent $P_qaa$ vertex changes by $[D(s)w]_q$ with $w\propto\kappa_2$. Hence class B is unchanged on shell, while the off-shell vertex changes.
- **Does the tadpole-ratio class (a), $2c^A_q(c^A_q-2c^P_q)A_0(m_a^2)/f_a^2$, survive? No, not as defined.**
  1. The raw coefficient is not invariant. The $\kappa$ frame built from the dressed frame has $c^A=\kappa$, $c^P=\kappa-c^{AP}$ and raw $S_2=\widehat S_2-4\kappa\widehat P_1-4\kappa^2m_q\neq0$. Its linear-spurion combination $-4m_qc^A(c^A-2c^P)=4\kappa\widehat P_1+4\kappa^2m_q$ depends on $\kappa$ and differs from $\widehat S_2$. The invariant is $\widehat S_2=S_2-4m_qc^A(c^A-2c^P)$ ([second-order record](../../../../chiral-transformations/source-invariants/dressed-sources/alp-dressed-mass-second-order-invariants.md)). It equals the linear-spurion combination only in frames with $S_2=0$.
  2. The off-shell entry $\delta M^2_{qq}$ is frame dependent (item above).
  3. Even its invariant replacement ($\widehat S_2$ alone) is not the full on-shell coefficient of $A_0(m_a^2)$ at order $\epsilon^2$. The terms with $Y$ in class A are of the same order.

  What survives is class A, the coefficient of $A_0(m_a^2)$ in the pole shift, with the boxed invariant coefficient. The old crossing $4\pi f_\pi/(\sqrt2\,c\,m_a)$ corresponds to the representative $\widehat S_2=-4c^2m_q$, $Y=0$.

## Result (c): numerical table (check only)

Neutral pion, $\lambda_3$ approximation (assumption 6). Then
- class A: $|\delta m^2/m^2|=\epsilon^2|A_0(m_a^2)|\,\mathcal I_A/(2f_\pi^2)$, with the dimensionless invariant combination $\mathcal I_A=|\sum_{q=u,d}(\widehat S_2+4Y\widehat P_1-4Y^2m_q)_{qq}|/(m_u+m_d)$;
- class B: $|\delta m^2/m^2|=\epsilon^4m_{\pi^0}^2\mathcal I_B^2|\bar B|/(32\pi^2f_\pi^2)$, with $\mathcal I_B=|\widehat P_{2,uu}-\widehat P_{2,dd}|/(m_u+m_d)$ ($\widehat\Theta_2$ drops because $\operatorname{tr}\lambda_3=0$).

The symbols $\mathcal I_A$, $\mathcal I_B$ are only labels for the table columns (**flagged for approval** if they are to be used elsewhere). The table lists $\epsilon_*\sqrt{\mathcal I}$, i.e. the crossing for unit invariant combination. L = $\overline{\rm MS}$ logarithms kept, $\mu=0.77$ GeV; S = size estimate ($|A_0|\to m_a^2/16\pi^2$, $|\bar B|\to1$).

| $m_a$ [GeV] | $A_0(m_a^2)$ [GeV$^2$] | A: $\epsilon_*\sqrt{\mathcal I_A}$ (L) | A (S) | $\bar B(m_{\pi^0}^2;m_a,m_a)$ | B: $\epsilon_*\sqrt{\mathcal I_B}$ (L) | B (S) |
|---|---|---|---|---|---|---|
| 0.1 | $+0.000322$ | 7.258 | 16.36 | $+4.4625$ | 2.395 | 3.482 |
| 0.3 | $+0.001644$ | 3.211 | 5.454 | $+1.9197$ | 2.958 | 3.482 |
| 0.5 | $+0.002950$ | 2.397 | 3.272 | $+0.8758$ | 3.599 | 3.482 |
| 1.0 | $+0.003022$ | 2.368 | 1.636 | $-0.5197$ | 4.101 | 3.482 |
| 2.0 | $-0.023026$ | 0.858 | 0.818 | $-1.9083$ | 2.962 | 3.482 |

- Inputs give $B_0=2.6405$ GeV. Class A does not depend on $B_0$ or the quark masses, because they cancel against $m_{\pi^0}^2=B_0(m_u+m_d)$.
- $A_0(m_a^2)=0$ at $m_a=\sqrt e\,\mu=1.27$ GeV, and $\bar B(m_{\pi^0}^2;m_a,m_a)=0$ near $m_a\approx0.75$ GeV. These zeros come from the subtraction and are not physical suppressions.
- Consistency with the [physical-input record](alp-loop-crossing-physical-inputs.md): there, class (i) used the representative $c^P=0$, $c^A=1$. That corresponds to $\mathcal I_A=4$ and $Y=0$, and dividing the A columns by 2 reproduces its entries (e.g. $0.818$ S and $1.18$ L at $m_a=1$ GeV). That record is representative dependent.

## Verification

[`test_derive_alp_invariant_loop_crossing.py`](../../../../../../scripts/test_derive_alp_invariant_loop_crossing.py) has 20 tests:
- the exact sparse builder equals the sympy builder;
- no linear terms, no derivative cubic vertices and no odd-derivative quartic terms in either frame;
- charged masses equal $B_0(m_q+m_r)$;
- pole masses agree between frames;
- on-shell coefficients agree between frames per loop function at $\epsilon=1/3$;
- exact $D X+X^{\mathsf T}D$ structure, off-shell entries differ, cubic vertices congruent;
- $\kappa_2$ vertex identity with $w\propto\kappa_2$;
- non-invariance of the linear-spurion combination;
- class A closed form against the exact coefficient (neutral and the three charged pairs), relative residual $\propto\epsilon^2$; symbolic reduction of the charged dressed-vertex result to the pole-aligned form, residual $0$;
- class B against exact, relative residual $\propto\epsilon^2$;
- class C against exact, residual $\propto\epsilon^4$;
- one-point terms scale at least as $\epsilon^4$;
- the bubble function against its closed form and its imaginary part above threshold;
- the zero of $A_0$;
- the table crossings solve ratio $=1$.

Test values of the sources are rationals used only for checks.

## Reproducible calculation

[`derive_alp_invariant_loop_crossing.py`](../../../../../../scripts/derive_alp_invariant_loop_crossing.py) (sympy for the Lagrangian, the vertices and the symbolic checks; numpy/scipy for the pole vectors and loop functions).

```sh
.venv/bin/python scripts/derive_alp_invariant_loop_crossing.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_invariant_loop_crossing.py
```
