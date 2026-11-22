---
title: Power counting of the dressed-frame vertices in the invariants, per loop class of an on-shell ALP amplitude
status: mixed
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd, Manohar:1983md, Gasser:1983yg]
provenance: analytical-derivation
---

# Power counting of the dressed-frame vertices in the invariants, per loop class

This record is task 3.1 of the OpenSpec change `estimate-loop-tree-crossing-finite-mixing`. It restates the vertex formula of section L1 of the [power-counting record](alp-loop-power-counting.md) with the invariants $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$ only, and gives, for each one-loop class of the observable class of the [selection record](alp-loop-observable-selection-a-to-3pi0.md) (on-shell amplitude with the ALP-like pole $A$ at first order in $1/f_a$ at tree level, e.g. $A\to P_{n_1}P_{n_2}P_{n_3}$ of the [tree record](alp-finite-mixing-three-meson-tree.md)), the power of $\epsilon$, the loop factor and the invariant contraction. No loop integral is evaluated. The choice of mass in the loop factor and the naturalness assumption are task 3.2; cancellations among classes are task 3.3.

## Conventions and units

- Units GeV. Notation of [@Alda:2026nrd]: $U=u^2=\exp(i\Phi/f)$, $\Phi=P_i\lambda_i$, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, $\langle\cdot\rangle$ the flavour trace, $N_f=3$, $x=a/f_a$, $\epsilon=f_\pi/f_a$, $f\equiv f_\pi$.
- Dressed frame ([derivative-coupling removal](alp-dressed-frame-derivative-coupling-removal.md)), $K=\mathbf 1$:
$$
\mathcal L=\frac{f^2}{4}\langle u_\mu u^\mu\rangle+\frac{f^2}{4}\langle\chi_+\rangle-\frac{f^2m_0^2}{4N_f}\Big(\frac{\langle\Phi\rangle}{f}+2\widehat\Theta(x)\Big)^2+\frac12(\partial a)^2-\frac12m_{a,0}^2a^2,
$$
  $\chi=2B_0\widehat M(x)$, $\widehat M(x)=\sum_n\frac{x^n}{n!}(\widehat S_n+i\widehat P_n)$, $\widehat\Theta(x)=\sum_n\frac{x^n}{n!}\widehat\Theta_n$.
- $\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1$, $\mathbb M_s(X)=\frac12\{\chi_0,X\}+\frac{m_0^2}{N_f}\langle X\rangle\mathbf 1-sX$, $\chi_0=2B_0\widehat S_0$, $Y=\mathbb M^{-1}_{m_A^2}(\mathcal C)$, $N(s)=1+2\epsilon^2\langle\mathbb M^{-1}_s(\mathcal C)^2\rangle$ ([operator-inverse record](alp-lo-finite-mixing-operator-inverse.md), [tree record](alp-finite-mixing-three-meson-tree.md)).
- Vertices are derivatives of $\mathcal L$ at zero fields, amputated, without the factor $i$, all momenta incoming. $\lambda_n$ is the generator of the unmixed meson $n$ ([invariant crossing record](alp-invariant-loop-crossing-meson-poles.md)).
- An *ALP end* is one end of a propagator of the ALP-like pole $A$ at a vertex, or the external $A$ leg.

## Assumptions

1. LO $U(3)$ Lagrangian above, one loop, stable real poles, $A$ not degenerate with any meson pole. No pure-ALP potential besides $\frac12m_{a,0}^2a^2$ (as in the tree record).
2. $\widehat S_0=m_q$ diagonal with three different entries; $\widehat P_0=0$, $\widehat\Theta_0=0$. All invariants are Hermitian matrices that commute with the quark charge (a $d$-$s$ entry allowed); this is an assumption. No value, representative or relation among invariants or quark masses is used; $\langle\lambda_3\widehat P_1\rangle$ is not set to zero.
3. $\widehat S_1=-i[c_q^{VS},m_q]$ ([first-order record](alp-dressed-mass-first-order-invariants.md)). This is not an extra assumption: since $m_q$ is diagonal, $(\widehat S_1)_{ij}=-i(c_q^{VS})_{ij}(m_j-m_i)$, so **$\widehat S_1$ has zero diagonal entries** for every source. With assumption 2 its only entries are $(\widehat S_1)_{ds}$ and $(\widehat S_1)_{sd}=(\widehat S_1)_{ds}^*$, i.e. it changes strangeness by one unit and conserves charge, and $\widehat S_1=0$ when $(c_q^{VS})_{ds}=0$. In particular $\langle\widehat S_1\rangle=0$ follows, so the earlier flagged assumption of task 2.1 is no longer needed. $\widehat S_1$ is otherwise general.
4. On-shell equality of amplitudes between the dressed frame and other frames: as in section L6 of the [power-counting record](alp-loop-power-counting.md) (assumption beyond tree level).

## Inputs

- General vertex for $n_a$ ALP and $n_P$ meson legs in $\widehat M_{n_a}$: [power-counting record](alp-loop-power-counting.md), section L1 (reused; restated below in the invariants).
- Class A contraction $\widehat S_2+4Y\widehat P_1-4Y^2m_q$ and the term $m_a^2\langle[Y,\lambda_n]^2\rangle$ of the meson pole shift: [invariant crossing record](alp-invariant-loop-crossing-meson-poles.md) (reused, generalized below to non-commuting matrices and $\widehat S_1\neq0$).
- Pole vectors at finite $\epsilon$: $A$ has meson part $2\epsilon Y/\sqrt{N(m_A^2)}$ and ALP part $1/\sqrt{N(m_A^2)}$; tree vertex list and leading tree term: [tree record](alp-finite-mixing-three-meson-tree.md) (reused).
- Loop factor $1/(16\pi^2)$ per loop and $|A_0(m^2)|\to m^2/(16\pi^2)$: [@Manohar:1983md], [power-counting record](alp-loop-power-counting.md) (reused).

## Result 1: vertices in the invariants (new form of L1, verified)

For $n_a$ ALP legs and $n_P$ meson legs with generators $\lambda_1,\dots,\lambda_{n_P}$, the momentum-independent vertex is
$$
V_{n_a;\lambda_1\cdots\lambda_{n_P}}=\epsilon^{n_a}f^{2-n_P-n_a}\Big[\frac{(-1)^{\lfloor n_P/2\rfloor}B_0}{n_P!}\sum_{\sigma\in S_{n_P}}\big\langle X_{n_a}\lambda_{\sigma_1}\cdots\lambda_{\sigma_{n_P}}\big\rangle+T_{n_a,n_P}\Big],\qquad
X_{n_a}=\begin{cases}\widehat S_{n_a}, & n_P\ \text{even},\\ \widehat P_{n_a}, & n_P\ \text{odd},\end{cases}
$$
with the topological part
$$
T_{n_a,0}=-\frac{m_0^2}{N_f}\sum_{k=1}^{n_a-1}\binom{n_a}{k}\widehat\Theta_k\widehat\Theta_{n_a-k},\qquad
T_{n_a,1}=-\frac{m_0^2}{N_f}\widehat\Theta_{n_a}\langle\lambda_1\rangle,\qquad
T_{0,2}=-\frac{m_0^2}{2N_f}\langle\lambda_1\rangle\langle\lambda_2\rangle,
$$
$T=0$ otherwise, plus $-m_{a,0}^2$ at $(n_a,n_P)=(2,0)$. The kinetic term gives vertices with $n_a=0$ only (momentum dependent). Consequences:

- Each ALP leg of a vertex costs $\epsilon/f$, each meson leg $1/f$; every vertex is $\epsilon^{n_a}$ times an $\epsilon$-independent contraction.
- One meson leg: $V_{n_a;\lambda}=\epsilon^{n_a}f^{1-n_a}\big\langle\lambda\big(B_0\widehat P_{n_a}-\frac{m_0^2}{N_f}\widehat\Theta_{n_a}\mathbf 1\big)\big\rangle$, which is $\epsilon\langle\lambda\,\mathcal C\rangle$ at $n_a=1$.
- Cubic vertices exist in general: $V_{1;\lambda\lambda'}\propto\widehat S_1$, $V_{2;\lambda}\propto B_0\widehat P_2-\frac{m_0^2}{N_f}\widehat\Theta_2\mathbf 1$, $V_{3}\propto B_0\langle\widehat S_3\rangle-6\frac{m_0^2}{N_f}\widehat\Theta_1\widehat\Theta_2$. The statement of L1 that only $n_a+n_P$ even survives and that there are no cubic vertices holds only in the representative with constant diagonal $c_q^A$, $c_q^P$ (there $\widehat S_{\rm odd}=0$, $\widehat P_{\rm even}=0$).

## Result 2: powers of $\epsilon$ from ALP ends (verified at the vertex level)

- **End on the ALP-like pole.** Contracting one leg of any vertex with the pole vector of $A$ gives
$$
\frac{\epsilon}{\sqrt{N(m_A^2)}}\Big[\frac{1}{\epsilon}V(\dots,a)+V(\dots,2Y)\Big],
$$
  where $V(\dots,2Y)$ is the vertex with that leg replaced by the meson matrix $2Y$; the bracket does not depend on $\epsilon$. The ALP component and the meson components of $A$ contribute at the same order. Two ends on $A$ give $\epsilon^2/N(m_A^2)$ times an $\epsilon$-independent contraction.
- **End on a meson-like pole.** Its ALP component is $v^{(n)}_a=\epsilon\langle\lambda_n\mathcal C\rangle/(m_A^2-m_{n,0}^2)+O(\epsilon^3)$ at fixed $m_A^2$, so an ALP leg on a meson-like line costs $\epsilon$ times $\langle\lambda_n\mathcal C\rangle/(m_A^2-m_{n,0}^2)$; this is not small near $m_A^2\simeq m_{n,0}^2$.
- **Counting rule.** The power of $\epsilon$ of a diagram is its number of ALP ends (each propagator of $A$, internal or tree, has two), times $N(m_A^2)^{-1/2}$ per end on $A$. The leading tree term has one end (the external $A$), so it is $O(\epsilon)$. A one-loop diagram relative to the tree carries $\epsilon^{2k}$, with $k$ the number of additional propagators of $A$ plus the number of additional ALP-meson mixing insertions on meson-like lines. At finite $\epsilon$ each extra propagator of $A$ gives $\epsilon^2/N(m_A^2)$, $N(m_A^2)=1+2\epsilon^2\langle Y^2\rangle$, which does not grow without bound at large $\epsilon$.

## Result 3: contraction of a closed ALP-like line (verified)

Closed $A$ line on a quartic vertex with external meson generators $a,b$ (momenta $p,-p$; $A$ momenta $k,-k$), coefficient of $\epsilon^2/N(m_A^2)$:

- Potential part:
$$
-\frac{B_0}{2f^2}\Big\langle\big(\widehat S_2+2\{Y,\widehat P_1\}-Y^2m_q-m_qY^2-2Ym_qY\big)\{a,b\}\Big\rangle
+\frac{B_0}{3f^2}\Big\langle\big(\widehat P_1-\{Y,m_q\}\big)\big([a,[b,Y]]+[b,[a,Y]]\big)\Big\rangle
+\frac{B_0}{6f^2}\Big\langle m_q\big(\{[Y,[Y,a]],b\}+\{[Y,[Y,b]],a\}\big)\Big\rangle .
$$
  For commuting $a,b,Y$ only the first term survives, and for $a=b=\lambda$ it is $-\frac{B_0}{f^2}\langle\lambda^2(\widehat S_2+4Y\widehat P_1-4Y^2m_q)\rangle$, the class A combination. The first term alone, which is the second-order scalar source after the axial rotation $g_R=g_L^\dagger=e^{-ixY}$, is **not** the full result for non-commuting matrices (the script finds a nonzero residual with a $d$-$s$ entry).
- Kinetic part: $\frac{1}{3f^2}(p^2+k^2)\langle[Y,a][Y,b]\rangle$, no $p\cdot k$ term.
- On shell ($k^2=m_A^2$, $p\cdot k=0$, using $\mathbb M_{m_A^2}(Y)=\mathcal C$) the sum is
$$
-\frac{B_0}{2f^2}\Big\langle\big(\widehat S_2+2\{Y,\widehat P_1\}-Y^2m_q-m_qY^2-2Ym_qY\big)\{a,b\}\Big\rangle+\frac{m_A^2}{f^2}\langle[Y,a][Y,b]\rangle+\frac{1}{6f^2}\Big(\big\langle[Y,a][Y,p^2b-B_0\{m_q,b\}]\big\rangle+(a\leftrightarrow b)\Big);
$$
  the last bracket vanishes when $\frac12\{\chi_0,b\}=p^2b$ (external meson on its unmixed shell), which reproduces class A, including the term $m_A^2\langle[Y,\lambda]^2\rangle$ (checked for $\lambda_1$, also with a $d$-$s$ entry in $Y$).
- Closed $A$ line on a cubic vertex (one meson generator $\lambda$): $\frac1f\big\langle\lambda\big(B_0(\widehat P_2-2\{Y,\widehat S_1\})-\frac{m_0^2}{N_f}\widehat\Theta_2\mathbf 1\big)\big\rangle$. It contains only CP-odd invariants ($\widehat P_2$, $\widehat\Theta_2$, $\widehat S_1$). Because $\widehat S_1$ has only $d$-$s$ entries (assumption 3), its term is, for a diagonal $\lambda$ with entries $\lambda_u,\lambda_d,\lambda_s$,
$$
-\frac{2B_0}{f}\langle\lambda\{Y,\widehat S_1\}\rangle=-\frac{2B_0}{f}(\lambda_d+\lambda_s)\big[Y_{ds}(\widehat S_1)_{sd}+Y_{sd}(\widehat S_1)_{ds}\big],\qquad
Y_{ds}=\frac{\mathcal C_{ds}}{B_0(m_d+m_s)-m_A^2}=\frac{B_0(\widehat P_1)_{ds}}{B_0(m_d+m_s)-m_A^2},
$$
  since $\mathbb M_s$ acts entry by entry on off-diagonal entries for diagonal $\chi_0$ (verified). It is nonzero only if $\widehat P_1$ also has a $d$-$s$ entry. For $\lambda_{6,7}$ the term is $-\frac{2B_0}{f}(Y_{dd}+Y_{ss})\big[\lambda_{ds}(\widehat S_1)_{sd}+\lambda_{sd}(\widehat S_1)_{ds}\big]$ (verified). For the charged generators it vanishes.

## Result 3b: cubic vertex with one ALP leg and zero-diagonal $\widehat S_1$ (verified)

$V_{1;ab}=-\frac{\epsilon B_0}{2f}\langle\widehat S_1\{a,b\}\rangle$ is nonzero only when $\{a,b\}$ has a $d$-$s$ entry: exactly the pairs (diagonal generator, $\lambda_6$ or $\lambda_7$) and ($\lambda_1$ or $\lambda_2$, $\lambda_4$ or $\lambda_5$), i.e. a meson pair of zero charge and strangeness $\pm1$ (neutral diagonal meson with a $K^0$-like meson, or $\pi^\pm$ with $K^\mp$). It vanishes for two diagonal generators and for $\{\lambda_6,\lambda_6\}$, $\{\lambda_6,\lambda_7\}$, $\{\lambda_7,\lambda_7\}$. The vertex with the ALP leg replaced by the meson part $2Y$ of $A$ is a three-meson vertex, which is zero ($\widehat P_0=0$, $T_{0,3}=0$, kinetic term even in $\Phi$).

## Result 3c: closed meson line with $\widehat S_1$ (verified)

The meson propagator at order $\epsilon^0$, $\Delta=(k^2-M^2_{\rm unmixed})^{-1}$, is invariant under $\Phi\to e^{i\alpha}\Phi e^{-i\alpha}$ with $\alpha$ diagonal, because $\chi_0$ is diagonal and the $m_0^2$ term depends only on $\langle\Phi\rangle$. Hence the closed meson line $\sum_{ij}\Delta_{ij}\{\lambda_i,\lambda_j\}$ is a diagonal matrix (equivalently $\langle\widehat S_1\lambda_n^2\rangle=0$ for each unmixed meson $n$), and its trace with the zero-diagonal $\widehat S_1$ vanishes. A closed meson line attached to $V_{1;ab}$ therefore gives no $A$ one-point function at this order. Insertions of ALP-meson mixing on the closed line cost a further $\epsilon^2$ and can give products such as $(\widehat S_1)_{sd}(\widehat P_1)_{ds}$ (not checked).

## Result 4: loop classes of the amplitude

All entries are relative to the leading tree term, whose contraction is the quartic vertex with one end on $A$ ([tree record](alp-finite-mixing-three-meson-tree.md), "Leading term"). The loop factor is written as $m^2/(16\pi^2f^2)$; the mass $m$ listed is the scale that appears in the loop function or the vertices, and which one is used in the estimate is fixed in task 3.2.

| class | extra ALP ends | power of $\epsilon$ (relative to tree) | loop factor | invariant contraction |
|---|---|---|---|---|
| 1. meson lines only | 0 | $\epsilon^0$ | $m_P^2/(16\pi^2f^2)$ | meson vertices ($m_q$, kinetic) and the tree end on $A$; the usual chiral correction, independent of $\epsilon$ |
| 2. closed $A$ line at one vertex with two or more remaining meson legs (tadpole, also inside external-leg self-energies) | 2 | $\epsilon^2/N(m_A^2)$ | $|A_0(m_A^2)|/f^2\to m_A^2/(16\pi^2f^2)$ | Result 3 (quartic); on the six-leg vertex of $A\to PPP$: Result 1 with two legs contracted as in Result 2. The closed line on a cubic vertex leaves one meson leg and is class 6 |
| 3. one $A$ line in a bubble or triangle | 2 | $\epsilon^2/N(m_A^2)$ | $1/(16\pi^2f^2)$ times $s\sim m_A^2$ or the vertex scale | products of vertices with one end on $A$: quartic as in the tree term (no restriction); cubic $-\frac{B_0}{2f}\langle\widehat S_1\{a,b\}\rangle$ only for the pairs of Result 3b, so for neutral diagonal external legs the meson line next to $A$ at a cubic vertex is a $K^0$-like line (or a $\pi^\pm K^\mp$ pair if both meson legs are internal). Since the $\epsilon^0$ meson propagator conserves strangeness (Result 3c), that line must end on a second vertex with a $d$-$s$ entry, so the cubic part is quadratic in $d$-$s$ entries ($(\widehat S_1)_{sd}$ times $(\widehat S_1)_{ds}$, $(\widehat P_1)_{ds}$ or $Y_{ds}$) |
| 4. two $A$ lines (bubble, triangle) | 4 | $\epsilon^4/N(m_A^2)^2$ | $m^2/(16\pi^2f^2)$, $m^2\sim s$ | cubic with two ends on $A$ (Result 3, CP odd; its $\widehat S_1$ term needs $Y_{ds}$ for diagonal legs), squared; quartic with two ends (Result 3); cubic with one end on $A$ restricted as in class 3 |
| 5. meson loop with one extra $A$ propagator or mixing insertion on a meson-like line | 2 | $\epsilon^2$ | $m_P^2/(16\pi^2f^2)$ | products $\langle\lambda_k\mathcal C\rangle\langle\lambda_l\mathcal C\rangle/[(m_A^2-m_k^2)(m_A^2-m_l^2)]$ with meson vertices; not a single invariant (as class C of the invariant crossing record) |
| 6. meson one-point functions from a closed $A$ line | 2 | $\epsilon^2/N(m_A^2)$ | $m_A^2/(16\pi^2f^2)$ | cubic closed-line contraction of Result 3: $\widehat P_2$ and $\widehat\Theta_2$ terms, and for neutral diagonal mesons the $\widehat S_1$ term $\propto(\lambda_d+\lambda_s)\,\mathrm{Re}\big[Y_{ds}(\widehat S_1)_{sd}\big]$, nonzero only if $\widehat P_1$ has a $d$-$s$ entry. The $A$ one-point from a closed meson line vanishes at this order (Result 3c) |

Remarks:
- The shift of the mass of the external $A$ is absorbed into $m_{a,0}$ (test (c) of design D1). The wave-function factors of the external legs and the shifts of the meson poles that enter $m_n$ and $Y$ follow the counting rule and belong to the classes above.
- Classes 2, 3, 5 and 6 are all of relative order $\epsilon^2$; they differ in the mass in the loop ($m_A^2$ versus $m_P^2$) and in the contraction. This is input for task 3.3.
- $\widehat S_1\neq0$ only if $(c_q^{VS})_{ds}\neq0$ (assumption 3). With $\widehat S_1=0$ the $\widehat S_1$ terms of classes 3, 4 and 6 vanish; class 6 keeps its $\widehat P_2$, $\widehat\Theta_2$ terms.
- **Class 6 for neutral external mesons.** Strangeness does not remove it: the $\widehat S_1$ term survives through $Y_{ds}(\widehat S_1)_{sd}$ because the $A$ pole is not a strangeness eigenstate when $\widehat P_1$ has a $d$-$s$ entry, and the $\widehat P_2$, $\widehat\Theta_2$ terms do not involve $\widehat S_1$. All its terms are CP odd, so the class would vanish only under an additional CP-conservation assumption, which is not made here. Only the $A$ one-point part of the former class 6 vanishes (Result 3c); no class is empty.
- **Open (flagged):** how the one-point functions of class 6 are treated (kept as a class, or removed by re-aligning the vacuum at one loop) is not decided here.

## Superseded content (representative dependent)

In the [power-counting record](alp-loop-power-counting.md): the L1 sentence "only $n_a+n_P$ even survives: there are no cubic vertices", the closed forms with $c_q^A$, $c_q^P$ in the "Result table" and in the L3 table, and the regime numbers of L4 hold only for constant diagonal $c_q^A$, $c_q^P$. Their invariant replacements are Results 1–4 of this record. The powers of $f_\pi/f_a$ in that table for "ALP line, meson entry", "$aP$ bubble" and "$aa$ bubble" agree with classes 2, 3 and 4.

## Reused versus new

- Reused: the vertex formula of L1 in $\widehat M_{n_a}$; the loop factor; class A of the meson pole shift; the pole vectors and vertex list of the tree record.
- New: the form of L1 in $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$ including the topological vertices; the existence of cubic vertices for general invariants; the counting rule by ALP ends with the factor $N(m_A^2)^{-1/2}$; the non-commuting closed-line contraction and its kinetic part; the list of classes for the amplitude.

## Verification

[`scripts/derive_alp_dressed_vertex_power_counting.py`](../../../scripts/derive_alp_dressed_vertex_power_counting.py), tests in [`scripts/test_derive_alp_dressed_vertex_power_counting.py`](../../../scripts/test_derive_alp_dressed_vertex_power_counting.py) (30 tests, all residuals exactly 0, exact rationals, test inputs only). Test inputs use $\widehat S_1=-i[c_q^{VS},m_q]$ built from a test $c_q^{VS}$, with zero diagonal asserted: $\widehat S_1=0$ for diagonal inputs, only $d$-$s$ entries otherwise (the tree record's test $\widehat S_1$ is diagonal and is overridden here).

- Result 1 against the fourth-order Lagrangian of [`derive_alp_invariant_loop_crossing.py`](../../../scripts/derive_alp_invariant_loop_crossing.py): all 1000 monomials with $n_a+n_P\le4$ in the ten fields (diagonal invariants, hence $\widehat S_1=0$ in this comparison; $\widehat S_1$ is tested by the $d$-$s$ comparisons below). With a $d$-$s$ entry: 144 vertices against the general vertex code `derive_alp_lo_joint_vertices.vertex`, and all 1000 monomials against a direct expansion. The $S$ and $P$ parts equal the L1 code of [`derive_alp_loop_power_counting.py`](../../../scripts/derive_alp_loop_power_counting.py). The one-meson combination equals $\epsilon\langle\lambda\mathcal C\rangle$ at $n_a=1$, for $n_a\le3$ and all generators.
- Homogeneity: every monomial to sixth order is exactly $\epsilon^{n_a}$ ($\epsilon$ kept symbolic).
- Kinetic term free of the ALP; topological term without vertices with $n_P\ge2$ except the mass term.
- Sixth order (needed for tadpoles in $A\to PPP$): 1715 monomials of the potential against a direct Taylor expansion (fields $P_u,P_d,P_s,a,P_1,P_6,P_7$, $d$-$s$ entry). The six-meson kinetic vertex is not covered by L1 and is not checked.
- Cubic vertices: $V_{1;ab}=0$ for $\widehat S_1=0$; with the $d$-$s$ entry, all 45 generator pairs give exactly the ten nonzero pairs of Result 3b, equal to the general vertex code. $V_{2;\lambda}$, $V_3$ nonzero; all absent in the $c_q^A$, $c_q^P$ representative.
- Result 3c: $\sum_{ij}\Delta_{ij}\{\lambda_i,\lambda_j\}$ diagonal for the $9\times9$ unmixed mass matrix (with the $m_0^2$ term) at a test $k^2$; its trace with $\widehat S_1$ and $\langle\widehat S_1\lambda_n^2\rangle$ for all nine unmixed mesons are 0.
- $\widehat S_1$ term of the cubic closed line: the formula of Result 3 for $\lambda_3,\lambda_8,\lambda_0$ and for $\lambda_{6,7}$; zero when $\widehat P_1$ is diagonal (diagonal $\lambda$), nonzero when both $\widehat P_1$ and $\widehat S_1$ have $d$-$s$ entries; $Y_{ds}=\mathcal C_{ds}/(B_0(m_d+m_s)-m_A^2)$.
- Result 2: one and two ends on $A$ on vertices with three and four legs; $1/\sqrt N$ factors out exactly; series of $v^{(n)}_a$ for the three neutral poles and the $K^0$-like pole with a $d$-$s$ entry.
- Result 3: commuting class A form; failure of the rotated-source term alone and validity of the full matrix form on all 45 generator pairs with a $d$-$s$ entry; kinetic part from the general vertex code; on-shell form and its reduction to class A for $\lambda_1$; cubic closed-line form for all nine field-basis generators and all nine Gell-Mann generators, with the zero-diagonal $\widehat S_1$.

Not verified: the class-level powers in Result 4 are obtained by counting ends with the verified per-end factors; no diagram is evaluated. Existing suites still pass: `test_derive_alp_loop_power_counting.py` (20), `test_derive_alp_invariant_loop_crossing.py` (20), `test_derive_alp_finite_mixing_three_meson_tree.py` (15).

```sh
.venv/bin/python scripts/derive_alp_dressed_vertex_power_counting.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_dressed_vertex_power_counting.py
```
