---
title: Second-order dressed scalar and pseudoscalar sources $\widehat S_2$, $\widehat P_2$
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Second-order dressed scalar and pseudoscalar sources $\widehat S_2$, $\widehat P_2$

## Convention

Taylor coefficients are derivatives at $x=0$, $x=a/f_a$, as in the [general record](alp-general-chiral-redefinition-invariants.md):
$$
F(x)=\sum_{n\ge0}\frac{x^n}{n!}F_n,\qquad F_n=\left.\frac{d^nF}{dx^n}\right|_{x=0}.
$$
For $n\le1$ this is also the coefficient of $x^n$; the operator multiplying $a^2/f_a^2$ has coefficient $F_2/2$.

## Inputs and assumptions

- Restricted partonic Lagrangian, transformation laws $M'=g_RMg_L^\dagger$, $C_\eta'=g_\eta C_\eta g_\eta^\dagger+ig_\eta\partial_xg_\eta^\dagger$, transport matrices $\partial_xW_\eta=iC_\eta W_\eta$, $W_\eta(0)=\mathbf1$, dressed source $\widehat M=W_R^\dagger MW_L$ with $M=S+iP$, and $\widehat S=(\widehat M+\widehat M^\dagger)/2$, $\widehat P=(\widehat M-\widehat M^\dagger)/(2i)$: [general record](alp-general-chiral-redefinition-invariants.md).
- Sources: $C_{R,n}=C_{V,n}+C_{A,n}$, $C_{L,n}=C_{V,n}-C_{A,n}$; $S_n$, $P_n$, $C_{V,n}$, $C_{A,n}$ Hermitian, noncommuting flavour matrices.
- First-order specialization (last section only): Alda:2026nrd parametrization $S_0=m_q$, $P_0=0$, $S_1=i[c_q^S,m_q]$, $P_1=\{c_q^P,m_q\}$, $C_{V,0}=c_q^V$, $C_{A,0}=c_q^A$ [@Alda:2026nrd], from the [first-order record](alp-dressed-mass-first-order-invariants.md).
- Same truncation as the general record (no dipoles, no two-ALP-derivative interactions); $g_\eta(0)=\mathbf1$. Under a residual constant rotation at $x=0$, $\widehat S_2,\widehat P_2$ transform covariantly.

## Derivation

1. The transport recursion of the general record gives
$$
W_{\eta,1}=iC_{\eta,0},\qquad W_{\eta,2}=iC_{\eta,1}-C_{\eta,0}^2 .
$$
2. The Leibniz rule for $W_R^\dagger(S+iP)W_L$ at second order gives $\widehat S_2$ and $\widehat P_2$ as the Hermitian and anti-Hermitian parts. Grouped with commutators and anticommutators, in terms of $S_n$, $P_n$, $C_{R,n}$, $C_{L,n}$:
$$
\boxed{
\begin{aligned}
\widehat S_2={}&S_2-\tfrac i4\Big[C_{R,0}+C_{L,0},\,4S_1-i[C_{R,0}+C_{L,0},S_0]+\{C_{R,0}-C_{L,0},P_0\}\Big]\\
&+\tfrac14\Big\{C_{R,0}-C_{L,0},\,4P_1-i[C_{R,0}+C_{L,0},P_0]-\{C_{R,0}-C_{L,0},S_0\}\Big\}
-\tfrac i2[C_{R,1}+C_{L,1},S_0]+\tfrac12\{C_{R,1}-C_{L,1},P_0\},\\[4pt]
\widehat P_2={}&P_2-\tfrac i4\Big[C_{R,0}+C_{L,0},\,4P_1-i[C_{R,0}+C_{L,0},P_0]-\{C_{R,0}-C_{L,0},S_0\}\Big]\\
&-\tfrac14\Big\{C_{R,0}-C_{L,0},\,4S_1-i[C_{R,0}+C_{L,0},S_0]+\{C_{R,0}-C_{L,0},P_0\}\Big\}
-\tfrac i2[C_{R,1}+C_{L,1},P_0]-\tfrac12\{C_{R,1}-C_{L,1},S_0\}.
\end{aligned}}
$$
3. The same result in $C_{V,n}$, $C_{A,n}$:
$$
\begin{aligned}
\widehat S_2={}&S_2-i\Big[C_{V,0},\,2S_1-i[C_{V,0},S_0]+\{C_{A,0},P_0\}\Big]
+\Big\{C_{A,0},\,2P_1-i[C_{V,0},P_0]-\{C_{A,0},S_0\}\Big\}
-i[C_{V,1},S_0]+\{C_{A,1},P_0\},\\
\widehat P_2={}&P_2-i\Big[C_{V,0},\,2P_1-i[C_{V,0},P_0]-\{C_{A,0},S_0\}\Big]
-\Big\{C_{A,0},\,2S_1-i[C_{V,0},S_0]+\{C_{A,0},P_0\}\Big\}
-i[C_{V,1},P_0]-\{C_{A,1},S_0\}.
\end{aligned}
$$
These nested forms have the fewest terms among the groupings tested by the script; this is not a proof that no shorter form exists. Both are Hermitian for Hermitian inputs. They contain products and nested commutators of the order-0 and order-1 coefficients; they are not same-order differences, as anticipated in the general record.

## Specialization to real quark masses

Assumption: $S_0=m_q$ Hermitian (in practice real diagonal) and $P_0=0$, with $S_1$, $P_1$, $S_2$, $P_2$, $C_{V,0}$, $C_{A,0}$, $C_{V,1}$, $C_{A,1}$ generic Hermitian. The general forms above reduce to
$$
\boxed{
\begin{aligned}
\widehat S_2={}&S_2-2i[C_{V,0},S_1]+2\{C_{A,0},P_1\}-i[C_{V,1},m_q]-[C_{V,0},[C_{V,0},m_q]]-\{C_{A,0},\{C_{A,0},m_q\}\},\\
\widehat P_2={}&P_2-2i[C_{V,0},P_1]-2\{C_{A,0},S_1\}-\{C_{A,1},m_q\}
+i\Big([C_{V,0},\{C_{A,0},m_q\}]+\{C_{A,0},[C_{V,0},m_q]\}\Big).
\end{aligned}}
$$
In $C_{R,n}$, $C_{L,n}$ (equivalent, but longer):
$$
\begin{aligned}
\widehat S_2={}&S_2-i[C_{R,0}+C_{L,0},S_1]+\{C_{R,0}-C_{L,0},P_1\}-\tfrac i2[C_{R,1}+C_{L,1},m_q]\\
&-\tfrac14[C_{R,0}+C_{L,0},[C_{R,0}+C_{L,0},m_q]]-\tfrac14\{C_{R,0}-C_{L,0},\{C_{R,0}-C_{L,0},m_q\}\},\\
\widehat P_2={}&P_2-i[C_{R,0}+C_{L,0},P_1]-\{C_{R,0}-C_{L,0},S_1\}-\tfrac12\{C_{R,1}-C_{L,1},m_q\}\\
&+\tfrac i2[C_{R,0}^2-C_{L,0}^2,m_q]-i\left(C_{R,0}m_qC_{L,0}-C_{L,0}m_qC_{R,0}\right).
\end{aligned}
$$
The cross term of $\widehat P_2$ can be rewritten with the Leibniz rule $[X,\{Y,Z\}]=\{[X,Y],Z\}+\{Y,[X,Z]\}$, e.g.
$$
[C_{V,0},\{C_{A,0},m_q\}]+\{C_{A,0},[C_{V,0},m_q]\}=\{[C_{V,0},C_{A,0}],m_q\}+2\{C_{A,0},[C_{V,0},m_q]\},
$$
but no rewriting tested is shorter. It contains six distinct matrix products, while a single nested bracket of $C_{V,0}$, $C_{A,0}$, $m_q$ contains at most four, so it cannot be written as one bracket. A diagonal $m_q$ does not shorten it, since $m_q$ still does not commute with $C_{V,0}$, $C_{A,0}$ in general.

## Specialization to the Alda:2026nrd first-order parametrization

Insert the first-order inputs listed above and keep $C_{V,1}$, $C_{A,1}$, $S_2$, $P_2$ generic:
$$
\begin{aligned}
\widehat S_2&=S_2-i[C_{V,1},m_q]+\big[c_q^V,[2c_q^S-c_q^V,m_q]\big]+\big\{c_q^A,\{2c_q^P-c_q^A,m_q\}\big\},\\
\widehat P_2&=P_2-\{C_{A,1},m_q\}-i\big[c_q^V,\{2c_q^P-c_q^A,m_q\}\big]-i\big\{c_q^A,[2c_q^S-c_q^V,m_q]\big\},
\end{aligned}
$$
with $2c_q^S-c_q^V=c_q^S-c_q^{VS}$ and $2c_q^P-c_q^A=c_q^P-c_q^{AP}$. At fixed $c_q^{VS}$, $c_q^{AP}$ the result still depends on $c_q^S$, $c_q^P$: the first-order data do not enter only through the first-order invariants. This dependence is compensated by the transformation of the raw second-order coefficients $S_2$, $P_2$, $C_{V,1}$, $C_{A,1}$, which are not invariant themselves.

## Reproducible calculation and checks

[`derive_alp_second_order_invariants.py`](../../../scripts/derive_alp_second_order_invariants.py) derives steps 1–3 and both specializations with noncommuting SymPy symbols and verifies every grouped form against the expanded Hermitian and anti-Hermitian parts, symbolically and with exact rational $3\times3$ matrices. Agreement of $\widehat M_2$ with the transport recursion of `derive_alp_invariant_coefficient_functions.py` is a consistency check, not an independent one. Independent checks:

- invariance of the closed forms using only the transformation laws of $S_n$, $P_n$, $C_{\eta,n}$ (not the transport derivation), with exact rational $3\times3$ matrices and $g_\eta$ built as Cayley transforms $(\mathbf1-iH)^{-1}(\mathbf1+iH)$, $H=xH_1+x^2H_2$, $[H_1,H_2]\neq0$, $\operatorname{tr}H_1\neq0$; the raw $S_2$ changes while $\widehat S_2$, $\widehat P_2$ do not;
- constant connections: comparison with $e^{-ixC_R}M(x)e^{ixC_L}$;
- numerical integration of the transport equation with `mpmath` (30 digits) for $C_\eta$ linear in $x$; relative deviation $1.4\times10^{-16}$ (floating-point check, tolerance $10^{-12}$);
- first-order regression: $\widehat S_1=-i[c_q^{VS},m_q]$, $\widehat P_1=-\{c_q^{AP},m_q\}$ of the [first-order record](alp-dressed-mass-first-order-invariants.md);
- vacuum-angle identity $\widehat\Theta+\tfrac12\arg\det\widehat M=\Theta+\tfrac12\arg\det M$ at $O(x^2)$, linking this record to the [second-order topological record](alp-dressed-topological-second-order.md);
- Hermiticity; unit tests confirm that perturbed compact forms are rejected.

These are series checks with fixed seeds, not all-orders proofs; the all-orders statement follows from $W_\eta'=g_\eta W_\eta$ in the general record.

```sh
.venv/bin/python scripts/derive_alp_second_order_invariants.py
.venv/bin/python -m unittest scripts.test_derive_alp_second_order_invariants
```
