---
title: Vertices of the leading Lagrangian and allowed loop topologies in the two-point functions
status: derived
created: 2026-11-06
updated: 2026-11-06
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds, Herrera-Siklody:1996tqr]
provenance: analytical-derivation
---

# Vertices and loop topologies

## Assumptions and inputs

- The leading Lagrangian is $\mathcal L^{(0)}=\frac{f_\pi^2}{4}\langle u_\mu u^\mu\rangle+\frac{f_\pi^2}{4}\langle\chi_+\rangle-\frac{f_\pi^2m_0^2}{4N_f}X_{TC}^2$ of [lagrangian.md](../lagrangian.md) and [the topological-source extract](../../extracts/Alda:2026nrd/lo-topological-action.md) [@Alda:2026nrd; @Kaiser:1998ds].
- The CCWZ block is $u=\exp(i\phi)$ with $\phi=\Phi/(2f_\pi)$, $\Phi=\sum_P\lambda_PP$, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$ and $u^2=U$. The vielbein $u_\mu$ and the spurions $\chi_\pm$ are the all-order expressions of [Appendix A of Alda:2026nrd](../../extracts/Alda:2026nrd/appendix-a-expansion.md): $u_\mu$ at linear order in the ALP with $c_q^A$ and $c_q^V$, and $\chi_S$, $\chi_P$ with $c_q^S$, $c_q^P$. The gluonic coupling $c_G$ enters through $X_{TC}$.
- **Assumption (project scope):** the ALP enters at first order in $1/f_a$. There is no CP-violating source and no constant $\theta$ term. The flavour matrices $c_q^{A,V,S,P}$, $\chi_0=2B_0m_q$ and $\Phi$ do not commute.
- The loop-class definitions are in [loop-classes.md](../counting/loop-classes.md) (classes G-LO and G-ALP).

The script is [`scripts/derive_nlo_loop_vertices.py`](../../../scripts/derive_nlo_loop_vertices.py):

```bash
.venv/bin/python scripts/derive_nlo_loop_vertices.py          # all checks
.venv/bin/python scripts/derive_nlo_loop_vertices.py --cubic  # vertices with three fields
```

## Method

1. Expand $\sin\mathcal C_\phi/\mathcal C_\phi$, $\cos\mathcal C_\phi$, $\sin\mathcal C_\phi$, $\cos\mathcal A_\phi$ and $\sin\mathcal A_\phi$ in the formulas of Appendix A to total field order four. Each map is applied to the noncommuting matrices.
2. Insert them in $\mathcal L^{(0)}$. Since $\log u=i\phi$ exactly, $X_{TC}=2i\,(\langle\log u\rangle+ic_Ga/f_a)$ is linear in $\Phi$ and $a$.
3. Keep terms by field content and at most one power of $a$. The coefficient of each class is a sum of traces of ordered products of $\Phi$, $\partial_\mu\Phi$, $\chi_0$ and the source matrices.
4. Cross-check: the quadratic terms $\Phi^2$ and $a\Phi$ reproduce [quadratic-lagrangian.md](../quadratic-lagrangian.md) (tests of the script).

## Result

| Source | $\Phi^3$ | $a\Phi^2$ | $\Phi^4$ | $a\Phi^3$ |
|---|---|---|---|---|
| none | $0$ | — | $\neq0$ | — |
| $c_q^A$ | — | $0$ | — | $\neq0$ |
| $c_q^V$ (generic) | — | $\neq0$ | — | $0$ |
| $c_q^S$ (generic) | — | $\neq0$ | — | $0$ |
| $c_q^P$ | — | $0$ | — | $\neq0$ |
| $c_G$ | — | $0$ | — | $0$ |

The $a\Phi^2$ term also vanishes for $c_q^V\propto\mathbb 1$ and for $[c_q^S,m_q]=0$. It is non-zero whenever the flavour matrix does not commute with $\Phi$ (for $c_q^V$) or with $m_q$ (for $c_q^S$), so it also exists for flavour-diagonal but non-universal $c_q^V$.

The vertices are:

$$
\mathcal L^{(0)}_{\Phi^4}=\frac{1}{f_\pi^2}\left[-\frac1{24}\langle\Phi^2\,\partial_\mu\Phi\,\partial^\mu\Phi\rangle+\frac1{24}\langle\Phi\,\partial_\mu\Phi\,\Phi\,\partial^\mu\Phi\rangle+\frac1{48}\langle\Phi^4\chi_0\rangle\right],
$$

$$
\mathcal L^{(0)}_{a\Phi^2}\Big|_{c_q^V}=\frac{i}{2f_a}\,\partial^\mu a\,\langle\Phi\,c_q^V\,\partial_\mu\Phi-\Phi\,\partial_\mu\Phi\,c_q^V\rangle,\qquad
\mathcal L^{(0)}_{a\Phi^2}\Big|_{c_q^S}=\frac{i\,a}{4f_a}\langle\Phi^2\,[\chi_0,c_q^S]\rangle,
$$

and

$$
\mathcal L^{(0)}_{a\Phi^3}\Big|_{c_q^A}=\frac{\partial^\mu a}{6f_\pi f_a}\langle\Phi^2\partial_\mu\Phi\,c_q^A+\Phi^2c_q^A\,\partial_\mu\Phi-2\,\Phi\,c_q^A\,\Phi\,\partial_\mu\Phi\rangle,
\qquad
\mathcal L^{(0)}_{a\Phi^3}\Big|_{c_q^P}=-\frac{a}{12f_\pi f_a}\langle\Phi^3\{\chi_0,c_q^P\}\rangle .
$$

## Loop topologies allowed in each two-point function

A one-loop two-point function has two kinds of graphs: a tadpole, with one vertex of four fields (two external, two joined in the loop), and a bubble, with two vertices of three fields each and two internal lines.

1. **$\Sigma_{ij}$ (meson–meson, order $a^0$).** The tadpole uses the $\Phi^4$ vertex. The bubble needs the $\Phi^3$ vertex, which is zero. $\Sigma_{ij}$ is tadpole-only.
2. **$\Sigma_{i\hat a}$ (meson–ALP, order $1/f_a$).** The tadpole uses the $a\Phi^3$ vertices of $c_q^A$ and $c_q^P$ (the vertices of $c_q^V$, $c_q^S$ and $c_G$ vanish). A bubble has one vertex with the external $\Phi_i$ and two internal lines, which must be mesons at first order in $1/f_a$. That vertex is $\Phi^3$, which is zero. $\Sigma_{i\hat a}$ is tadpole-only.
3. **The $a\Phi^2$ vertices of $c_q^V$ and $c_q^S$** therefore do not enter $\Sigma_{ij}$ or $\Sigma_{i\hat a}$ at first order in $1/f_a$. They enter a closed meson loop attached to the ALP alone (a term linear in $a$), and the bubble with two such vertices belongs to the ALP self-energy at order $1/f_a^2$. Both are outside the scope of this change.

Consequence: in $\Sigma_{ij}$ and $\Sigma_{i\hat a}$ the loop momentum integral does not depend on the external momentum, so the dependence on $s$ is polynomial (at most linear, from the derivative couplings) and no unitarity cut appears at this order. This is the structure used in the following records.

## Verification

`--test` and the default run exit with code 0. The checks are: the $\Phi^3$ coefficient is zero; the $a\Phi^2$ coefficient is zero for $c_q^A$, $c_q^P$, $c_G$, $c_q^V\propto\mathbb 1$ and $[c_q^S,m_q]=0$, and non-zero for generic $c_q^V$ and $c_q^S$; the quadratic terms agree with the quadratic-Lagrangian record.
