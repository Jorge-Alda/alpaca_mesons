---
title: Vertex rules up to five legs for generic sources in the one-loop sigma to pi+ pi- estimate
status: derived
created: 2026-11-26
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Vertex rules up to five legs for generic $\widehat S_{r,s}$, $\widehat P_{r,s}$

Task 2.2 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama`. No loop is evaluated.

## Assumptions and inputs

- The Lagrangian, the vertex convention, and the assumptions are those of the [tree and legs record](murayama-sigma-pipi-loop-tree-and-legs.md). Units are GeV.
- The closed formulas are compared with [the estimate record](../../../couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md). The radial mass is taken from [the general-source radial record](../../../couplings-and-decays/finite-mixing/radial-decay-polar-cartesian-general-sources.md), Step 2: $m_\sigma^2=2\lambda f_a^2+B_0\epsilon^2(\langle\widehat S_{1,0}\rangle-2\langle\widehat S_{2,0}\rangle)$.
- $f_\pi$ does not depend on $\sigma$ (assumption). The anomalous singlet term gives no vertex, and $\widehat\Theta_1=0$.

## Result

For legs $x_1,\dots,x_n$, let $R$ be the radial legs, $A$ the legs that contribute their ALP component $\xi_x$, and $P$ the legs that contribute their meson matrix $\lambda_x$. Let $r=|R|$, $s=|A|$, $m=|P|$, and $X_{r,s}=\widehat S_{r,s}$ ($m$ even) or $\widehat P_{r,s}$ ($m$ odd). Then:
$$
V^{(n)}_{\rm pot}=B_0f_\pi^{2-n}\sum_{(R,A,P)}r!\,s!\,\epsilon^{r+s}\frac{(-1)^{m(m-1)/2}}{m!}\prod_{x\in A}\xi_x\sum_{\pi\in S_m}\big\langle X_{r,s}\lambda_{\pi_1}\cdots\lambda_{\pi_m}\big\rangle ,
$$
$$
V^{(4)}_{\rm kin}=\frac{1}{48f_\pi^2}\sum_{\pi\in S_4}\big(-k_{\pi_2}\cdot k_{\pi_4}\big)\big\langle[\lambda_{\pi_1},\lambda_{\pi_2}][\lambda_{\pi_3},\lambda_{\pi_4}]\big\rangle\quad(\text{four legs with meson matrices}),
$$
$$
V_{{\rm kin},a}(\sigma;x,y)=-\frac{2\epsilon}{f_\pi}\xi_x\xi_y\,k_x\cdot k_y,\qquad V_{{\rm kin},a}(\sigma,\sigma;x,y)=-\frac{2\epsilon^2}{f_\pi^2}\xi_x\xi_y\,k_x\cdot k_y ,
$$
$$
V(\sigma,\sigma,\sigma)=-6\lambda f_a+\frac{6B_0\epsilon^3}{f_\pi}\langle\widehat S_{3,0}\rangle
=-\frac{3\epsilon}{f_\pi}\Big[m_\sigma^2-\epsilon^2B_0\big(\langle\widehat S_{1,0}\rangle-2\langle\widehat S_{2,0}\rangle\big)\Big]+\frac{6B_0\epsilon^3}{f_\pi}\langle\widehat S_{3,0}\rangle ,
$$
$$
V(\sigma^4)=-6\lambda+\frac{24B_0\epsilon^4}{f_\pi^2}\langle\widehat S_{4,0}\rangle,\qquad V(\sigma^5)=\frac{120B_0\epsilon^5}{f_\pi^3}\langle\widehat S_{5,0}\rangle .
$$
The $\langle\widehat S_{r,0}\rangle$ terms are the $m=s=0$ case of $V_{\rm pot}$. $\lambda$ is eliminated through $m_\sigma^2$ with $\lambda f_a=\epsilon[m_\sigma^2-\epsilon^2B_0(\langle\widehat S_{1,0}\rangle-2\langle\widehat S_{2,0}\rangle)]/(2f_\pi)$. There is no cubic or quintic meson kinetic vertex, and no vertex coupling $\sigma$ to the meson kinetic term. The six-meson kinetic vertex is not needed, because no graph of $\sigma\to\pi^+\pi^-$ at one loop has a vertex with more than five legs. A vertex is the sum of all the parts that apply to its legs. For the generic source letters, $\widehat P_{0,0}=\widehat P_{1,0}=0$ are applied as matrix identities, and $\langle\widehat S_{0,1}\rangle=\langle\widehat S_{1,1}\rangle=0$ as trace identities.

## Comparison with the estimate record

The closed formulas of [the estimate record](../../../couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md), with the overall factor $f_\pi^{2-n}$ of $V^{(n)}_{\rm pot}$, the factor $1/f_\pi^2$ of $V^{(4)}_{\rm kin}$ and $C(\sigma)=2\epsilon/f_\pi$, $C(\sigma,\sigma)=2\epsilon^2/f_\pi^2$ in $V_{{\rm kin},a}$, agree with the direct expansion for every leg pattern with $n=3,4,5$. The vertex function of [`derive_radial_loop_trace_contractions.py`](../../../../../../scripts/derive_radial_loop_trace_contractions.py) includes $f_\pi^{2-n}$.

The powers of $\epsilon$ of each vertex by field end (radial, ALP field, meson) at large $\epsilon$, with $\lambda$ fixed, are in the [large-$\epsilon$ vertex record](../asymptotic-behavior/murayama-sigma-pipi-loop-large-epsilon-vertices-tree.md), and the leg factors that attach to the field ends in the [large-$\epsilon$ input record](../asymptotic-behavior/murayama-sigma-pipi-loop-large-epsilon-inputs.md).

## Verification

- **Construction 1** is the direct noncommutative expansion of $\mathcal L$, with $\Phi=\sum_j\phi_jx_j$, $a=\sum_j\phi_j\xi_j$ and $\sigma=\sum\phi_j$ over the radial legs, and extraction of the coefficient of $\phi_1\cdots\phi_n$. **Construction 2** is the closed rules above, both as words and as compact trees. They agree as exact identities of cyclic words for every pattern $n=3,4,5$ with $r=0,\dots,n$ radial legs and generic pseudoscalar legs, where each pseudoscalar leg has both $\lambda_x$ and $\xi_x$.
- After the up-only substitution $\widehat S_{r,s}+i\widehat P_{r,s}=m_u(-i)^s P_u/s!$ for $r\le1$ (zero for $r\ge2$), the cubic and quartic rules agree with the record formulas. $V(\sigma;x,y)$ of the record agrees for generic sources, and $V(\sigma,\sigma;x,y)$ for the up-only sources. The record [`murayama-loop-generator-trace-products.md`](../../loop-estimates/diagram-traces/murayama-loop-generator-trace-products.md) does not write the vertex rules out. It uses the same vertex function as the estimate script, so it is checked through the formulas of the estimate record.
- Dimensions: under $X\to c^{d_X}X$, with $d=1$ for $B_0$, $f_\pi$, $m_\sigma$ and each source letter, and $d=2$ for each $k_i\cdot k_j$, every vertex scales as $c^{4-n}$.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_vertices.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_vertices
```

[Script](../../../../../../scripts/derive_sigma_pipi_loop_vertices.py), [tests](../../../../../../scripts/test_derive_sigma_pipi_loop_vertices.py), [shared algebra](../../../../../../scripts/sigma_pipi_loop_algebra.py).
