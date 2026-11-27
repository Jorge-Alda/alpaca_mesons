---
title: Cartesian frame - up-only sources, leading large-epsilon size estimate of the one-loop sigma to pi+ pi- amplitude
status: derived
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Cartesian frame: up-only sources, $M_1^{\rm est}/M_0$ at large $\epsilon$

Task 4.3 (up-only part) of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama-cartesian`. **No loop integral, scalar loop function, loop coefficient or loop amplitude is computed.** $R$ is a size ratio, not a computed correction. Units: GeV for $B_0$, $m_u$, $f_\pi$, $m_1$, $m_P$, $m_{\pi^\pm}$; $\epsilon=f_\pi/f_a$ is dimensionless.

**Caveat.** The estimate is an order of magnitude per graph class. Graphs with the same $\epsilon$ scaling can cancel in the exact loop, and these cancellations can differ between the Cartesian and the polar parametrization. Agreement with the polar estimate, including its up-only prefactor $-21/(64\pi^2)$, is not expected, and neither estimate is authoritative. The comparison with the polar estimate is in the [polar comparison record](cartesian-sigma-pipi-loop-polar-comparison.md) (task 4.4).

## Inputs

- Up-only sources (decision D-A4) and $M_0^{\text{up-only}}=-B_0m_u\epsilon/f_\pi$: [Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md), Step 6.
- Up-only pole data ($A=\kappa_\sigma=B_0m_u$, $m_a^2-m_\sigma^2\to2B_0m_u$, invariance of the span of $P_u$ and $\mathbf 1$ under $L$): [Cartesian large-$\epsilon$ input record](cartesian-sigma-pipi-loop-large-epsilon-inputs.md).
- $X_\Gamma$ and $\delta Z$ of every kept graph in every routing: [Cartesian 1PI trace record](cartesian-sigma-pipi-loop-traces-1pi.md), [Cartesian external-leg record](cartesian-sigma-pipi-loop-traces-external-legs.md).
- Large-$\epsilon$ routing choices, regions of $r=\kappa_\sigma/A$ and generic leading terms: [Cartesian large-$\epsilon$ graph record](cartesian-sigma-pipi-loop-large-epsilon-graphs.md).

## Assumptions

1. **Up-only sources (model assumption).** $\widehat M=m_uP_u\varphi^*/f_a$: $\widehat S^C_{0,0}=\widehat S^C_{1,0}=m_uP_u$, $\widehat P^C_{0,1}=-m_uP_u$, every other $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$ zero. Then $A=\kappa_\sigma=B_0m_u$ and $\mathcal C=-B_0m_uP_u$.
2. **Matrix relations (no matrix components).** $P_u^2=P_u$, $P_d^2=P_d$, $P_uP_d=P_dP_u=0$, $\langle P_u\rangle=\langle P_d\rangle=1$, $\langle\mathbf 1\rangle=3$. For $\pi^+=u\bar d$: $\lambda_{\pi^+}\lambda_{\pi^-}=2P_u$, $\lambda_{\pi^-}\lambda_{\pi^+}=2P_d$, $P_u\lambda_{\pi^+}=\lambda_{\pi^+}P_d=\lambda_{\pi^+}$, $\lambda_{\pi^+}P_u=P_d\lambda_{\pi^+}=0$, the same with $u\leftrightarrow d$ for $\lambda_{\pi^-}$, and $\lambda_{\pi^\pm}^2=0$. A trace that these relations do not reduce stops the script; none occurs.
3. All assumptions of the [large-$\epsilon$ graph record](cartesian-sigma-pipi-loop-large-epsilon-graphs.md) (fixed sources, $\lambda$, $f_\pi$, $B_0$, $m_1$, $m_P$, $m_{\pi^\pm}$; $m_P$ below $m_a$ and $m_\sigma$).

## Method

1. **$\mathbb M^{-1}_{m_a^2}(\mathcal C)$ exactly.** $L(X)=B_0\{\widehat S^C_{0,0},X\}+\frac{m_1^2}{3}\langle X\rangle\mathbf 1$ maps the span of $P_u$ and $\mathbf 1$ into itself: $L(P_u)=2B_0m_uP_u+\frac{m_1^2}{3}\mathbf 1$, $L(\mathbf 1)=2B_0m_uP_u+m_1^2\mathbf 1$. Solving $L(X)-sX=\mathcal C$ with $X$ in this span gives
$$
\mathbb M^{-1}_s(\mathcal C)=\frac{3B_0m_u(s-m_1^2)\,P_u+B_0m_1^2m_u\,\mathbf 1}{3s^2-3m_1^2s-6B_0m_us+4B_0m_1^2m_u},
$$
and $N(s)=1+2\epsilon^2\langle(\mathbb M^{-1}_s(\mathcal C))^2\rangle$ follows from $P_u^2=P_u$, $\langle P_u\rangle=1$, $\langle\mathbf 1\rangle=3$: with the coefficient of $P_u$ written first and the coefficient of $\mathbf 1$ second, $\langle(\mathbb M^{-1}_s(\mathcal C))^2\rangle$ is the square of the first, plus twice their product, plus three times the square of the second. All exact.
2. **Pole data.** $m_a^2$ is the root of the exact $D_a(s)=s-\epsilon^2A+2\epsilon^2\langle\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\rangle$, $\langle\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\rangle=-B_0m_u$ times the sum of the coefficients of $P_u$ and $\mathbf 1$ in $\mathbb M^{-1}_s(\mathcal C)$, obtained as a series in $u^2$ through $u^6$ ($K=3$; the residual of $u^2D_a$ starts at $u^{10}$). $m_\sigma^2=2\lambda f_\pi^2/\epsilon^2+\epsilon^2B_0m_u$ is exact.
3. **Graphs.** Each kept graph is taken with the large-$\epsilon$ routing choice of the generic sources (D12 at large $\epsilon$), the up-only substitution is made and the traces are reduced with assumption 2. $\max(m_P,\dots)$ is replaced by the larger of $m_a$, $m_\sigma$, and $\max(m_a,m_\sigma)=m_a$, because $m_a^2-m_\sigma^2\to2B_0m_u>0$ (below).
4. **Leading terms.** The exact up-only expressions are written as series in $u$ only to read the leading power. Each building block ($m_a^2$, the coefficients of $P_u$ and $\mathbf 1$ in $\mathbb M^{-1}_{m_a^2}(\mathcal C)$, $N$) carries the first power at which it stops being exact, and products track it. The result is recomputed with the root through $u^{10}$ ($K=5$) and must agree.
5. **Boundary.** The up-only sources sit on the boundary $r=\kappa_\sigma/A=1$ of the two regions of the routing choice. The computation is repeated with the routing choices of the neighbouring regions $r<1$ and $r>1$; for $r>1$ also with $\max(m_a,m_\sigma)\to m_\sigma$, the branch of that region.

## Result

Shares of $M_1^{\rm est}/M_0$ that do not vanish (1PI graph: $X_\Gamma/M_0$; $s$ leg: $\delta Z/2$; $\pi^\pm$ legs: $\delta Z$, both legs), with the routing choices at $r=1$:

| graph | share of $M_1^{\rm est}/M_0$: leading power | leading coefficient |
|---|---|---|
| $C2$ | $\epsilon^{0}$ | $- \frac{B_0 m_u}{16 \pi^{2} f_\pi^{2}}$ |
| $C3$ | $\epsilon^{0}$ | $\frac{B_0 m_u}{16 \pi^{2} f_\pi^{2}}$ |
| $C12$ | $\epsilon^{-2}$ | $\frac{2 B_0 m_u + f_\pi^{2} \lambda}{12 \pi^{2} f_\pi^{2}}$ |
| $C13$ | $\epsilon^{-2}$ | $- \frac{7 B_0 m_u + \left(m_{\pi^\pm}\right)^{2}}{12 \pi^{2} f_\pi^{2}}$ |
| $C15$ | $\epsilon^{-2}$ | $\frac{5 B_0 m_u - f_\pi^{2} \lambda + \left(m_{\pi^\pm}\right)^{2}}{12 \pi^{2} f_\pi^{2}}$ |
| $C30$ | $\epsilon^{0}$ | $\frac{B_0 m_u}{16 \pi^{2} f_\pi^{2}}$ |
| $C40$ | $\epsilon^{-2}$ | $\frac{3 \lambda}{8 \pi^{2}}$ |
| $C\sigma4$ | $\epsilon^{-4}$ | $- \frac{9 f_\pi^{2} \lambda^{2}}{16 \pi^{2} B_0 m_u}$ |
| $C\sigma6$ | $\epsilon^{-4}$ | $- \frac{\left(2 B_0 m_u + f_\pi^{2} \lambda\right)^{2}}{16 \pi^{2} B_0 f_\pi^{2} m_u}$ |
| $C\sigma7$ | $\epsilon^{-2}$ | $- \frac{B_0 m_u}{4 \pi^{2} f_\pi^{2}}$ |
| $C\sigma9$ | $\epsilon^{2}$ | $- \frac{B_0^{2} m_u^{2}}{8 \pi^{2} f_\pi^{2} m_P^{2}}$ |
| $C\pi2$ | $\epsilon^{0}$ | $\frac{B_0 m_u}{48 \pi^{2} f_\pi^{2}}$ |
| $C\pi3$ | $\epsilon^{0}$ | $- \frac{B_0 m_u}{48 \pi^{2} f_\pi^{2}}$ |
| $C\pi8$ | $\epsilon^{0}$ | $- \frac{B_0 m_u}{16 \pi^{2} f_\pi^{2}}$ |

Exactly zero: $C1$, $C10$, $C11$, $C14$, $C19$, $C20$, $C21$, $C39$, $C41$, $C42$, $C43$, $C44$, $C45$, $C46$, $C47$, $C48$, $C7$, $C8$, $C9$, $C\pi7$, $C\pi9$, $C\sigma5$, $C\sigma8$.

$$
\frac{M_1^{\rm est}}{M_0}\to - \frac{B_0^{2} m_u^{2}}{8 \pi^{2} f_\pi^{2} m_P^{2}}\,\epsilon^{2},\qquad R\to \frac{B_0^{2} m_u^{2}}{8 \pi^{2} f_\pi^{2} m_P^{2}}\,\epsilon^{2}
$$

Graphs at the leading power: $C\sigma9$.

$m_a^2-m_\sigma^2\to2 B_0 m_u\,\epsilon^{0}$, so $\max(m_a,m_\sigma)=m_a$.

| routing choices and branch | leading power | leading coefficient | graphs | generic leading power of the region | its coefficient for up-only sources | graphs where the up-only and generic leading terms coincide |
|---|---|---|---|---|---|---|
| $r < 1$ | $\epsilon^{2}$ | $- \frac{B_0^{2} m_u^{2}}{8 \pi^{2} f_\pi^{2} m_P^{2}}$ | $C\sigma9$ | $\epsilon^{4}$ | $0$ | $C30$, $C\sigma9$, $C\pi2$, $C\pi3$, $C\pi8$ |
| $r = 1$ | $\epsilon^{2}$ | $- \frac{B_0^{2} m_u^{2}}{8 \pi^{2} f_\pi^{2} m_P^{2}}$ | $C\sigma9$ | $\epsilon^{4}$ | $0$ | $C30$, $C\sigma9$, $C\pi2$, $C\pi3$, $C\pi8$ |
| $r > 1$ | $\epsilon^{2}$ | $- \frac{B_0^{2} m_u^{2}}{8 \pi^{2} f_\pi^{2} m_P^{2}}$ | $C\sigma9$ | $\epsilon^{4}$ | $0$ | $C30$, $C\sigma9$, $C\pi2$, $C\pi3$, $C\pi8$ |

**Leading term.** For the up-only sources
$$
\frac{M_1^{\rm est}}{M_0}\to-\frac{B_0^2m_u^2}{8\pi^2f_\pi^2m_P^2}\,\epsilon^2,
$$
from the graph $C\sigma9$ alone: its all-generator-sum term (two meson lines taking the generator sum, two $s$ ends, scale $m_P$), which D4 does not exclude. Its generic leading term, $\delta Z^{\max}\to-\frac{\epsilon^2}{16\pi^2f_\pi^2m_P^2}\big(\kappa_\sigma^2+3B_0^2\langle\widehat S^C_{1,0}\widehat S^C_{1,0}\rangle+4B_0\kappa_\sigma\langle\widehat S^C_{2,0}\rangle+4B_0^2\langle\widehat S^C_{2,0}\rangle^2\big)$ ([large-$\epsilon$ graph record](cartesian-sigma-pipi-loop-large-epsilon-graphs.md)), gives the same value for the up-only sources.

**Why the generic $\epsilon^4$ term is absent.** Every product of traces in the generic leading numerator contains $\widehat S^C_{3,0}$, $\widehat S^C_{1,2}$, $\widehat S^C_{2,1}$, $\widehat S^C_{0,2}$, $\widehat S^C_{1,1}$ or $\widehat S^C_{2,0}$. All of them vanish for the up-only sources, because $\widehat M=m_uP_u\varphi^*/f_a$ is linear in $s$ and $p$. The script checks that the generic $\epsilon^4$ coefficient is zero for the up-only sources in each case.

**Lower powers.** At $\epsilon^0$ the leading terms cancel in pairs ($C2$ with $C3$, $C30$ with $C\pi8$, $C\pi2$ with $C\pi3$); the remaining graphs are of order $\epsilon^{-2}$ or $\epsilon^{-4}$ or vanish exactly. These lower terms do not affect the leading power.

**Boundary and neighbouring regions.** The up-only sources lie on the boundary $r=\kappa_\sigma/A=1$, where $m_a^2-m_\sigma^2$ is of order $\epsilon^0$ instead of $\epsilon^2$; its first nonvanishing term $2B_0m_u>0$ puts the ALP-like pole above the radial pole. The leading term is the same with the routing choices of $r<1$, of $r=1$ and of $r>1$ (the last also with $\max(m_a,m_\sigma)=m_\sigma$), so the limit agrees from both neighbouring regions.

**Status of the leading graph.** The result rests on a single graph, so no cancellation between graphs with the same leading power can occur within the estimate. It depends on the symbolic scale $m_P$, the largest pole mass of the nine meson-like states, which rule D3 assigns to the bare generator sum, with no lower limit on the loop scale. At large $\epsilon$, $m_P<m_\sigma$, so $(m_P^2)^{-1}$ from the residue derivative is larger than it would be with a scale $m_\sigma$. This is a consequence of the D3 assumption and is listed in the [Cartesian limits record](cartesian-sigma-pipi-loop-limits.md).

## Verification

Exact symbolic checks in the script and the tests, no numerical matrix: $\langle\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\rangle$ and $\langle(\mathbb M_s^{-1}(\mathcal C))^2\rangle$ from the exact $\mathbb M^{-1}_s(\mathcal C)$ agree with the large-$s$ series $-\sum_kc_k/s^{k+1}$ and $\sum_n(n+1)c_n/s^{n+2}$ built from the generic trace sums ($k,n\le2$); $D_a$ equals the up-only $D_a$ of the input script; the root residual and its agreement with the generic root through $u^4$; $m_a^2-m_\sigma^2\to2B_0m_u$; the trace reductions (including an unreduced trace that must stop the script); the tree $-B_0m_u\epsilon/f_\pi$; the leading term, its graph, the exactly vanishing graphs, the agreement between $K=3$ and $K=5$, and the sum of the graphs in the three cases; the series arithmetic against the lowest powers of the numerator and denominator of the fully substituted exact expression for $C2$ and $C\sigma9$; the generic $\epsilon^4$ coefficient vanishes for up-only sources, and for every graph whose up-only and generic leading powers coincide ($C30$, $C\sigma9$, $C\pi2$, $C\pi3$, $C\pi8$) the generic coefficient gives the up-only one. The result block is the `--markdown` output of the script (tested).

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_up_only_estimate.py
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_up_only_estimate.py --markdown
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_up_only_estimate
```

[Script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_up_only_estimate.py) (about five minutes), [tests](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_up_only_estimate.py), built on [`derive_sigma_pipi_cartesian_loop_large_epsilon_graphs.py`](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_large_epsilon_graphs.py); the generic word-rewriting helper and the projector rules come from [`derive_sigma_pipi_loop_tree.py`](../../../../../../scripts/derive_sigma_pipi_loop_tree.py) (no polar vertex is used).
