---
title: Cartesian frame - assembled size estimate M1_est of the one-loop sigma to pi+ pi- amplitude and the ratio R
status: derived
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Cartesian frame: $M_1^{\rm est}$ and $R$ for $\sigma\to\pi^+\pi^-$

Task 4.1 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama-cartesian`. **No loop integral, scalar loop function, loop coefficient or loop amplitude is computed.** Each loop enters only through the size factor $S_\Gamma(m_{\rm term}^2)^{2-r_\Gamma}/(16\pi^2)$ per term (D1, D2) and through the residue prescription for the external legs.

**Caveat.** $M_1^{\rm est}$ is an order-of-magnitude estimate per graph class. Graphs with the same $\epsilon$ scaling can cancel in the exact loop, and these cancellations can differ between the Cartesian and the polar parametrization. Agreement with the polar estimate, including its up-only prefactor $-21/(64\pi^2)$, is not expected, and neither estimate is authoritative. No comparison with the polar estimate is made here (task 4.4).

## Conventions

Units GeV. Cartesian fields $\varphi=f_a+s+ip$, $\epsilon=f_\pi/f_a$, all momenta incoming, $M=V$ for a contact graph, as in the [Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md). $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$ are generic noncommuting Hermitian matrices of mass dimension one; $B_0$, $f_\pi$ and all masses have mass dimension one, $m_1^2$ two; $\epsilon$, $N(m_a^2)$, $\mathbb M^{-1}_{m_a^2}(\mathcal C)$ and $\lambda_{\pi^\pm}$ are dimensionless.

## Inputs

- Tree $M_0=-\frac{B_0\epsilon}{2f_\pi}\big\langle\widehat S^C_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\big\rangle$: [Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md), Step 6.
- $X_\Gamma^{\max}$ of the 26 kept 1PI graphs (termwise routing maximum, D12, with the per-term scales of D3): [Cartesian 1PI trace record](cartesian-sigma-pipi-loop-traces-1pi.md).
- $\delta Z^{\max}$ of the 6 kept $s$ self-energies and of the 5 kept self-energies of each charged pion, and the residue prescription: [Cartesian external-leg record](cartesian-sigma-pipi-loop-traces-external-legs.md).
- Rules R1–R10: [Cartesian selection record](cartesian-sigma-pipi-loop-graph-selection.md); graph list: [Cartesian graph-list record](cartesian-sigma-pipi-loop-graph-list.md).

## Method

**Step 1 (assembly).**
$$
M_1^{\rm est}=\sum_\Gamma X_\Gamma^{\max}+\frac{M_0}{2}\big(\delta Z_s^{\max}+\delta Z_{\pi^+}^{\max}+\delta Z_{\pi^-}^{\max}\big),
$$
where the sum runs over the kept 1PI graphs and each $\delta Z^{\max}$ is the sum of the termwise routing maxima of the kept self-energies of that leg. $\delta Z_{\pi^+}^{\max}$ is taken from the trace script; $\delta Z_{\pi^-}^{\max}$ is computed again with the roles of the two legs exchanged ($\pi^-$ in, $\pi^+$ out). The $\pi^\pm$ self-energies do not depend on the routing.

**Step 2 (ratio).** $R=|M_1^{\rm est}|/|M_0|$ with the exact $M_0$. When the tree vanishes, $R$ is undefined and only $M_1^{\rm est}$ is reported. The script does not divide in that case; it is checked with $\widehat S^C_{1,0}=0$, for which $M_0=0$ while $M_1^{\rm est}$ keeps 57 products of traces.

**Step 3 (checks).** Exact symbolic identities, listed under Verification.

## Result

$M_1^{\rm est}$ is a sum of products of traces with exact coefficients in $B_0$, $f_\pi$, $\epsilon$, $N(m_a^2)$, $m_\sigma$, $m_a$, $m_{\pi^\pm}$, $m_P$ and $\max(\dots)$ of these masses. 34 products have coefficients that depend on the masses through Piecewise branches (routing choice, D12; branches in homogeneous inequalities between squared masses). The trace content is (script output):

| quantity | value |
|---|---|
| kept 1PI graphs | 26 |
| kept s self-energies | 6 |
| kept pi self-energies (each charged leg) | 5 |
| products in the 1PI sum | 116 |
| products in dZ_s^max | 18 |
| products in dZ_pi^max | 6 |
| products in the external-leg term | 48 |
| products in M_1^est | 164 |
| Piecewise products in M_1^est | 34 |
| products with S^C_00 in the 1PI sum | 0 |
| products with S^C_00 in the external-leg term | 0 |
| products with S^C_00 in M_1^est | 0 |
| source letters in $M_1^{\rm est}$ | $\widehat P^C_{0,1}$, $\widehat P^C_{1,1}$, $\widehat P^C_{2,0}$, $\widehat S^C_{0,1}$, $\widehat S^C_{0,2}$, $\widehat S^C_{1,0}$, $\widehat S^C_{1,1}$, $\widehat S^C_{1,2}$, $\widehat S^C_{2,0}$, $\widehat S^C_{2,1}$, $\widehat S^C_{3,0}$ |
| graphs with $\widehat S^C_{0,0}$ | $C12$, $C13$, $C15$ |
| products with $\widehat S^C_{0,0}$ at $m_{\rm term}=m_a$ | $C12$: 72, $C13$: 72, $C15$: 72 |
| products with $\widehat S^C_{0,0}$ at $m_{\rm term}=\max\left(m_P, m_a\right)$ | $C13$: 66, $C15$: 66 |
| $\widehat S^C_{0,0}$ part cancels at each $m_{\rm term}$ | True |
| $\widehat S^C_{0,0}$ absent from the routing-dependent branches | True |

Besides the source letters, the products of traces contain $\lambda_{\pi^\pm}$ and $\mathbb M^{-1}_{m_a^2}(\mathcal C)$, which depends on $\widehat P^C_{0,1}$ and $\widehat S^C_{0,0}$ ([frame record](cartesian-sigma-pipi-loop-frame-and-tree.md), Step 3). $\widehat S^C_{0,1}$ and $\widehat S^C_{1,1}$ appear only inside longer traces; their single traces vanish by assumption.

**$\widehat S^C_{0,0}$ as an explicit letter.** It enters the graphs $C12$, $C13$, $C15$ (the bubbles with $\{\pi^+,\pi^-\}$ at one vertex and $s$ at the other, with two ALP-like lines, one ALP-like and one meson line, and two meson lines), through the four-meson source vertex $\widehat S^C_{0,0}$ at $\{\pi^+,\pi^-\}$. The script checks that its products cancel in the sum of the three graphs, separately at $m_{\rm term}=m_a$ and at $m_{\rm term}=\max(m_P,m_a)$, and that they never enter a routing-dependent branch. So $\widehat S^C_{0,0}$ appears in $M_1^{\rm est}$ only through $\mathbb M^{-1}_{m_a^2}(\mathcal C)$, $N(m_a^2)$ and the masses. This cancellation agrees with rule R3 (an ALP-like line plus pieces (ii)+(iii)+(iv) of a meson line is the bare $p$ line at the common scale $m_a$, and the bare $p$ line has no meson matrix for $\widehat S^C_{0,0}$ to attach to). It is a property of the estimate with the D3 scales, not a statement about the exact loop.

**Ratio.** For generic sources $M_0\neq0$ and $R=|M_1^{\rm est}|/|M_0|$ is the quotient of the two exact expressions. Its leading large-$\epsilon$ behaviour is in the [Cartesian large-$\epsilon$ graph record](cartesian-sigma-pipi-loop-large-epsilon-graphs.md), and the up-only limit in the [Cartesian up-only record](cartesian-sigma-pipi-loop-up-only-estimate.md). The limits of the estimate are in the [Cartesian limits record](cartesian-sigma-pipi-loop-limits.md).

## Verification

Exact symbolic checks in the script and the tests, no random or numerical matrix:
- the sum of the parts: $M_1^{\rm est}-\sum_\Gamma X_\Gamma^{\max}-M_0(\delta Z_s^{\max}+\delta Z_{\pi^+}^{\max}+\delta Z_{\pi^-}^{\max})/2=0$;
- $\delta Z_{\pi^+}^{\max}=\delta Z_{\pi^-}^{\max}$;
- the tree equals $-\frac{B_0\epsilon}{2f_\pi}\langle\widehat S^C_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\rangle$;
- no floating-point number in any part;
- dimensions: under $m\to cm$ for every mass, $B_0\to cB_0$, $f_\pi\to cf_\pi$, $m_1^2\to c^2m_1^2$ and every source letter $\to c\times$letter, each $X_\Gamma^{\max}$, $M_0$ and $M_1^{\rm est}$ scale as $c$ and each $\delta Z^{\max}$ as $c^0$, in every combination of Piecewise branches. This uses $\max(cx_1,\dots,cx_n)=c\max(x_1,\dots,x_n)$ for $c>0$, checked for two and three arguments in every ordering, and the homogeneity of every Piecewise condition (its two sides scale with the same power of $c$). Negative control: $M_1^{\rm est}$ does not scale as $c^0$ or $c^2$;
- $R$ is returned undefined for $\widehat S^C_{1,0}=0$, where $M_1^{\rm est}\neq0$;
- the letter content and the $\widehat S^C_{0,0}$ cancellation of the table above; the table is the `--markdown` output of the script.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_estimate.py
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_estimate.py --markdown
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_estimate
```

[Script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_estimate.py), [tests](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_estimate.py); inputs from `estimate_inputs()` and `entries()` of [`derive_sigma_pipi_cartesian_loop_traces.py`](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_traces.py).
