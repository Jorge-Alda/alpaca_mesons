---
title: Behaviour at large epsilon of the ALP-like pole data in the one-loop sigma to pi+ pi- estimate
status: mixed
created: 2026-11-26
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Behaviour at large $\epsilon$ of the ALP-like pole data

Task 2.3 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama` (decision D10). The regime of interest is $\epsilon=f_\pi/f_a\gg1$. No step below assumes $\epsilon\ll1$. No loop integral, loop function or loop amplitude is evaluated.

## Conventions, units and inputs

- Natural units: $s$, $m_a^2$, $m_1^2$ in GeV$^2$; $B_0$, $f_\pi$ and the sources in GeV. Notation of Alda:2026nrd and DiLuzio:2026poh.
- Lagrangian, polar fields, connection-free frame, the normalization of the nine $U(3)$ generators (trace of the product of two generators equal to $2$ for equal generators and $0$ otherwise) and the assumptions: [tree and legs record](../one-loop-estimate/murayama-sigma-pipi-loop-tree-and-legs.md). Vertex rules: [vertex-rule record](../one-loop-estimate/murayama-sigma-pipi-loop-vertex-rules.md).
- $\widehat\Theta_1=0$, so $\mathcal C=B_0\widehat P_{0,1}$ ([implicit crossing record](../../loop-estimates/power-counting/murayama-radial-crossing-implicit.md)).
- Definitions, as in [the estimate record for $a\to\pi^+\pi^-\pi^0$](../../../couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md):
$$
\mathbb M_s(X)=B_0\{\widehat S_{0,0},X\}+\frac{m_1^2}{3}\langle X\rangle\mathbf 1-sX,\qquad N(s)=1+2\epsilon^2\big\langle(\mathbb M_s^{-1}(\mathcal C))^2\big\rangle,
$$
$$
v_a=\Big(\frac{2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)}{\sqrt{N(m_a^2)}},\ \frac{1}{\sqrt{N(m_a^2)}}\Big),\qquad \xi_a=\frac{1}{\sqrt{N(m_a^2)}} .
$$
- $\mathbb M_{m_a^2}$ is invertible (assumption 5 of the tree and legs record). Matrices are noncommuting symbols inside traces; no components are used.
- The script writes the trace $\langle(\mathbb M^{-1}_{m_a^2}(\mathcal C))^2\rangle$ as the symbol `q`; this record writes the trace in full.

## Step 1. Quadratic Lagrangian from the vertex rules

The two-leg case of the potential vertex rule, plus the two-point kinetic terms of $\frac{f_\pi^2}{4}\langle u_\mu u^\mu\rangle$ and $\frac12(\partial a)^2$ and the anomalous term $-\frac{m_1^2}{12}\langle\Phi\rangle^2$ of the Lagrangian, gives the second derivatives at momentum $p$, $p^2=s$:
$$
D_{PP}(X,Y)=\frac s2\langle XY\rangle-\frac{B_0}{2}\langle\widehat S_{0,0}\{X,Y\}\rangle-\frac{m_1^2}{6}\langle X\rangle\langle Y\rangle=-\frac12\langle X\,\mathbb M_s(Y)\rangle ,
$$
$$
D_{aP}(X)=\epsilon\langle\mathcal C X\rangle,\qquad D_{aa}=s+2B_0\epsilon^2\langle\widehat S_{0,2}\rangle .
$$
So the ALP mass entry is $-2B_0\epsilon^2\langle\widehat S_{0,2}\rangle$. There is no kinetic mixing in the connection-free frame. $\mathbb M_s$ is self-adjoint for the trace, $\langle X\mathbb M_s(Y)\rangle=\langle\mathbb M_s(X)Y\rangle$. Neither $\mathbb M_s$ nor $\mathcal C$ contains $\epsilon$.

## Step 2. Elimination of the mesons: pole condition and normalization

Stationarity of the quadratic Lagrangian in the meson field along any matrix $X$, $D_{PP}(X,\Phi)+a\,D_{aP}(X)=0$, is solved by $\Phi=2\epsilon a\,\mathbb M_s^{-1}(\mathcal C)$. The script checks this as a trace identity for a generic probe $X$, using only the defining relation $\mathbb M_s(\mathbb M_s^{-1}(\mathcal C))=\mathcal C$. Inserting the solution gives the ALP inverse propagator after the mesons are integrated out:
$$
D_a(s)=s+2B_0\epsilon^2\langle\widehat S_{0,2}\rangle+2\epsilon^2\langle\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\rangle .
$$
Since $d\mathbb M_s/ds=-1$, $d\mathbb M_s^{-1}/ds=\mathbb M_s^{-2}$, and by self-adjointness $D_a'(s)=N(s)$. $N(s)$ is also the kinetic norm $\frac12\langle ZZ\rangle+1$ of the pole vector $(Z,1)$, $Z=2\epsilon\,\mathbb M_s^{-1}(\mathcal C)$. Dividing by $\sqrt N$ gives $v_a$. This agrees with the general-source form of [the operator-inverse record](../../../../meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-finite-mixing-operator-inverse.md) after the replacements of the Murayama model.

## Step 3. Exact identities of the meson-line residue (any $\epsilon$)

The internal meson line carries $\mathbf 1_{10}-v_av_a^{\rm T}$ on the space (nine meson matrices) $\oplus$ (ALP field). Its entries are:

- ALP-field entry: $1-1/N(m_a^2)$;
- mixed entries: $-2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)/N(m_a^2)$;
- meson block: $\mathbf 1_9-\frac{4\epsilon^2}{N}\,\mathbb M^{-1}_{m_a^2}(\mathcal C)\otimes\mathbb M^{-1}_{m_a^2}(\mathcal C)$, i.e. $\mathbf 1_9$ minus $(1-1/N)$ times the projector on the meson direction of the ALP-like state. Acting on a meson matrix $X$, this projector gives $\mathbb M^{-1}_{m_a^2}(\mathcal C)\,\langle\mathbb M^{-1}_{m_a^2}(\mathcal C)\,X\rangle/\langle(\mathbb M^{-1}_{m_a^2}(\mathcal C))^2\rangle$.

The script checks exactly, for any $\epsilon$, that $v_a^{\rm T}v_a=1$, so the residue is a projector of rank nine, that $4\epsilon^2/N=2(1-1/N)/\langle(\mathbb M^{-1}_{m_a^2}(\mathcal C))^2\rangle$, that the projector on the meson direction is idempotent, and that its trace over the nine mesons is one (completeness). The graph scripts use these exact forms.

## Step 4. Large $\epsilon$ at fixed sources (assumption)

**Assumption.** The sources are fixed by the model and do not depend on $\epsilon$; the ALP mass entry is $\epsilon^2A$, with $A=-2B_0\langle\widehat S_{0,2}\rangle$ fixed and nonzero. Then
$$
D_a(s)=s-\epsilon^2F(s),\qquad F(s)=A-2\langle\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\rangle .
$$

**Large-$s$ expansion.** Write $\mathbb M_s=L-s$ with $L(X)=B_0\{\widehat S_{0,0},X\}+\frac{m_1^2}{3}\langle X\rangle\mathbf 1$. The exact operator identity $(L-s)\big(-\frac1s\sum_{k=0}^{K}(L/s)^k\big)=1-(L/s)^{K+1}$ (checked for $K=0,\dots,4$ with a noncommuting symbol) gives the expansion $\mathbb M_s^{-1}=-\frac1s\sum_k(L/s)^k$. It holds for $s$ above the spectrum of $L$, i.e. above the unmixed meson masses squared. With $c_k=\langle\mathcal C\,L^k(\mathcal C)\rangle$ and $L$ self-adjoint:
$$
F(s)=A+\frac{2c_0}{s}+\frac{2c_1}{s^2}+O(s^{-3}),\qquad \big\langle(\mathbb M_s^{-1}(\mathcal C))^2\big\rangle=\sum_n\frac{(n+1)c_n}{s^{n+2}} .
$$

**Root.** The root of $D_a$, as a series in $1/\epsilon^2$, is
$$
m_a^2=\epsilon^2A+\frac{2c_0}{A}+\frac{1}{\epsilon^2}\Big(\frac{2c_1}{A^2}-\frac{4c_0^2}{A^3}\Big)+O(\epsilon^{-4}) .
$$
Its existence for large $\epsilon$ follows from the intermediate value theorem. $D_a(\epsilon^2A/2)=\epsilon^2(A/2-F)<0$ and $D_a(2\epsilon^2A)=\epsilon^2(2A-F)>0$, because $F\to A$ at large $s$. At the root:
$$
\mathbb M^{-1}_{m_a^2}(\mathcal C)=-\frac{\mathcal C}{\epsilon^2A}+O(\epsilon^{-4}),\qquad N(m_a^2)=1+\frac{2c_0}{\epsilon^2A^2}+O(\epsilon^{-4})\to1,\qquad \xi_a\to1 .
$$
The meson part of $v_a$ is of order $\epsilon^{-1}$. $N$ was computed both as $1+2\epsilon^2\langle(\mathbb M^{-1})^2\rangle$ and as $D_a'$, and the two agree. The ALP-like state is mostly the ALP field, and its mass grows as $\epsilon\sqrt A$.

**Assumptions:** $A>0$; $F$ has no zero at large $s$ (it tends to $A$); the root found is the continuation of the ALP-like pole; $\epsilon^2A$ lies above the spectrum of $L$, so that the expansion holds. The use of leading-order chiral perturbation theory at $m_a^2\sim\epsilon^2A$ is a limit of the estimate ([limits record](../one-loop-estimate/murayama-sigma-pipi-loop-limits.md)).

**Up-only check.** The sources $\widehat S_{r,s}+i\widehat P_{r,s}=m_u(-i)^sP_u/s!$ give $\widehat S_{0,2}=-m_uP_u/2$ and $\widehat P_{0,1}=-m_uP_u$. With only the projector relations $P_u^2=P_u$, $\widehat S_{0,0}P_u=P_u\widehat S_{0,0}=m_uP_u$, $\langle P_u\rangle=1$, this gives $A=B_0m_u>0$, $c_0=B_0^2m_u^2$, $c_1=B_0^2m_u^2(2B_0m_u+m_1^2/3)$, and
$$
m_a^2=\epsilon^2B_0m_u+2B_0m_u+\frac{2m_1^2}{3\epsilon^2}+\dots,\qquad N=1+\frac{2}{\epsilon^2}+\dots,
$$
with meson matrix $2P_u/\epsilon+\dots$ in $v_a$.

The large-$\epsilon$ behaviour of the radial mass (with $\lambda$ fixed), the residue pieces and the leg factors is in the [large-$\epsilon$ input record](murayama-sigma-pipi-loop-large-epsilon-inputs.md); the vertices and the tree in the [large-$\epsilon$ vertex record](murayama-sigma-pipi-loop-large-epsilon-vertices-tree.md); the graphs, $M_1^{\rm est}$ and $R\sim\epsilon^4$ in the [large-$\epsilon$ graph record](murayama-sigma-pipi-loop-large-epsilon-graphs.md).

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_large_epsilon.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_large_epsilon
```

[Script](../../../../../../scripts/derive_sigma_pipi_loop_large_epsilon.py) (sections 1 to 4), [tests](../../../../../../scripts/test_derive_sigma_pipi_loop_large_epsilon.py). All checks are exact symbolic identities: quadratic Lagrangian, stationarity, $D_a$, $N=D_a'$ and the norm, the residue identities, the resolvent series, the root and normalization at fixed sources, and the up-only values.
