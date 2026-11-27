---
title: Cartesian frame - assumptions, quadratic Lagrangian, pole data and tree amplitude of sigma to pi+ pi-
status: mixed
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Cartesian frame: assumptions, quadratic Lagrangian, pole data and tree of $\sigma\to\pi^+\pi^-$

Task 1.1 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama-cartesian`. The Cartesian estimate is an independent order-of-magnitude estimate of the one-loop $\sigma\to\pi^+\pi^-$ amplitude; neither it nor the polar estimate is authoritative. No loop integral, scalar loop function, loop coefficient or loop amplitude is evaluated here. Everything below is derived from the Cartesian Lagrangian; no polar vertex, graph or field map is used.

## Conventions

Natural units, GeV. Conventions and notation of Alda:2026nrd and DiLuzio:2026poh: $\epsilon=f_\pi/f_a$, $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$, $\lambda_0=\sqrt{2/3}\,\mathbf 1$, $U=e^{i\Phi/f_\pi}$, $\chi=2B_0\widehat M$. Cartesian fields $\varphi=f_a+s+ip$; $s$ is the Cartesian radial field, $p$ the Cartesian ALP field, and the external $\sigma$ is the pole of $s$. The Lagrangian is
$$
\mathcal L=\frac{f_\pi^2}{4}\langle u_\mu u^\mu\rangle+B_0f_\pi^2\sum_{r,s}\Big(\frac{s}{f_a}\Big)^r\Big(\frac{p}{f_a}\Big)^s\Big[\big\langle\widehat S^C_{r,s}\cos\tfrac{\Phi}{f_\pi}\big\rangle+\big\langle\widehat P^C_{r,s}\sin\tfrac{\Phi}{f_\pi}\big\rangle\Big]-\frac{m_1^2}{12}\langle\Phi\rangle^2+\frac12(\partial s)^2+\frac12(\partial p)^2-V_\varphi ,
$$
$$
V_\varphi=\frac{m_\varphi^2}{2}|\varphi|^2+\frac{\lambda}{4}|\varphi|^4,\qquad |\varphi|^2=(f_a+s)^2+p^2 .
$$
The source part is $\frac{f_\pi^2}{4}\langle\chi_+\rangle$ with $\widehat M(s,p)=\sum_{r,s}(\widehat S^C_{r,s}+i\widehat P^C_{r,s})(s/f_a)^r(p/f_a)^s$. The second index of $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$ counts the powers of $p/f_a$ (notation of the OpenSpec change); it is not the field $s$. The coefficients are generic noncommuting Hermitian $3\times3$ matrices of mass dimension one.

A vertex $V(x_1,\dots,x_n)$ is the $n$-th derivative of $\mathcal L$ with respect to the fields of its legs, with $\partial_\mu\to-ik_\mu$ and all momenta incoming, so that $M=V$ for a contact graph (same convention as the [polar vertex record](../one-loop-estimate/murayama-sigma-pipi-loop-vertex-rules.md)). A leg is a radial leg $s$ or a pseudoscalar leg $(\lambda_x,\xi_x)$, where $\lambda_x$ is its meson matrix and $\xi_x$ its $p$ component. The external legs are $s$ (momentum $p_\sigma$), $\pi^+=(\lambda_{\pi^+},0)$ and $\pi^-=(\lambda_{\pi^-},0)$; the charged legs have no $p$ component because $p$ is neutral (charge conservation, a fact). In the two-point functions below, $\mathbb M_{k^2}$ is evaluated at the squared momentum $k^2$ of the line.

## Inputs

- Cartesian Lagrangian above: Alda:2026nrd, DiLuzio:2026poh ([model extract](../../../../../extracts/DiLuzio:2026poh/model-lagrangian-and-field-normalization.md)).
- Polar assumptions 1–10, used only for their wording and numbering: [polar tree and legs record](../one-loop-estimate/murayama-sigma-pipi-loop-tree-and-legs.md), [general-source radial record](../../../couplings-and-decays/finite-mixing/radial-decay-polar-cartesian-general-sources.md), [polar selection record](../one-loop-estimate/murayama-sigma-pipi-loop-graph-selection.md).
- Charge rule (assumption B1): [charged decoupling record](../../../../meson-mixing/finite-mixing/spectrum-and-poles/quadratic-mixing/alp-lo-charged-meson-decoupling.md).

## Method

The script expands $\mathcal L$ directly with formal legs: $s=\sum\phi_j$ over the radial legs, $p=\sum_j\xi_j\phi_j$, $\Phi=\sum_j\phi_j\lambda_j$, $\partial_\mu\phi_j\to-ik_{j,\mu}\phi_j$. It extracts the coefficient of $\phi_1\cdots\phi_n$ as exact cyclic words. The source part is expanded in the exponential form $\frac{f_\pi^2}{4}2B_0\langle\widehat MU^\dagger+U\widehat M^\dagger\rangle$. The one-point ($n=1$) and two-point ($n=2$, momenta $k$ and $-k$) functions give the conditions and the quadratic data. The mesons are then eliminated at the stationary point of the quadratic Lagrangian.

## Result

### Step 1: one-point functions

$$
V(s)=B_0\epsilon f_\pi\langle\widehat S^C_{1,0}\rangle-\frac{f_\pi}{\epsilon}\big(m_\varphi^2+\lambda f_a^2\big),\qquad
V(\lambda_x,\xi_x)=B_0f_\pi\langle\widehat P^C_{0,0}\lambda_x\rangle+B_0\epsilon f_\pi\,\xi_x\langle\widehat S^C_{0,1}\rangle .
$$
The meson matrices span all Hermitian matrices, so no tadpole of the mesons or of $p$ requires $\widehat P^C_{0,0}=0$ (matrix identity) and $\langle\widehat S^C_{0,1}\rangle=0$. The $s$ tadpole is removed by the stationarity condition
$$
m_\varphi^2+\lambda f_a^2=\epsilon^2B_0\langle\widehat S^C_{1,0}\rangle .
$$

### Step 2: two-point functions after stationarity

$$
V(s,s)=k^2-m_\sigma^2,\qquad m_\sigma^2=2\lambda f_a^2+\epsilon^2\kappa_\sigma,\qquad \kappa_\sigma=B_0\big(\langle\widehat S^C_{1,0}\rangle-2\langle\widehat S^C_{2,0}\rangle\big),
$$
$$
V(s;\lambda_x,\xi_x)=\epsilon B_0\langle\widehat P^C_{1,0}\lambda_x\rangle+\epsilon^2B_0\,\xi_x\langle\widehat S^C_{1,1}\rangle ,
$$
$$
V\big((\lambda_x,\xi_x),(\lambda_y,\xi_y)\big)=\xi_x\xi_y\big(k^2-\epsilon^2A\big)+\epsilon\,\xi_x\langle\mathcal C\lambda_y\rangle+\epsilon\,\xi_y\langle\mathcal C\lambda_x\rangle-\frac12\big\langle\lambda_x\,\mathbb M_{k^2}(\lambda_y)\big\rangle ,
$$
with
$$
A=B_0\big(\langle\widehat S^C_{1,0}\rangle-2\langle\widehat S^C_{0,2}\rangle\big),\qquad \mathcal C=B_0\widehat P^C_{0,1},\qquad \mathbb M_{k^2}(X)=B_0\{\widehat S^C_{0,0},X\}+\frac{m_1^2}{3}\langle X\rangle\mathbf 1-k^2X .
$$
$m_\sigma^2$ is $\partial^2V/\partial s^2$ at the stationary point, where $V$ is $V_\varphi$ minus the source term at $U=\mathbf 1$, $p=0$. The $p$ mass entry is $\epsilon^2A$; $A>0$ and $\kappa_\sigma>0$ are assumptions. $V_\varphi$ has no $sp$ term, so no $s$–$p$ mixing requires $\langle\widehat S^C_{1,1}\rangle=0$, and no $s$–meson mixing requires $\widehat P^C_{1,0}=0$ (matrix identity). $\mathbb M_{k^2}$ is self-adjoint: $\langle X\mathbb M_{k^2}(Y)\rangle=\langle\mathbb M_{k^2}(X)Y\rangle$.

### Step 3: elimination of the mesons, pole equation and pole vector

Stationarity of the quadratic Lagrangian in $\Phi$ at fixed $p$ gives $\Phi=2\epsilon p\,\mathbb M^{-1}_{k^2}(\mathcal C)$. Inserting it gives the inverse propagator of the ALP-like state,
$$
D_a(k^2)=k^2-\epsilon^2A+2\epsilon^2\big\langle\mathcal C\,\mathbb M^{-1}_{k^2}(\mathcal C)\big\rangle ,
$$
and the ALP-like pole equation is $D_a(m_a^2)=0$. Differentiating the defining relation $\mathbb M_{k^2}(\mathbb M^{-1}_{k^2}(\mathcal C))=\mathcal C$ and using self-adjointness gives
$$
N(k^2)=D_a'(k^2)=1+2\epsilon^2\big\langle\big(\mathbb M^{-1}_{k^2}(\mathcal C)\big)^2\big\rangle ,
$$
which is also the kinetic norm $\frac12\langle Z^2\rangle+1$ of the vector $(Z,1)$ on (mesons $\oplus$ $p$), $Z=2\epsilon\,\mathbb M^{-1}_{k^2}(\mathcal C)$. The pole vector is
$$
v_a=\Big(\frac{2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)}{\sqrt{N(m_a^2)}},\ \frac{1}{\sqrt{N(m_a^2)}}\Big)\quad\text{on (mesons}\oplus p).
$$
The symbols $A$, $\kappa_\sigma$, $\mathcal C$, $N$, $\mathbb M_s$ have the same names as in the polar records; their Cartesian definitions are the ones above. They are not compared with the polar ones here.

### Step 4: CP

CP acts as $s\to s$, $p\to-p$, $\Phi\to-\Phi^{\rm T}$. The script computes the vertex of $\mathcal L(s,-p,-\Phi^{\rm T})$ from the direct expansion, using $\langle W\rangle=\langle W^{\rm T}\rangle$, and finds that it equals the vertex of $\mathcal L$ with every source $X^C_{r,s}$ replaced by $\eta X^{C\,{\rm T}}_{r,s}$, with $\eta=(-1)^s$ for $\widehat S^C$ and $(-1)^{s+1}$ for $\widehat P^C$ (all patterns with $n\le4$). The meson matrices span all Hermitian matrices, so CP invariance is equivalent to
$$
\widehat S^C_{r,s}=(-1)^s\big(\widehat S^C_{r,s}\big)^{\rm T},\qquad \widehat P^C_{r,s}=(-1)^{s+1}\big(\widehat P^C_{r,s}\big)^{\rm T},
$$
hence $\langle\widehat S^C_{r,s}\rangle=0$ for odd $s$, in particular $\langle\widehat S^C_{0,1}\rangle=\langle\widehat S^C_{1,1}\rangle=0$.

### Step 5: charge rule

The generic proof is reused: for letters with $[Q,X_k]=q_kX_k$, $(q_1+\dots+q_n)\langle X_1\cdots X_n\rangle=\langle[Q,X_1\cdots X_n]\rangle=0$. With $[Q,\widehat S^C_{0,0}]=0$, $[Q,\mathbb M_{k^2}(X)]=\mathbb M_{k^2}([Q,X])$, so $\mathbb M^{-1}_{m_a^2}(\mathcal C)$ is neutral when $[Q,\widehat P^C_{0,1}]=0$.

### Step 6: tree amplitude

$$
M_0=V(s;\pi^+,\pi^-)=-\frac{B_0\epsilon}{2f_\pi}\big\langle\widehat S^C_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\big\rangle .
$$
It has no momentum dependence. The $V_\varphi$ term $-2\lambda f_a\xi_x\xi_y$ of $V(s;x,y)$ and the $\widehat P^C_{1,1}$, $\widehat S^C_{1,2}$ terms need a $p$ component, which the charged legs do not have. The external $s$ has no meson or $p$ component (Step 2 with assumption 2) and a canonical kinetic term, so $M_0$ is the on-shell amplitude.

**Up-only limit (D-A4).** With the Yukawa coupling to $\varphi^*$ (DiLuzio:2026poh), $\widehat M=m_uP_u\varphi^*/f_a=m_uP_u\big(1+(s-ip)/f_a\big)$, i.e. $\widehat S^C_{0,0}=\widehat S^C_{1,0}=m_uP_u$, $\widehat P^C_{0,1}=-m_uP_u$, all other coefficients zero, with $P_u$ the up projector. Then $A=\kappa_\sigma=B_0m_u$ and $\mathcal C=-B_0m_uP_u$. With $\lambda_{\pi^+}\lambda_{\pi^-}=2P_u$, $\lambda_{\pi^-}\lambda_{\pi^+}=2P_d$, $P_u^2=P_u$, $P_d^2=P_d$, $P_uP_d=0$, $\langle P_u\rangle=\langle P_d\rangle=1$ (no explicit matrix):
$$
M_0^{\text{up-only}}=-\frac{B_0m_u\epsilon}{f_\pi}=-\frac{B_0m_u}{f_a} .
$$
The script obtains the same value from the substituted tree and from the direct expansion with $\widehat M=m_uP_u\varphi^*/f_a$.

## Assumptions in the Cartesian invariants

| polar no. | Cartesian statement | status |
|---|---|---|
| 1 | No tree tadpoles: $\widehat P^C_{0,0}=0$, $\langle\widehat S^C_{0,1}\rangle=0$; $s$ tadpole removed by $m_\varphi^2+\lambda f_a^2=\epsilon^2B_0\langle\widehat S^C_{1,0}\rangle$ | derived, Step 1 |
| 2 | No quadratic $s$–$p$ mixing: $\langle\widehat S^C_{1,1}\rangle=0$ ($V_\varphi$ has no $sp$ term); no $s$–meson mixing: $\widehat P^C_{1,0}=0$ | derived, Step 2 |
| 3 | Anomalous singlet term $-\frac{m_1^2}{12}\langle\Phi\rangle^2$ independent of $s$ and $p$ | assumption |
| 4 | Leading $U(3)$ chiral order; Cartesian frame with no derivative coupling of $s$ or $p$ to the mesons | assumption; see the remark below |
| 5 | Regular simple poles; $\mathbb M_{m_a^2}$ invertible | assumption |
| 6 | The $p$–meson quadratic entry comes only from the source: $\mathcal C=B_0\widehat P^C_{0,1}$ (Cartesian form of $\widehat\Theta_1=0$) | derived from 3 and 4, Step 2 |
| 7 | All $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$ commute with $Q$ | assumption; charge rule derived, Step 5 |
| 8 | $f_\pi$ independent of $s$ and $p$: no $s$–meson or $p$–meson kinetic coupling | assumption; $M_0$ has no momentum dependence |
| 9 | One-loop tadpoles: those of $\pi^0$, $\eta$, $\eta'$ and $p$ vanish under CP; the one-loop $s$ tadpole is set to zero (vacuum condition $\langle{\rm Re}\,\varphi\rangle=f_a$ at one loop) | assumption |
| 10 | CP with $s$ even, $p$ and $\Phi$ odd: $\widehat S^C_{r,s}=(-1)^s\widehat S^{C\,{\rm T}}_{r,s}$, $\widehat P^C_{r,s}=(-1)^{s+1}\widehat P^{C\,{\rm T}}_{r,s}$; flavour-diagonal sources and quark masses | CP conditions derived, Step 4; flavour-diagonality is an assumption |

The numbers refer to the assumption table of the [polar tree and legs record](../one-loop-estimate/murayama-sigma-pipi-loop-tree-and-legs.md).

**Remark on assumptions 4 and 6.** In the polar records $\widehat\Theta=\Theta+\int_0^x\langle(C_R-C_L)/2\rangle$ ([frame record](../one-loop-estimate/murayama-sigma-pipi-loop-frame-invariance.md)) collects the ALP dependence of the anomalous singlet term and the trace of the axial derivative coupling. The Cartesian Lagrangian above has neither: the anomalous term does not depend on $s$ or $p$ (assumption 3), and $s$, $p$ have no derivative coupling to the mesons (assumption 4). So the only $p$–meson quadratic entry is $\epsilon\langle\mathcal C\lambda_x\rangle$ with $\mathcal C=B_0\widehat P^C_{0,1}$ (Step 2).

**Remark on assumption 9.** This is the Cartesian form of the polar condition "$\sigma$ tadpole $=0$". The two vacuum conditions are not the same choice of parameters at one loop: they differ by a term $\langle p^2\rangle/(2f_a)$ (statement of decision D-A1 of this change, not derived here). The size of this difference belongs to the limits task. Nothing about it is computed here.

**Use in the trace algebra.** As in the polar estimate, assumption 10 is not imposed as matrix identities: the sources stay generic and noncommuting. $\widehat P^C_{0,0}=\widehat P^C_{1,0}=0$ are applied as matrix identities, and $\langle\widehat S^C_{0,1}\rangle=\langle\widehat S^C_{1,1}\rangle=0$ as trace identities.

## Verification

All checks are exact symbolic identities in the script and the tests. No random or numerical matrix is used.
- The stationarity condition is solved from the $s$ one-point function. The two-point functions after stationarity equal the forms of Step 2, term by term: the $\xi_x\xi_y$ term, the $p$–meson terms, and the meson–meson term $-\frac12\langle\lambda_x\mathbb M_{k^2}(\lambda_y)\rangle$.
- Meson stationarity, $D_a$, $N=D_a'$ through the derivative of the defining relation and self-adjointness, and $N$ equal to the kinetic norm.
- CP: the vertex of the transformed Lagrangian equals the vertex with $X\to\eta X^{\rm T}$ for every pattern with $n\le4$. The test also checks that the wrong sign for $\widehat P^C$ fails.
- Charge rule for words of length 1 to 4, and $[Q,\mathbb M_{k^2}(X)]=\mathbb M_{k^2}([Q,X])$.
- $M_0$ from the closed vertex rule equals the direct expansion and the compact form. The up-only value is obtained in two ways.
- A test checks once, with explicit $U(3)$ generators and symbolic $3\times3$ matrices (never printed), the normalization $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$ and the completeness identities $\sum_i\langle A\lambda_i\rangle\langle\lambda_iB\rangle=2\langle AB\rangle$ and $\sum_i\langle A\lambda_iB\lambda_i\rangle=2\langle A\rangle\langle B\rangle$.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_tree.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_tree
```

[Script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_tree.py), [tests](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_tree.py), [Cartesian algebra](../../../../../../scripts/sigma_pipi_cartesian_loop_algebra.py) (generic trace tools reused from [`sigma_pipi_loop_algebra.py`](../../../../../../scripts/sigma_pipi_loop_algebra.py); the generic charge-rule proof and projector reduction come from [`derive_sigma_pipi_loop_tree.py`](../../../../../../scripts/derive_sigma_pipi_loop_tree.py), whose polar vertex functions are not used).
