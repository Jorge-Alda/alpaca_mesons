---
title: External legs, kinematics, charge rule and tree amplitude for the one-loop sigma to pi+ pi- estimate
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# External legs, kinematics, charge rule and tree amplitude of $\sigma\to\pi^+\pi^-$

Tasks 1.1 and 1.2 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama`. No loop is evaluated.

## Conventions

Natural units, GeV. Conventions and notation of Alda:2026nrd and DiLuzio:2026poh, as in [the general-source radial record](radial-decay-polar-cartesian-general-sources.md): polar fields in the connection-free frame, $\epsilon=f_\pi/f_a$, $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$, $\chi=2B_0\widehat M$. The Lagrangian is
$$
\mathcal L=\frac{f_\pi^2}{4}\langle u_\mu u^\mu\rangle+B_0f_\pi^2\sum_{r,s}\Big(\frac{\sigma}{f_a}\Big)^r\Big(\frac{a}{f_a}\Big)^s\Big[\big\langle\widehat S_{r,s}\cos\tfrac{\Phi}{f_\pi}\big\rangle+\big\langle\widehat P_{r,s}\sin\tfrac{\Phi}{f_\pi}\big\rangle\Big]-\frac{m_1^2}{12}\langle\Phi\rangle^2+\frac12(\partial\sigma)^2+\frac12\Big(1+\frac{\sigma}{f_a}\Big)^2(\partial a)^2-V_\varphi ,
$$
with $V_\varphi=\frac{m_\varphi^2}{2}(f_a+\sigma)^2+\frac{\lambda}{4}(f_a+\sigma)^4$. The source part is $\frac{f_\pi^2}{4}\langle\chi_+\rangle$ written with $U=e^{i\Phi/f_\pi}$. A vertex $V(x_1,\dots,x_n)$ is the $n$-th derivative of $\mathcal L$ with respect to the fields of its legs, with $\partial_\mu\to-ik_\mu$ and all momenta incoming, so that $M=V$ for a contact graph.

## Legs and kinematics

The legs follow [the estimate record](murayama-radial-up-only-estimate-result.md): a leg is a pair $(\lambda_x,\xi_x)$ or a radial leg $\sigma$.

- $\sigma$: radial leg, momentum $p_\sigma$. It has no meson component, because $\widehat P_{1,0}=0$ and $\langle\widehat S_{1,1}\rangle=0$ (assumption 2 below).
- $\pi^+$: $(\lambda_{\pi^+},0)$ with momentum $p_1$; $\pi^-$: $(\lambda_{\pi^-},0)$ with momentum $p_2$. The matrices $\lambda_{\pi^\pm}$ are treated as opaque charge eigen-matrices, $[Q,\lambda_{\pi^\pm}]=\pm\lambda_{\pi^\pm}$. They may contain $K^\pm$ components. A nondegenerate charged pole is assumed.

Momentum conservation $p_\sigma+p_1+p_2=0$ with $p_\sigma^2=m_\sigma^2$ and $p_1^2=p_2^2=m_{\pi^\pm}^2$ gives exactly
$$
p_1\cdot p_2=\frac{m_\sigma^2-2m_{\pi^\pm}^2}{2},\qquad p_\sigma\cdot p_1=p_\sigma\cdot p_2=-\frac{m_\sigma^2}{2}.
$$
The physical decay has outgoing momenta $-p_1$, $-p_2$, and the dot products are the same.

## Charge rule

Assumption B1: all $\widehat S_{r,s}$, $\widehat P_{r,s}$ commute with the quark charge matrix $Q$ ([charged decoupling record](alp-lo-charged-meson-decoupling.md)). Then:

1. For letters with $[Q,X_k]=q_kX_k$ and $W=X_1\cdots X_n$, $[Q,W]=(q_1+\dots+q_n)W$ and $\langle[Q,W]\rangle=0$. Hence $\langle W\rangle=0$ whenever $q_1+\dots+q_n\neq0$. The script verifies both identities for words of length one to four with symbolic charges.
2. $\mathbb M_s(X)=B_0\{\widehat S_{0,0},X\}+\tfrac{m_1^2}{3}\langle X\rangle\mathbf 1-sX$ commutes with $[Q,\cdot\,]$ when $[Q,\widehat S_{0,0}]=0$. The script verifies this for symbolic noncommuting $X$. So $\mathbb M_s$ maps $Q$-commuting matrices onto themselves, and so does its inverse. Therefore $\mathbb M^{-1}_{m_a^2}(\mathcal C)$, with $\mathcal C=B_0\widehat P_{0,1}$, is neutral.

The rule is applied with charge $+1$ for $\lambda_{\pi^+}$, $-1$ for $\lambda_{\pi^-}$, and $0$ for the sources and for $\mathbb M^{-1}_{m_a^2}(\mathcal C)$. Traces with nonzero total charge vanish. The rule is applied only after the generator sums of the internal meson lines have been carried out.

## Tree amplitude

From the general vertex rule ([vertex-rule record](murayama-sigma-pipi-loop-vertex-rules.md)), with the assumptions and the charge rule applied:
$$
M_0=V(\sigma;\pi^+,\pi^-)=-\frac{B_0\epsilon}{2f_\pi}\big\langle\widehat S_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\big\rangle .
$$
This is the $\widehat S_{1,0}$ term of the polar cubic Lagrangian of [the general-source radial record](radial-decay-polar-cartesian-general-sources.md), Step 3. The $\widehat P_{1,1}$ and $\widehat S_{1,2}$ terms need an ALP component, which the charged legs do not have.

**Up-only check.** Insert $\widehat S_{1,0}=m_uP_u$ and the relations $\lambda_{\pi^+}\lambda_{\pi^-}=2P_u$, $\lambda_{\pi^-}\lambda_{\pi^+}=2P_d$, $P_u^2=P_u$, $P_d^2=P_d$, $P_uP_d=P_dP_u=0$, $\langle P_u\rangle=\langle P_d\rangle=1$. Here $P_u$ and $P_d$ are the up and down projectors, and no explicit matrix is used. This gives $M_0=-B_0\epsilon m_u/f_\pi=-B_0m_u/f_a$. The result equals minus the charged-pion vertex of [`derive_murayama_sigma_pp.py`](../../../scripts/derive_murayama_sigma_pp.py), which uses the convention $\mathcal L\supset-g\,\sigma\pi^+\pi^-$. The comparison substitutes $f_\varphi=f_\pi/\epsilon$ and is an exact symbolic identity.

## Assumptions (task 1.2)

| assumption | stated in |
|---|---|
| 1. No tadpoles: $\widehat P_{0,0}=0$, $\langle\widehat S_{0,1}\rangle=0$ | [general-source radial record](radial-decay-polar-cartesian-general-sources.md), assumption 1 |
| 2. No quadratic radial mixing: $\widehat P_{1,0}=0$, $\langle\widehat S_{1,1}\rangle=0$ ($\langle\widehat S_{1,1}\rangle=0$ also follows from assumption 10) | same record, assumption 2 |
| 3. Anomalous singlet term independent of $\sigma$ and $a$ | same record, assumption 3 |
| 4. Connection-free polar frame, leading $U(3)$ chiral order | same record, assumption 4 |
| 5. Regular simple poles, $\mathbb M_{m_a^2}$ invertible | same record, assumption 5 |
| 6. $\widehat\Theta_1=0$, so $\mathcal C=B_0\widehat P_{0,1}$ | [implicit crossing record](murayama-radial-crossing-implicit.md), Model; frame independence in [the frame record](murayama-sigma-pipi-loop-frame-invariance.md) |
| 7. Sources commute with $Q$ (B1) | [charged decoupling record](alp-lo-charged-meson-decoupling.md) |
| 8. $f_\pi$ does not depend on $\sigma$: no $\sigma$-meson kinetic coupling (B4) | this record (assumption) |
| 9. One-loop tadpoles: those of $\pi^0$, $\eta$, $\eta'$ and of the ALP field vanish exactly under assumption 10; the $\sigma$ tadpole is set to zero by assumption (a vacuum condition on $f_a$, consistent with eliminating $\lambda$ through $m_\sigma^2$) | this record (assumption); [selection record](murayama-sigma-pipi-loop-graph-selection.md), rule 5 |
| 10. D11: the source terms are CP-invariant, with $\sigma$ CP-even and $a$, $\Phi$ CP-odd ($\Phi\to-\Phi^{\rm T}$), so $\widehat S_{r,s}=(-1)^s\widehat S_{r,s}^{\rm T}$, $\widehat P_{r,s}=(-1)^{s+1}\widehat P_{r,s}^{\rm T}$ and $\langle\widehat S_{r,s}\rangle=0$ for odd $s$; $\sigma$ does not mix with the ALP or any neutral pseudoscalar at any loop order. The sources and the quark mass matrix are flavour-diagonal, so $\pi^\pm$ does not mix with another charged pole (this implies assumption 7). Not imposed as matrix identities in the trace algebra; used only to remove the one-particle-reducible mixing graphs of the external legs | [selection record](murayama-sigma-pipi-loop-graph-selection.md), rule 4 (assumption) |

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_tree.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_tree
```

[Script](../../../scripts/derive_sigma_pipi_loop_tree.py), [tests](../../../scripts/test_derive_sigma_pipi_loop_tree.py), [shared algebra](../../../scripts/sigma_pipi_loop_algebra.py). All tests are exact symbolic identities.
