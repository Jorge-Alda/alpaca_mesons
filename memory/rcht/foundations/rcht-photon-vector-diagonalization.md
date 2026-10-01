---
title: Quadratic photon-vector diagonalization in antisymmetric-tensor RChT
status: derived
created: 2026-09-23
updated: 2026-09-23
bibtex_keys: [Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# Quadratic Photon--Vector Diagonalization in RChT

## Scope and assumptions

This record diagonalizes the tree-level quadratic sector containing the photon and the three neutral members of a vector nonet. The vectors use the antisymmetric-tensor RChT convention and the $F_V$ operator recorded in [the RChT summary](rcht-basics.md) [@Ecker:1988te]. The calculation is performed by [`scripts/derive_rcht_photon_vector_diagonalization.py`](../../../scripts/derive_rcht_photon_vector_diagonalization.py).

Assumptions are:

- natural units and the metric/sign conventions implicit in the recorded Ecker Lagrangian;
- a canonical bare photon term $-F_{\mu\nu}F^{\mu\nu}/4$;
- three real, diagonal, orthonormal $U(3)$ generators $T_i$ for neutral vector mass eigenstates, with distinct input masses $M_i$ allowed;
- a common real $F_V$ coupling; replacing it by $F_{V_i}$ only changes $F_Vc_i\to F_{V_i}c_i$ below;
- no ideal-mixing assignment, loops, widths, or additional electromagnetic counterterms;
- exact algebra in the mixing parameters within this truncated quadratic tree-level Lagrangian. Consequently, the resummed $O(e^2)$ mass and charge shifts are not a complete electromagnetic calculation at that order.

## Electromagnetic source and generic diagonal generators

Take

$$
D_\mu q=(\partial_\mu+ieQA_\mu)q,
\qquad
Q=\operatorname{diag}\!\left(\frac23,-\frac13,-\frac13\right),
$$

so the partonic interaction is $-eA_\mu J_Q^\mu$, with $J_Q^\mu=\bar q\gamma^\mu Qq$. The corresponding external RChT sources are

$$
\ell_\mu=r_\mu=-eQA_\mu.
$$

In the Goldstone vacuum $u=1$,

$$
F_L^{\mu\nu}=F_R^{\mu\nu}=-eQF^{\mu\nu},
\qquad
f_+^{\mu\nu}=-2eQF^{\mu\nu}.
$$

These source conventions and the $F_V$ interaction are literature inputs; the following diagonalization is derived. [@Ecker:1988te]

Write the neutral vectors as

$$
V_{\mu\nu}=\sum_{i=1}^3T_iV^i_{\mu\nu},
\qquad
T_i=\operatorname{diag}(t_{ui},t_{di},t_{si}),
\qquad
\langle T_iT_j\rangle=\delta_{ij}.
$$

No $T_i$ is identified with $\rho^0$, $\omega$, or $\phi$. All dependence on the neutral flavor basis enters through

$$
c_i\equiv\langle QT_i\rangle
=\frac{2t_{ui}-t_{di}-t_{si}}{3}.
$$

For any complete orthonormal basis of diagonal generators,

$$
\sum_i c_i^2=\langle Q^2\rangle=\frac23.
$$

The normalized singlet has zero overlap, $\langle Q\mathbf1/\sqrt3\rangle=0$, because $Q$ is traceless. A generic mass eigenstate may nevertheless couple to the photon when it mixes singlet and traceless directions; this statement does not assume ideal mixing.

## Tensor mixing term

The RChT operator gives

$$
\begin{aligned}
\mathcal L_{F_V}
&=\frac{F_V}{2\sqrt2}\left\langle V_{\mu\nu}f_+^{\mu\nu}\right\rangle\\
&=-\frac{eF_V}{\sqrt2}\sum_i c_iV^i_{\mu\nu}F^{\mu\nu}
=-\frac12\sum_i g_iV^i_{\mu\nu}F^{\mu\nu},
\end{aligned}
$$

where

$$
g_i=\sqrt2\,eF_Vc_i,
\qquad
\epsilon_i\equiv\frac{g_i}{M_i}
=\frac{\sqrt2\,eF_Vc_i}{M_i}.
$$

The $g_i$ have mass dimension one and the $\epsilon_i$ are dimensionless.

## First-order reduction to Proca variables

For each vector introduce an auxiliary field $P^i_\mu$. Up to an integration by parts, the tensor kinetic term has the equivalent first-order form

$$
\mathcal L_i^{(1)}
=\frac12P^i_\mu P_i^\mu
-\frac12P^i_{\mu\nu}V_i^{\mu\nu}
+\frac{M_i^2}{4}V^i_{\mu\nu}V_i^{\mu\nu}
-\frac{g_i}{2}V^i_{\mu\nu}F^{\mu\nu},
$$

with $P^i_{\mu\nu}=\partial_\mu P^i_\nu-\partial_\nu P^i_\mu$. Eliminating $P^i_\mu$ returns the antisymmetric-tensor kinetic term. Instead eliminating $V^i_{\mu\nu}$ gives

$$
V^i_{\mu\nu}=\frac{P^i_{\mu\nu}+g_iF_{\mu\nu}}{M_i^2},
$$

and therefore

$$
\mathcal L_i^{(1)}
=\frac12P^i_\mu P_i^\mu
-\frac{(P^i_{\mu\nu}+g_iF_{\mu\nu})^2}{4M_i^2}.
$$

Defining $P^i_\mu=M_iR^i_\mu$ yields the exact quadratic Proca form

$$
\mathcal L_i
=-\frac14(R^i_{\mu\nu}+\epsilon_iF_{\mu\nu})^2
+\frac12M_i^2R^i_\mu R_i^\mu.
$$

The local term $-\epsilon_i^2F^2/4$ is required by the field transformation; dropping it changes the exact quadratic theory. This illustrates why representation changes require their accompanying local operators. [@Ecker:1989yg]

## Exact photon decoupling

Let $\boldsymbol\epsilon=(\epsilon_1,\epsilon_2,\epsilon_3)^T$ and

$$
Z=1+\boldsymbol\epsilon^T\boldsymbol\epsilon.
$$

Including the bare photon term, the quadratic Lagrangian is

$$
\mathcal L_{\rm quad}
=-\frac14ZF_{\mu\nu}F^{\mu\nu}
-\frac14\sum_iR^i_{\mu\nu}R_i^{\mu\nu}
-\frac12\sum_i\epsilon_iF_{\mu\nu}R_i^{\mu\nu}
+\frac12\sum_iM_i^2R^i_\mu R_i^\mu.
$$

The exact triangular field transformation

$$
A_{\gamma\mu}=\sqrt Z\left(A_\mu+\frac{\boldsymbol\epsilon^T\mathbf R_\mu}{Z}\right),
\qquad
A_\mu=\frac{A_{\gamma\mu}}{\sqrt Z}
-\frac{\boldsymbol\epsilon^T\mathbf R_\mu}{Z},
$$

removes every photon--vector kinetic cross term. The result is

$$
\mathcal L_{\rm quad}
=-\frac14F_{\gamma\mu\nu}F_\gamma^{\mu\nu}
-\frac14\mathbf R_{\mu\nu}^T K_V\mathbf R^{\mu\nu}
+\frac12\mathbf R_\mu^T D\mathbf R^\mu,
$$

where

$$
K_V=\mathbf1-\frac{\boldsymbol\epsilon\boldsymbol\epsilon^T}{Z},
\qquad
D=\operatorname{diag}(M_1^2,M_2^2,M_3^2).
$$

$K_V$ has eigenvalue $1/Z$ along $\boldsymbol\epsilon$ and eigenvalue $1$ in the two orthogonal directions, so it is positive definite. The photon is exactly canonical and massless.

The partonic current interaction becomes

$$
-eA_\mu J_Q^\mu
=-\frac{e}{\sqrt Z}A_{\gamma\mu}J_Q^\mu
+\frac{e}{Z}\boldsymbol\epsilon^T\mathbf R_\mu J_Q^\mu.
$$

Thus $e$ in the starting Lagrangian is a bare parameter and the canonically normalized tree-level photon coupling is

$$
e_\gamma=\frac{e}{\sqrt Z}.
$$

If $e$ is instead defined from the measured electric charge, its counterterm or matching condition must absorb this normalization; one must not simultaneously identify the starting $e$ with the physical charge and retain the finite factor as an additional prediction.

## Canonical massive vectors

An explicit positive inverse square root is

$$
K_V^{-1/2}=\mathbf1
+\frac{\boldsymbol\epsilon\boldsymbol\epsilon^T}{\sqrt Z+1}.
$$

Let $O$ be a real orthogonal matrix satisfying

$$
O^TK_V^{-1/2}DK_V^{-1/2}O
=\operatorname{diag}(m_1^2,m_2^2,m_3^2).
$$

Then

$$
\mathbf R_\mu=K_V^{-1/2}O\,\mathbf W_\mu
$$

produces canonical, fully diagonal massive fields $W^a_\mu$. This construction leaves the generator overlaps $c_i$ generic and determines the final electromagnetic coupling of each $W^a$ from

$$
\frac{e}{Z}\boldsymbol\epsilon^TK_V^{-1/2}O.
$$

For nondegenerate $M_i$ and nonzero $\epsilon_i$, the massive poles can be written as the roots of

$$
\boxed{
1-\sum_{i=1}^3\frac{\epsilon_i^2M_i^2}{s-M_i^2}=0
}
$$

The universally valid equation is

$$
\boxed{\det(sK_V-D)=0.}
$$

The rational form was obtained by dividing by $\prod_i(s-M_i^2)$. If some $\epsilon_i=0$, the corresponding uncoupled root $s=M_i^2$ must be restored; at degenerate masses, modes orthogonal to $\boldsymbol\epsilon$ can likewise remain at the common unmixed mass. The orthogonal matrix within a degenerate subspace is not unique, although physical projectors and current couplings are well defined.

In the physical transverse spin-one sector, after gauge fixing and omission of irrelevant nonzero kinematic factors, the characteristic polynomial is

$$
s\det(sK_V-D).
$$

Its factor $s$ exhibits the massless photon mode. This transverse quadratic check is consistent with gauge invariance but does not replace a Ward-identity check for a complete interacting amplitude.

## Leading electromagnetic order

Since $\epsilon_i=O(e)$,

$$
Z=1+O(e^2),
\qquad
K_V=\mathbf1+O(e^2),
$$

and

$$
A_{\gamma\mu}=A_\mu+\sum_i\epsilon_iR^i_\mu+O(e^2).
$$

Therefore the vector masses first shift at $O(e^2)$. If the calculation retains only one electromagnetic insertion, this linear field redefinition is sufficient. The exact formulas above resum repeated mixing insertions in the truncated tree-level quadratic model; they do not supply the other loops and counterterms required for a complete $O(e^2)$ RChT analysis.

## Reproduction

Run

```text
.venv/bin/python scripts/derive_rcht_photon_vector_diagonalization.py
```

The script derives the tensor elimination, verifies the photon-decoupling congruence transformation, constructs $K_V^{-1/2}$, derives the generic rational and universal determinant pole equations, checks an electromagnetically decoupled vector direction, checks the massless transverse photon mode, and verifies $\sum_i c_i^2=2/3$ for a generic three-angle orthogonal basis of diagonal generators.
