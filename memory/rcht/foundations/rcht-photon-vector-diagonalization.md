---
title: Quadratic photon-vector diagonalization in antisymmetric-tensor RChT
status: mixed
created: 2026-09-23
updated: 2026-10-19
bibtex_keys: [Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# Quadratic Photon--Vector Diagonalization in RChT

## Scope and assumptions

This record diagonalizes the tree-level quadratic sector containing the photon and the three neutral members of a vector nonet, including the vacuum-spurion source terms used by the focused fit. The vectors use the antisymmetric-tensor RChT convention and the universal $F_V$ operator recorded in [the RChT summary](rcht-basics.md) [@Ecker:1988te]. The spurion-dependent source extension is an explicit fit hypothesis, not a claim that it is the complete RChT source action. The calculation is performed by [`scripts/derive_rcht_photon_vector_diagonalization.py`](../../../scripts/derive_rcht_photon_vector_diagonalization.py).

Assumptions are:

- natural units and the metric/sign conventions implicit in the recorded Ecker Lagrangian;
- a canonical bare photon term $-F_{\mu\nu}F^{\mu\nu}/4$;
- three real, diagonal, orthonormal $U(3)$ generators $T_i$ for neutral vector mass eigenstates, with distinct input masses $M_i$ allowed;
- a real diagonal vacuum spurion $\chi_+=\operatorname{diag}(x_l,x_l,x_s)$ and the source extension with real $z_V,y_V$ defined below; $[\chi_+,Q]=0$ in this vacuum;
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

The source conventions and universal $F_V$ term are literature inputs [@Ecker:1988te]. The spurion extension is the explicit source hypothesis used in the [focused preliminary fit](../odd-parity/odd-preliminary-fit.md). For that hypothesis, extend the universal term by

$$
\mathcal L_{\rm source}=
\frac{F_V}{2\sqrt2}\langle V_{\mu\nu}f_+^{\mu\nu}\rangle
+\frac{z_V}{2\sqrt2}\langle V_{\mu\nu}\{\chi_+,f_+^{\mu\nu}\}\rangle
+\frac{y_V}{2\sqrt2}\langle V_{\mu\nu}\rangle\langle\chi_+f_+^{\mu\nu}\rangle.
$$

Here $F_V$ has units GeV and $z_V,y_V$ have units GeV$^{-1}$. The single-trace fit sets $y_V=0$; the double-trace fit retains all three terms. These extra terms are the declared fit extension. Evaluating $f_+^{\mu\nu}=-2eQF^{\mu\nu}$ gives the source matrix

$$
S_\gamma=F_VQ+z_V\{\chi_+,Q\}+y_V\langle\chi_+Q\rangle\mathbf1.
$$

For $\chi_+=\operatorname{diag}(x_l,x_l,x_s)$, this is

$$
S_\gamma=\operatorname{diag}\!\left(
\frac{2F_V+4x_lz_V+y_V(x_l-x_s)}{3},
\frac{-F_V-2x_lz_V+y_V(x_l-x_s)}{3},
\frac{-F_V-2x_sz_V+y_V(x_l-x_s)}{3}
\right).
$$

The ensuing diagonalization is derived for this stated quadratic source sector.

Write the neutral vectors as

$$
V_{\mu\nu}=\sum_{i=1}^3T_iV^i_{\mu\nu},
\qquad
T_i=\operatorname{diag}(t_{ui},t_{di},t_{si}),
\qquad
\langle T_iT_j\rangle=\delta_{ij}.
$$

No $T_i$ is identified with $\rho^0$, $\omega$, or $\phi$. For the universal term, the overlap is

$$
c_i\equiv\langle QT_i\rangle
=\frac{2t_{ui}-t_{di}-t_{si}}{3},
\qquad \sum_i c_i^2=\langle Q^2\rangle=\frac23.
$$

With the spurion source, the relevant residue is instead

$$
R_i\equiv\langle T_iS_\gamma\rangle
=F_Vc_i+z_V\langle T_i\{\chi_+,Q\}\rangle
+y_V\langle\chi_+Q\rangle\langle T_i\rangle.
$$

For a complete orthonormal diagonal-generator basis, $\sum_iR_i^2=\langle S_\gamma^2\rangle$. The spurion terms give the normalized-singlet residue $\langle S_\gamma\rangle/\sqrt3=(x_l-x_s)(2z_V+3y_V)/(3\sqrt3)$, even though $\langle Q\rangle=0$.

## Tensor mixing term

Using the source matrix above, the neutral mixing term is

$$
\mathcal L_{V\gamma}
=-\frac{e}{\sqrt2}\sum_iR_iV^i_{\mu\nu}F^{\mu\nu}
=-\frac12\sum_i g_iV^i_{\mu\nu}F^{\mu\nu},
$$

where

$$
g_i=\sqrt2\,eR_i,
\qquad
\epsilon_i\equiv\frac{g_i}{M_i}=\frac{\sqrt2\,eR_i}{M_i}.
$$

The universal result follows at $z_V=y_V=0$, when $R_i=F_Vc_i$.

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

produces canonical, fully diagonal massive fields $W^a_\mu$. This construction leaves the generalized residues $R_i$ generic and determines the final electromagnetic coupling of each $W^a$ from

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

The script derives the spurion-dependent residues from $S_\gamma$, checks their universal limit and basis-invariant norm, then derives the tensor elimination and photon-decoupling congruence transformation with $g_i=\sqrt2 eR_i$. It constructs $K_V^{-1/2}$, derives the generic rational and determinant pole equations, checks an electromagnetically decoupled vector direction and the massless transverse photon mode, and verifies $\sum_i c_i^2=2/3$ in the universal limit.
