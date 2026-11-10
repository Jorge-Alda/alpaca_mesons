---
title: One-loop derivatives and pole residues
status: derived
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, Gasser:1983yg, Kaiser:1998ds]
provenance: analytical-derivation
---

# One-loop derivatives and pole residues

## Inputs and assumptions

The finite tadpoles and flavour conventions are those of the [meson](meson-self-energy.md), [mixed](mixed-alp-meson-self-energy.md), and [subtraction](loop-functions-and-subtraction.md) records [@Gasser:1983yg; @Kaiser:1998ds]. The inverse-propagator sign and LSZ blocks follow the [LSZ loop record](../loop-corrections-to-alp-meson-lsz-mixing.md) [@Alda:2026nrd]: $Q=sK-M^2+\Sigma$, $A=Q_{PP}$, $b=Q_{P\hat a}$. Masses and decay constants are in GeV, $s$ and $\Sigma$ in GeV$^2$, and $\partial_s\Sigma$ is dimensionless.

Work at one loop and first order in $1/f_a$, with the retained tadpole topologies and no resonance loops; see [vertices and topologies](vertices-and-topologies.md) and [loop classes](../counting/loop-classes.md). The expressions below assume an invertible LO meson block at the ALP pole and nondegenerate LO meson poles. The mixed finite contraction has been independently compared with the vertex expansion for one non-commuting component, rather than for every finite coefficient.

## Derivatives

The external-momentum term in the meson tadpole Lagrangian gives, using $\Phi=\sum_i\lambda_iP_i$ and the finite spectral contractions of the [subtraction record](loop-functions-and-subtraction.md),

$$
\begin{aligned}
\partial_s\Sigma_{ij}={1\over24f_\pi^2}\sum_{c,d}\langle\phi_c\phi_d\rangle_{\rm fin}\,
\operatorname{tr}\!\left(&\lambda_c\lambda_i\lambda_d\lambda_j
+\lambda_c\lambda_j\lambda_d\lambda_i\\
&-\lambda_c\lambda_d\lambda_i\lambda_j
-\lambda_c\lambda_d\lambda_j\lambda_i\right).
\end{aligned}
$$

This is the Hessian of the kinetic term in the [meson self-energy record](meson-self-energy.md). It is constant in $s$. Its singlet row and column vanish. In particular, with the charged-field contractions printed by the script,

$$
\partial_s\Sigma_{38}={\sqrt3\over6f_\pi^2}\left(\langle K^0\bar K^0\rangle_{\rm fin}-\langle K^+K^-\rangle_{\rm fin}\right),\qquad
\partial_s\Sigma_{00}=\partial_s\Sigma_{30}=\partial_s\Sigma_{80}=0.
$$

Differentiating the mixed expression gives

$$
\partial_s\Sigma_{i\hat a}={1\over6f_\pi f_a}\sum_{c,d}\langle\phi_c\phi_d\rangle_{\rm fin}
\operatorname{tr}\!\left(c_q^A[\lambda_c\lambda_d\lambda_i+\lambda_i\lambda_c\lambda_d-2\lambda_d\lambda_i\lambda_c]\right).
$$

Only $c_q^A$ contributes to this derivative: the $c_q^P$ tadpole is constant in $s$, and the $c_q^V$, $c_q^S$ and $c_G$ mixed tadpoles vanish under these assumptions. Both derivatives are independent of $s$; their full flavour entries are printed by the [script](../../../scripts/derive_nlo_loop_self_energies.py).

## LSZ components and pole shifts

Write $A=A_0+\delta A$, $b=b_0+\delta b$ with the local NLO and loop corrections included in $\delta A,\delta b$. At fixed $s$, the LSZ equation gives

$$
w_0^{(a)}=-A_0^{-1}b_0,\qquad
\delta w^{(a)}=-A_0^{-1}[\delta b+\delta A\,w_0^{(a)}].
$$

For the loop contribution alone, $\delta A=\Sigma_{ij}$ and $\delta b=\Sigma_{i\hat a}$. If an independent ALP pole shift is supplied, evaluation at the shifted pole adds

$$
\delta s_a\,\partial_s w_0^{(a)}(s_a^{(0)})=-\delta s_a\,A_0^{-1}(s_a^{(0)})\left(K_{PP}^{(0)}w_0^{(a)}(s_a^{(0)})+K_{P\hat a}^{(0)}\right).
$$

These are the fixed-$s$ and pole-shift terms of the [LSZ loop record](../loop-corrections-to-alp-meson-lsz-mixing.md). The ALP self-energy and mixing-induced ALP pole shift are $\mathcal O(1/f_a^2)$ and are excluded here.

For a unit-normalized LO meson eigenvector, a nondegenerate LO pole and a canonical LO kinetic matrix, the loop contribution to the pole shift is minus the expectation value of $\Sigma$ in that eigenvector, evaluated at the LO pole. To first order, the inverse residue is one plus the expectation value of $\partial_s\Sigma$ in the same eigenvector and at the same pole. Local NLO contributions enter the corresponding expressions with the full $\delta A$ in place of $\Sigma$. For an ALP-like pole, use the full eigenvector $u^{(a)}=(w^{(a)},1)^T$ and the normalization $Z_a^{-1}=u^{(a)T}Q'(s_a)u^{(a)}$; this is distinct from the meson-only expression above. Amplitudes contain both the corrected components and the residue.

## Verification

Run `.venv/bin/python scripts/derive_nlo_loop_self_energies.py --poles`. The script prints the derivatives, verifies the tree-level LSZ limit and the first-order pole shift against a two-meson determinant. The derivation and its focused unit tests are described in [the pole script record](pole-script-derivation.md).
