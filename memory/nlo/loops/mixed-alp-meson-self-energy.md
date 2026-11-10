---
title: One-loop mixed ALP–meson self-energy of the U(3) nonet
status: derived
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds, Gasser:1983yg]
provenance: analytical-derivation
---

# Mixed ALP–meson self-energy

## Inputs and assumptions

Use the vertices and topologies of [vertices and topologies](vertices-and-topologies.md), the spectral tadpole contraction and subtraction of [loop functions and subtraction](loop-functions-and-subtraction.md), the leading mass matrix of [neutral meson matrices](../neutral-meson-matrices.md), and the local quadratic terms of [quadratic Lagrangian](../quadratic-lagrangian.md). The conventions are those of Alda:2026nrd: $Q=sK-M^2+\Sigma$, $\Phi=\sum_i\lambda_iP_i$, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$ and $\chi_0=2B_0\operatorname{diag}(m_u,m_d,m_s)$. Masses and decay constants are in GeV and $s$ in GeV$^2$. Work to one loop and first order in $1/f_a$, with arbitrary flavour matrices, no CP-violating sources or constant $\theta$, and no resonance loops ([loop classes](../counting/loop-classes.md)).

## Tadpole contraction

For each external generator $\lambda_i$, the following expression transcribes the three placements of the external meson in the $a\Phi^3$ vertices. Its coincident propagator is understood as the finite spectral sum given in [loop functions and subtraction](loop-functions-and-subtraction.md), including the full LO neutral mixing:

$$
\begin{aligned}
\Sigma_{i\hat a}(s)={}&\sum_{c,d}\langle\phi_c\phi_d\rangle_{\rm fin}\,\bigg[\frac{s}{6f_\pi f_a}\operatorname{tr}\!\left(c_q^A[\lambda_c\lambda_d\lambda_i+\lambda_i\lambda_c\lambda_d-2\lambda_d\lambda_i\lambda_c]\right)\\
&-\frac{1}{12f_\pi f_a}\operatorname{tr}\!\left(c_q^P\left\{\chi_0,\lambda_i\lambda_c\lambda_d+\lambda_c\lambda_i\lambda_d+\lambda_c\lambda_d\lambda_i\right\}\right)\bigg].
\end{aligned}
$$

The sum runs over all nine real nonet generators. This expression retains non-commuting source matrices and gives the axial and pseudoscalar contributions separately. The $c_q^V$, $c_q^S$, and $c_G$ contributions to this one-loop mixed function vanish under the stated assumptions: their $a\Phi^3$ vertices vanish, and a mixed bubble would require the absent $\Phi^3$ vertex. The flavour-violating cubic $a\Phi^2$ vertices and their role are listed in [vertices and topologies](vertices-and-topologies.md). The $c_G$ leading mixed term itself need not vanish.

The loop expression is combined with the renormalized local terms and the divergent meson field redefinition in [counterterm cancellation](counterterm-cancellation.md). The latter changes the leading mixed block as described in [mixed ultraviolet field renormalization](mixed-uv-field-renormalization.md). The finite contraction is a full-nonet prescription; its comparison with the vertex expansion for arbitrary noncommuting matrices is given in [finite mixed contraction verification](finite-mixed-contraction-verification.md).

## Verification

Run `.venv/bin/python scripts/derive_nlo_loop_self_energies.py --mixed`. It prints the source coefficients for all nine external fields, verifies the absent vertices and the calculated pole cancellation, and independently checks the local mixed poles. The full finite contraction is checked by `scripts/verify_mixed_finite_contractions.py`. These checks are conditional on the stated input operator basis and subtraction.
