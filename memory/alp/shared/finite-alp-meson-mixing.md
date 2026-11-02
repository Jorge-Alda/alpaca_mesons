---
title: Finite ALP-meson mixing and one-ALP amplitudes
status: derived
created: 2026-11-02
updated: 2026-11-02
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Finite ALP-meson mixing and one-ALP amplitudes

## Scope and inputs

This record specializes to the leading-order chiral Lagrangian of Alda:2026nrd, with a canonical bare ALP quadratic term

$$
\mathcal L_{a,0}=\frac12(\partial_\mu a)^2-\frac12m_{a,0}^2a^2.
$$

Use a CP-conserving stationary vacuum, a complete normalized Hermitian $U(N_f)$ generator basis with $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, and Hermitian $c_q^A,c_q^P$. The quark spurion is linear in $a$ in the chosen source coordinates; no independent $a^2$ quark-spurion term is added. Thus the LO kinetic and anomaly squares are kept exactly, while $c_q^S,c_q^V$ do not contribute to the quadratic blocks. This is an exact diagonalization of the quadratic theory at LO, not a resummation of higher chiral orders. The cited paper [@DiLuzio:2026poh] remains a distinct example with a dynamical PQ scalar inside the quark-mass spurion.

Use the CCWZ building blocks and spurion transformations in the [Appendix-A source record](../../extracts/Alda:2026nrd/appendix-a-expansion.md), the topological term and singlet convention in the [LO topological source record](../../extracts/Alda:2026nrd/lo-topological-action.md), and the retained LO/NLO operators in the [NLO Lagrangian](../../../nlo/lagrangian.md). The full local NLO quadratic expansion is in [the quadratic Lagrangian record](../../../nlo/quadratic-lagrangian.md). Loops and operators beyond that retained basis are not included here.

## Chiral covariance

The starting point is the invariant chiral Lagrangian, not a component mixing angle. Under $U(3)_L\times U(3)_R$, take $U\to LUR^\dagger$ and the mass spurion $\chi\to L\chi R^\dagger$; external sources and the CCWZ blocks transform as specified in the Appendix-A record. If the ALP phase appears in a dynamical spurion, its transformation must be included in the spurion rule. Thus $\langle u_\mu u^\mu\rangle$, $\langle\chi_+\rangle$, and the topological invariant remain covariant/invariant before any expansion. A linear field redefinition among $a$ and pseudoscalar coordinates is a change of coordinates on this same theory and does not break chiral covariance.

## Exact quadratic theory

Expand the invariant LO Lagrangian plus $\mathcal L_{a,0}$ about the stationary vacuum, without expanding in $f_\pi/f_a$. With $P_i$ the Goldstone coordinates, the quadratic terms are

$$
\mathcal L_2=\frac12\partial_\mu\Phi_I K_{IJ}\partial^\mu\Phi_J-\frac12\Phi_I M^2_{IJ}\Phi_J,\qquad \Phi_I=(P_i,a).
$$

For the neutral sector, the blocks in this generator basis are

$$
K_{ij}=\delta_{ij},\qquad
K_{ia}=-v_i,\qquad
K_{aa}=1+v_i v_i,
\qquad
v_i=\frac{f_\pi}{f_a}\operatorname{tr}(c_q^A\lambda_i),
$$

$$
M^2_{ij}=\frac14\operatorname{tr}\!\left(\{\lambda_i,\lambda_j\}\chi_0\right)
+\frac{m_0^2}{2N_f}\operatorname{tr}(\lambda_i)\operatorname{tr}(\lambda_j),
\qquad \chi_0=2B_0m_q,
$$

$$
M^2_{ia}=-\frac{f_\pi}{2f_a}\operatorname{tr}\!\left(\lambda_i\{c_q^P,\chi_0\}\right)
+\frac{f_\pi m_0^2c_G}{N_ff_a}\operatorname{tr}(\lambda_i),
\qquad
M^2_{aa}=m_{a,0}^2+\frac{2f_\pi^2m_0^2c_G^2}{N_ff_a^2}.
$$

The second term in $M^2_{ij}$ and the $c_G$ terms come from the LO topological invariant. The last term in $M^2_{aa}$ is induced by its square; the bare ALP mass is only one contribution to the total diagonal entry. Likewise the chiral kinetic term adds $v_iv_i$ to $K_{aa}$, so the full ALP coordinate is not canonically normalized before diagonalization. Charged Goldstones have no bilinear mixing with a flavor-neutral ALP source. These formulas use the project’s Appendix-A/topological-source conventions [@Alda:2026nrd].

The full kinetic matrix is positive definite for real $v_i$, since

$$
\mathcal L_{\rm kin}^{(2)}=\frac12(\partial_\mu P_i-v_i\partial_\mu a)^2+\frac12(\partial_\mu a)^2.
$$

Set $P_i=Q_i+v_i a$. This preserves chiral covariance as a change of field coordinates, makes the quadratic kinetic terms canonical, and gives

$$
M_{Q Q}^2=M_{PP}^2,\qquad
M_{Qa}^2=b_i=M^2_{ia}+M^2_{ij}v_j,\qquad
M_{aa,\rm can}^2=d=M^2_{aa}+2v_iM^2_{ia}+v_iM^2_{ij}v_j.
$$

Thus the exact canonical mass matrix is $\begin{pmatrix}M^2_{PP}&b\\b^T&d\end{pmatrix}$. It must be diagonalized as a whole. In particular, $m_{a,0}^2$ alone is not the physical ALP mass, and the off-diagonal block $b$ need not be small when $f_a\simeq f_\pi$. Stability requires the full mass Hessian to be nonnegative (strictly positive for massive, nondegenerate poles); a positive bare $m_{a,0}^2$ alone does not ensure this for arbitrary source choices.

For an isolated single meson coordinate with diagonal mass $m_P^2$, the two exact eigenvalues are

$$
m_\pm^2=\frac{m_P^2+d}{2}\pm\frac12\sqrt{(d-m_P^2)^2+4b^2}.
$$

At degeneracy the splitting is $2|b|$, so the small-mixing expansion in $b/(d-m_P^2)$ fails precisely where the levels approach. In the full neutral sector this is replaced by diagonalizing the complete $(P_i,a)$ matrix; a particular pion can receive a mass shift only through the components and source structure that mix with it.

Equivalently, before the kinetic-coordinate change, solve the generalized eigenproblem

$$
M^2 v^{(n)}=m_n^2 K v^{(n)},\qquad v^{(n)T}K v^{(m)}=\delta_{nm}.
$$

If $O$ is the orthogonal matrix that diagonalizes the canonical matrix, the pole masses are its eigenvalues and the canonical fields are $(Q_i,a)^T=O_{In}\varphi_n$. Equivalently, for a simple pole $s=m_n^2$, the meson components in the original coordinates relative to an $a$ component equal to one are

$$
P_i=(O_{in}+v_iO_{an})\varphi_n,\qquad a=O_{an}\varphi_n.
$$

This relation displays both the meson content and the overlap of the physical pole with the original ALP field. It is exact for the stated LO quadratic theory.

## One-ALP amplitudes

Let $v^{(a)}$ be the $K$-normalized eigenvector for the physical pole identified continuously with the ALP as the mixing is turned off. For any amputated interaction vertex with one external scalar leg and arbitrary other on-shell legs, the physical amplitude is the coherent projection

$$
\mathcal M(a_{\rm phys},X)=v_I^{(a)}\,\Gamma_I(X),
$$

where $\Gamma_I$ is the set of vertices with one original-coordinate field $\Phi_I$ on that leg, evaluated at the physical pole with the same field coordinates and full quadratic prescription. In the canonical coordinates the expression is

$$
\mathcal M(\varphi_n,X)=O_{an}\Gamma_a(X)+(O_{in}+v_iO_{an})\Gamma_i(X).
$$

For the pole continuously connected to the bare ALP, the second term is its meson component and can be order one. This applies to initial or final one-ALP states; crossing changes momenta and conjugation conventions, not the pole eigenvector. If eigenvalues are degenerate or nearly degenerate, use the full pole subspace instead of assigning a unique ALP label by continuity. All other external legs must also be projected onto their physical eigenvectors.

The result is invariant under a nonsingular linear field-coordinate change when $K,M^2$, vertices, and external eigenvectors are all transformed together. In particular, removing $K_{i\hat a}$ and then keeping only the shifted direct vertex is incorrect: the mass mixing and the transformed interaction vertices must also be retained. At first order in mixing, the general expression reduces to the direct-plus-meson-pole LSZ organization in [the invariant LSZ split](../../../nlo/lsz-invariant-split.md); that record explicitly assumes first order and is not itself a finite-mixing result.

## What changes relative to the small-mixing regime

- All neutral pseudoscalar masses and eigenstates in the mixed sector must be recomputed from the full generalized eigenproblem. Charged states that do not share the relevant conserved quantum numbers have no bilinear mixing, although their parameters can still be affected if the dynamical spurion changes the invariant potential.
- Each one-ALP decay or production amplitude receives coherent contributions from every meson component with the same quantum numbers. Interference can enhance, suppress, or cancel the direct ALP contribution.
- Couplings inferred as linear in $1/f_a$ from an unmixed ALP field need not describe the physical pole when $f_a$ is comparable to $f_\pi$; both the pole residue and the state composition change.
- The chiral expansion is still an expansion in momenta and spurion order. Resumming a quadratic mixing matrix does not resum omitted chiral operators or loops, and large mixing does not guarantee chiral convergence.

The illustrative paper [@DiLuzio:2026poh] shows this distinction concretely: the phase of a dynamical PQ scalar is part of $M_q$; its CP-odd Hessian mixes with $\pi^0,\eta,\eta'$ and changes the neutral-pion mass. That leading mass shift follows from the model-specific spurion and potential. It is not a universal consequence of merely setting $f_a\sim f_\pi$ in the generic external-source Lagrangian.

## Verification status

The LO Hessian blocks and the canonical coordinate change are derived symbolically by [`scripts/derive_alp_lo_finite_mixing.py`](../../../scripts/derive_alp_lo_finite_mixing.py). The general pole and LSZ identities remain in [`scripts/derive_alp_finite_mixing.py`](../../../scripts/derive_alp_finite_mixing.py). No numerical spectrum or channel-specific width is supplied; those require numerical source matrices, quark masses, and $m_{a,0}$.
