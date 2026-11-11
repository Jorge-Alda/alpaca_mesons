---
title: Finite ALP-meson mixing and one-ALP amplitudes
status: derived
created: 2026-11-02
updated: 2026-11-11
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

Use the CCWZ building blocks and spurion transformations in the [Appendix-A source record](../../extracts/Alda:2026nrd/appendix-a-expansion.md), the topological term and singlet convention in the [LO topological source record](../../extracts/Alda:2026nrd/lo-topological-action.md), and the retained LO/NLO operators in the [NLO Lagrangian](../../nlo/lagrangian.md). The full local NLO quadratic expansion is in [the quadratic Lagrangian record](../../nlo/quadratic-lagrangian.md). Loops and operators beyond that retained basis are not included here.

## Chiral covariance

The LO theory is expressed through the CCWZ building blocks $u_\mu$, $\chi_\pm$, and $X_{TC}$. Under a chiral transformation they transform as $u_\mu\to h u_\mu h^\dagger$, $\chi_\pm\to h\chi_\pm h^\dagger$, and $X_{TC}\to X_{TC}$, with external sources and the anomalous singlet source transformed consistently as recorded in the [Appendix-A extract](../../extracts/Alda:2026nrd/appendix-a-expansion.md) and [topological-source extract](../../extracts/Alda:2026nrd/lo-topological-action.md). Consequently, $\langle u_\mu u^\mu\rangle$, $\langle\chi_+\rangle$, and $X_{TC}^2$ are unchanged. The field redefinition below is treated in these same CCWZ variables.

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

Set $P_i=P'_i+v_i a$. This makes the quadratic kinetic terms canonical and gives

$$
M_{P'P'}^2=M_{PP}^2,\qquad
M_{P'\hat a}^2=M^2_{i\hat a}+M^2_{ij}v_j,\qquad
M_{\hat a\hat a,\rm can}^2=M^2_{\hat a\hat a}+2v_iM^2_{i\hat a}+v_iM^2_{ij}v_j.
$$

The exact canonical mass matrix is the full matrix with these blocks, and must be diagonalized as a whole. In particular, $m_{a,0}^2$ alone is not the physical ALP mass, and $M^2_{i\hat a}+M^2_{ij}v_j$ need not be small when $f_a\simeq f_\pi$. Stability requires the full mass Hessian to be nonnegative (strictly positive for massive, nondegenerate poles); a positive bare $m_{a,0}^2$ alone does not ensure this for arbitrary source choices.

For an isolated single meson coordinate with diagonal mass $m_P^2$, the two exact eigenvalues are

$$
m_\pm^2=\frac{m_P^2+M_{\hat a\hat a,\rm can}^2}{2}
\pm\frac12\sqrt{(M_{\hat a\hat a,\rm can}^2-m_P^2)^2+4(M_{P\hat a,\rm can}^2)^2}.
$$

At degeneracy the splitting is $2|M_{P\hat a,\rm can}^2|$, so the small-mixing expansion fails precisely where the levels approach. In the full neutral sector this is replaced by diagonalizing the complete $(P_i,a)$ matrix; a particular pion can receive a mass shift only through the components and source structure that mix with it.

For the physical pole $s=m_a^2$, write the eigenvector in Alda:2026nrd's normalization as $(w_i^{(a)}(s),1)^T$. It obeys the exact mixed-field equation

$$
\begin{aligned}
(sK_{ij}-M^2_{ij})w_j^{(a)}(s)&+(sK_{i\hat a}-M^2_{i\hat a})=0,\\
w_i^{(a)}(s)&=[(M^2-sK)^{-1}]_{ij}(sK_{j\hat a}-M^2_{j\hat a}).
\end{aligned}
$$

When the inverse exists, this gives the meson overlap exactly; at exceptional points use the full generalized eigenproblem. The pole normalization is

$$
K_{\hat a\hat a}+2K_{i\hat a}w_i^{(a)}(s)+K_{ij}w_i^{(a)}(s)w_j^{(a)}(s).
$$

For a simple stable pole, its external-leg factor is the inverse square root of this expression, evaluated at $s=m_a^2$.

## One-ALP amplitudes

At finite mixing, use the exact ALP-pole eigenvector $w_i^{(a)}(m_a^2)$ from the preceding section. The general one-leg amplitude, with all other external legs already projected onto their physical states, is

$$
\mathcal M(a_{\rm phys},X)=
\frac{\mathcal V_{\hat a}(m_a^2;X)+\mathcal V_i(m_a^2;X)w_i^{(a)}(m_a^2)}
{\sqrt{K_{\hat a\hat a}+2K_{i\hat a}w_i^{(a)}+K_{ij}w_i^{(a)}w_j^{(a)}}}.
$$

Here $\mathcal V_{\hat a}$ and $\mathcal V_i$ are amputated interaction vertices with one $a$ or $P_i$ leg. This formula applies to an incoming or outgoing ALP, with the usual crossing of momenta and complex conjugation.

Now make a finite ALP-dependent axial change of variables with group element

$$
L_a=\exp\!\left[-\frac{ia}{2f_\pi}w_i^{(a)}(m_a^2)\lambda_i\right],
\qquad R_a=L_a^\dagger.
$$

In CCWZ language, the transformed building blocks obey

$$
u_\mu' = h u_\mu h^\dagger,\qquad
\chi_\pm'=h\chi_\pm h^\dagger,\qquad
X_{TC}'=X_{TC},
$$

where $h$ is the compensator associated with the composed chiral transformation. The external left and right sources transform with the local chiral connection terms, and the anomalous singlet source includes the axial Jacobian. Therefore the LO invariants are unchanged:

$$
\langle u_\mu'u'^{\mu}\rangle=\langle u_\mu u^{\mu}\rangle,\qquad
\langle\chi_+'\rangle=\langle\chi_+\rangle,\qquad
(X_{TC}')^2=X_{TC}^2.
$$

This establishes covariance without using $U$ or $\chi$ as working variables. The mixing vector $w_i^{(a)}(m_a^2)$ is kept exact. Only the coordinate map is expanded about the vacuum to determine the quadratic Hessian: with $\phi=\Phi/(2f_\pi)$,

$$
\phi'=\phi-\frac{a}{2f_\pi}w_i^{(a)}(m_a^2)\lambda_i
+O(\text{fields}^2),
\qquad P'_i=P_i-w_i^{(a)}(m_a^2)a+O(\text{fields}^2).
$$

This truncation is sufficient for the quadratic blocks because their Hessian at the vacuum depends only on the linear part of a field redefinition. It is not an expansion in the mixing or in $f_\pi/f_a$. For interaction amplitudes with several mesons, higher terms in the CCWZ coordinate map generate additional vertices. The complete transformed direct vertex must include them; retaining only the tangent-space shift would generally omit contact contributions. If the other external legs are also expressed in primed coordinates, their projections must be transformed as well.

In the following block formulas, $w_i^{(a)}$ means $w_i^{(a)}(m_a^2)$. The finite linear map gives the following quadratic-block transformation (the displayed blocks are exact for the quadratic Lagrangian):

$$
\begin{aligned}
K'_{ij}&=K_{ij}, & K'_{i\hat a}&=K_{i\hat a}+K_{ij}w_j^{(a)},
&K'_{\hat a\hat a}&=K_{\hat a\hat a}+2w_i^{(a)}K_{i\hat a}+w_i^{(a)}K_{ij}w_j^{(a)},\\
M_{ij}^{2\prime}&=M^2_{ij}, & M_{i\hat a}^{2\prime}&=M^2_{i\hat a}+M^2_{ij}w_j^{(a)},
&M_{\hat a\hat a}^{2\prime}&=M^2_{\hat a\hat a}+2w_i^{(a)}M^2_{i\hat a}+w_i^{(a)}M^2_{ij}w_j^{(a)}.
\end{aligned}
$$

Consequently, the mixed inverse-propagator block transforms as

$$
sK'_{i\hat a}-M_{i\hat a}^{2\prime}
=(sK_{i\hat a}-M^2_{i\hat a})+(sK_{ij}-M^2_{ij})w_j^{(a)}.
$$

The ALP-pole equation is $(m_a^2K_{ij}-M^2_{ij})w_j^{(a)}+(m_a^2K_{i\hat a}-M^2_{i\hat a})=0$. Therefore

$$
m_a^2K'_{i\hat a}-M_{i\hat a}^{2\prime}=0,\qquad
w_i^{(a)\prime}(m_a^2)=w_i^{(a)}(m_a^2)-w_i^{(a)}(m_a^2)=0.
$$

Thus the mixing vanishes on the chosen pole in the primed coordinates. The kinetic and mass mixed blocks need not vanish separately, nor does this choice diagonalize the propagator at every $s$; it aligns the fields with this one exact pole.

For the interaction source, write the amputated term linear in the external leg as

$$
\mathcal V^{(1)}=a\mathcal V_{\hat a}+P_i\mathcal V_i.
$$

Substituting $P_i=P'_i+w_i^{(a)}a$ gives

$$
\mathcal V'_{\hat a}=\mathcal V_{\hat a}+w_i^{(a)}\mathcal V_i,
\qquad \mathcal V'_i=\mathcal V_i.
$$

For a vertex linear in the external fields, the linear part of the redefinition gives the direct primed source

$$
\mathcal V'_{\hat a}(m_a^2;X)
=\mathcal V_{\hat a}(m_a^2;X)+\mathcal V_i(m_a^2;X)w_i^{(a)}(m_a^2).
$$

The normalization transforms into the ALP diagonal kinetic entry:

$$
K'_{\hat a\hat a}
=K_{\hat a\hat a}+2K_{i\hat a}w_i^{(a)}+K_{ij}w_i^{(a)}w_j^{(a)}.
$$

For the full chiral transformation, define $\mathcal V_{\hat a,\mathrm{full}}'$ by substituting the exact exponential map and transformed spurions into the complete LO interaction Lagrangian, then taking the amputated vertex with one $a$ leg. This includes the nonlinear contact vertices mentioned above. The physical amplitude in primed variables is

$$
\boxed{\displaystyle
\mathcal M(a_{\rm phys},X)=
\frac{\mathcal V'_{\hat a,\mathrm{full}}(m_a^2;X)}
{\sqrt{K'_{\hat a\hat a}}}.
}
$$

In CCWZ terms, construct $\mathcal L'_{\rm LO}=\mathcal L_{\rm LO}[u_\mu',\chi_+',X_{TC}']$ with the transformed sources, then differentiate with respect to the external $a$ and expand in the meson coordinates needed for $X$. The chain rule differentiates the complete building blocks, including their nonlinear dependence on $\phi$; this generates the contact terms. Thus $\mathcal V_{\hat a}+w_i^{(a)}\mathcal V_i$ is the tangent-space one-leg contribution, while the full CCWZ prescription supplies the rest for multi-meson channels.

The quadratic pole alignment and normalization are exact. In the transformed interaction, all meson contributions are now part of the direct physical-ALP vertex; this includes the explicit linear contribution $\mathcal V_{\hat a}+w_i^{(a)}\mathcal V_i$ and, for multi-meson channels, the nonlinear contact terms fixed by the exponential map. Substituting the invertible chiral field redefinition into the full covariant LO Lagrangian leaves the generating functional unchanged (with the external sources and anomalous singlet source transformed as stated), so its on-shell S-matrix is unchanged. This is the finite-mixing counterpart of Alda:2026nrd's no-mixing organization. It does not permit retaining the original direct ALP vertices while discarding meson vertices. At first order, the transformed source reduces to the local prescription organized in [the LSZ invariant-split record](../../nlo/lsz-invariant-split.md); at finite mixing it depends on the full pole vector and the complete transformed interactions.

## What changes relative to the small-mixing regime

- All neutral pseudoscalar masses and eigenstates in the mixed sector must be recomputed from the full generalized eigenproblem. Charged states that do not share the relevant conserved quantum numbers have no bilinear mixing, although their parameters can still be affected if the dynamical spurion changes the invariant potential.
- Each one-ALP decay or production amplitude receives coherent contributions from every meson component with the same quantum numbers. Interference can enhance, suppress, or cancel the direct ALP contribution.
- Couplings inferred as linear in $1/f_a$ from an unmixed ALP field need not describe the physical pole when $f_a$ is comparable to $f_\pi$; both the pole residue and the state composition change.
- The chiral expansion is still an expansion in momenta and spurion order. Resumming a quadratic mixing matrix does not resum omitted chiral operators or loops, and large mixing does not guarantee chiral convergence.

The illustrative paper [@DiLuzio:2026poh] shows this distinction concretely: the phase of a dynamical PQ scalar is part of $M_q$; its CP-odd Hessian mixes with $\pi^0,\eta,\eta'$ and changes the neutral-pion mass. That leading mass shift follows from the model-specific spurion and potential. It is not a universal consequence of merely setting $f_a\sim f_\pi$ in the generic external-source Lagrangian.

## Verification status

The congruence $D'(s)=T^{\mathsf T}D(s)T$ and the invariance of poles, residues and amplitudes are in the [congruence record](alp-inverse-propagator-congruence.md). The pole-aligned shift, its regularity condition and the exceptional poles are in the [single-pole record](alp-single-pole-alignment.md); the alignment of several poles in the [several-pole record](alp-multi-pole-alignment.md); and the identification of $L_a$ with $g_R$ for constant generators $G_R=-G_L=-(f_a/2f_\pi)w_i\lambda_i$ in the [chiral-generator record](alp-pole-alignment-chiral-generator.md). The newer records write $L_a$ and $R_a$ of this record as $\tilde g_R$ and $\tilde g_L$, respectively.

The LO Hessian blocks and the canonical coordinate change are derived symbolically by [`scripts/derive_alp_lo_finite_mixing.py`](../../../scripts/derive_alp_lo_finite_mixing.py). The general pole and LSZ identities remain in [`scripts/derive_alp_finite_mixing.py`](../../../scripts/derive_alp_finite_mixing.py). No numerical spectrum or channel-specific width is supplied; those require numerical source matrices, quark masses, and $m_{a,0}$.
