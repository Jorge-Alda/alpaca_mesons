---
title: Bai vector-only WZW vertices and matching to canonical RChT vectors
status: mixed
created: 2026-09-28
updated: 2026-09-28
bibtex_keys: [Bai:2025fvl, Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg, Alda:2026bar]
provenance: literature
---

# Vector-only background fields

## Bai conventions and exact printed formulae

Bai defines the left/right background one-forms through their vector and axial combinations. With axial-vector fields switched off, the relevant equation is

$$
 \mathbb B_V=\mathbb B_L+\mathbb B_R
 =g\begin{pmatrix}
 \rho_0+\omega&\sqrt2\rho^+&\sqrt2K^*_+\\
 \sqrt2\rho^-&-\rho_0+\omega&\sqrt2K^*_0\\
 \sqrt2K^*_-&\sqrt2\bar K^*_0&\sqrt2\phi
 \end{pmatrix}
 +(\boldsymbol k_L+\boldsymbol k_R)\frac{da}{f_a},
 \qquad
 \mathbb B_A=(\boldsymbol k_L-\boldsymbol k_R)\frac{da}{f_a}.
$$

This is Bai Eq. (2.8), source `main.tex`, lines 257--272. Their normalization is $g\simeq\sqrt{12\pi}$; the fields are fictitious gauge backgrounds above confinement and are promoted to dynamical massive mesons below it. The ideal flavor matrices implicit in the displayed matrix are

$$
 T_{\rho^0}=\frac{\operatorname{diag}(1,-1,0)}{\sqrt2},\qquad
 T_\omega=\frac{\operatorname{diag}(1,1,0)}{\sqrt2},\qquad
 T_\phi=\operatorname{diag}(0,0,1),
$$

with off-diagonal entries normalized by the displayed factors of $\sqrt2$ (Bai Eq. (2.8), lines 258--272). The Goldstone convention is $U=\exp[(\sqrt2 i/f_\pi)\Phi]$ with $\operatorname{Tr}(t^at^b)=2\delta^{ab}$ (Bai Eq. (2.12), lines 290--298).

The complete anomalous functional is Bai Eq. (3.1), source lines 401--416:

$$
 \Gamma_{\rm WZW}^{\rm full}=\Gamma_0(U)+\mathcal C\int\operatorname{Tr}\{\text{all gauged }\mathcal A_{L,R},\alpha=dUU^\dagger,\beta=U^\dagger dU\text{ terms}\}+\Gamma_c,
 \qquad \mathcal C=\frac{N_c}{48\pi^2},
$$

where the full displayed polynomial is in Bai Eq. (3.1), and the Bardeen counterterm is Bai Eq. (2.5), source lines 233--243. The counterterm is part of the definition, not an optional contact term.

For a spin-one pair Bai prints the exact vector-only Lorentz form (their Eq. (4.4), source lines 682--689):

$$
 \mathcal L_{a\mathcal A_1\mathcal A_2}
 =S^{\rm vx}[\mathcal A_1\mathcal A_2]C_{\mathcal A_1\mathcal A_2}
 \mathcal F(m_a)\,\epsilon_{\mu\nu\rho\sigma}a\,\partial^\mu\mathcal A_1^\nu\partial^\rho\mathcal A_2^\sigma,
$$

with $S^{\rm vx}=2$ for identical fields and $1$ otherwise. In the vector-only sector the ALP enters through $\mathbb B_A\propto(k_L-k_R)da/f_a$; hence $aVV$ and $aV\gamma$ depend on $c_q^L-c_q^R$ after the full pole-plus-contact sum. Bai states this at source lines 607--608. The associated width is their Eq. (4.5), source lines 684--689.

The explicitly worked vector example is $a\to\omega\gamma$ (Bai Eqs. (3.20)--(3.22), source lines 610--651):

$$
 c_{\omega\gamma}^{\rm eff}=c_{\rm WZW}+c_{\rm ano}+\sum_{p=\pi^0,\eta,\eta'}c_p\theta_{pa},
$$

$$
 c_{\rm WZW}=\frac{egN_c[c_d^LQ_d-c_d^RQ_d+c_u^LQ_u-c_u^RQ_u-2Q_dc_{gg}\kappa_d-2Q_uc_{gg}\kappa_u]}{16\pi^2f_a},
 \quad
 c_{\rm ano}=-\frac{egN_c(Q_d\kappa_dc_{gg}+Q_u\kappa_uc_{gg})}{8\pi^2f_a},
$$

and $c_p$ are the pseudoscalar-pole coefficients printed in Bai Eq. (3.21). The $\kappa_q$ terms cancel in the displayed effective result, which is given in Bai Eq. (3.22). The same decomposition—direct WZW, anomalous/rotation term, and pseudoscalar pole—is the only exact pole-cancellation formula established in the source for a vector final state.

**Evidence limitation.** Bai does not print a separate component-level $P\to VV$ formula, nor separate $a\to\rho\rho$, $a\to\rho\omega$, or $a\to\phi\phi$ coefficients. The exact source-level statement for those channels is the functional Eq. (3.1), the vector-only background substitution Eq. (2.8), and the generic $aVV$ Lorentz form above. Any claimed component coefficient requires an explicit expansion of Eq. (3.1), and is not literature evidence until that expansion is recorded.

## Rotation cancellation and homogeneous ambiguity

Bai's field redefinition is

$$
 U\to U_LUU_R^\dagger,
 \qquad \mathbb B_{L/R}\to\mathbb B_{L/R}+c_{gg}\boldsymbol\theta_{L/R}\frac{da}{f_a},
 \qquad \boldsymbol\theta_{L/R}=\boldsymbol\delta_q\mp\boldsymbol\kappa_q,
$$

(Bai Eq. (3.10), source lines 437--440). The transformed functional contains

$$
 \Gamma_{\rm WZW}+\Gamma_c\to\Gamma_{\rm WZW}+\Gamma_c-
 \delta[\Gamma_{\rm WZW}+\Gamma_c](c_{gg}\boldsymbol\theta_L,c_{gg}\boldsymbol\theta_R),
$$

(Bai Eq. (3.11), source lines 442--458). Thus cancellation requires the full WZW term, counterterm, anomalous variation, and pseudoscalar pole/mixing contribution together. Bai explicitly says the vector-only WZW/HLS results agree, but also notes that a generalized HLS description has fourteen homogeneous counterterms (source line 460). Their coefficients are not fixed by anomaly matching; omitting or choosing them is a model/basis assumption. The printed Bai calculation therefore does not establish a unique homogeneous contribution to $P VV$ or $aVV$.

## Can Bai one-forms be identified with repository $W^A$?

Not as a literal field identity. Bai's $\mathbb B_V$ is a flavor-matrix connection with coupling $g$ absorbed into the one-form and with ideal quark-flavor assignments. The repository's $W^A_\mu$ are canonical massive eigenfields obtained from antisymmetric tensors after photon diagonalization. For each neutral tensor direction, the exact quadratic conversion already recorded in `memory/rcht-photon-vector-diagonalization.md` is

$$
 V^i_{\mu\nu}=\frac{P^i_{\mu\nu}+g_iF_{\mu\nu}}{M_i^2},
 \qquad P^i_\mu=M_iR^i_\mu,
 \qquad
 R_\mu^i=(K_V^{-1/2})_{ij}O_{jA}W^A_\mu.
$$

Therefore a consistent interface is to insert the *same flavor generators* into Bai's background functional and replace the one-form coefficients by the matched linear combination of $W^A_\mu$, including the $K_V^{-1/2}O$ rotation and all induced local terms. Componentwise $\mathbb B_V=gT_iV^i$ may equal a chosen $g_{iA}T_iW^A$ only after this matching; generally it does not equal $gT_AW^A$ when masses are nondegenerate or photon mixing is retained. Also, Bai's connection transformation law and the RChT CCWZ matter-field law are different: the identification is a phenomenological source matching, not a proof that antisymmetric RChT vectors are gauge bosons. Dropping the local terms generated by the tensor-to-Proca redefinition breaks representation equivalence (Ecker 1989, and the repository derivation).

## Ruiz--Femenia VVP basis and the $aVV$ insertion

Ruiz--Femenia et al. define

$$
 \mathcal L_{VVP}=\sum_{i=1}^4d_i\mathcal O_{VVP}^i,
$$

with (their Eq. (10), arXiv source HTML, Sec. 2):

$$
 \begin{aligned}
 \mathcal O^1_{VVP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\nu},V^{\rho\alpha}\}\nabla_\alpha u^\sigma\rangle,\\
 \mathcal O^2_{VVP}&=i\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\nu},V^{\rho\sigma}\}\chi_-\rangle,\\
 \mathcal O^3_{VVP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla_\alpha V^{\mu\nu},V^{\rho\alpha}\}u^\sigma\rangle,\\
 \mathcal O^4_{VVP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla^\sigma V^{\mu\nu},V^{\rho\alpha}\}u_\alpha\rangle.
 \end{aligned}
$$

With Alda's ALP inserted only through the axial derivative source and with Goldstones set to their vacuum for a direct $aVV$ vertex, $u_\mu$ is linear in $\partial_\mu a/f_a$, whereas $\chi_-$ has no such derivative-source insertion. Consequently the direct vertex receives $d_1$, $d_3$, and $d_4$ structures; $d_2$ contributes only if a pseudoscalar/scalar source or an ALP-induced pseudoscalar insertion is also retained. This is a source insertion analysis, not a published Alda matching result. The commonly quoted combinations $d_1+8d_2-d_3$ and $d_3$ belong to Ruiz--Femenia's $\omega\to\pi\gamma$ and $\omega\to3\pi$ amplitudes (Eqs. (27), (30)); they should not be transplanted to direct $aVV$ without expanding the Alda source and specifying external states.

**Status.** Literature evidence: all Bai equations and Ruiz--Femenia operator definitions cited above. Derived interface: the $W^A$ matching paragraph and the $d_1,d_3,d_4$ source-insertion classification. Assumptions: no axial-vector background, direct vacuum $aVV$ insertion, and no extra homogeneous odd-parity terms. Unresolved: explicit Bai component-level $P\to VV$ coefficients and the numerical matching of $g$ to the RChT odd-parity couplings.
