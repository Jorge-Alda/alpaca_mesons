---
title: Resonance saturation of NLO constants and explicit pseudoscalar-resonance mixing
status: mixed
created: 2026-11-02
updated: 2026-11-05
bibtex_keys: [Ecker:1988te, Ecker:1989yg, Kaiser:1998ds, Kaiser:2000gs, Alda:2026nrd]
provenance: analytical-derivation
---

# Resonance contributions to the NLO constants and explicit resonance mixing

## Scope and conventions

Use the $U(3)$ operator basis in [lagrangian.md](lagrangian.md), the neutral-field conventions in [neutral-meson-matrices.md](neutral-meson-matrices.md), and the quadratic Hessians in [quadratic-lagrangian.md](quadratic-lagrangian.md). In particular, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$ and $\Phi=\lambda_iP_i$. The tree-level matching discussed below follows Ecker:1988te's $SU(3)$ resonance theory. This is enough to identify the established resonance contributions to the ordinary $L_4$--$L_8$ structures, but it does not uniquely saturate the extra $U(3)$ singlet constants $\Lambda_1,\Lambda_2,L_{18},L_{25}$. The latter need an explicitly specified singlet/topological resonance Lagrangian and matching prescription. The $SU(3)$ formulas must therefore not be presented as predictions for all constants in this project's $U(3)$ basis [@Ecker:1988te; @Kaiser:1998ds; @Kaiser:2000gs].

## Which resonance contributes to which standard NLO constant?

For a resonance with mass $M_R$, quadratic operator $\mathcal D_R$, and linear Goldstone current $J_R$, tree-level elimination gives

$$
\Delta\mathcal L_{\rm loc}=-\frac12J_R\mathcal D_R^{-1}J_R.
$$

Expanding $\mathcal D_R^{-1}$ in derivatives and spurions gives local $p^4$ operators. The following are the leading single-resonance contributions in the usual $SU(3)$ basis, in terms of the standard resonance couplings $c_d,c_m$ (scalar octet), $\widetilde c_d,\widetilde c_m$ (scalar singlet), $d_m$ (pseudoscalar octet), $G_V,F_V$ (vector octet), and $F_A$ (axial-vector octet). A constant not listed has zero contribution from the stated exchange at this order.

| NLO constant | Resonance contribution at tree level | Origin |
|---|---|---|
| $L_4$ | $-c_dc_m/(3M_{S_8}^2)+\widetilde c_d\widetilde c_m/M_{S_1}^2$ | scalar octet and scalar singlet |
| $L_5$ | $c_dc_m/M_{S_8}^2$ | scalar octet |
| $L_6$ | $-c_m^2/(6M_{S_8}^2)+\widetilde c_m^2/(2M_{S_1}^2)$ | scalar octet and scalar singlet |
| $L_7$ | $+d_m^2/(6M_{P_8}^2)$ from the octet projector; a pseudoscalar singlet contributes $-\widetilde d_m^2/(2M_{P_1}^2)$ | pseudoscalar octet and singlet |
| $L_8$ | $c_m^2/(2M_{S_8}^2)-d_m^2/(2M_{P_8}^2)$ | scalar octet and pseudoscalar octet |
| $L_1$ | $G_V^2/(8M_V^2)-c_d^2/(6M_{S_8}^2)+\widetilde c_d^2/(2M_{S_1}^2)$ | vector octet and scalar octet/singlet |
| $L_2$ | $G_V^2/(4M_V^2)$ | vector octet |
| $L_3$ | $-3G_V^2/(4M_V^2)+c_d^2/(2M_{S_8}^2)$ | vector and scalar octets |
| $L_9$ | $F_VG_V/(2M_V^2)$ | vector octet |
| $L_{10}$ | $-F_V^2/(4M_V^2)+F_A^2/(4M_A^2)$ | vector and axial-vector octets |

The scalar singlet pieces cancel the corresponding octet pieces in $L_1,L_4,L_6$ in the degenerate nonet limit with the large-$N_c$ coupling relations. This is why the classic resonance estimate has $L_4,L_6$ suppressed, while $L_5,L_8$ are leading. Ecker:1988te finds that vector and axial-vector exchange feed $L_1,L_2,L_3,L_9,L_{10}$, scalar exchange feeds $L_1,L_3$--$L_6,L_8$ (with $L_5$ dominated by the scalar octet), and the light pseudoscalar singlet controls the important $L_7$ contribution. The detailed contribution depends on retained multiplets, representation, and short-distance assumptions [@Ecker:1988te; @Ecker:1989yg].

In this project's NLO quadratic matrices, $L_4,L_5$ enter the kinetic Hessian and $L_6,L_7,L_8$ enter the mass Hessian. Thus the scalar estimates for $L_4,L_5$ change kinetic normalization, while scalar and pseudoscalar estimates for $L_6,L_7,L_8$ change mass entries. The vector and axial-vector exchanges in the classic basis do not generate these neutral mass/kinetic constants at this order; their main effects are in four-derivative and external-source operators. Crucially, the standard $SU(3)$ saturation table does **not** fix $\Lambda_1,\Lambda_2,L_{18},L_{25}$. In the project's convention, $\Lambda_1$ contributes to the singlet kinetic term; $\Lambda_2$ and $L_{25}$ to singlet-related mass entries; $L_{18}$ to singlet-trace kinetic structures. Assigning these to a particular resonance requires extending the resonance fields and couplings to the $U(3)$ topological-source sector and matching those operators explicitly. Resonance saturation is an estimate of selected renormalized constants at a scale, not a symmetry theorem [@Kaiser:1998ds; @Kaiser:2000gs].

## Why a scalar resonance shifts the mass matrix without a bilinear mixing angle

Parity forbids a scalar resonance from mixing bilinearly with one pseudoscalar around the parity-even vacuum. Its effects on the pseudoscalar quadratic Lagrangian arise through scalar exchange between two-Goldstone vertices and through the induced scalar vacuum/tadpole shift. Integrating that scalar out generates the local $L_5,L_8$ (and singlet-dependent $L_4,L_6$) terms above. If the scalar is kept explicit, these same exchange contributions must be retained as propagator graphs or as the corresponding coupled-field elimination; they must not also be included in the local LECs.

The pseudoscalar resonance does have a bilinear mixing with Goldstones, because its linear coupling is to $\chi_-$, whose expansion contains one Goldstone field proportional to $\chi_0$. This gives a concrete dynamical derivation of the pseudoscalar-exchange contribution to the NLO masses.

## Explicit pseudoscalar resonance: coupled quadratic system

For an explicit neutral pseudoscalar resonance nonet $R=R_A\lambda_A/\sqrt2$, take canonical kinetic terms and, initially, a common mass $M_R$. Here $A=3,8$ describes the neutral members of the octet and $A=0$ the singlet; in a realistic matching the octet and singlet masses and couplings are independent. Use the conventional linear interaction $i d_m\langle R\chi_-\rangle$. With the field normalization of this repository, the linear part of $\chi_-$ is $-2i\{\Phi,\chi_0\}/(2f_\pi)=-i\{\Phi,\chi_0\}/f_\pi$. The mixed Hessian of this interaction defines the bilinear matrix below, with resonance index $A$ and Goldstone index $a$:

$$
B_{Aa}=\frac{d_m}{\sqrt2f_\pi}\operatorname{tr}\!\left(\lambda_A\{\lambda_a,\chi_0\}\right).
$$

The sign of $B$ depends on the sign convention for $P_R$ and has no effect on the second-order pole shifts. The full quadratic Lagrangian is

$$
\mathcal L^{(2)}=\frac12\partial P^T\partial P+\frac12\partial R^T\partial R
-\frac12P^TM_P^2P-\frac12R^TM_R^2R+R^TBP,
$$

so the inverse propagator at $s=p^2$ is

$$
\mathcal Q(s)=\begin{pmatrix}sI-M_P^2&B^T\\B&sI-M_R^2\end{pmatrix}.
$$

The exact pole equation is $\det\mathcal Q(s)=0$. Away from a resonance pole it is equivalently

$$
\det\!\left[sI-M_P^2-B^T(sI-M_R^2)^{-1}B\right]=0.
$$

This Schur complement retains the full resonance pole and momentum dependence. For a single nondegenerate Goldstone mass $m_P$ and resonance mass $M_P$, with small bilinear entry $B$, the two eigenvalues and the Goldstone-like eigenvector are

$$
m_{P,\rm phys}^2=m_P^2+\frac{B^2}{m_P^2-M_R^2}+O(B^4),\qquad
M_{R,\rm phys}^2=M_R^2+\frac{B^2}{M_R^2-m_P^2}+O(B^4),
$$

$$
\frac{R}{P}\bigg|_{s\simeq m_P^2}=\frac{B}{M_R^2-m_P^2}+O(B^3).
$$

For a multiplet, first diagonalize the Goldstone mass matrix, then use the same formulas between every nondegenerate Goldstone eigenstate and resonance eigenstate. Near degeneracy, perturbation theory fails and the corresponding $2\times2$ (or larger) block must be diagonalized exactly. If the Goldstone sector has NLO kinetic and mass matrices, canonicalize it first or use the generalized eigenproblem $M^2v=sKv$ within the full coupled system.

Expanding the resonance propagator for $|s|\ll M_R^2$ gives

$$
-B^T(sI-M_R^2)^{-1}B
=B^TM_R^{-2}B+s B^TM_R^{-4}B+O(s^2/M_R^6).
$$

Thus the Goldstone-only theory receives the induced local corrections

$$
\delta M_P^2=-B^TM_R^{-2}B,\qquad
\delta K_P=B^TM_R^{-4}B.
$$

These are the low-energy expansion of the **same** explicit exchange. In the pseudoscalar octet matching, the leading mass term matches the $-d_m^2/(2M_P^2)$ part of $L_8$ together with the octet-projector contribution $+d_m^2/(6M_P^2)$ to $L_7$. Since $B=O(p^2)$, the induced mass term is $O(p^4)$, while the displayed kinetic correction is $O(p^6)$. Once the resonance is explicit, retain $B$ and $M_R$ and omit these same induced pieces from the local LECs; any independent short-distance remainder remains in the local Lagrangian.

## NLO masses and mixings in the Goldstone block

Even before adding explicit resonance fields, the project's raw NLO neutral quadratic terms are $\tfrac12\partial P^T\delta K\partial P-\tfrac12P^T\delta M^2P$. Let $M_0^2$ denote the LO neutral mass matrix including the $U(1)_A$ singlet mass. To first NLO order, canonical normalization $P=(I-\tfrac12\delta K)P_c$ gives

$$
M_{c}^2=M_0^2+\Delta_c,
\qquad \Delta_c=\delta M^2-\frac12\{\delta K,M_0^2\}.
$$

Let $U_0^TM_0^2U_0=\operatorname{diag}(\mu_1^2,\mu_2^2,\mu_3^2)$, and set $\widehat\Delta=U_0^T\Delta_cU_0$. For nondegenerate LO eigenvalues,

$$
m_i^2=\mu_i^2+\widehat\Delta_{ii}+O(\mathrm{NLO}^2),\qquad
U=U_0(I+\epsilon),\qquad
\epsilon_{ij}=\frac{\widehat\Delta_{ij}}{\mu_j^2-\mu_i^2}\quad(i\ne j),\quad \epsilon_{ii}=0.
$$

The antisymmetric $\epsilon$ rotates the eigenvectors; for a two-state block the NLO angle shift is $\delta\theta=\widehat\Delta_{12}/(\mu_2^2-\mu_1^2)$. This gives the explicit NLO mass and mixing correction using the entries already listed in [neutral-meson-matrices.md](neutral-meson-matrices.md).

More generally, in the LO mass basis with residual local $\delta K,\delta M^2$ and explicit pseudoscalar resonances $A$ coupled by $B_{Ai}$, the first nonzero corrections are

$$
\delta m_i^2=\delta M^2_{ii}-m_i^2\delta K_{ii}-\sum_A\frac{B_{Ai}^2}{M_A^2-m_i^2},
\qquad
v_A^{(i)}=\frac{B_{Ai}}{M_A^2-m_i^2},
$$

$$
v_j^{(i)}=\frac{\delta M^2_{ji}-m_i^2\delta K_{ji}-\displaystyle\sum_A B_{Aj}B_{Ai}/(M_A^2-m_i^2)}{m_i^2-m_j^2},\qquad j\ne i.
$$

Here $v_A^{(i)}$ is the resonance component and $v_j^{(i)}$ the other-Goldstone component in the original field coordinates. Goldstone--resonance mixing starts at first order in $B$; the pole mass shift and induced mixing between Goldstones start at second order. Reversing the sign of a resonance field reverses the sign of its eigenvector component and $B$, but not the pole shifts. These expressions require separated poles and small bilinear mixing.

With the pseudoscalar resonance kept dynamical, replace the Goldstone inverse propagator by its full coupled block $\mathcal Q(s)$. The NLO local matrix $M_0^2+\delta M^2_{\rm residual}$ and kinetic matrix $I+\delta K_{\rm residual}$ remain in the Goldstone block, and the extra resonance shifts/mixes the poles by the formulas above. Equivalently, solve the Schur-complement pole equation and take its eigenvectors. In the low-energy expansion this reproduces the resonance contribution to $\delta M^2$ and higher-order $\delta K$; using both the explicit block and the same saturated LEC pieces would double count. The $\eta_8$--$\eta_0$ and isospin-breaking $\pi_3$ mixings are therefore retained, with additional mixing with the explicit pseudoscalar resonance whenever $B_{Aa}\ne0$.

For an ALP, append the ALP to the Goldstone-sector field vector and retain its mixed kinetic and mass entries as in [lsz-invariant-split.md](lsz-invariant-split.md). If the explicit resonance also has an ALP bilinear coupling, add the corresponding resonance row/column to $\mathcal Q(s)$; the LSZ pole eigenvector then contains meson and resonance components. The invariant source construction still applies to the full block system [@Alda:2026nrd].

## Inputs, limitations, and verification

The resonance matching in the table is the leading single-resonance tree approximation of Ecker:1988te, expressed with standard $SU(3)$ resonance couplings; the source's equations for vector, axial-vector, scalar, and pseudoscalar exchange give these assignments. The neutral NLO matrices and normalization are the repository derivations linked above. The explicit bilinear formulas follow by expanding $i d_m\langle P_R\chi_-\rangle$ and block-inverting the quadratic form. They apply for a parity-even vacuum and small mixing; exactly or nearly degenerate states require exact finite-dimensional diagonalization. Loop effects, energy-dependent widths and self-energies, flavor-dependent resonance masses, and OZI-violating $U(3)$ resonance couplings are not included in these tree formulas. The executable derivation is [`scripts/derive_dynamical_pseudoscalar_masses.py`](../../scripts/derive_dynamical_pseudoscalar_masses.py); run `.venv/bin/python scripts/derive_dynamical_pseudoscalar_masses.py` to print the matrices, poles, and mixing coefficients.

See also [the separation of explicit resonance exchange from the residual constants $L_i^r(\mu)$](counting/resonance-exchange-and-residual-constants.md) and the declared power counting in [counting/INDEX.md](counting/INDEX.md).
