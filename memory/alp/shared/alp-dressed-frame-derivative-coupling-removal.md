---
title: Removal of the ALP derivative couplings in the LO meson theory and on-shell equivalence at finite $f_\pi/f_a$
status: derived
created: 2026-11-21
updated: 2026-11-21
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Removal of the ALP derivative couplings in the LO meson theory and on-shell equivalence at finite $f_\pi/f_a$

## Statement

Because the ALP enters the quark currents only through the single scalar $x=a/f_a$, the external connections $r_\mu=(\partial_\mu x)C_R(x)$ and $\ell_\mu=(\partial_\mu x)C_L(x)$ are pure gauge. The ALP-dependent chiral rotation by the transport matrices $W_{R,L}(x)$ removes them completely. In the resulting frame (called the dressed frame below) the LO Lagrangian has

- no ALP-quark current sources, $r_\mu=\ell_\mu=0$, and vanishing field strengths $F^{R,L}_{\mu\nu}=0$;
- the mass source $\chi=2B_0\widehat M(x)$ with $\widehat M=W_R^\dagger MW_L$;
- the topological coefficient $\widehat\Theta(x)=\Theta(x)+\int_0^x\operatorname{tr}C_A$;
- canonical, unmixed kinetic terms for $(P_i,a)$ ($K=\mathbf 1$); ALP-meson mixing occurs only through $M^2$.

The on-shell $2\to2$ tree amplitudes among the physical poles are the same in the original (derivative) frame and the dressed frame, exactly in $f_\pi/f_a$. Both statements were tested and hold.

## Conventions and units

Natural units; masses in GeV, $s,t$ and squared masses in GeV$^2$. $U=\exp(i\Phi/f_\pi)$, $\Phi=P_i\lambda_i$, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$ ([Appendix-A extract](../../extracts/Alda:2026nrd/appendix-a-expansion.md)). $C_R=C_V+C_A$, $C_L=C_V-C_A$; $\partial_xW_\eta=iC_\eta W_\eta$, $W_\eta(0)=\mathbf 1$ ([general-redefinition record](alp-general-chiral-redefinition-invariants.md)). Vertices are derivatives of the LO Lagrangian with all momenta incoming and the overall factor $i$ omitted ([joint-vertex record](alp-lo-joint-vertex-expansion.md)). Conventions of Alda:2026nrd [@Alda:2026nrd] throughout.

## Inputs

- Connection law $C'_\eta=g_\eta C_\eta g_\eta^\dagger+ig_\eta\partial_xg_\eta^\dagger$, transport matrices, $\widehat M$ and $\partial_\mu\widehat U=W_R^\dagger(D_\mu U)W_L$: [general-redefinition record](alp-general-chiral-redefinition-invariants.md).
- $\widehat\Theta$ and $X_{TC}=2i(\operatorname{tr}\log\hat u+i\widehat\Theta)$ in the frame without connections: [dressed topological record](alp-dressed-topological-invariant.md). Replacing Alda:2026nrd's linear $c_Ga/f_a$ by a general $\Theta(x)$ is an extension (assumption stated there); here $\Theta(x)=c_Gx$, so only Alda:2026nrd's form is needed.
- $\widehat M_1=-i[c_q^{VS},m_q]-i\{c_q^{AP},m_q\}$: [first-order record](alp-dressed-mass-first-order-invariants.md).
- CCWZ vertex generator and quadratic blocks: [joint-vertex record](alp-lo-joint-vertex-expansion.md). Pole-projected amplitudes with exact mixed propagator $D(s)^{-1}$: [pole-projected record](alp-lo-pole-projected-tree-amplitudes.md). Quadratic blocks of the original frame: [finite-mixing record](finite-alp-meson-mixing.md).

## Method

1. **Connection removal (H1 i).** Choose $g_\eta=W_\eta^\dagger$. Then $C'_\eta=W_\eta^\dagger C_\eta W_\eta+iW_\eta^\dagger\partial_xW_\eta=W_\eta^\dagger C_\eta W_\eta-W_\eta^\dagger C_\eta W_\eta=0$. The script builds $W_\eta$ from the Taylor recursion $W_{n+1}=i\sum_r\binom nrC_rW_{n-r}$ and evaluates $C'$ through $O(x^3)$ for generic $x$-dependent, non-commuting Hermitian $C(x)=\sum_{n\le3}C_nx^n/n!$: symbolic $2\times2$ (16 real parameters) and random rational complex $3\times3$ matrices.
2. **Field strength (H1 ii).** With $x=x(y)$ a function of two coordinates, $F_{01}=\partial_0r_1-\partial_1r_0-i[r_0,r_1]$ is evaluated symbolically: the second-derivative and $\partial_0x\,\partial_1x\,C'(x)$ terms cancel by antisymmetry and $[C(x),C(x)]=0$.
3. **Dressed sources.** $\widehat M_n=\sum_{r+s+t=n}\frac{n!}{r!s!t!}W^\dagger_{R,r}M_tW_{L,s}$ and $\widehat\Theta_n=\Theta_n+\operatorname{tr}C_{A,n-1}$. For the linear spurion of Alda:2026nrd, $M(x)=m_q+ix\{c_q^P,m_q\}$, constant $c_q^A$, $c_q^V=c_q^S=0$, $\Theta=c_Gx$:
$$\widehat M(x)=e^{-ixc_q^A}\big(m_q+ix\{c_q^P,m_q\}\big)e^{-ixc_q^A},\qquad \widehat\Theta(x)=(c_G+\operatorname{tr}c_q^A)\,x=c_{GA}\,x .$$
The script checks the recursion against this exponential form for a non-diagonal $c_q^A$, and against the diagonal closed form $\widehat M_n=m_q(-2i)^n(c^A_q)^{n-1}(c^A_q-n\,c^P_q)$.
4. **Frame comparison (H2).** For each frame, `quadratic_blocks` gives $(K,M^2)$, the poles and $K$-orthonormal pole vectors are found, and the $2\to2$ amplitude (four-point contact plus $s$-, $t$-, $u$-channel exchanges with the full $D(q^2)^{-1}$) is evaluated with `projected_amplitude`. The original frame uses Taylor arrays $\chi=(2B_0m_q,\,2iB_0\{c_q^P,m_q\},0,\dots)$, $C_R=-C_L=c_q^A$, $\Theta_1=c_G$; the dressed frame uses $2B_0\widehat M_n$ ($n\le4$), $C_R=C_L=0$, $\widehat\Theta_1=c_{GA}$. Through four external legs these Taylor orders are complete, so no expansion in $f_\pi/f_a$ is made. The external states are identified as the ALP-like pole (largest $|v_a|$) and the pion-like pole (largest $|v_{P_3}|$ among the others).

**Numerical precision (stated, not an approximation of the physics).** The source arrays, all basis-field vertices and $s,t$ are exact; the pole masses are roots of quartic or quintic polynomials and are computed, with the pole vectors and the on-shell Gram matrix, with 50 significant digits. Equality is tested to relative precision $10^{-40}$.

**Illustrative inputs (GeV units; not a phenomenological choice).** $m_q=\operatorname{diag}(1/500,1/200,1/10)$, $B_0=27/10$, $f_\pi=23/250$, $m_0^2=4/5$, $m_{a,0}^2=4$, $c_q^A=\operatorname{diag}(1,1/2,-1/3)$, $c_q^P=\operatorname{diag}(1/5,-1/4,1/2)$, $c_G=1/3$; $f_\pi/f_a\in\{1/10,1,3\}$. A second case adds the real symmetric $d$-$s$ entry $(c_q^A)_{ds}=(c_q^A)_{sd}=1/4$ (CP conserving), tested at $f_\pi/f_a\in\{1,3\}$. Kinematic points $(s,t)=(3,-3/5)$ and $(7,-7/4)$ GeV$^2$, with $u=\sum_im_i^2-s-t$, i.e. two Gram matrices per channel and ratio, away from all exchanged poles. With $m_{a,0}^2=4$ GeV$^2$ the full Hessian is positive for all ratios: e.g. diagonal case poles $\{0.0187,0.280,1.077,4.019\}$, $\{0.0187,0.186,1.044,5.850\}$, $\{0.0176,0.0261,1.012,19.80\}$ GeV$^2$ for $f_\pi/f_a=1/10,1,3$ (with $m_{a,0}^2=1/2$ GeV$^2$ the vacuum is unstable at $f_\pi/f_a=3$). Channels: $(A,\pi)\to(A,\pi)$, $(\pi,\pi)\to(A,A)$, $(A,A)\to(A,A)$, $A$ the ALP-like pole.

**Basis in the $d$-$s$ case.** With real symmetric $c_q^{A,P}$ and real diagonal $m_q$, both frames are invariant under $U\to U^{\mathsf T}$ with $a$ fixed ($c_q^A$ symmetric gives $D_\mu(U^{\mathsf T})=(D_\mu U)^{\mathsf T}$, and $\widehat M(x)$ is a symmetric matrix). $P_7$ is odd under this symmetry and $P_0,P_3,P_8,P_6,a$ are even, so $P_7$ neither mixes with nor is exchanged between the states used; the basis $(P_3,P_8,P_0,P_6,a)$ is complete for these amplitudes. Charged mesons do not couple to these states since the sources commute with the quark charge matrix.

## Result

- H1: $C'=0$ through $O(x^3)$ and $F_{01}=0$ for generic non-commuting $C(x)$; the dressed kinetic matrix is exactly $K=\mathbf 1$ (the original-frame $K$ has $K_{ia}\ne0$). Consequently all ALP couplings in the dressed frame come from $\frac{f_\pi^2}{4}\langle\chi_+\rangle$ with $\chi=2B_0\widehat M(x)$ and from the topological term with $\widehat\Theta(x)$. In particular, apart from the free ALP kinetic term, no vertex containing an ALP depends on momenta (checked for three- and four-leg vertices).
- H2: for all ratios, both source choices, both kinematic points and the three channels, the pole masses, the ALP components of the pole vectors (the ALP coordinate is not redefined) and the on-shell amplitudes coincide between the two frames to better than $10^{-40}$ relative precision. Example ($f_\pi/f_a=1$, diagonal sources, $s=3$, $t=-3/5$ GeV$^2$): $\mathcal M(A\pi\to A\pi)=5.6177855852551$, $\mathcal M(AA\to AA)=53.516771032960$ in both frames.
- Observation: in the diagonal case the dressed-frame amplitudes do not depend on $s,t$: the neutral diagonal generators commute, so $\langle u_\mu u^\mu\rangle$ has no four-meson vertex among them, and all cubic vertices vanish with CP-conserving sources (see the [partial-wave record](alp-lo-partial-wave-unitarity.md)). The $d$-$s$ case has genuine kinematic dependence through $[\lambda_3,\lambda_6]\ne0$.

## Verification

Unit tests in `scripts/test_derive_alp_dressed_frame_unitarity.py`, classes `H1DressedFrame` and `H2OnShellEquivalence`, plus a test that the cached multilinear vertex contraction equals `ccwz_vertex` and a test of the $P_7$ decoupling (zero quadratic mixing and zero three- and four-leg vertices with one $P_7$, in both frames). All pass. The existing suites `test_derive_alp_lo_joint_vertices.py` and `test_derive_alp_pole_projected_amplitudes.py` also pass. Limitations: LO only, tree level, stable real poles; photons and weak sources are not included; the agreement is a numerical test at the stated points, while the general statement follows from the invertible field redefinition $\widehat U=W_R^\dagger UW_L$.

## Invocation

```sh
.venv/bin/python scripts/derive_alp_dressed_frame_unitarity.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_dressed_frame_unitarity.py
```

Script: [`derive_alp_dressed_frame_unitarity.py`](../../../scripts/derive_alp_dressed_frame_unitarity.py).
