---
title: Double-trace photon source term in the vector mass spectrum and the KK pi current
status: derived
created: 2026-11-17
updated: 2026-11-17
bibtex_keys: [Alda:2026nrd, Ecker:1988te]
provenance: analytical-derivation
---

# Double-trace photon source term: vector masses and the $KK\pi$ current

This record settles task 4.1 of the change `extend-joint-vector-model-diagnostics`, which is item 3 of the [candidate causes](joint-vector-fit-inadequacy-candidate-causes.md). The term had been deferred in the [likelihood variation options](joint-vector-likelihood-variation-options.md) on 2026-11-07.

**Verdict.** The term does not contribute to the vector mass spectrum at the order used in the joint fit, so it cannot move the $\omega,\phi$ masses. It does contribute to the $KK\pi$ current, but only through the photon couplings of $\omega$ and $\phi$. The extension is therefore **applicable to the $KK\pi$ current and to the $\omega,\phi$ photon couplings, and not applicable to the mass spectrum**.

## Operator

The term written $y\,\mathrm{tr}(\chi_+Q)\mathbb 1$ in the change proposal is the contribution of the operator

$$
\mathcal L_y=\frac{y_V}{2\sqrt2}\langle V_{\mu\nu}\rangle\langle\chi_+f_+^{\mu\nu}\rangle,
\qquad [y_V]=\mathrm{GeV}^{-1},
$$

to the photon source matrix. With $f_+^{\mu\nu}=-2eQF^{\mu\nu}$ in the Goldstone vacuum, it gives

$$
S_\gamma=F_VQ+z_V\{\chi_+,Q\}+y_V\langle\chi_+Q\rangle\mathbf 1 .
$$

This is the definition in the [photon–vector diagonalization record](../foundations/rcht-photon-vector-diagonalization.md) and the [leptonic-width record](../foundations/vector-leptonic-width.md). The same term is implemented as `y` in `source_residues` of [`scripts/implementation/vector_sector.py`](../../../scripts/implementation/vector_sector.py), and as the `double_trace` hypothesis of `action_sector` in [`scripts/derive_odd_low_energy_jacobians.py`](../../../scripts/derive_odd_low_energy_jacobians.py). It has two flavor traces: the vector field appears only through its trace $\langle V\rangle$, and the source appears only through $\langle\chi_+Q\rangle$.

**Assumption (identification of the operator):** the proposal's $y$ is taken to be this $y_V$. A different double trace, $\frac{\kappa_V}{4}\langle V_{\mu\nu}\rangle\langle V^{\mu\nu}\chi_+\rangle$, does change the $\omega,\phi$ masses ([mass extension](../foundations/vector-mass-extension.md)). That term is already a coordinate of the joint fit, together with $\Delta_V$, so the four ground-state masses are already described by four mass parameters. The user confirmed this identification: $y$ is $y_V$.

## Assumptions

- Tree level, one vector nonet, antisymmetric-tensor fields, canonical kinetic terms ([RChT basics](../foundations/rcht-basics.md)) [@Ecker:1988te].
- Mass terms are counted at order $e^0$, as in the [nonideal mass record](../foundations/vector-nonideal-mass-action.md) and in the joint fit. Electromagnetic mass shifts are omitted for every source term ($F_V$, $z_V$ and $y_V$). This is an approximation stated explicitly here, not a result.
- The vacuum spurion is real and diagonal, $\chi_+=\mathrm{diag}(x_u,x_d,x_s)$. The $KK\pi$ substitution below uses exact strong isospin, $x_u=x_d=x_l$, as in the fit.
- The term acts on the ground-state nonet only. The excited-vector poles of the BaBar current have their own fixed couplings and are not changed by $y_V$.

## Inputs

- Chiral building blocks and Goldstone matrix: Alda:2026nrd conventions, $\Phi=\sum_P\lambda_PP$ with $\mathrm{tr}(\lambda_i\lambda_j)=2\delta_{ij}$ including the singlet ($\eta'$) direction, $u=\exp(i\Phi/(2f))$, $\chi=2B_0\,\mathrm{diag}(m_u,m_d,m_s)$ ([Alda:2026nrd generators extract](../../extracts/Alda:2026nrd/neutral-meson-generators.md)) [@Alda:2026nrd]. Definitions of $\chi_\pm$ and $f_\pm$: [RChT basics](../foundations/rcht-basics.md).
- Photon sources $\ell_\mu=r_\mu=-eQA_\mu$, the source matrix $S_\gamma$, the residues $R_i=\langle T_iS_\gamma\rangle$ and the $O(e^2)$ pole equation $\det(sK_V-D)=0$: [photon–vector diagonalization](../foundations/rcht-photon-vector-diagonalization.md).
- Ideal generators $T_\rho=\mathrm{diag}(1,-1,0)/\sqrt2$, $T_{\omega_I}=\mathrm{diag}(1,1,0)/\sqrt2$, $T_{\phi_I}=\mathrm{diag}(0,0,1)$ and the definitions $F_l=F_V+2z_Vx_l$, $F_s=F_V+2z_Vx_s$: [nonideal mass record](../foundations/vector-nonideal-mass-action.md), [leptonic-width record](../foundations/vector-leptonic-width.md).
- Structure of the photon insertions in the BaBar $KK\pi$ current (`source_rho`, `source_light`, `source_strange`): [`scripts/fit_babar_kkpi_global_vector_sector.py`](../../../scripts/fit_babar_kkpi_global_vector_sector.py), described in the [joint-likelihood baseline](babar-kkpi-global-vector-sector-fit-baseline.md).

## Method (step by step)

The script [`scripts/derive_double_trace_vector_mass_and_kkpi.py`](../../../scripts/derive_double_trace_vector_mass_and_kkpi.py) performs each step with SymPy.

1. **Goldstone expansion of the source trace.** For a general Hermitian $\Phi$ (nine real fields), expand $u$ to third order and build $\chi_+=u^\dagger\chi u^\dagger+u\chi u$ and $f_+=u F_Lu^\dagger+u^\dagger F_Ru$ with $F_L=F_R=-eQF$. The trace is
   $$
   \langle\chi_+f_+\rangle=-4eF\langle\chi Q\rangle+O(\Phi^2),
   $$
   that is $-2eF\langle\chi_+Q\rangle$ at $\Phi=0$. The $\Phi^1$ and $\Phi^3$ coefficients vanish, and the $\Phi^2$ coefficient does not. The reason is that $\Phi\to-\Phi$ exchanges $u$ and $u^\dagger$, which leaves $\chi_+$ and $f_+$ unchanged when $\chi$ is real and $F_L=F_R$.
2. **No vector mass term.** Insert a general nonet $V$ (nine components, one Lorentz structure). $\mathcal L_y$ has degree one in the vector components, so its second derivative with respect to them vanishes. Every term is proportional to $eF^{\mu\nu}$, so $\mathcal L_y=0$ when the photon is switched off. At order $e^0$ the quadratic vector terms, and therefore the $\rho,\omega,\phi,K^*$ masses and the $\omega$–$\phi$ angle, are unchanged.
3. **Which vectors couple.** The derivative of $\mathcal L_y$ vanishes for all off-diagonal components ($\rho^\pm$, $K^{*}$) and for the traceless neutral directions ($\rho^0$, $\omega_8$). Only the trace $\langle V\rangle$ couples. In the physical basis this is $\langle T_\omega\rangle=\sqrt2\cos\theta_V-\sin\theta_V$ and $\langle T_\phi\rangle=\sqrt2\sin\theta_V+\cos\theta_V$.
4. **Photon residues.** With $\chi_+=\mathrm{diag}(x_l,x_l,x_s)$,
   $$
   R_\rho=\frac{F_l}{\sqrt2},\qquad
   R_{\omega_I}=\frac{F_l+2y_V(x_l-x_s)}{3\sqrt2},\qquad
   R_{\phi_I}=-\frac{F_s-y_V(x_l-x_s)}{3}.
   $$
   The shifts $\delta R_{\omega_I}=\sqrt2y_V(x_l-x_s)/3$ and $\delta R_{\phi_I}=y_V(x_l-x_s)/3$ agree with the [leptonic-width record](../foundations/vector-leptonic-width.md).
5. **$KK\pi$ current.** At tree level, the photon enters the $e^+e^-\to KK\pi$ current either directly (WZW and local terms) or through a photon–vector transition followed by the vector's decay into $KK\pi$. $\mathcal L_y$ always contains one vector, so it has no direct photon–$KK\pi$ term. Its $\Phi^0$ part is a photon–$\omega/\phi$ transition. Its $\Phi^2$ part is a $\gamma VPP$ vertex. That vertex cannot contribute to a three-pseudoscalar final state at tree level, because the vector would then have to turn into a single pseudoscalar. The only Lorentz structure for that is $\langle V_{\mu\nu}\rangle p^\mu p^\nu$, which vanishes by antisymmetry (checked in the script). There are no odd-$\Phi$ terms. Hence $y_V$ enters the current only through $R_{\omega_I},R_{\phi_I}$. The script multiplies these residues by the ideal-basis neutral propagator matrix and recovers the fit code's `source_light` and `source_strange` combinations exactly, with the replacements
   $$
   F_l\to F_l+2y_V(x_l-x_s)\ \text{in the }\omega_I\text{ insertion},\qquad
   F_s\to F_s-y_V(x_l-x_s)\ \text{in the }\phi_I\text{ insertion},
   $$
   while the $\rho^0$ insertion keeps $F_l$. In the current code, `source_rho` and `source_light` share the factor `light`. Task 4.2 must therefore separate the isovector factor from the isoscalar one.
6. **Order $e^2$ (for completeness).** For a vector of mass $M_i$ mixing with the photon through $g_i=\sqrt2eR_i$, the pole equation $\det(sK_V-D)=0$ gives $m_i^2=M_i^2+2e^2R_i^2+O(e^4)$; the script checks this for two vectors. Through this shift $y_V$ changes the $\omega,\phi$ masses by $2e^2(R_i^2|_{y_V}-R_i^2|_{y_V=0})$. The joint fit omits this shift for every source term, so including it for $y_V$ alone would be inconsistent. This record does not include it.
7. **Weak charged current.** For an off-diagonal current generator, $\langle\chi_+T\rangle=0$ in the diagonal vacuum. $\mathcal L_y$ therefore gives nothing to $\tau\to KK\pi\nu$-type charged currents.

## Result

- **(a) Vector mass spectrum: not applicable.** $\mathcal L_y$ is linear in $V$ and proportional to $e$, so it gives no order-$e^0$ quadratic vector term. The $\rho,\omega,\phi,K^*$ masses and the $\omega$–$\phi$ angle are unchanged. The expectation in candidate cause 3, that this term moves the $\omega,\phi$ masses independently, does not hold at the order of the fit. Independent $\omega,\phi$ mass freedom comes from $\Delta_V$ and $\kappa_V$, which are already fitted.
- **(b) $KK\pi$ current: applicable.** It enters only through the photon residues of the isoscalar states: $\delta R_{\omega_I}=\sqrt2y_V(x_l-x_s)/3$ and $\delta R_{\phi_I}=y_V(x_l-x_s)/3$. Equivalently, $F_l\to F_l+2y_V(x_l-x_s)$ in the $\omega_I$ insertion and $F_s\to F_s-y_V(x_l-x_s)$ in the $\phi_I$ insertion. The $\rho^0$, $\rho^\pm$ and $K^*$ couplings are unchanged.
- **Other blocks of the joint likelihood.** The same residue shifts enter the $\omega,\phi$ leptonic widths and the photon-insertion parts of the odd radiative rows. Both blocks already take the residues from `source_residues`, which accepts `y`. Together with $F_V$ and $z_V$, $y_V$ makes $(R_\rho,R_{\omega_I},R_{\phi_I})$ independent whenever $x_l\ne x_s$ (determinant $-2(x_l-x_s)^2/9$, [leptonic-width record](../foundations/vector-leptonic-width.md)).

## Verification

```text
.venv/bin/python scripts/derive_double_trace_vector_mass_and_kkpi.py
.venv/bin/python -m unittest scripts.test_derive_double_trace_vector_mass_and_kkpi
```

The script asserts each step above and prints "All assertions passed." The test module [`scripts/test_derive_double_trace_vector_mass_and_kkpi.py`](../../../scripts/test_derive_double_trace_vector_mass_and_kkpi.py) runs 9 tests, all passing. Besides the symbolic assertions, it includes two independent numerical checks. The symbolic residues agree with `vector_sector.source_residues` at a random point. The source trace computed with an exact matrix exponential (`scipy.linalg.expm`) is even in $\Phi$ to all orders.

## Open points

- Task 4.2 must add a $y_V$ coordinate to the joint likelihood and separate the isovector and isoscalar light factors in the BaBar current. The user approved the notation $y_V$ as the same coefficient as the proposal's $y$, with the descriptive name "double-trace singlet-source structure" (the coordinate is the coefficient of the double-trace singlet source $\langle\chi_+Q\rangle\mathbf 1$).
- The $O(e^2)$ electromagnetic mass shifts remain omitted for all source terms.
