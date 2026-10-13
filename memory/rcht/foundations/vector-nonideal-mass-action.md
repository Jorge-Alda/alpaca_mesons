---
title: Minimal tensor mass action for nonideal omega-phi mixing
status: derived
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Ecker:1988te, Ecker:1989yg, Moussallam:2021dpk]
provenance: analytical-derivation
---

# Minimal tensor mass action for nonideal omega-phi mixing

## Inputs and assumptions

The canonical antisymmetric-tensor kinetic term, generator normalization $\langle T_iT_j\rangle=\delta_{ij}$, positive-strange ideal $\phi$ phase and universal $F_V$ interaction are the inputs from [RChT basics](rcht-basics.md) [@Ecker:1988te]. Photon source signs and overlaps follow [quadratic photon-vector diagonalization](rcht-photon-vector-diagonalization.md). The [Moussallam extract](../../extracts/Moussallam:2021dpk/vector-mixing.md) motivates flavor and nonet breaking [@Moussallam:2021dpk]. This record derives an explicit *assumed tensor extension*, rather than asserting a complete translation of that paper's vector-field action or odd-parity interactions.

Assume tree level, exact strong isospin, a real vacuum spurion $\chi_+=\operatorname{diag}(x_l,x_l,x_s)$, and canonical flavor-independent kinetic terms. In the standard vacuum, $x_l=4B_0\hat m$, $x_s=4B_0m_s$; those values must eventually come from the selected Alda arXiv conventions/inputs, not a stale local manuscript. Natural units apply: $x_l,x_s,\Delta_V$ have mass dimension two and $e_m^V$ is dimensionless. Retain only one single-trace spurion mass operator and one singlet mass shift. Other double-trace spurion operators, kinetic breaking, loops, widths, higher spurion powers and electromagnetic mass shifts are omitted explicitly. This is a starting truncation, not the most general symmetry-breaking action.

## Action and executable derivation

Run:

```text
.venv/bin/python scripts/derive_vector_nonideal_mixing.py
```

The [script](../../../scripts/derive_vector_nonideal_mixing.py) constructs the explicit flavor matrices already defined in the input record, performs traces, differentiates the quadratic action, solves its determinant polynomial, rotates generators and sources, and computes the flavor response. It suppresses the common Lorentz contraction only when extracting mass coefficients.

Add to the canonical kinetic baseline:

$$
\mathcal L_{\mathrm{mass}}=
\frac{M_V^2}{4}\langle V_{\mu\nu}V^{\mu\nu}\rangle
-e_m^V\langle V_{\mu\nu}V^{\mu\nu}\chi_+\rangle
+\frac{\Delta_V}{12}\langle V_{\mu\nu}\rangle\langle V^{\mu\nu}\rangle.
$$

The last coefficient gives precisely a $\Delta_V$ mass-squared shift to the normalized singlet $T_0=\mathbf1/\sqrt3$. In the ideal basis $T_{\omega_I}=\operatorname{diag}(1,1,0)/\sqrt2$, $T_{\phi_I}=\operatorname{diag}(0,0,1)$, the script extracts

$$
\mathcal M_V^2=\begin{pmatrix} A&B\\B&D\end{pmatrix}
=\begin{pmatrix}
M_V^2-4e_m^Vx_l+2\Delta_V/3&\sqrt2\Delta_V/3\\
\sqrt2\Delta_V/3&M_V^2-4e_m^Vx_s+\Delta_V/3
\end{pmatrix}.
$$

The off-diagonal sign differs from Moussallam's displayed ideal-basis matrix because his ideal phi has negative strange phase; the positive-strange convention requires the corresponding phase change at all vertices.

The charged-generator calculation gives

$$
m_\rho^2=M_V^2-4e_m^Vx_l,\qquad
m_{K^*}^2=M_V^2-2e_m^V(x_l+x_s).
$$

Consequently $D=2m_{K^*}^2-m_\rho^2+\Delta_V/3$. These are conditional tree-level mass relations; measured masses do not guarantee exact satisfaction by the truncation.

## Physical generators and poles

Define $\mathbf V_I=R\mathbf V_{\rm phys}$ with

$$
R=\begin{pmatrix}\cos\theta_V&\sin\theta_V\\-\sin\theta_V&\cos\theta_V\end{pmatrix}.
$$

The derived off-diagonal entry is $B\cos2\theta_V+(A-D)\sin2\theta_V/2$. Substituting $t=\tan\theta_V$ in this expression, SymPy obtains $[(A-D)t+B(1-t^2)]/(1+t^2)=0$ and solves the numerator. The root connected to zero for $B\to0$ at $D>A$ is $t=[A-D+\sqrt{(D-A)^2+4B^2}]/(2B)$. The equivalent atan2 expression remains usable at $B=0$ when the eigenvalues are distinct. On the branch continuously connected to a lighter ideal nonstrange state when $D>A$,

$$
\theta_V=\frac12\operatorname{atan2}(2B,D-A),\qquad
m_{\omega,\phi}^2=\frac{A+D\mp\sqrt{(D-A)^2+4B^2}}2.
$$

The degenerate case $B=0,A=D$ has no unique angle; physical projectors or complete state sums must be used. Positive eigenvalues are required for a stable chosen parameter point. The physical matrices are

$$
T_\omega=\operatorname{diag}\left(\frac{\cos\theta_V}{\sqrt2},\frac{\cos\theta_V}{\sqrt2},-\sin\theta_V\right),\qquad
T_\phi=\operatorname{diag}\left(\frac{\sin\theta_V}{\sqrt2},\frac{\sin\theta_V}{\sqrt2},\cos\theta_V\right).
$$

## Photon overlaps and propagation

With $Q=\operatorname{diag}(2/3,-1/3,-1/3)$, the script obtains

$$
c_\omega=\langle QT_\omega\rangle=\frac{\cos\theta_V}{3\sqrt2}+\frac{\sin\theta_V}3,\qquad
c_\phi=\frac{\sin\theta_V}{3\sqrt2}-\frac{\cos\theta_V}3.
$$

The tensor photon couplings are $g_i=\sqrt2 eF_Vc_i$ as in the input diagonalization record. The isoscalar norm is $c_\omega^2+c_\phi^2=1/6$; including $c_\rho=1/\sqrt2$ restores $\sum_i c_i^2=2/3$.

For $\mathbf c=(1/(3\sqrt2),-1/3)^T$, the scalar flavor response, before $F_V$ factors and tensor kinematic/contact structures, is

$$
\mathbf c^T(\mathcal M_V^2-s\mathbf1)^{-1}\mathbf c
=\frac{2A+D+2\sqrt2B-3s}{18[(A-s)(D-s)-B^2]}
=\frac{c_\omega^2}{m_\omega^2-s}+\frac{c_\phi^2}{m_\phi^2-s}.
$$

To derive the vector response itself, restrict the quadratic action to the transverse electromagnetic tensor subspace $V_i^{\mu\nu}=r_iF^{\mu\nu}$ at $s=q^2$. The tensor kinetic term contributes $-s\mathbf1$ in this subspace. Suppressing the common $F_{\mu\nu}F^{\mu\nu}$ factor gives

$$
\mathcal L_{\mathrm{EM}}=\frac14\mathbf r^T(\mathcal M_V^2-s\mathbf1)\mathbf r-\frac12\mathbf g^T\mathbf r,
\qquad \mathbf g=\sqrt2 eF_V\mathbf c.
$$

The script differentiates this action and solves the two equations of motion to obtain the exact response within the stated quadratic truncation:

$$
\mathbf r=(\mathcal M_V^2-s\mathbf1)^{-1}\mathbf g
=\frac{eF_V}{3[(A-s)(D-s)-B^2]}
\begin{pmatrix}D-s+\sqrt2 B\\-\sqrt2(A-s)-B\end{pmatrix}.
$$

This electromagnetic-subspace response describes an external photon insertion at leading order in $e$; it does not add an electromagnetic mass shift or a second photon-vector insertion. The full tensor propagator also carries Lorentz projectors and contact structures absent from this restricted calculation. No complete propagator identity beyond this subspace is asserted. For a generic consistently rotated interaction vector $\mathbf j$, the script separately proves

$$
\sum_i\frac{(R^T\mathbf c)_i(R^T\mathbf j)_i}{M_{\rm common}^2-s}
=\frac{\mathbf c^T\mathbf j}{M_{\rm common}^2-s}.
$$

Thus rotating every vertex while retaining a common internal pole cannot change a complete exchange sum. Split masses and/or declared coupling breaking must enter an amplitude calculation explicitly.

## Verification and remaining work

The script runs successfully and prints each derived expression with SymPy's LaTeX renderer. Assertions verify orthogonality, positive-strange singlet-mixing sign, the selected diagonal eigenvalues, the electromagnetic-subspace equations of motion, ideal zero-mixing response, common-mass response and common-pole exchange invariance. No phenomenological mass or decay fit has been performed here.

Next: fit or constrain the declared mass parameters with documented mass inputs and truncation uncertainties; allow declared residue-breaking terms when calibrating leptonic/radiative data; rotate every odd vertex and redo the complete split-pole TFF/OPE analysis. Proca/tensor interaction matching requires induced local terms [@Ecker:1989yg]; this quadratic derivation supplies no shortcut for that uncomputed matching. A precise dynamical prescription additionally needs self-energies and momentum-dependent mixing rather than this constant real angle.
