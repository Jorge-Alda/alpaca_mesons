---
title: Tensor-vector leptonic widths and minimal source breaking
status: derived
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# Tensor-vector leptonic widths and minimal source breaking

## Inputs and assumptions

Use the canonical tensor action and generator normalization in [RChT basics](rcht-basics.md) [@Ecker:1988te], the leading electromagnetic reduction in [photon-vector diagonalization](rcht-photon-vector-diagonalization.md), and the positive-strange physical rotation in [nonideal mass action](vector-nonideal-mass-action.md). No additional literature input or local Alda manuscript is used. The extension below is an explicitly assumed invariant truncation, not a literature-established determination of its coefficient. Natural units, metric $(+---)$, tree level, real couplings, exact strong isospin and canonical vector kinetic terms apply. Electromagnetic mass shifts, loops and radiative corrections are omitted. The fine-structure constant is the physical Thomson-limit input at this order; finite photon normalization shifts are not added again.

## Current and decay amplitude

Define a signed dimension-one residue $R_i$ by the tensor transition

$$
\mathcal L_{V\gamma}=-\frac{e R_i}{\sqrt2}V^i_{\mu\nu}F^{\mu\nu}.
$$

The input first-order tensor reduction gives $\epsilon_i=\sqrt2 eR_i/M_i$ and the leading current coupling $e\epsilon_i$. Equivalently the electromagnetic-current matrix element is

$$
\langle0|J_Q^\mu|V_i(q,\varepsilon)\rangle
=\sqrt2R_iM_i\varepsilon^\mu,
\qquad
\mathcal M(V_i\to\ell^-\ell^+)
=\frac{\sqrt2e^2R_i}{M_i}\varepsilon_\mu\bar u\gamma^\mu v.
$$

An overall state phase has no effect on the width. This definition avoids confusing $R_i$ with a dimensionless photon-vector mixing parameter.

## Executable spin sum and phase space

Run `.venv/bin/python scripts/derive_vector_leptonic_width.py`. The [script](../../../scripts/derive_vector_leptonic_width.py) uses SymPy's explicit Dirac matrices to evaluate the lepton trace at $q=(M,0,0,0)$ and $p_\mp=(M/2,0,0,\pm p)$, then substitutes $p^2=M^2/4-m_\ell^2$. The three physical polarization sum is the spatial trace; dividing by three averages the initial spin. The result generated from that trace is

$$
\overline{|\mathcal M|^2}=
\frac{8e^4R_i^2}{3M_i^2}(M_i^2+2m_\ell^2).
$$

Integrating $d\Phi_2=\beta_\ell\,d\cos\vartheta/(16\pi)$ after azimuthal integration gives $\Phi_2=\beta_\ell/(8\pi)$, where $\beta_\ell=\sqrt{1-4m_\ell^2/M_i^2}$. Applying the flux $1/(2M_i)$ and substituting $e^2=4\pi\alpha_{\rm em}$ yields

$$
\Gamma(V_i\to\ell^+\ell^-)=
\frac{8\pi\alpha_{\rm em}^2R_i^2}{3M_i}
\left(1+\frac{2m_\ell^2}{M_i^2}\right)
\sqrt{1-\frac{4m_\ell^2}{M_i^2}}.
$$

The channel is closed for $M_i<2m_\ell$. Widths determine $|R_i|$, never its relative sign. A fit must retain a declared sign branch or test discrete branches with other amplitudes.

## Minimal one-spurion source interaction

Choose and declare

$$
\mathcal L_{\rm source}=
\frac{F_V}{2\sqrt2}\langle V_{\mu\nu}f_+^{\mu\nu}\rangle
+\frac{z_V}{2\sqrt2}\langle V_{\mu\nu}\{\chi_+,f_+^{\mu\nu}\}\rangle.
$$

Here $[F_V]=1$, $[z_V]=-1$, $[\chi_+]=2$, $[V_{\mu\nu}]=1$ and $[f_+]=2$ in mass units. It is a first-power spurion extension; higher powers and additional trace structures are omitted. In the real vacuum $\chi_+=\operatorname{diag}(x_l,x_l,x_s)$ and $f_+=-2eQF$, the script takes explicit traces of the already specified flavor generators and derives

$$
F_l=F_V+2z_Vx_l,\qquad F_s=F_V+2z_Vx_s,
$$

$$
R_\rho=\frac{F_l}{\sqrt2},\qquad
R_{\omega_I}=\frac{F_l}{3\sqrt2},\qquad
R_{\phi_I}=-\frac{F_s}{3},
$$

$$
R_\omega=\frac{F_l\cos\theta_V}{3\sqrt2}
+\frac{F_s\sin\theta_V}{3},\qquad
R_\phi=\frac{F_l\sin\theta_V}{3\sqrt2}
-\frac{F_s\cos\theta_V}{3}.
$$

Thus the same rotation used for mass eigenstates rotates every source overlap. Universal $F_V$ is the limit $z_V=0$. Distinct $F_l,F_s$ are only equivalent to independent $F_V,z_V$ when $x_l\ne x_s$; at fixed physical spurions they may be fitted first, then converted. This is source breaking, not kinetic breaking or an independent change of generators.

## Singlet trace limitation

The unbroken double trace $\langle V\rangle\langle f_+\rangle$ gives no vacuum electromagnetic transition because $\operatorname{Tr}Q=0$. Therefore an independent coefficient for that operator cannot repair the leptonic residues. A spurion double trace such as

$$
\frac{y_V}{2\sqrt2}\langle V_{\mu\nu}\rangle
\langle\chi_+f_+^{\mu\nu}\rangle,
\qquad [y_V]=-1,
$$

is nonzero: $\operatorname{Tr}(\chi_+Q)=(x_l-x_s)/3$. It adds
$\delta R_i=y_V\operatorname{Tr}T_i\,(x_l-x_s)/3$, with ideal-light/strange weights $\sqrt2$ and $1$. It leaves $\rho$ unchanged in exact isospin and supplies another isoscalar residue direction when $x_l\ne x_s$. This possible extension must be separately declared and calibrated; three leptonic widths alone would not determine four continuous parameters $(F_l,F_s,\theta_V,y_V)$.

The script differentiates $(R_\rho,R_{\omega_I},R_{\phi_I})$ with respect to $(F_V,z_V,y_V)$ and obtains determinant $-2(x_l-x_s)^2/9$. Thus at fixed mixing and nondegenerate spurions the three coefficients span all three signed neutral residues. Once this full rank extension is allowed, the three widths supply no remaining test of the mixing angle by themselves.

## Verification

The derivation script executes successfully. It generates the spin trace, spin average, integrated phase space, finite-mass width and all flavor residues using SymPy; assertions separately verify the massless width and rotation invariance of the isoscalar residue norm. No decay-data fit or numerical parameter selection is performed in this record.
