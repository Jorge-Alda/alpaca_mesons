---
title: Mixed ultraviolet pole and meson field redefinition
status: derived
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds]
provenance: analytical-derivation
---

# Mixed ultraviolet pole and meson field redefinition

## Inputs and assumptions

Use the incomplete mixed pole calculation in [mixed tadpole pole check](mixed-tadpole-pole-check.md), the leading mass matrix and field conventions in [neutral meson matrices](../neutral-meson-matrices.md), and the inverse-propagator blocks $A,b,d$ and LSZ components $w^{(a)}$ in [loop corrections to ALP-meson LSZ mixing](../loop-corrections-to-alp-meson-lsz-mixing.md). The flavour-basis neutral fields are $\pi_3,\eta_8,\eta_0$. The calculation is at first order in $1/f_a$, with no CP-violating sources and no constant $\theta$ term. Flavour couplings are arbitrary matrices, and $\chi_0=2B_0\operatorname{diag}(m_u,m_d,m_s)$. Use $f_\pi=93\,$MeV as in Alda:2026nrd; masses and decay constants are in GeV, and $s$ is in GeV$^2$.

## Mixed pole for generic couplings

The script adds the calculated loop and local pole coefficients from the linked diagnostic for every external generator, and leaves the source matrices arbitrary inside the flavour traces. For example, the $\pi_3$ entry is $\lambda(\mu)$ times

$$
\begin{aligned}
&\frac{2 B_0 s}{3f_af_\pi}\operatorname{tr}\!\left[c_q^A\operatorname{diag}(-3m_u-2m_d-m_s,\,2m_u+3m_d+m_s,\,m_u-m_d)\right]\\
&+\frac{4 B_0^2}{3f_af_\pi}\operatorname{tr}\!\left[c_q^P\operatorname{diag}\!\left(m_u(3m_u+2m_d+m_s),\,-m_d(2m_u+3m_d+m_s),\,m_s(m_d-m_u)\right)\right].
\end{aligned}
$$

Only diagonal entries of the source matrices enter this particular flavour-basis entry. Flavour-off-diagonal couplings contribute to other entries, as printed by the script. For the currently implemented terms, the $\eta_0$ residual vanishes and the $c_q^V,c_q^S,c_G$ UV coefficients vanish. These statements concern the current one-loop and local pole diagnostic, not the full leading ALP–meson mixing or a completed renormalization calculation. The source coefficients are based on the existing implementation; a full independent check of all local terms and induced field transformations remains open.

## Change of fields

Let $\delta\Phi_i$ be a local, linear, ALP-dependent shift of the old meson field, so that the old field equals the new field plus $\delta\Phi_i$, and let $\delta\Phi_i/a$ be independent of the fields. At first order in the shift, substituting into the quadratic form with blocks $A,b,d$ gives

$$
b_{i,\mathrm{new}}(s)=b_{i,\mathrm{old}}(s)+A_{ij}(s)\frac{\delta\Phi_j}{a},\qquad
 d_{\mathrm{new}}(s)=d_{\mathrm{old}}(s)+2b_{i,\mathrm{old}}(s)\frac{\delta\Phi_i}{a}.
$$

For several coupled mesons, this follows from multiplication of the full symmetric inverse propagator by the triangular change-of-fields matrix and its transpose; the script performs the one-meson case explicitly. At first order in ALP–meson mixing, the second term in $d_{\mathrm{new}}$ is of second order in that mixing. For an invertible meson block, $w^{(a)}_{i,\mathrm{new}}=w^{(a)}_{i,\mathrm{old}}-\delta\Phi_i/a$. Direct ALP vertices also change when the old fields are replaced, so the mixing contribution alone is basis dependent; physical amplitudes require both terms. Pole positions are invariant under an invertible field transformation.

An $s$-dependent residual cannot be declared removed solely because a local shift changes $b$: $A_{ij}^{(0)}=s\delta_{ij}-(M_0^2)_{ij}$ also generates a mass-dependent term. Its cancellation requires checking all transformed terms and field renormalizations. A divergent shift also needs a consistent treatment of renormalized fields. No such cancellation is claimed here.

## Verification

Run `.venv/bin/python scripts/derive_mixed_field_redefinition.py`. SymPy prints the coefficient matrices and the quadratic-form transformation; unit tests reproduce the $c_q^A=\lambda_3$ diagnostic, verify the calculated relation between pseudoscalar and axial coefficient matrices, and check the change of the mixed block. This tests the algebra conditional on the inputs; it does not establish completeness of those inputs.
