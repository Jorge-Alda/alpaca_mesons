---
title: Mixed ultraviolet poles and meson field renormalization
status: derived
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds]
provenance: analytical-derivation
---

# Mixed ultraviolet poles and meson field renormalization

## Inputs, conventions and assumptions

Use the loop and local coefficients in [mixed tadpole pole check](mixed-tadpole-pole-check.md), the leading mixed couplings from [quadratic Lagrangian](../quadratic-lagrangian.md), the leading mass matrix and the meson field redefinition used in [meson self-energy](meson-self-energy.md), and the source and field conventions of Alda:2026nrd and Kaiser:1998ds. Set $\chi_0=2B_0\operatorname{diag}(m_u,m_d,m_s)$ and $f_\pi=93\,\mathrm{MeV}$. The neutral flavour fields are $(\pi_3,\eta_8,\eta_0)$. All masses and decay constants are in GeV, $s$ in GeV$^2$, and $Q=sK-M^2+\Sigma$. The calculation is at first order in $1/f_a$ and one-loop order. The flavour coupling matrices are arbitrary and are not assumed to commute with $\chi_0$. No CP-violating sources or constant $\theta$ term are included. Resonance loops are excluded as in [loop classes](../counting/loop-classes.md).

## Field transformation and induced mixed terms

For the full nonet, write $\Phi_{\rm old}=\Phi+\lambda(\mu)\delta\Phi$ with the same field redefinition that cancels the meson poles:

$$
\delta\Phi=\frac1{f_\pi^2}\left[-\frac12\{\chi_0,\Phi\}-\frac13\operatorname{tr}(\chi_0)\Phi+\frac13\operatorname{tr}(\chi_0\Phi)\mathbf1+\frac13\operatorname{tr}(\Phi)\chi_0\right].
$$

Writing its component matrix as $\delta\Phi_i=R_{ij}\Phi_j$, the mixed block transforms at one-loop order as

$$
b_{i,\rm new}=b_{i,\rm old}+\lambda(\mu)R_{ji}b_{j}^{(0)}.
$$

The leading mixed blocks used to obtain the induced terms are

$$
b_i^{(0)}(s)=-\frac{f_\pi s}{f_a}\operatorname{tr}(c_q^A\lambda_i)
+\frac{f_\pi}{2f_a}\operatorname{tr}\!\left(c_q^P\{\chi_0,\lambda_i\}\right)
-\frac{f_\pi m_0^2}{3f_a}c_G\operatorname{tr}(\lambda_i).
$$

The $R^Tb^{(0)}$ terms are therefore ALP dependent even though $R$ is not. The script derives them with the nine generators without imposing commutation on the sources. For example, for external $\pi_3$, $c_q^A=\lambda_3$ and $c_q^P=c_G=0$, the induced pole is $+\lambda(\mu)s\,2B_0(5m_u+5m_d+2m_s)/(3f_af_\pi)$, exactly opposite to the residual found before the field transformation.

For arbitrary $c_q^A$ and $c_q^P$, every calculated mixed pole cancels entry by entry after adding $R^Tb^{(0)}$. The singlet row $R_{0j}$ vanishes, so the leading $c_G$ mixed block does not generate an induced pole. The implemented $c_q^V$ and $c_q^S$ quadratic mixed blocks and their loop and local poles vanish. An additional ALP-dependent linear meson shift is not needed to cancel the calculated quadratic poles. A general such shift would also change the mixed block by $(s\mathbf1-M_0^2)_{ij}\delta\Phi_j/a$; its induced kinetic and mass terms must be retained together. ALP self-energy changes from a shift are of second order in $1/f_a$ and outside the present calculation.

This check establishes cancellation for the implemented quadratic loop and local coefficients. It does not independently establish that their input operator list is complete, nor does it derive the transformations of interaction vertices or a renormalized physical amplitude. Those terms must also be transformed in an amplitude calculation.

## Verification

Run

```bash
.venv/bin/python scripts/derive_mixed_uv_field_renormalization.py
.venv/bin/python scripts/derive_nlo_loop_self_energies.py --mixed
```

The first script prints initial and induced source coefficient matrices and tests the meson cancellation, all nine mixed entries, the absence of an induced $c_G$ pole, and the quadratic change-of-fields algebra. The second prints the mixed tadpoles, local and induced poles and asserts their cancellation. The checks are exact symbolic checks, conditional on the linked input coefficients.
