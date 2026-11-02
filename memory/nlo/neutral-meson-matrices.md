---
title: NLO kinetic and mass matrices in the neutral Gell-Mann basis
status: derived
created: 2026-11-02
updated: 2026-11-02
bibtex_keys: [Kaiser:1998ds, Alda:2026nrd]
provenance: analytical-derivation
---

# Neutral pseudoscalar matrices

## Assumptions and inputs

This record projects the two-meson NLO Lagrangian in [quadratic-lagrangian.md](quadratic-lagrangian.md), based on the operator definition in [lagrangian.md](lagrangian.md) and the field normalization in [Alda:2026nrd](../extracts/Alda:2026nrd/neutral-meson-generators.md). The literature inputs are [@Kaiser:1998ds; @Alda:2026nrd].

There is no ALP. Take the real quark mass matrix to be diagonal, $\chi_0=2B_0\operatorname{diag}(m_u,m_d,m_s)$, and use

$$
\Phi=\lambda_3\pi_3+\lambda_8\eta_8+\lambda_0\eta_0,\quad
\lambda_3=\operatorname{diag}(1,-1,0),\quad
\lambda_8=\frac{\operatorname{diag}(1,1,-2)}{\sqrt3},\quad
\lambda_0=\sqrt{\frac23}\mathbf 1.
$$

Thus $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$ for $i,j\in\{3,8,0\}$. Only the displayed local $\mathcal{L}^{(1)}$ and $\mathcal{L}^{(2)}$ terms are included; loops, higher chiral orders, and the leading mass matrix are not included.

## Definition and result

Set $P=(\pi_3,\eta_8,\eta_0)^T$ and define the raw NLO quadratic terms by

$$
\mathcal{L}_{\mathrm{NLO,quad}}=
\frac12\partial_\mu P^T\delta K\,\partial^\mu P
-\frac12P^T\delta M^2P.
$$

Both matrices are symmetric. Their independent entries are

$$
\begin{aligned}
\delta K_{33}&=\frac{8B_0}{f_\pi^2}\left[2L_4(m_u+m_d+m_s)+L_5(m_u+m_d)\right],\\
\delta K_{38}&=\frac{8B_0L_5}{\sqrt3f_\pi^2}(m_u-m_d),\\
\delta K_{30}&=\frac{4\sqrt6B_0}{3f_\pi^2}(2L_5+3L_{18})(m_u-m_d),\\
\delta K_{88}&=\frac{8B_0}{3f_\pi^2}\left[6L_4(m_u+m_d+m_s)+L_5(m_u+m_d+4m_s)\right],\\
\delta K_{80}&=\frac{4\sqrt2B_0}{3f_\pi^2}(2L_5+3L_{18})(m_u+m_d-2m_s),\\
\delta K_{00}&=\Lambda_1+\frac{16B_0}{3f_\pi^2}(3L_4+L_5+3L_{18})(m_u+m_d+m_s).
\end{aligned}
$$

$$
\begin{aligned}
\delta M^2_{33}&=\frac{32B_0^2}{f_\pi^2}\left[L_6(m_u+m_d+m_s)(m_u+m_d)+L_7(m_d-m_u)^2+L_8(m_u^2+m_d^2)\right],\\
\delta M^2_{38}&=-\frac{32\sqrt3B_0^2(m_d-m_u)}{3f_\pi^2}\left[L_6(m_u+m_d+m_s)+L_7(m_u+m_d-2m_s)+L_8(m_u+m_d)\right],\\
\delta M^2_{30}&=\frac{\sqrt6B_0(m_d-m_u)}{3f_\pi^2}\left[48B_0L_{25}(m_u+m_d)-32B_0(m_u+m_d+m_s)(L_6+L_7+L_8)-\Lambda_2f_\pi^2\right],\\
\delta M^2_{88}&=\frac{32B_0^2}{3f_\pi^2}\left[L_6(m_u+m_d+m_s)(m_u+m_d+4m_s)+L_7(m_u+m_d-2m_s)^2+L_8(m_u^2+m_d^2+4m_s^2)\right],\\
\delta M^2_{80}&=\frac{\sqrt2B_0}{3f_\pi^2}\left[32B_0(L_6+L_7)(m_u+m_d+m_s)(m_u+m_d-2m_s)\right.\\
&\left.\qquad\qquad\quad +(32B_0L_8-48B_0L_{25})(m_u^2+m_d^2-2m_s^2)+\Lambda_2f_\pi^2(m_u+m_d-2m_s)\right],\\
\delta M^2_{00}&=\frac{64B_0^2}{3f_\pi^2}\left[(L_6+L_7)(m_u+m_d+m_s)^2+(L_8-3L_{25})(m_u^2+m_d^2+m_s^2)\right]+\frac{4B_0\Lambda_2}{3}(m_u+m_d+m_s).
\end{aligned}
$$

The matrix ordering is $(3,8,0)$; obtain lower-triangular entries by symmetry. These are raw coefficients before canonical normalization. In particular, the local NLO $\delta M^2$ contains no $m_0^2$. If one subsequently canonically normalizes the fields, the NLO mass correction receives the separate contribution $-\tfrac12\{\delta K,M^2_{\mathrm{LO}}\}$.

## Method and reproducibility

The executable SymPy derivation substitutes the three diagonal generators into each operator of the two-meson Lagrangian and takes field or derivative Hessians. For the above convention, the kinetic Hessian is $\delta K$ and minus the field Hessian of the nondifferentiated terms is $\delta M^2$. Run:

    .venv/bin/python scripts/derive_nlo_neutral_matrices.py

The script prints the per-coupling matrices and the six independent entries of each summed matrix. Its internal checks cover generator normalization, decoupling of $\pi_3$ from $(\eta_8,\eta_0)$ at $m_u=m_d$, and the chiral limit.
