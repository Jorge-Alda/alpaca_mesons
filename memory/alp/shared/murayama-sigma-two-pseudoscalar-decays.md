---
title: Radial-mode decays to two pseudoscalars in the only-up model
status: derived
created: 2026-11-18
updated: 2026-11-18
bibtex_keys: [DiLuzio:2026poh, Alda:2026nrd]
provenance: analytical-derivation
---

# Radial-mode decays to two pseudoscalars

## Inputs and assumptions

Use the model and canonical polar normalization of [DiLuzio:2026poh](../../extracts/DiLuzio:2026poh/model-lagrangian-and-field-normalization.md), its [neutral mass matrix](../../extracts/DiLuzio:2026poh/angular-mode-couplings.md), and the $\mathbb M_s$ operator and exact pole construction of the [finite-mixing record](alp-lo-finite-mixing-operator-inverse.md). The [radial and phase coefficient convention](alp-murayama-radial-phase-mass-coefficients.md) uses ordinary coefficients of dimensionful fields; the [dressed-source construction](alp-general-chiral-redefinition-invariants.md) gives their meaning under allowed chiral redefinitions. The [pole-aligned transformation](alp-pole-alignment-chiral-generator.md) specifies the corresponding ALP-dependent chiral change of fields.

Set $(n_u,n_d,n_s)=(1,0,0)$, $f_a=f_\varphi$, $\theta=0$ and use a CP-conserving vacuum. Work at leading $U(3)$ chiral order and tree level, with no expansion in $f_\pi/f_a$, electromagnetic corrections, higher chiral orders or scalar resonances. Fields, masses, $B_0$, $f_\pi$ and $f_a$ are in GeV; widths and couplings are in GeV. The calculation assumes stable simple neutral poles. The model takes $m_\varphi^2>0$ and $\lambda>0$; a stable vacuum suffices for the local expansion.

## Dressed source and field expansion

In the partonic convention $-\bar q_R\widehat M q_L-\bar q_L\widehat M^\dagger q_R$, $\widehat M=\widehat S+i\widehat P$. Expand the invariant sources in *dimensionful* $\sigma$ and $a$ as in the radial-phase record:
$$
\widehat S(a,\sigma)=\sum_{r,s\ge0}\widehat S_{rs}\sigma^r a^s,\qquad
\widehat P(a,\sigma)=\sum_{r,s\ge0}\widehat P_{rs}\sigma^r a^s.
$$
These are ordinary monomial coefficients, so $\widehat S_{rs}=\left.\partial_\sigma^r\partial_a^s\widehat S\right|_0/(r!s!)$, and likewise for $\widehat P_{rs}$. Every coefficient is a constant flavour matrix. In particular $\widehat S_{10}$ is an independent radial coupling, not a derivative of the constant $\widehat S_{00}$.

The paper uses the conjugate source $M_q$ in the $\bar q_Lq_R$ term. Its $P_{rs}$ has the opposite sign to $\widehat P_{rs}$ in the partonic convention used here. In the original connection-free representative, dressing leaves the source unchanged and
$$
\widehat M=\operatorname{diag}\!\left[m_u\left(1+\frac\sigma{f_a}\right)e^{-ia/f_a},m_d,m_s\right].
$$
Direct expansion gives
$$
\begin{aligned}
\widehat S_{00}&=\operatorname{diag}(m_u,m_d,m_s),&
\widehat S_{10}&=\frac1{f_a}\operatorname{diag}(m_u,0,0),\\
\widehat P_{01}&=-\frac1{f_a}\operatorname{diag}(m_u,0,0),&
\widehat P_{11}&=-\frac1{f_a^2}\operatorname{diag}(m_u,0,0),\\
\widehat S_{02}&=-\frac1{2f_a^2}\operatorname{diag}(m_u,0,0),&
\widehat S_{12}&=-\frac1{2f_a^3}\operatorname{diag}(m_u,0,0).
\end{aligned}
$$
The first index counts $\sigma$, the second counts $a$. The factor $1/2$ in $\widehat S_{02}$ and $\widehat S_{12}$ is required by ordinary Taylor coefficients. The earlier derivative convention with $x=a/f_a$ is related by $\widehat P_1=f_a\widehat P_{01}$, $\widehat S_2=2f_a^2\widehat S_{02}$, $\partial_\sigma\widehat P_1|_0=f_a\widehat P_{11}$ and $\partial_\sigma\widehat S_2|_0=2f_a^2\widehat S_{12}$.

The values $c_q^{AP}=\operatorname{diag}(1/2,0,0)$ and $c_{GA}=0$ describe the angular field in this model; they do not define a scalar invariant. For example replacing the up-quark source by $(m_u+y_u\sigma)e^{-ia/f_a}$ leaves the angular coupling fixed while changing $\widehat S_{10}$ and hence the radial decay. The scalar potential is further independent input.

## Pole equation and ALP-dependent field change

The operator on Hermitian neutral flavour matrices is
$$
\mathbb M_s(X)=B_0\{\widehat S_{00},X\}+\frac{m_1^2}{3}\operatorname{tr}(X)\mathbf1-sX.
$$
For a regular simple pole $s=m_i^2$ with nonzero angular component, the meson field shift is obtained from $-2B_0f_\pi\mathbb M_{m_i^2}^{-1}(\widehat P_{01})$; its generator implements the ALP-dependent chiral rotation of the linked pole-alignment record. It eliminates the mixed inverse propagator at the selected pole, not at all momenta. All four physical poles with nonzero angular component solve
$$
s+2B_0f_\pi^2\operatorname{tr}(\widehat S_{02})+2B_0^2f_\pi^2\operatorname{tr}\!\left[\widehat P_{01}\mathbb M_s^{-1}(\widehat P_{01})\right]=0.
$$
At each pole the normalization factor for the unnormalized angular component equal to one is
$$
1+2B_0^2f_\pi^2\operatorname{tr}\!\left[\mathbb M_{m_i^2}^{-1}(\widehat P_{01})^2\right].
$$
An orthogonal transformation of the complete neutral fields also diagonalizes the quadratic terms, but no separate meson flavour rotation is needed. Each physical field's sign is arbitrary; below its angular component is chosen positive.

## Neutral coupling

With $\mathcal L\supset-\frac12\sigma\sum_{ij}g_{ij}P_iP_j$, substitute both external pole directions and their normalization into the cubic source and kinetic vertices. This gives, without implicit eigenvector components,
$$
\begin{aligned}
g_{ij}={}&\left\{
2B_0^3f_\pi^2\operatorname{tr}\!\left[\widehat S_{10}\left\{\mathbb M_{m_i^2}^{-1}(\widehat P_{01}),\mathbb M_{m_j^2}^{-1}(\widehat P_{01})\right\}\right]\right.\\
&\left. -2B_0^2f_\pi^2\operatorname{tr}\!\left[\widehat P_{11}\left(\mathbb M_{m_i^2}^{-1}(\widehat P_{01})+\mathbb M_{m_j^2}^{-1}(\widehat P_{01})\right)\right]
-2B_0f_\pi^2\operatorname{tr}(\widehat S_{12})
+\frac{m_\sigma^2-m_i^2-m_j^2}{f_a}\right\}\\
&\times\prod_{r=i,j}\left\{1+2B_0^2f_\pi^2\operatorname{tr}\!\left[\mathbb M_{m_r^2}^{-1}(\widehat P_{01})^2\right]\right\}^{-1/2}.
\end{aligned}
$$
The first three terms are generated by $\widehat S_{10}$, $\widehat P_{11}$ and $\widehat S_{12}$; the last is the on-shell contribution of $\sigma(\partial a)^2/f_a$ in the polar kinetic term. The formula includes the axion-like neutral state among the daughters. This expression is specialized to the stated polar scalar kinetic term and model; the first-order angular invariants alone do not fix the scalar couplings.

The [Cartesian and polar field-coordinate record](alp-polar-cartesian-amplitude-equivalence.md) gives the general on-shell equivalence. In Cartesian coordinates the cubic terms independently obtained for this model are
$$
\mathcal L_{\sigma PP}=-\frac{B_0m_u}{f_a}\sigma\left[\frac12\left(\pi^0+\frac{\eta_8}{\sqrt3}+\sqrt{\frac23}\eta_0\right)^2+\pi^+\pi^-+K^+K^-\right]-\lambda f_a\sigma a^2.
$$
The vacuum condition yields
$$
m_\sigma^2=\frac{B_0m_uf_\pi^2}{f_a^2}+2\lambda f_a^2,\qquad
m_\varphi^2=\frac{B_0m_uf_\pi^2}{f_a^2}-\lambda f_a^2.
$$
Using this relation and the pole equations gives the same on-shell coupling in both field parametrizations. The quartic $\lambda$ is not fixed by the pseudoscalar spectrum.

## Charged states, width and limits

At this order
$$
g_{\pi^+\pi^-}=g_{K^+K^-}=B_0\operatorname{tr}(\widehat S_{10})=\frac{B_0m_u}{f_a},\qquad
g_{K^0\bar K^0}=0.
$$
Other two-meson charge combinations are excluded by electric charge and strangeness conservation. For $m_\sigma>m_i+m_j$,
$$
\Gamma(\sigma\to P_iP_j)=\frac{|g_{ij}|^2}{16\pi m_\sigma^3(1+\delta_{ij})}
\sqrt{[m_\sigma^2-(m_i+m_j)^2][m_\sigma^2-(m_i-m_j)^2]}.
$$
The width vanishes at or below threshold; particle-antiparticle pairs use $\delta_{ij}=0$. The operator expression requires an invertible $\mathbb M_{m_i^2}$ and a nonzero angular component. A decoupled state or an exceptional meson pole needs a separate projection or a controlled limit; degenerate poles require an orthonormal basis of their full eigenspace. Broad or overlapping states, rescattering, higher chiral orders and scalar resonances require additional treatment. The [source estimate](../../extracts/DiLuzio:2026poh/radial-mode-pion-decay-estimate.md) is approximate and does not resolve the charge channels.

For redefinitions depending on $\sigma$, radial connections and all induced derivative terms must be retained; the displayed kinetic term uses the original connection-free polar representation. The on-shell physical amplitude is invariant under the corresponding complete field changes.

## Verification

The [operator derivation](../../../scripts/derive_murayama_sigma_operator_amplitude.py) calculates the pole directions, normalization and polar vertex and directly differentiates the dressed source using the ordinary $(\sigma,a)$ coefficients above. The [Cartesian derivation](../../../scripts/derive_murayama_sigma_pp.py) checks the polar and Cartesian vertices on shell. Their accompanying unit tests verify the coefficients, source contractions, mass-operator equation, charged vertices and on-shell equality. Additional independent representations and checks remain in `scripts/derive_murayama_sigma_explicit_mixing.py`, `scripts/derive_murayama_sigma_normal_basis.py`, `scripts/derive_murayama_sigma_invariant_couplings.py`, `scripts/derive_murayama_sigma_matrix_invariants.py` and `scripts/derive_murayama_sigma_source_coefficients.py`, each paired with its unit test. Their determinant and adjugate forms are alternative algebraic residue representations, not additional physical assumptions.

```sh
.venv/bin/python scripts/derive_murayama_sigma_operator_amplitude.py
.venv/bin/python -m unittest scripts.test_derive_murayama_sigma_operator_amplitude scripts.test_derive_murayama_sigma_pp
```
