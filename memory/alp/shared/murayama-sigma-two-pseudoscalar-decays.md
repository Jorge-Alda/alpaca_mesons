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

Use $\Phi=\sum_{k=3,8,0}P_k\lambda_k$, $\operatorname{tr}(\lambda_k\lambda_l)=2\delta_{kl}$ and $U=\exp(i\Phi/f_\pi)$. In the partonic source convention above, the source and anomaly terms at quadratic order are
$$
\begin{aligned}
\mathcal L_2={}&\frac14\operatorname{tr}(\partial_\mu\Phi\partial^\mu\Phi)+\frac12(\partial a)^2
-\frac{B_0}{2}\operatorname{tr}(\widehat S_{00}\Phi^2)
-B_0f_\pi a\operatorname{tr}(\widehat P_{01}\Phi)\\
&+B_0f_\pi^2a^2\operatorname{tr}(\widehat S_{02})
-\frac{m_1^2}{12}(\operatorname{tr}\Phi)^2.
\end{aligned}
$$
For a momentum with $p^2=s$, variation of $\Phi$ and $a$ gives
$$
\begin{aligned}
\mathbb M_s(\Phi)&=-2B_0f_\pi\widehat P_{01}\,a,\\
\left[s+2B_0f_\pi^2\operatorname{tr}(\widehat S_{02})\right]a
&=B_0f_\pi\operatorname{tr}(\widehat P_{01}\Phi),
\end{aligned}
$$
where the operator on Hermitian neutral flavour matrices is
$$
\mathbb M_s(X)=B_0\{\widehat S_{00},X\}+\frac{m_1^2}{3}\operatorname{tr}(X)\mathbf1-sX.
$$
For an invertible $\mathbb M_s$, the first equation gives
$$
\Phi=-2B_0f_\pi\mathbb M_s^{-1}(\widehat P_{01})\,a.
$$
Substitution into the second equation gives the pole equation for $a\ne0$:
$$
s+2B_0f_\pi^2\operatorname{tr}(\widehat S_{02})+2B_0^2f_\pi^2\operatorname{tr}\!\left[\widehat P_{01}\mathbb M_s^{-1}(\widehat P_{01})\right]=0.
$$
At a regular simple pole $s=m_i^2$, this is the linear meson shift of the ALP-dependent chiral rotation in the linked pole-alignment record. The shift eliminates the mixed inverse propagator at the selected pole, not at all momenta. For its normalization use the existing pole-shift coefficients $w_k^{(i)}$:
$$
\sum_k w_k^{(i)}\lambda_k=-2B_0f_\pi\mathbb M_{m_i^2}^{-1}(\widehat P_{01}),\qquad
\frac12\operatorname{tr}\left[\left(\sum_k w_k^{(i)}\lambda_k\right)^2\right]=\sum_k\left(w_k^{(i)}\right)^2.
$$
Thus the squared norm of the pole direction with angular component one is
$$
1+2B_0^2f_\pi^2\operatorname{tr}\!\left[\mathbb M_{m_i^2}^{-1}(\widehat P_{01})^2\right].
$$
### Basis with no ALP–meson mixing

Select the physical ALP-like pole $s=m_a^2$ and its shift coefficients $w_i^{(a)}$ above. Their squared norm is the established $N(m_a^2)=1+\sum_i(w_i^{(a)})^2$. The ALP-dependent chiral transformation is
$$
e^{i\Phi/f_\pi}=e^{ia'\sum_iw_i^{(a)}\lambda_i/(2f_\pi)}e^{i\Phi'/f_\pi}e^{ia'\sum_iw_i^{(a)}\lambda_i/(2f_\pi)},\qquad a=a',
$$
whose linear part is $P_i=P_i'+w_i^{(a)}a'$. At the order needed for $\sigma PP$, this linear part suffices: the next term in the symmetric chiral map is cubic in the pseudoscalar fields.

The pole shift alone leaves a mixed kinetic term. Complete the change of variables by
$$
a'=a''-\frac{\sum_iw_i^{(a)}P_i'}{N(m_a^2)}.
$$
Using the pole equation, the quadratic Lagrangian in $(P_i',a'')$ is
$$
\begin{aligned}
\mathcal L_2={}&\frac12\sum_{ij}\left(\delta_{ij}-\frac{w_i^{(a)}w_j^{(a)}}{N(m_a^2)}\right)\partial_\mu P_i'\partial^\mu P_j'
-\frac12\sum_{ij}P_i'P_j'\left[\frac{\operatorname{tr}(\lambda_i\mathbb M_0(\lambda_j))}{2}-\frac{m_a^2w_i^{(a)}w_j^{(a)}}{N(m_a^2)}\right]\\
&+\frac{N(m_a^2)}2\left[(\partial a'')^2-m_a^2(a'')^2\right].
\end{aligned}
$$
There is no mixed ALP–meson kinetic or mass term at any momentum. Define $a_{\rm phys}=\sqrt{N(m_a^2)}\,a''$; canonically normalize the positive meson kinetic matrix and diagonalize its remaining symmetric mass matrix with a meson-only rotation. This does not reintroduce ALP mixing. The completion is a general linear field change, in addition to the chiral shift; a single shift does not achieve the displayed decoupling at every momentum. In a degenerate eigenspace choose an orthonormal basis.

The full transformation before the final meson-only mass rotation is
$$
\begin{pmatrix}P_i\\a\end{pmatrix}
=\begin{pmatrix}
\displaystyle\delta_{ij}-\frac{w_i^{(a)}w_j^{(a)}}{\sqrt{N(m_a^2)}[\sqrt{N(m_a^2)}+1]}&\displaystyle\frac{w_i^{(a)}}{\sqrt{N(m_a^2)}}\\[2mm]
\displaystyle-\frac{w_j^{(a)}}{\sqrt{N(m_a^2)}}&\displaystyle\frac1{\sqrt{N(m_a^2)}}
\end{pmatrix}
\begin{pmatrix}P_j''\\a_{\rm phys}\end{pmatrix}.
$$
The fields $P_j''$ have canonical kinetic terms. The matrix is orthogonal and its final column is the normalized ALP-like pole direction. Each physical field's sign is arbitrary; below its angular component is chosen positive.

## Neutral coupling

Expand the dressed mass term $\tfrac12B_0f_\pi^2\operatorname{tr}(\widehat M U+\widehat M^\dagger U^\dagger)$ to one $\sigma$ and two outgoing fields. After dropping its linear term,
$$
\begin{aligned}
\mathcal L_{\sigma PP}={}&-\frac{B_0\sigma}{2}\operatorname{tr}(\widehat S_{10}\Phi^2)
-B_0f_\pi\sigma a\operatorname{tr}(\widehat P_{11}\Phi)
+B_0f_\pi^2\sigma a^2\operatorname{tr}(\widehat S_{12})
+\frac{\sigma}{f_a}(\partial a)^2.
\end{aligned}
$$
The last term follows from the polar kinetic term:
$$
\frac12\left(1+\frac{\sigma}{f_a}\right)^2(\partial a)^2
=\frac12(\partial a)^2+\frac{\sigma}{f_a}(\partial a)^2+O(\sigma^2).
$$
Substitute the chiral shift $\Phi=\Phi'+a'\sum_iw_i^{(a)}\lambda_i$ into the cubic terms:
$$
\begin{aligned}
\mathcal L_{\sigma PP}={}&-\frac{B_0\sigma}{2}\operatorname{tr}(\widehat S_{10}\Phi'^2)
-B_0\sigma a'\left[\frac12\operatorname{tr}\!\left(\widehat S_{10}\left\{\Phi',\sum_iw_i^{(a)}\lambda_i\right\}\right)+f_\pi\operatorname{tr}(\widehat P_{11}\Phi')\right]\\
&+\sigma(a')^2\left[-\frac{B_0}{2}\operatorname{tr}\!\left(\widehat S_{10}\left(\sum_iw_i^{(a)}\lambda_i\right)^2\right)-B_0f_\pi\operatorname{tr}\!\left(\widehat P_{11}\sum_iw_i^{(a)}\lambda_i\right)+B_0f_\pi^2\operatorname{tr}(\widehat S_{12})\right]
+\frac{\sigma}{f_a}(\partial a')^2.
\end{aligned}
$$
Then use $a'=(a_{\rm phys}-\sum_iw_i^{(a)}P_i'')/\sqrt{N(m_a^2)}$ and $\Phi'=\sum_{ij}\lambda_i[(\mathbf1-ww^{\mathsf T}/N(m_a^2))^{-1/2}]_{ij}P_j''$, followed by the meson-only mass rotation. For the resulting complete orthogonal field transformation $T$, with the same notation as the [multi-pole record](alp-multi-pole-alignment.md), the transformed on-shell vertices are explicitly
$$
\begin{aligned}
g_{ij}={}&\frac{B_0}{2}\operatorname{tr}\left[\widehat S_{10}\left\{\sum_kT_{ki}\lambda_k,\sum_lT_{lj}\lambda_l\right\}\right]
+B_0f_\pi\operatorname{tr}\left[\widehat P_{11}\left(T_{ai}\sum_kT_{kj}\lambda_k+T_{aj}\sum_kT_{ki}\lambda_k\right)\right]\\
&+\left[-2B_0f_\pi^2\operatorname{tr}(\widehat S_{12})+\frac{m_\sigma^2-m_i^2-m_j^2}{f_a}\right]T_{ai}T_{aj}.
\end{aligned}
$$
These include all $\sigma a_{\rm phys}a_{\rm phys}$, $\sigma a_{\rm phys}P_i$ and $\sigma P_iP_j$ vertices, with their derivative contribution. The transformed coefficients depend on the selected pole, while the on-shell $g_{ij}$ do not.

To write the physical result without displaying the auxiliary meson-only rotation, insert the explicit pole directions derived above:
$$
\Phi=\sum_{r=i,j}\frac{P_r}{\sqrt{1+\sum_k(w_k^{(r)})^2}}\sum_k w_k^{(r)}\lambda_k,
\qquad
a=\sum_{r=i,j}\frac{P_r}{\sqrt{1+\sum_k(w_k^{(r)})^2}}.
$$
For $i=j$, the displayed expression uses one field $P_i$ rather than adding it twice. Write the interaction as $\mathcal L\supset-\frac12\sigma\sum_{ij}g_{ij}P_iP_j$. The derivatives of the first three terms give, before external normalization,
$$
\begin{aligned}
&\frac{B_0}{2}\operatorname{tr}\!\left[\widehat S_{10}\left\{\sum_kw_k^{(i)}\lambda_k,\sum_lw_l^{(j)}\lambda_l\right\}\right]
+B_0f_\pi\operatorname{tr}\!\left[\widehat P_{11}\left(\sum_kw_k^{(i)}\lambda_k+\sum_lw_l^{(j)}\lambda_l\right)\right]\\
&\hspace{70mm}-2B_0f_\pi^2\operatorname{tr}(\widehat S_{12}).
\end{aligned}
$$
The two derivatives of $a^2$ give the factor $2$ in the final term. With outgoing $p_i,p_j$, the kinetic vertex and on-shell relation give
$$
\frac{2p_i\cdot p_j}{f_a}
=\frac{m_\sigma^2-m_i^2-m_j^2}{f_a},\qquad
(p_i+p_j)^2=m_\sigma^2.
$$
Substituting the explicit pole shifts and their normalization yields
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

## Cartesian derivation

Use the coordinates $s,p$ of the [polar and Cartesian record](alp-polar-cartesian-amplitude-equivalence.md), with $f_a=f_\varphi$. The Cartesian radial field $s$ in this section is not the momentum-squared argument of $\mathbb M_s$. In natural units,
$$
f_a\varphi=\frac{f_a+\sigma}{\sqrt2}e^{ia/f_a}
=\frac{f_a+s+ip}{\sqrt2},\qquad
f_a^2|\partial_\mu\varphi|^2=\frac12(\partial s)^2+\frac12(\partial p)^2.
$$
Near the vacuum the coordinate change and its inverse begin as
$$
\begin{aligned}
s&=\sigma-\frac{a^2}{2f_a}+O(\text{fields}^3),&
p&=a+\frac{\sigma a}{f_a}+O(\text{fields}^3),\\
\sigma&=s+\frac{p^2}{2f_a}+O(\text{fields}^3),&
a&=p-\frac{sp}{f_a}+O(\text{fields}^3).
\end{aligned}
$$
The identity linear transformation implies the same quadratic neutral mixing and the same no-mixing field transformation derived above. The nonlinear terms modify the cubic vertices, so $s$ and $\sigma$ must be distinguished there.

At $U=\mathbf1$, including the vacuum quark source, the potential is
$$
V=\frac{m_\varphi^2}{2}\big[(f_a+s)^2+p^2\big]
+\frac\lambda4\big[(f_a+s)^2+p^2\big]^2
-B_0m_uf_\pi^2\left(1+\frac{s}{f_a}\right)+\text{constant}.
$$
The stationarity condition, scalar masses and cubic potential derivative are
$$
\begin{aligned}
\left.\frac{\partial V}{\partial s}\right|_0&=f_am_\varphi^2+\lambda f_a^3-\frac{B_0m_uf_\pi^2}{f_a}=0,\\
m_\varphi^2&=\frac{B_0m_uf_\pi^2}{f_a^2}-\lambda f_a^2,\\
\left.\frac{\partial^2V}{\partial s^2}\right|_0&=m_\sigma^2=\frac{B_0m_uf_\pi^2}{f_a^2}+2\lambda f_a^2,\\
\left.\frac{\partial^2V}{\partial p^2}\right|_0&=\frac{B_0m_uf_\pi^2}{f_a^2},&
\left.\frac{\partial^3V}{\partial s\,\partial p^2}\right|_0&=2\lambda f_a.
\end{aligned}
$$
The last two equations concern the unmixed Cartesian field $p$, not the physical ALP mass.

The Cartesian mass spurion and the expansion of the Goldstone field are
$$
\widehat M=\widehat S_{00}+s\widehat S_{10}+ip\widehat P_{01},\qquad
U=\mathbf1+\frac{i\Phi}{f_\pi}-\frac{\Phi^2}{2f_\pi^2}+O(\Phi^3),
\qquad \Phi=\sum_{k=3,8,0}P_k\lambda_k+\text{charged fields}.
$$
In this equation $\widehat S_{00},\widehat S_{10},\widehat P_{01}$ have the values in the polar source expansion above; the higher Cartesian coefficients $\widehat P_{11}$ and $\widehat S_{12}$ vanish. This is the original connection-free representative, in which the source is linear in the Cartesian fields. The source term and potential yield
$$
\begin{aligned}
\mathcal L_{sPP}
&=-\frac{B_0s}{2}\operatorname{tr}(\widehat S_{10}\Phi^2)-\lambda f_a s p^2\\
&=-\frac{B_0m_u}{f_a}s\left[\frac12\left(\pi^0+\frac{\eta_8}{\sqrt3}+\sqrt{\frac23}\eta_0\right)^2+\pi^+\pi^-+K^+K^-\right]-\lambda f_a s p^2.
\end{aligned}
$$
In particular the source has no $sp\Phi$ vertex: the polar $\widehat P_{11}$ and $\widehat S_{12}$ vertices are replaced by the Cartesian scalar-potential vertex.

Apply the same complete orthogonal no-mixing transformation $T$ displayed above, with $p=\sum_iT_{ai}P_i$ and $\Phi=\sum_iP_i\sum_kT_{ki}\lambda_k$. Differentiating the Cartesian cubic Lagrangian gives
$$
g_{ij}=\frac{B_0}{2}\operatorname{tr}\!\left[\widehat S_{10}\left\{\sum_kT_{ki}\lambda_k,\sum_lT_{lj}\lambda_l\right\}\right]+2\lambda f_aT_{ai}T_{aj}.
$$
For regular poles insert the explicit directions from the operator equation; this eliminates the entries of $T$ and gives
$$
\begin{aligned}
g_{ij}={}&\left\{2\lambda f_a+
2B_0^3 f_\pi^2\operatorname{tr}\!\left[\widehat S_{10}\left\{\mathbb M_{m_i^2}^{-1}(\widehat P_{01}),\mathbb M_{m_j^2}^{-1}(\widehat P_{01})\right\}\right]\right\}\\
&\times\prod_{r=i,j}\left\{1+2B_0^2f_\pi^2\operatorname{tr}\!\left[\mathbb M_{m_r^2}^{-1}(\widehat P_{01})^2\right]\right\}^{-1/2}.
\end{aligned}
$$
To compare with the polar vertex, the angular row of the pole equation gives, for each $r=i,j$,
$$
\left(\frac{B_0m_uf_\pi^2}{f_a^2}-m_r^2\right)
=\frac{B_0m_uf_\pi}{f_a}
\left(w_3^{(r)}+\frac{w_8^{(r)}}{\sqrt3}+\sqrt{\frac23}w_0^{(r)}\right).
$$
Before the common external normalization, substituting the model's source coefficients into the polar expression and using the pole equations gives
$$
\begin{aligned}
&\frac{B_0m_uf_\pi^2}{f_a}
\prod_{r=i,j}\left[\frac1{f_a}-\frac1{f_\pi}\left(w_3^{(r)}+\frac{w_8^{(r)}}{\sqrt3}+\sqrt{\frac23}w_0^{(r)}\right)\right]
+\frac{m_\sigma^2-m_i^2-m_j^2}{f_a}\\
&\hspace{8mm}=2\lambda f_a+\frac{B_0m_u}{f_a}
\prod_{r=i,j}\left(w_3^{(r)}+\frac{w_8^{(r)}}{\sqrt3}+\sqrt{\frac23}w_0^{(r)}\right).
\end{aligned}
$$
The left-hand side is the polar vertex and the right-hand side the Cartesian vertex. Their common normalization gives $g_{ij}^{\rm polar}=g_{ij}^{\rm Cartesian}$ on shell. The comparison without the normalized shift does not require dividing by an angular component. The quartic $\lambda$ is independent scalar-potential input and is not fixed by the pseudoscalar spectrum.

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

The operator script constructs the shifted and kinetic-orthogonalized matrices for two generic meson fields and verifies that both ALP–meson blocks vanish. This checks the algebra for a general symmetric meson block; the displayed three-field expression follows by the same matrix multiplication.

The [operator derivation](../../../scripts/derive_murayama_sigma_operator_amplitude.py) calculates the pole directions, normalization and polar vertex, directly differentiates the dressed source using the ordinary $(\sigma,a)$ coefficients above, and expands the diagonal source and neutral meson fields to the displayed cubic order. The [Cartesian derivation](../../../scripts/derive_murayama_sigma_pp.py) checks the polar and Cartesian vertices on shell. Their accompanying unit tests verify the coefficients, cubic source expansion, source contractions, mass-operator equation, charged vertices and on-shell equality. Additional independent representations and checks remain in `scripts/derive_murayama_sigma_explicit_mixing.py`, `scripts/derive_murayama_sigma_normal_basis.py`, `scripts/derive_murayama_sigma_invariant_couplings.py`, `scripts/derive_murayama_sigma_matrix_invariants.py` and `scripts/derive_murayama_sigma_source_coefficients.py`, each paired with its unit test. Their determinant and adjugate forms are alternative algebraic residue representations, not additional physical assumptions.

```sh
.venv/bin/python scripts/derive_murayama_sigma_operator_amplitude.py
.venv/bin/python -m unittest scripts.test_derive_murayama_sigma_operator_amplitude scripts.test_derive_murayama_sigma_pp
```
