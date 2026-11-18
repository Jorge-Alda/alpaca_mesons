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

Use the model and canonical polar normalization of [DiLuzio:2026poh](../../extracts/DiLuzio:2026poh/model-lagrangian-and-field-normalization.md), its [neutral mass matrix](../../extracts/DiLuzio:2026poh/angular-mode-couplings.md), and the $\mathbb M_s$ operator and exact pole construction of the [finite-mixing record](alp-lo-finite-mixing-operator-inverse.md). The [radial and phase coefficient record](alp-murayama-radial-phase-mass-coefficients.md) gives the matrix Taylor coefficients of the source; the [dressed-source construction](alp-general-chiral-redefinition-invariants.md) gives their meaning under allowed chiral redefinitions. The [pole-aligned transformation](alp-pole-alignment-chiral-generator.md) specifies the corresponding ALP-dependent chiral change of fields.

Set $(n_u,n_d,n_s)=(1,0,0)$, $f_a=f_\varphi$, $\theta=0$ and use a CP-conserving vacuum. Work at leading $U(3)$ chiral order and tree level, with no expansion in $f_\pi/f_a$, electromagnetic corrections, higher chiral orders or scalar resonances. Fields, masses, $B_0$, $f_\pi$ and $f_a$ are in GeV; widths and couplings are in GeV. The calculation assumes stable simple neutral poles. The model takes $m_\varphi^2>0$ and $\lambda>0$; a stable vacuum suffices for the local expansion.

The source-coefficient convention below is the one of the draft (section on radial-mode decays): coefficients multiply the dimensionless ratios $\sigma/f_a$ and $a/f_a$ and have mass dimension one. An earlier version of this record used coefficients of the dimensionful fields; the two are related by $\widehat S_{rs}^{\rm here}=f_a^{r+s}\widehat S_{rs}^{\rm earlier}$, and likewise for $\widehat P_{rs}$. The earlier version also wrote both linear $\widehat P$ terms of the Lagrangian with the opposite sign; the two sign changes compensate in $g_{ij}$. The signs here follow from $\langle\chi_+\rangle$ with $\chi=2B_0\widehat M$ and are checked by the general script below.

The general-source version of the polar–Cartesian comparison is in [the general-source record](radial-decay-polar-cartesian-general-sources.md).

## Dressed source and field expansion

In the partonic convention $-\bar q_R\widehat M q_L-\bar q_L\widehat M^\dagger q_R$, $\widehat M=\widehat S+i\widehat P$. Expand the invariant sources in $\sigma/f_a$ and $a/f_a$:
$$
\widehat S(a,\sigma)=\sum_{r,s\ge0}\widehat S_{rs}\left(\frac{\sigma}{f_a}\right)^r\left(\frac{a}{f_a}\right)^s,\qquad
\widehat P(a,\sigma)=\sum_{r,s\ge0}\widehat P_{rs}\left(\frac{\sigma}{f_a}\right)^r\left(\frac{a}{f_a}\right)^s,
$$
so that
$$
\widehat S_{rs}=\frac{f_a^{r+s}}{r!\,s!}\left.\frac{\partial^{r+s}\widehat S}{\partial\sigma^r\partial a^s}\right|_{a=\sigma=0},
$$
and likewise for $\widehat P_{rs}$. Every coefficient is a constant flavour matrix of mass dimension one. In particular $\widehat S_{10}$ is an independent radial coupling, not a derivative of the constant $\widehat S_{00}$.

The paper uses the conjugate source $M_q$ in the $\bar q_Lq_R$ term. Its $P_{rs}$ has the opposite sign to $\widehat P_{rs}$ in the partonic convention used here. In the original connection-free representative, dressing leaves the source unchanged and
$$
\widehat M=\operatorname{diag}\!\left[m_u\left(1+\frac\sigma{f_a}\right)e^{-ia/f_a},m_d,m_s\right].
$$
Direct expansion gives, with $e_1=\operatorname{diag}(1,0,0)$,
$$
\begin{aligned}
\widehat S_{00}&=\operatorname{diag}(m_u,m_d,m_s),&
\widehat S_{10}&=m_ue_1,\\
\widehat P_{01}&=-m_ue_1,&
\widehat P_{11}&=-m_ue_1,\\
\widehat S_{02}&=-\frac{m_u}{2}e_1,&
\widehat S_{12}&=-\frac{m_u}{2}e_1.
\end{aligned}
$$
All other coefficients with $r\le1$ and $s\le2$ vanish. The first index counts $\sigma$, the second counts $a$. The factor $1/2$ in $\widehat S_{02}$ and $\widehat S_{12}$ is required by ordinary Taylor coefficients.

The values $c_q^{AP}=\operatorname{diag}(1/2,0,0)$ and $c_{GA}=0$ describe the angular field in this model; they do not define a scalar invariant. For example replacing the up-quark source by $(m_u+y_u\sigma)e^{-ia/f_a}$ leaves the angular coupling fixed while changing $\widehat S_{10}$ to $y_uf_ae_1$, and hence changes the radial decay. The scalar potential is further independent input.

## Pole equation and ALP-dependent field change

Use $\Phi=\sum_{k=3,8,0}P_k\lambda_k$, $\langle\lambda_k\lambda_l\rangle=2\delta_{kl}$, $u=\exp[i\Phi/(2f_\pi)]$, $U=u^2$ and the meson Lagrangian $\tfrac{f_\pi^2}{4}\langle u_\mu u^\mu+\chi_+\rangle$ with $\chi=2B_0\widehat M$. Then
$$
\langle\chi_+\rangle=4B_0\left\langle\widehat S+\frac{\widehat P\Phi}{f_\pi}-\frac{\widehat S\Phi^2}{2f_\pi^2}\right\rangle+\mathcal O(\Phi^3),
$$
and the quadratic Lagrangian is
$$
\begin{aligned}
\mathcal L_2={}&\frac14\langle\partial_\mu\Phi\partial^\mu\Phi\rangle+\frac12(\partial a)^2
-\frac{B_0}{2}\langle\widehat S_{00}\Phi^2\rangle
+\frac{B_0f_\pi}{f_a}a\langle\widehat P_{01}\Phi\rangle\\
&+\frac{B_0f_\pi^2}{f_a^2}a^2\langle\widehat S_{02}\rangle
-\frac{m_1^2}{12}\langle\Phi\rangle^2.
\end{aligned}
$$
For a momentum with $p^2=s$, variation of $\Phi$ and $a$ gives
$$
\begin{aligned}
\mathbb M_s(\Phi)&=\frac{2B_0f_\pi}{f_a}\widehat P_{01}\,a,\\
\left[s+\frac{2B_0f_\pi^2}{f_a^2}\langle\widehat S_{02}\rangle\right]a
&=-\frac{B_0f_\pi}{f_a}\langle\widehat P_{01}\Phi\rangle,
\end{aligned}
$$
where the operator on Hermitian neutral flavour matrices is
$$
\mathbb M_s(X)=B_0\{\widehat S_{00},X\}+\frac{m_1^2}{3}\langle X\rangle\mathbf1-sX.
$$
For an invertible $\mathbb M_s$, the first equation gives $\Phi=(2B_0f_\pi/f_a)\mathbb M_s^{-1}(\widehat P_{01})\,a$. Substitution into the second equation gives the pole equation for $a\ne0$:
$$
s+\frac{2B_0f_\pi^2}{f_a^2}\left[\langle\widehat S_{02}\rangle+B_0\left\langle\widehat P_{01}\,\mathbb M_s^{-1}(\widehat P_{01})\right\rangle\right]=0.
$$
At a regular simple pole $s=m_i^2$, this is the linear meson shift of the ALP-dependent chiral rotation in the linked pole-alignment record. The shift eliminates the mixed inverse propagator at the selected pole, not at all momenta. For its normalization use the pole-shift coefficients $w_k^{(i)}$:
$$
\sum_k w_k^{(i)}\lambda_k=\frac{2B_0f_\pi}{f_a}\mathbb M_{m_i^2}^{-1}(\widehat P_{01}),\qquad
\frac12\left\langle\left(\sum_k w_k^{(i)}\lambda_k\right)^2\right\rangle=\sum_k\left(w_k^{(i)}\right)^2.
$$
Thus the squared norm of the pole direction with angular component one is
$$
N(m_i^2)=1+\sum_k\left(w_k^{(i)}\right)^2=1+\frac{2B_0^2f_\pi^2}{f_a^2}\left\langle\mathbb M_{m_i^2}^{-1}(\widehat P_{01})^2\right\rangle.
$$

### Basis with no ALP–meson mixing

Select the physical ALP-like pole $s=m_a^2$ and its shift coefficients $w_i^{(a)}$ above. Their squared norm is $N(m_a^2)=1+\sum_i(w_i^{(a)})^2$. The ALP-dependent chiral transformation is
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
-\frac12\sum_{ij}P_i'P_j'\left[\frac{\langle\lambda_i\mathbb M_0(\lambda_j)\rangle}{2}-\frac{m_a^2w_i^{(a)}w_j^{(a)}}{N(m_a^2)}\right]\\
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
The fields $P_j''$ have canonical kinetic terms. The matrix is orthogonal and its final column is the normalized ALP-like pole direction. After the meson-only rotation every column of the complete orthogonal transformation is a normalized pole direction $(w_k^{(r)},1)/\sqrt{N(m_r^2)}$. Each physical field's sign is arbitrary; below its angular component is chosen positive. Keeping the two physical fields of a vertex,
$$
\Phi=\sum_{r=i,j}\frac{P_r}{\sqrt{N(m_r^2)}}\sum_k w_k^{(r)}\lambda_k,
\qquad
a=\sum_{r=i,j}\frac{P_r}{\sqrt{N(m_r^2)}}.
$$
For $i=j$, the displayed expression uses one field $P_i$ rather than adding it twice.

## Polar vertex

Expanding $\tfrac{f_\pi^2}{4}\langle\chi_+\rangle$ and the polar kinetic term to one $\sigma$ and two outgoing fields,
$$
\mathcal L_{\sigma PP}=\frac{\sigma}{f_a}\left[-\frac{B_0}{2}\langle\widehat S_{10}\Phi^2\rangle
+\frac{B_0f_\pi}{f_a}a\langle\widehat P_{11}\Phi\rangle
+\frac{B_0f_\pi^2}{f_a^2}a^2\langle\widehat S_{12}\rangle
+(\partial a)^2\right].
$$
The last term follows from
$$
\frac12\left(1+\frac{\sigma}{f_a}\right)^2(\partial a)^2
=\frac12(\partial a)^2+\frac{\sigma}{f_a}(\partial a)^2+O(\sigma^2).
$$
The kinetic term $\tfrac{f_\pi^2}{4}\langle u_\mu u^\mu\rangle$ does not depend on $\sigma$ in this representative. Write the interaction as $\mathcal L\supset-\frac12\sigma\sum_{ij}g_{ij}P_iP_j$. With outgoing $p_i,p_j$, the kinetic vertex and on-shell relation give
$$
\frac{2p_i\cdot p_j}{f_a}
=\frac{m_\sigma^2-m_i^2-m_j^2}{f_a},\qquad
(p_i+p_j)^2=m_\sigma^2.
$$
Inserting the physical fields,
$$
\begin{aligned}
g_{ij}={}&\frac1{f_a}\Big\{\frac{B_0}{2}\Big\langle\widehat S_{10}\Big\{\sum_kw_k^{(i)}\lambda_k,\sum_lw_l^{(j)}\lambda_l\Big\}\Big\rangle
-\frac{B_0f_\pi}{f_a}\Big\langle\widehat P_{11}\Big(\sum_kw_k^{(i)}\lambda_k+\sum_kw_k^{(j)}\lambda_k\Big)\Big\rangle\\
&-\frac{2B_0f_\pi^2}{f_a^2}\langle\widehat S_{12}\rangle+m_\sigma^2-m_i^2-m_j^2\Big\}\frac{1}{\sqrt{N(m_i^2)\,N(m_j^2)}}.
\end{aligned}
$$
These include all $\sigma a_{\rm phys}a_{\rm phys}$, $\sigma a_{\rm phys}P_i$ and $\sigma P_iP_j$ vertices, with their derivative contribution. The coefficients after the field change depend on the pole selected for the shift, while the on-shell $g_{ij}$ do not. In terms of the sources alone,
$$
\begin{aligned}
g_{ij}={}&\frac1{f_a}\Big\{
\frac{2B_0f_\pi^2}{f_a^2}\Big[B_0^2\left\langle\widehat S_{10}\left\{\mathbb M_{m_i^2}^{-1}(\widehat P_{01}),\mathbb M_{m_j^2}^{-1}(\widehat P_{01})\right\}\right\rangle
-B_0\left\langle\widehat P_{11}\left(\mathbb M_{m_i^2}^{-1}(\widehat P_{01})+\mathbb M_{m_j^2}^{-1}(\widehat P_{01})\right)\right\rangle\\
&-\langle\widehat S_{12}\rangle\Big]
+m_\sigma^2-m_i^2-m_j^2\Big\}
\prod_{r=i,j}\left\{1+\frac{2B_0^2f_\pi^2}{f_a^2}\left\langle\mathbb M_{m_r^2}^{-1}(\widehat P_{01})^2\right\rangle\right\}^{-1/2}.
\end{aligned}
$$
The three terms in square brackets are generated by $\widehat S_{10}$, $\widehat P_{11}$ and $\widehat S_{12}$; the last is the on-shell contribution of $\sigma(\partial a)^2/f_a$. The formula includes the axion-like neutral state among the daughters. It holds for general sources under the assumptions listed in the [general-source record](radial-decay-polar-cartesian-general-sources.md); the first-order angular invariants alone do not fix the scalar couplings.

## Cartesian derivation

Use the coordinates $s,p$ of the [polar and Cartesian record](alp-polar-cartesian-amplitude-equivalence.md), with $f_a=f_\varphi$. The Cartesian radial field $s$ is not the momentum-squared argument of $\mathbb M_s$. In natural units,
$$
f_a\varphi=\frac{f_a+\sigma}{\sqrt2}e^{ia/f_a}
=\frac{f_a+s+ip}{\sqrt2},\qquad
f_a^2|\partial_\mu\varphi|^2=\frac12(\partial s)^2+\frac12(\partial p)^2.
$$
Near the vacuum
$$
\sigma=s+\frac{p^2}{2f_a}-\frac{sp^2}{2f_a^2}+\cdots,\qquad
a=p-\frac{sp}{f_a}+\cdots.
$$
The identity linear transformation implies the same quadratic neutral mixing and the same no-mixing field transformation derived above. The nonlinear terms modify the cubic vertices.

At $U=\mathbf1$, including the vacuum quark source, the potential of the only-up model is
$$
V=\frac{m_\varphi^2}{2}\big[(f_a+s)^2+p^2\big]
+\frac\lambda4\big[(f_a+s)^2+p^2\big]^2
-B_0m_uf_\pi^2\left(1+\frac{s}{f_a}\right)+\text{constant}.
$$
Stationarity and the second and third derivatives give
$$
\begin{aligned}
m_\varphi^2&=\frac{B_0m_uf_\pi^2}{f_a^2}-\lambda f_a^2,&
m_\sigma^2&=\left.\frac{\partial^2V}{\partial s^2}\right|_0=\frac{B_0m_uf_\pi^2}{f_a^2}+2\lambda f_a^2,\\
\left.\frac{\partial^2V}{\partial p^2}\right|_0&=\frac{B_0m_uf_\pi^2}{f_a^2},&
\left.\frac{\partial^3V}{\partial s\,\partial p^2}\right|_0&=2\lambda f_a.
\end{aligned}
$$
The last two concern the unmixed Cartesian field $p$, not the physical ALP mass. In this model the source is linear in the Cartesian fields,
$$
\widehat M=\widehat S_{00}+\frac{s}{f_a}\widehat S_{10}+i\frac{p}{f_a}\widehat P_{01},
$$
with the coefficients above, and
$$
\begin{aligned}
\mathcal L_{sPP}
&=-\frac{B_0s}{2f_a}\langle\widehat S_{10}\Phi^2\rangle-\lambda f_a s p^2\\
&=-\frac{B_0m_u}{f_a}s\left[\frac12\left(\pi^0+\frac{\eta_8}{\sqrt3}+\sqrt{\frac23}\eta_0\right)^2+\pi^+\pi^-+K^+K^-\right]-\lambda f_a s p^2.
\end{aligned}
$$
The source gives no $sp\Phi$ vertex in this model: the polar $\widehat P_{11}$ and $\widehat S_{12}$ vertices are replaced by the scalar-potential vertex. For general sources an $sp\Phi$ vertex proportional to $\widehat P_{11}-\widehat P_{01}$ and further $sp^2$ trace terms appear; see the [general-source record](radial-decay-polar-cartesian-general-sources.md). Here
$$
g_{ij}=\Big\{2\lambda f_a+\frac{2B_0^3f_\pi^2}{f_a^3}\left\langle\widehat S_{10}\left\{\mathbb M_{m_i^2}^{-1}(\widehat P_{01}),\mathbb M_{m_j^2}^{-1}(\widehat P_{01})\right\}\right\rangle\Big\}
\prod_{r=i,j}\left\{1+\frac{2B_0^2f_\pi^2}{f_a^2}\left\langle\mathbb M_{m_r^2}^{-1}(\widehat P_{01})^2\right\rangle\right\}^{-1/2}.
$$

### On-shell equality in the only-up model

Every $\sigma$- or $a$-dependent coefficient is proportional to $e_1$. Since $\widehat S_{00}$ is diagonal, $\mathbb M_{m_r^2}^{-1}(e_1)$ is diagonal and $\langle e_1XY\rangle=\langle e_1X\rangle\langle e_1Y\rangle$ for diagonal $X,Y$. The pole equation at $s=m_r^2$ becomes
$$
m_r^2-\frac{B_0m_uf_\pi^2}{f_a^2}=-\frac{2B_0^2m_u^2f_\pi^2}{f_a^2}\left\langle e_1\,\mathbb M_{m_r^2}^{-1}(e_1)\right\rangle,
$$
and before the common normalization the equality to be shown reads
$$
\frac{B_0m_uf_\pi^2}{f_a^3}\prod_{r=i,j}\left[1-2B_0m_u\left\langle e_1\,\mathbb M_{m_r^2}^{-1}(e_1)\right\rangle\right]+\frac{m_\sigma^2-m_i^2-m_j^2}{f_a}
=2\lambda f_a+\frac{4B_0^3m_u^3f_\pi^2}{f_a^3}\prod_{r=i,j}\left\langle e_1\,\mathbb M_{m_r^2}^{-1}(e_1)\right\rangle.
$$
Inserting $m_\sigma^2$ and $m_i^2,m_j^2$, the terms of zeroth and first order in $\langle e_1\mathbb M_{m_r^2}^{-1}(e_1)\rangle$ cancel on the left and the rest coincides with the right. Hence $g_{ij}^{\rm polar}=g_{ij}^{\rm Cartesian}$ on shell. The comparison does not require dividing by an angular component. The [general-source record](radial-decay-polar-cartesian-general-sources.md) gives the same equality without the $e_1$ factorization, in trace form.

## Charged states, width and limits

At this order
$$
g_{\pi^+\pi^-}=g_{K^+K^-}=\frac{B_0\langle\widehat S_{10}\rangle}{f_a}=\frac{B_0m_u}{f_a},\qquad
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

[The general-source script](../../../scripts/derive_radial_decay_polar_cartesian_general.py) uses the convention of this record. It checks the only-up coefficients above, $\widehat P_{11}-\widehat P_{01}=0$, and the vanishing of the extra Cartesian trace terms. It also compares directly expanded polar and Cartesian Lagrangians with both trace formulas. See its record for residuals.

The older scripts use the earlier dimensionful-field coefficients and the opposite sign of both linear $\widehat P$ terms; their $g_{ij}$ agree after the conversion stated in the inputs. [The operator derivation](../../../scripts/derive_murayama_sigma_operator_amplitude.py) calculates the pole directions, normalization and polar vertex. [The Cartesian derivation](../../../scripts/derive_murayama_sigma_pp.py) checks the polar and Cartesian vertices of the model on shell. Further independent representations and checks are in `scripts/derive_murayama_sigma_explicit_mixing.py`, `scripts/derive_murayama_sigma_normal_basis.py`, `scripts/derive_murayama_sigma_invariant_couplings.py`, `scripts/derive_murayama_sigma_matrix_invariants.py` and `scripts/derive_murayama_sigma_source_coefficients.py`, each with its own unit test. Their determinant and adjugate forms are alternative algebraic representations of the residues, not additional physical assumptions.

```sh
.venv/bin/python -m unittest scripts.test_derive_radial_decay_polar_cartesian_general
.venv/bin/python -m unittest scripts.test_derive_murayama_sigma_operator_amplitude scripts.test_derive_murayama_sigma_pp
```
