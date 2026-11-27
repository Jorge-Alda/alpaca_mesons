---
title: Quartic vertex for a_phys -> P_i P_j P_k in invariant form, with the on-shell check against the basis without ALP-meson mixing
status: derived
created: 2026-11-19
updated: 2026-11-19
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Quartic vertex in invariant form

**Result.**

1. *Invariant form.* In the original polar fields, with the original source coefficients, the quartic vertex for four legs with an ALP component is a sum of traces of $\widehat S_{00}$, $\widehat P_{01}$, $\widehat S_{02}$, $\widehat P_{03}$, $\widehat S_{04}$ and of $\mathbb M^{-1}_{m_x^2}(\widehat P_{01})$ at the four pole masses, times $\prod_x N(m_x^2)^{-1/2}$ (Step 3). Each leg enters only through its invariant leg matrix $\tfrac{2B_0f_\pi}{f_a}\mathbb M^{-1}_{m_x^2}(\widehat P_{01})/\sqrt{N(m_x^2)}$ and its ALP component $1/\sqrt{N(m_x^2)}$. Charged and $K^0$-type legs enter through their meson generators.
2. *Diagonal PQ-charge sources* (model of DiLuzio:2026poh [@DiLuzio:2026poh]; only-up: $n=(1,0,0)$). The vertex reduces to the neutral form (4a), the closed form (4b) for one charged flavour and the charged-pair form (4c).
3. *Check: independence of the field basis on shell.* The vertex contracted with the four physical poles is the same in the basis without ALP–meson mixing of the draft (`eq:murayama-chiral-shift`); off shell the two vertices differ (section "Check" below). The cubic vertices and the $\sigma P_xP_y$ couplings $g_{xy}(q^2)$ are the same in both bases, also off shell.

## Conventions and units

Conventions of Alda:2026nrd [@Alda:2026nrd] and of the [radial-exchange record](alp-three-pseudoscalar-radial-exchange.md): polar coordinates at $\sigma=0$, $\boldsymbol U=e^{i\boldsymbol\Phi/f_\pi}$, $\boldsymbol\chi=2B_0\widehat{\boldsymbol M}$, $\widehat{\boldsymbol M}(a,0)=\sum_n(\widehat S_{0n}+i\widehat P_{0n})(a/f_a)^n$ (`eq:murayama-source-expansion`), $\mathbb M_s$ of `eq:murayama-operator`, $N(m^2)$ of `eq:murayama-pole-norm`. Vertices are derivatives of the Lagrangian without the factor $i$; all momenta outgoing, the parent carries $-p$; $s_{uv}=(p_u+p_v)^2$. GeV and GeV$^2$; the four-point vertex is dimensionless.

## Assumptions

Those of the [radial-exchange record](alp-three-pseudoscalar-radial-exchange.md): tree level, LO $U(3)$, no expansion in $f_\pi/f_a$, $\widehat P_{00}=0$, the anomalous singlet term $-\tfrac{m_1^2}{12}\langle\boldsymbol\Phi\rangle^2$ independent of $a$ and $\sigma$, regular simple poles with invertible $\mathbb M_{m_x^2}$, no bare ALP mass.

## Inputs

- Pole vectors in the original fields, `eq:murayama-physical-fields`, `eq:murayama-pole-norm`: a pole $x$ with an ALP component is the solution of the linear equations of motion at $p^2=m_x^2$ with meson matrix $\tfrac{2B_0f_\pi}{f_a}\mathbb M^{-1}_{m_x^2}(\widehat P_{01})/\sqrt{N(m_x^2)}$ and ALP component $1/\sqrt{N(m_x^2)}$ (checked against the Hessian eigenvectors, `PoleVectorTest`).
- Chiral shift and linear completion (only for the check): `eq:murayama-chiral-shift`, `eq:murayama-T`.
- Polar Lagrangian and $g_{xy}(q^2)$: [radial-exchange record](alp-three-pseudoscalar-radial-exchange.md); only-up coefficients and pole equation: [only-up radial-decay record](../../radial-mode/couplings-and-decays/finite-mixing/murayama-sigma-two-pseudoscalar-decays.md).

## Step 1: quartic Lagrangian in the original fields

SymPy (noncommuting $\boldsymbol\Phi$, $\partial\boldsymbol\Phi$, $\widehat S_{0n}$, $\widehat P_{0n}$; words compared up to cyclic order) gives, at $\sigma=0$,
$$
\mathcal L_4=\frac{1}{48f_\pi^2}\big\langle[\boldsymbol\Phi,\partial_\mu\boldsymbol\Phi][\boldsymbol\Phi,\partial^\mu\boldsymbol\Phi]\big\rangle
+B_0f_\pi^2\Big[\frac{\langle\widehat S_{00}\boldsymbol\Phi^4\rangle}{24f_\pi^4}-\frac{a\,\langle\widehat P_{01}\boldsymbol\Phi^3\rangle}{6f_af_\pi^3}-\frac{a^2\langle\widehat S_{02}\boldsymbol\Phi^2\rangle}{2f_a^2f_\pi^2}+\frac{a^3\langle\widehat P_{03}\boldsymbol\Phi\rangle}{f_a^3f_\pi}+\frac{a^4\langle\widehat S_{04}\rangle}{f_a^4}\Big].
$$
The ALP has no quartic kinetic term at $\sigma=0$, and the anomalous singlet term is quadratic.

## Step 2: the vertex on general legs

The vertex is the coefficient of $P_aP_iP_jP_k$ after inserting $\boldsymbol\Phi=\sum_xP_x\cdot(\text{meson matrix of }x)$, $a=\sum_xP_x\cdot(\text{ALP component of }x)$ and $\partial_\mu\to ip_{x\mu}$ (SymPy: brute-force fourth derivative of $\mathcal L_4$, `derive_quartic_vertex`).

The kinetic term of $\mathcal L_4$ gives, with outgoing momenta ($-p$ for the parent), a sum over the six pairs $\{x,y\}$ of $-\tfrac{1}{24f_\pi^2}$ times the product of the two momenta, times the traces $\langle[u,x][v,y]\rangle+\langle[v,x][u,y]\rangle$ of the meson matrices of the four legs, with $\{u,v\}$ the complementary pair. The trace combination of the pair $\{a,r\}$ equals that of $\{s,t\}$, and $p_s\!\cdot p_t-p\cdot p_r=s_{st}-\tfrac12\sum_xp_x^2$. The sum of the three trace combinations vanishes identically (SymPy), so only $s_{st}$ remains. This holds for any momenta with $p=p_i+p_j+p_k$.

## Step 3: the vertex on pole legs

Inserting the pole legs of the Inputs (SymPy, `derive_pole_leg_form`), for four legs with an ALP component,
$$
\begin{aligned}
\mathcal V^{(4)}_{aijk}={}&\frac{2B_0f_\pi^2}{3f_a^4}\prod_{x=a,i,j,k}\frac{1}{\sqrt{N(m_x^2)}}\;\Big\{
B_0^4\sum_{(x,y,z,w)}\big\langle\widehat S_{00}\,\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\,\mathbb M^{-1}_{m_y^2}(\widehat P_{01})\,\mathbb M^{-1}_{m_z^2}(\widehat P_{01})\,\mathbb M^{-1}_{m_w^2}(\widehat P_{01})\big\rangle\\
&-2B_0^3\sum_{x}\sum_{(y,z,w)}\big\langle\widehat P_{01}\,\mathbb M^{-1}_{m_y^2}(\widehat P_{01})\,\mathbb M^{-1}_{m_z^2}(\widehat P_{01})\,\mathbb M^{-1}_{m_w^2}(\widehat P_{01})\big\rangle
-6B_0^2\sum_{\{x,y\}}\big\langle\widehat S_{02}\big\{\mathbb M^{-1}_{m_u^2}(\widehat P_{01}),\mathbb M^{-1}_{m_v^2}(\widehat P_{01})\big\}\big\rangle\\
&+18B_0\sum_x\big\langle\widehat P_{03}\,\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\big\rangle+36\,\langle\widehat S_{04}\rangle\\
&-B_0^3\sum_{(r|st)}s_{st}\Big(\big\langle\big[\mathbb M^{-1}_{m_s^2}(\widehat P_{01}),\mathbb M^{-1}_{m_a^2}(\widehat P_{01})\big]\big[\mathbb M^{-1}_{m_t^2}(\widehat P_{01}),\mathbb M^{-1}_{m_r^2}(\widehat P_{01})\big]\big\rangle+(s\leftrightarrow t)\Big)\Big\}.
\end{aligned}
$$
Here $(x,y,z,w)$ runs over the 24 orderings of $(a,i,j,k)$; in the second term $x$ runs over the four legs and $(y,z,w)$ over the six orderings of the other three; in the third $\{x,y\}$ runs over the six pairs of legs and $\{u,v\}$ is the complementary pair; $(r|st)\in\{(i|jk),(j|ik),(k|ij)\}$. The second, third, fourth and fifth terms are those in which one, two, three or four legs enter through their ALP component.

**Charged and $K^0$-type legs.** For a leg $x$ without ALP component, replace $\tfrac{2B_0f_\pi}{f_a}\mathbb M^{-1}_{m_x^2}(\widehat P_{01})$ by its meson generator (for example $\sqrt2|f\rangle\langle f'|$), replace $N(m_x^2)$ by $1$, and drop every term in which $x$ enters through its ALP component (the ALP component is $0$). Checked for $\pi^\pm$, $K^\pm$, $K^0\bar K^0$ legs.

## Step 4: diagonal PQ-charge sources

For $\widehat{\boldsymbol M}=\operatorname{diag}[m_q(1+\sigma/f_a)^{n_q}e^{-in_qa/f_a}]$: $\widehat S_{00}=\operatorname{diag}(m_q)$, $\widehat P_{01}=-\operatorname{diag}(m_qn_q)$, $\widehat S_{02}=-\tfrac12\operatorname{diag}(m_qn_q^2)$, $\widehat P_{03}=\tfrac16\operatorname{diag}(m_qn_q^3)$, $\widehat S_{04}=\tfrac1{24}\operatorname{diag}(m_qn_q^4)$, and $\mathbb M^{-1}_{m_x^2}(\widehat P_{01})$ is diagonal.

**(4a) Four legs with an ALP component** ($a$ and three of $\pi^0,\eta,\eta'$, repetitions allowed). SymPy, symbolic $n_q$:
$$
\mathcal V^{(4)}_{aijk}=\frac{B_0f_\pi^2}{f_a^4}\prod_{x=a,i,j,k}\frac{1}{\sqrt{N(m_x^2)}}\;\sum_{q=u,d,s}m_q\prod_{x=a,i,j,k}\Big[2B_0\big(\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\big)_{qq}+n_q\Big].
$$
The bracket divided by $f_a\sqrt{N(m_x^2)}$ is $\tfrac{2B_0}{f_a}\big(\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\big)_{qq}/\sqrt{N(m_x^2)}+n_q/(f_a\sqrt{N(m_x^2)})$, i.e. $1/f_\pi$ times the $qq$ entry of the leg matrix plus $n_q/f_a$ times the ALP component: the argument of the cosine of flavour $q$ evaluated on pole $x$. For the only-up model the ALP equation of motion gives $2B_0\big(\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\big)_{uu}+1=\dfrac{m_x^2f_a^2}{B_0m_uf_\pi^2}$.

**(4b) One charged flavour.** With the neutral equations of motion ([cancellation record](murayama-distinct-neutral-contact-cancellation.md)), (4a) becomes
$$
\mathcal V^{(4)}_{aijk}=\frac{B_0f_\pi^2}{f_a^4}\Big(\frac{2B_0m_1^2}{3}\Big)^4\prod_{x=a,i,j,k}\frac{\big\langle\mathbb M^{-1}_{m_x^2}(\widehat P_{01})\big\rangle}{\sqrt{N(m_x^2)}}\;\sum_{q=u,d,s}\frac{m_q}{\prod_{x=a,i,j,k}\big[m_x^2-2B_0m_q\big(1+\tfrac{n_q^2f_\pi^2}{2f_a^2}\big)\big]},
$$
which vanishes for four distinct poles.

**(4c) One neutral leg $i$ and a pair of flavour-off-diagonal mesons** $\sqrt2|f\rangle\langle f'|$, $\sqrt2|f'\rangle\langle f|$ as $j,k$ ($\pi^+\pi^-$: $(f,f')=(u,d)$; $K^+K^-$: $(u,s)$; $K^0\bar K^0$: $(d,s)$). SymPy, symbolic $n_q$:
$$
\begin{aligned}
\mathcal V^{(4)}_{aijk}={}&\frac{B_0}{f_a^2\sqrt{N(m_a^2)N(m_i^2)}}\Big\{\frac{2B_0^2}{3}(m_f+m_{f'})\Big[\big(\mathbb M^{-1}_{m_a^2}(\widehat P_{01})\big)_{ff}+\big(\mathbb M^{-1}_{m_a^2}(\widehat P_{01})\big)_{f'f'}\Big]\Big[\big(\mathbb M^{-1}_{m_i^2}(\widehat P_{01})\big)_{ff}+\big(\mathbb M^{-1}_{m_i^2}(\widehat P_{01})\big)_{f'f'}\Big]\\
&+\frac{4B_0^2}{3}\sum_{q=f,f'}m_q\big(\mathbb M^{-1}_{m_a^2}(\widehat P_{01})\big)_{qq}\big(\mathbb M^{-1}_{m_i^2}(\widehat P_{01})\big)_{qq}
+\frac{2B_0}{3}\sum_{(q,q')}m_qn_q\Big[2\big(\mathbb M^{-1}_{m_a^2}(\widehat P_{01})+\mathbb M^{-1}_{m_i^2}(\widehat P_{01})\big)_{qq}+\big(\mathbb M^{-1}_{m_a^2}(\widehat P_{01})+\mathbb M^{-1}_{m_i^2}(\widehat P_{01})\big)_{q'q'}\Big]\\
&+\sum_{q=f,f'}m_qn_q^2\Big\}
+\frac{B_0^2\,(2s_{jk}-s_{ij}-s_{ik})}{3f_a^2\sqrt{N(m_a^2)N(m_i^2)}}\Big[\big(\mathbb M^{-1}_{m_a^2}(\widehat P_{01})\big)_{ff}-\big(\mathbb M^{-1}_{m_a^2}(\widehat P_{01})\big)_{f'f'}\Big]\Big[\big(\mathbb M^{-1}_{m_i^2}(\widehat P_{01})\big)_{ff}-\big(\mathbb M^{-1}_{m_i^2}(\widehat P_{01})\big)_{f'f'}\Big],
\end{aligned}
$$
where $(q,q')$ runs over $(f,f')$ and $(f',f)$. It is symmetric under $j\leftrightarrow k$ ($s_{ij}\leftrightarrow s_{ik}$, $f\leftrightarrow f'$). The last term comes from the kinetic term; $2s_{jk}-s_{ij}-s_{ik}=3s_{jk}-\sum_xm_x^2$ on shell.

## Check: independence of the field basis on shell

The draft also uses the basis without ALP–meson mixing (`eq:murayama-chiral-shift`, `eq:murayama-T`). The vertex above does not need it; this section shows that it gives the same on-shell vertex.

With $X=\tfrac{iaB_0}{f_a}\mathbb M^{-1}_{m_a^2}(\widehat P_{01})$ and $Y=i\boldsymbol\Phi'/f_\pi$, the chiral shift is $\boldsymbol U=e^Xe^{Y}e^X$. SymPy (noncommuting $X$, $Y$) gives
$$
\log\big(e^Xe^Ye^X\big)=2X+Y-\tfrac16\big([X,[X,Y]]+[Y,[X,Y]]\big)+\mathcal O(5),
$$
with vanishing second and fourth orders. Hence
$$
\begin{aligned}
\boldsymbol\Phi={}&\boldsymbol\Phi'+a\,\frac{2B_0f_\pi}{f_a}\mathbb M^{-1}_{m_a^2}(\widehat P_{01})
+\frac{a^2B_0^2}{6f_a^2}\Big[\mathbb M^{-1}_{m_a^2}(\widehat P_{01}),\Big[\mathbb M^{-1}_{m_a^2}(\widehat P_{01}),\boldsymbol\Phi'\Big]\Big]\\
&+\frac{aB_0}{6f_\pi f_a}\Big[\boldsymbol\Phi',\Big[\mathbb M^{-1}_{m_a^2}(\widehat P_{01}),\boldsymbol\Phi'\Big]\Big]+\mathcal O(5),
\end{aligned}
$$
(cubic term checked numerically: residual of order five), and $\langle\boldsymbol\Phi\rangle=\langle\boldsymbol\Phi'\rangle+a\tfrac{2B_0f_\pi}{f_a}\langle\mathbb M^{-1}_{m_a^2}(\widehat P_{01})\rangle$ exactly. The further steps of `eq:murayama-T` are linear.

The change of fields is therefore a linear part plus a cubic part. Inserting it in the Lagrangian, the cubic terms in the new fields are the old cubic terms with the linear part inserted. The quartic terms are the old quartic terms with the linear part inserted, plus the variation of the quadratic Lagrangian (the linear equations of motion) multiplied by the cubic part of the change of fields. In the four-point vertex this extra term places one leg $x$ in its linear equation of motion at $p_x^2$, which vanishes when the leg is an on-shell pole. Hence:

- on shell, $\mathcal V^{(4)}$ in the basis without ALP–meson mixing equals $\mathcal V^{(4)}$ of the original fields of Step 3;
- off shell they differ when the meson matrices of the legs and $\mathbb M^{-1}_{m_a^2}(\widehat P_{01})$ do not commute (for commuting matrices the cubic part of the change of fields vanishes);
- the cubic vertices are related by the linear part alone, so the selection rule of the [parity record](../murayama-cubic-vertex-parity-rule.md) holds in both bases;
- $\sigma$ is not transformed, so the $\sigma P_xP_y$ vertex, and hence $g_{xy}(q^2)$ of the [radial-exchange record](alp-three-pseudoscalar-radial-exchange.md), is the same in both bases for any $q^2$ and any leg momenta.

## Verification

[Script](../../../../scripts/derive_murayama_only_up_three_pseudoscalar_amplitude.py) (`derive_symmetric_chiral_map`, `derive_quartic_lagrangian`, `derive_quartic_vertex`, `derive_kinetic_channel_identity`, `derive_pole_leg_form`, `derive_diagonal_forms`, `quartic_vertex_from_legs`, `original_basis_direct_quartic`, `unmixed_basis_direct_quartic`, `alp_shift_matrix`, `inverse_mass_operator_of_p01`) and [unit tests](../../../../scripts/test_derive_murayama_only_up_three_pseudoscalar_amplitude.py) (`SymbolicTest`, `PoleVectorTest`, `QuarticVertexBasisTest`, `OnlyUpAmplitudeTest`, `DraftQuarticVertexTest`, `DraftDiagonalQuarticTest`). "Direct" expands the exact eleven-field polar Lagrangian to quartic order with all derivative terms; the basis without mixing is the same Lagrangian with $\boldsymbol U=\boldsymbol g_R\boldsymbol U'\boldsymbol g_R$, $\boldsymbol g_R=\exp\big(\tfrac{iaB_0}{f_a}\mathbb M^{-1}_{m_a^2}(\widehat P_{01})\big)$ (internal cross-check only). Random non-commuting sources: seed 2026 of the radial-exchange script. Parameters otherwise as in the [parity record](../murayama-cubic-vertex-parity-rule.md); Dalitz point $s_{12}=0.3$, $s_{13}=0.5$ GeV$^2$ (algebraic check, outside the physical region); off-shell momenta: four random momenta with $p=p_i+p_j+p_k$ (`off_shell_gram`). `Draft*` tests implement the formulas of Steps 3 and 4 literally, with $\mathbb M^{-1}_{m_x^2}(\widehat P_{01})$ and $N(m_x^2)$ only.

| check | residual |
|---|---|
| SymPy: even orders of $\log(e^Xe^Ye^X)$; third order as displayed | $0$; exact |
| SymPy: derived $\mathcal L_4$ vs Step 1; brute-force fourth derivative vs vertex on general legs; sum of the three kinetic trace combinations; pole-leg form of Step 3 | exact |
| SymPy: (4a) and (4c) for pairs $(u,d)$, $(d,s)$ from the general vertex | exact |
| Hessian pole vectors vs $\big(\tfrac{2B_0f_\pi}{f_a}\mathbb M^{-1}_{m^2}(\widehat P_{01}),1\big)/\sqrt N$ (only-up; seed 2026) | $<10^{-10}$ |
| on shell: vertex without mixing vs original fields, relative (seed 2026, four daughter sets; charge models $n=(1,0,0)$, $(2,0,1)$, seven final states each) | $\le1.4\times10^{-13}$; $\le2.9\times10^{-14}$ |
| off shell: same comparison, relative (seed 2026; charge models with charged legs) | $0.29$ to $1.15$; up to $2.0$ (all-neutral diagonal legs: $9\times10^{-15}$, since the cubic part of the change of fields vanishes) |
| Step 3 formula vs direct original-field vertex, on and off shell (seed 2026; charge models) | $\le4.3\times10^{-11}$ (precision of $\mathbb M^{-1}_{m^2}$); $\le1.6\times10^{-13}$ |
| Step 3 written literally (`draft_quartic_vertex`, incl. the charged-leg rule) vs direct | $\le4\times10^{-11}$; $\le1.6\times10^{-13}$ |
| (4a), (4c) literally vs direct amplitude minus $\sigma$ exchange, $n=(1,0,0)$, $(2,0,1)$, $(1,1,0)$, $(-1,0,\tfrac12)$, eight final states (relative to $|\mathcal M|$) | $\le3.9\times10^{-12}$; $\le5.8\times10^{-14}$ |
| (4b) vs direct contact, only-up $\pi^0\pi^0\eta$, $3\pi^0$, $\pi^0\eta\eta'$ | $<10^{-10}$ relative to $|\mathcal M|$ |
| $2B_0(\mathbb M^{-1}_{m_x^2}(\widehat P_{01}))_{uu}+1=m_x^2f_a^2/(B_0m_uf_\pi^2)$, only-up neutral poles | $<10^{-9}$ relative |

```sh
.venv/bin/python scripts/derive_murayama_only_up_three_pseudoscalar_amplitude.py
.venv/bin/python -m unittest scripts.test_derive_murayama_only_up_three_pseudoscalar_amplitude
```
