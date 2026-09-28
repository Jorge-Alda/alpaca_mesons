---
title: Vector resonances in the odd-intrinsic-parity sector
status: literature
created: 2026-09-23
updated: 2026-09-28
bibtex_keys: [Ruiz-Femenia:2003jdx]
provenance: literature
---

# Vector Resonances in the Odd-Intrinsic-Parity Sector

The auditable source is `memory/sources/Ruiz-Femenia:2003jdx/arxiv-source.tar`; the
equation locations below refer to `vvp.tex` in that archive.

The paper formulates the light vector octet as antisymmetric tensors $V_{\mu\nu}$ transforming covariantly under the unbroken chiral group. Its leading even-parity interactions are

$$
\mathcal L_2(V)=\frac{F_V}{2\sqrt{2}}\langle V_{\mu\nu}f_+^{\mu\nu}\rangle
+\frac{iG_V}{\sqrt{2}}\langle V_{\mu\nu}u^\mu u^\nu\rangle.
$$

The odd-parity action contains independent $VJP$ and $VVP$ operator bases in addition to the Wess--Zumino action. The paper uses $\epsilon_{0123}=+1$ (source `vvp.tex`, lines 352--354), and defines [@Ruiz-Femenia:2003jdx, Eq. (10)]:

$$
\begin{aligned}
\mathcal O^1_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\nu},f_+^{\rho\alpha}\}\nabla_\alpha u^\sigma\rangle,\\
\mathcal O^2_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\alpha},f_+^{\rho\sigma}\}\nabla_\alpha u^\nu\rangle,\\
\mathcal O^3_{VJP}&=i\epsilon_{\mu\nu\rho\sigma}\langle\{V^{\mu\nu},f_+^{\rho\sigma}\}\chi_-\rangle,\\
\mathcal O^4_{VJP}&=i\epsilon_{\mu\nu\rho\sigma}\langle V^{\mu\nu}[f_-^{\rho\sigma},\chi_+]\rangle,\\
\mathcal O^5_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla_\alpha V^{\mu\nu},f_+^{\rho\alpha}\}u^\sigma\rangle,\\
\mathcal O^6_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla_\alpha V^{\mu\alpha},f_+^{\rho\sigma}\}u^\nu\rangle,\\
\mathcal O^7_{VJP}&=\epsilon_{\mu\nu\rho\sigma}\langle\{\nabla^\sigma V^{\mu\nu},f_+^{\rho\alpha}\}u_\alpha\rangle.
\end{aligned}
$$

Their normalization is [@Ruiz-Femenia:2003jdx, Eq. (12)]:

$$
\mathcal L_{VJP}=\sum_{a=1}^{7}\frac{c_a}{M_V}\mathcal O^a_{VJP},\qquad
\mathcal L_{VVP}=\sum_{a=1}^{4}d_a\mathcal O^a_{VVP}.
$$

Here $M_V$ is inserted to make the $c_a$ dimensionless. The four $VVP$
operators are given in the companion extract below and in [@Ruiz-Femenia:2003jdx,
Eq. (11)].

The building blocks use

$$
u_\mu=i\left[u^\dagger(\partial_\mu-ir_\mu)u-u(\partial_\mu-i\ell_\mu)u^\dagger\right],
\quad \chi_\pm=u^\dagger\chi u^\dagger\pm u\chi^\dagger u,
\quad \chi=2B_0(s+ip),
$$

and

$$
f_\pm^{\mu\nu}=uF_L^{\mu\nu}u^\dagger\pm u^\dagger F_R^{\mu\nu}u
$$

[@Ruiz-Femenia:2003jdx, source `vvp.tex`, lines 198--205 and 286--289].

## Direct $aV\gamma$ classification (derived source insertion)

This paragraph is a repository derivation, not a formula asserted in the paper.
At the Goldstone vacuum, specialize the electromagnetic source to
$f_+^{\mu\nu}=-2eQF^{\mu\nu}$ and $f_-^{\mu\nu}=0$, and insert the physical-ALP
relations already derived in `two-source-alp-rcht-odd-parity-extension.md`,

$$
u_\mu\big|_{a_{\rm physical}}=-2\frac{\partial_\mu a}{f_a}X,
\qquad
\chi_-\big|_{a_{\rm physical}}=-4iB_0\frac{a}{f_a}
\{c_q^{AP}+X,m_q\}.
$$

Therefore the published VJP basis permits a direct $aV\gamma$ vertex from
$\mathcal O^1$, $\mathcal O^2$, $\mathcal O^3$, $\mathcal O^5$, $\mathcal O^6$,
and $\mathcal O^7$, subject to a nonzero flavor trace for the chosen vector and
ALP sources. $\mathcal O^4$ vanishes under $f_-=0$. The first, second, fifth,
sixth, and seventh structures are derivative-ALP terms through $u_\mu$; the
third is the mass/source insertion through $\chi_-$. This conclusion is before
integration by parts, on-shell reduction, vector-photon diagonalization, or
short-distance substitutions.

The short-distance constraints imply $F_VG_V=F^2$ and $2F_VG_V=F_V^2$, hence $F_V=\sqrt{2}F$ and $G_V=F/\sqrt{2}$ in the single-resonance approximation. The authors report that these values agree well with phenomenology. [@Ruiz-Femenia:2003jdx, Eqs. (15)--(16)]

For the odd-parity couplings, matching the single-resonance $\langle VVP\rangle$
function to the QCD OPE gives [@Ruiz-Femenia:2003jdx, Eq. (19)]:

$$
\begin{aligned}
4c_3+c_1&=0, & c_1-c_2+c_5&=0,\\
c_5-c_6&=\frac{N_C}{64\pi^2}\frac{M_V}{\sqrt2F_V},\\
d_1+8d_2&=-\frac{N_C}{64\pi^2}\frac{M_V^2}{F_V^2}+\frac{F^2}{4F_V^2},\\
d_3&=-\frac{N_C}{64\pi^2}\frac{M_V^2}{F_V^2}+\frac{F^2}{8F_V^2}.
\end{aligned}
$$

These constraints apply beyond the chiral limit according to the authors, because
the effective couplings do not depend on Goldstone masses. They do **not** fix
all $c_i,d_i$ individually: $c_4,c_7,d_4$ do not enter this matching, while the
equations fix only combinations of the remaining coefficients (notably
$d_1+8d_2$). Phenomenological channel calculations test combinations of these
couplings but are not a complete determination of the unfixed parameters.

For $\omega\to\pi\gamma$, the constrained theory gives $0.524$--$0.703\,\mathrm{MeV}$ as $F$ is varied from $F_\pi$ to its chiral-limit value, compared in the paper with $0.734\pm0.035\,\mathrm{MeV}$. The discrepancy is $5$--$30\%$, consistent with the expected size of omitted subleading $1/N_c$ effects. [@Ruiz-Femenia:2003jdx, Sec. 4.1, Eq. (29)]

Strict vector-meson dominance is not supported in the same truncation: vector exchange alone gives about $1.4\,\mathrm{MeV}$ for $\omega\to3\pi$, versus the quoted $7.52\pm0.06\,\mathrm{MeV}$, requiring a sizeable direct contribution. [@Ruiz-Femenia:2003jdx, Sec. 4.2, Eqs. (30)--(32)]

The paper explicitly warns against adding local $O(p^6)$ odd-parity terms already saturated by the retained vectors: doing so double counts and can spoil the short-distance behavior unless the resonance operator basis is enlarged consistently. [@Ruiz-Femenia:2003jdx, end of Sec. 3]
