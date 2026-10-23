---
title: Pion vector form factor and ALP pion production partial wave
status: mixed
created: 2026-10-24
updated: 2026-10-24
bibtex_keys: [Messerli:2025rnv, Alda:2026nrd, GomezDumm:2013sib, Watson:1952ji, Omnes:1958hv]
provenance: analytical-derivation
---

# Scope and inputs

This record identifies the two hadronic amplitudes needed in the low-energy
isovector spectral density of the ALP transition form factor. The spectral
normalization and scalar TFF convention follow the
[spectral-construction record](alp-isovector-tff-spectral-construction.md)
[@Messerli:2025rnv]. ALP sources are in Alda:2026nrd's no-mixing convention
[@Alda:2026nrd]. The elastic pion vector form factor is represented
dispersively, with its phase and continuation treated as inputs
[@GomezDumm:2013sib; @Watson:1952ji; @Omnes:1958hv].

# Charged-pion vector form factor

Define the form factor through the electromagnetic-current matrix element
$$
\langle \pi^+(p_+)\pi^-(p_-)|J^\mu_{\rm em}(0)|0\rangle
=(p_+-p_-)^\mu F_\pi^V(s),\qquad s=(p_++p_-)^2.
$$
Current conservation fixes the vector structure, and charge normalization
gives $F_\pi^V(0)=1$. In the elastic isovector region, Watson's theorem fixes
$$
\arg F_\pi^V(s+i0)=\delta_1^1(s),
$$
where $\delta_1^1$ is the isospin-one, angular-momentum-one $\pi\pi$ phase
shift [@Watson:1952ji]. A once-subtracted Omnes representation is
$$
F_\pi^V(s)=P(s)\,\Omega_1^1(s),\qquad
\Omega_1^1(s)=\exp\!\left[\frac{s}{\pi}\int_{4M_\pi^2}^{\infty}
\frac{\delta_1^1(x)}{x(x-s-i0)}\,dx\right],
$$
with $P(0)=1$ and no elastic right-hand-cut phase in $P$. Further subtraction
constants, inelastic effects, and the phase continuation above the elastic
domain are model inputs, not fixed by charge normalization
[@Omnes:1958hv; @GomezDumm:2013sib]. A tree-level rho pole is an alternative
approximation; it must not be combined with the dispersive rho strength
without a matching prescription.

# ALP production P wave

Start from the crossed current amplitude
$$
a(p_a)+\gamma^*(q)\longrightarrow\pi^+(p_+)\pi^-(p_-),
$$
with $s=(p_++p_-)^2$, $q^2=y$, and $z$ the pion scattering angle in the
dipion rest frame. Its transverse odd-parity part may be parameterized as
$$
\mathcal M^\mu=i\epsilon^{\mu\nu\rho\sigma}
p_{+\nu}p_{-\rho}q_\sigma\,\mathcal A_a(s,y,z),
$$
with the overall phase fixed by the selected S-matrix convention. Project
the coefficient onto angular momentum one,
$$
\mathcal F_{a\pi\pi\gamma^*}(s,y)=\mathcal N_P
\int_{-1}^{1} dz\,w_P(z)\,\mathcal A_a(s,y,z).
$$
The weight $w_P$ and normalization $\mathcal N_P$ must be derived in a
declared helicity or partial-wave convention and matched to
Messerli:2025rnv's $\mathcal F_{P\pi\pi\gamma^*}$ normalization by deriving
the two-pion unitarity relation. They are not fixed by an integrated
real-photon width.

The following is the tree odd-resonance contribution from the retained
seven-VJP/four-VVP operators and leading $G_V$ $VPP$ vertex. It is not the
complete amplitude: independent local anomalous terms, loops, rescattering,
and higher operators are not included. Use $k=-q$ as the outgoing-current
momentum in the crossed decay routing, $k^2=y$, and
$$
s=(p_++p_-)^2,\qquad t=(k+p_-)^2,\qquad u=(k+p_+)^2,
\qquad p_a=k+p_++p_-.
$$
Define $\mathcal A^{\rm odd}_{a\pi\pi\gamma^*}$ by
$$
\mathcal M^{\rm odd}=\epsilon_{\mu\nu\rho\sigma}
\varepsilon^{*\mu}p_+^\nu p_-^\rho k^\sigma\,
\mathcal A^{\rm odd}_{a\pi\pi\gamma^*}(s,t,u;y).
$$
For the charged pion pair and strong-isospin limit, the exchanged vector
propagators reduce to the neutral rho pole. The virtual photon response
remains the full neutral-vector matrix
$$
W_\gamma(y)=\sum_i\frac{R_iT_i}{M_i^2-y-i0}.
$$
Use the no-mixing source matrices and pion generators
$$
\mathcal U_a=\frac{2c_q^A}{f_a},\qquad
\mathcal Y_a=-\frac{2}{f_a}\{c_q^P,\chi\},\qquad
\mathcal U_{\pi^\pm}=-\frac{\lambda_\pm^\dagger}{f_\pi},\qquad
\mathcal Y_{\pi^\pm}=\frac{\{\lambda_\pm^\dagger,\chi\}}{f_\pi},
\qquad \chi=2B_0m_q.
$$
For any flavor matrix $C$ and source matrices $\mathcal U,\mathcal Y$, define
the trace contractions
$$
\begin{aligned}
D_W(C,\mathcal U;y)&=\frac12\operatorname{tr}
  (\{C,W_\gamma(y)\}\mathcal U),&
P_W(C,\mathcal Y;y)&=\frac12\operatorname{tr}
  (\{C,W_\gamma(y)\}\mathcal Y),\\
D_Q(C,\mathcal U)&=\frac12\operatorname{tr}
  (\{C,Q\}\mathcal U),&
P_Q(C,\mathcal Y)&=\frac12\operatorname{tr}
  (\{C,Q\}\mathcal Y).
\end{aligned}
$$
The retained amplitude from the original seven VJP and four VVP operators
and leading $G_V$ interaction follows the kernel decomposition in the
[odd $a\to\gamma PP$ derivation](../photon-two-pseudoscalar/alp-photon-pp-amplitude.md).
For the complete orthonormal Hermitian $U(3)$ vector basis,
$$
\sum_R\operatorname{tr}(T_RA)\operatorname{tr}(T_RB)
=\operatorname{tr}(AB),\qquad \operatorname{tr}(T_RT_S)=\delta_{RS}.
$$
Applying this completeness relation to the flavor sums gives the matrix
contractions above. For channel label $j=s,t,u$, set
$$
\begin{array}{c|c|c|c|c|c}
j&p_j&r_j&m_j^2&C_j&(\mathcal U_j,\mathcal Y_j)\\ \hline
s&-\dfrac{G_V}{\sqrt2 f_\pi^2}&s&m_a^2&
[\lambda_+^\dagger,\lambda_-^\dagger]&(\mathcal U_a,\mathcal Y_a)\\[1mm]
t&\dfrac{\sqrt2G_V}{f_\pi f_a}&t&m_\pi^2&
[\lambda_+^\dagger,c_q^A]&(\mathcal U_{\pi^-},\mathcal Y_{\pi^-})\\[1mm]
u&-\dfrac{\sqrt2G_V}{f_\pi f_a}&u&m_\pi^2&
[\lambda_-^\dagger,c_q^A]&(\mathcal U_{\pi^+},\mathcal Y_{\pi^+})
\end{array}
$$
and define each matrix numerator by
$$
\begin{aligned}
N_j(r_j,y)=p_j\Bigg\{&4\sqrt2e\Big(
  [(d_1-d_3)m_j^2+d_3(r_j+y)]D_W(C_j,\mathcal U_j;y)
  -4d_2P_W(C_j,\mathcal Y_j;y)\Big)\\
&-\frac{4e}{M_V}\Big(
  [(c_1+c_2-c_5)m_j^2+(-c_1+c_2+c_5-2c_6)r_j
  +(c_1-c_2+c_5)y]D_Q(C_j,\mathcal U_j)
  -4c_3P_Q(C_j,\mathcal Y_j)\Big)\Bigg\}.
\end{aligned}
$$
Then
$$
\boxed{\mathcal A^{\rm odd}=
\frac{N_s(s,y)}{s-M_\rho^2+i0}
+\frac{N_t(t,y)}{t-M_\rho^2+i0}
+\frac{N_u(u,y)}{u-M_\rho^2+i0}.}
$$
The matrix form keeps the virtual photon response as the full neutral-vector
matrix $W_\gamma(y)$ and does not introduce its diagonal entries. Under the
stated strong-isospin and diagonal-source assumptions, evaluation of the
traces reproduces $N_u(u,y)=N_t(u,y)$ and the component expressions in the
executable derivation. The complete normalized vector nonet is summed in
the trace contractions, while the propagators connecting the odd source to
the pion current reduce to the neutral rho pole. The general source-kernel
construction is given in the [odd $a\to\gamma PP$ derivation](../photon-two-pseudoscalar/alp-photon-pp-amplitude.md).
The charged-pion $c_4$ trace vanishes for $m_u=m_d$ because it is
proportional to $\chi_d-\chi_u$. The equations are the output of the
executable symbolic derivation
[`derive_alp_virtual_photon_pion_amplitude.py`](../../../scripts/derive_alp_virtual_photon_pion_amplitude.py).

For the partial-wave convention used in the existing fit-null audit, expand
the scalar coefficient in derivatives of Legendre polynomials and set
$$
f_1(s,y)=\frac34\int_{-1}^1dz\,(1-z^2)
\mathcal A^{\rm odd}_{a\pi\pi\gamma^*}(s,t(z),u(z);y).
$$
With $\lambda(s,m_a^2,y)=s^2+m_a^4+y^2-2sm_a^2-2sy-2m_a^2y$,
$$
\begin{aligned}
t(z)&=M_\pi^2+\frac{m_a^2-s+y}{2}
-\frac{\sqrt{\lambda(s,m_a^2,y)}}{2}\,\sigma_\pi(s)z,\\
u(z)&=M_\pi^2+\frac{m_a^2-s+y}{2}
+\frac{\sqrt{\lambda(s,m_a^2,y)}}{2}\,\sigma_\pi(s)z.
\end{aligned}
$$
The $s$-channel term is independent of $z$ and projects unchanged. The
crossed numerator $N_t(r,y)$ is affine in $r$, so its derivative with
respect to $r$ is constant. Define the crossed-pole center and half-span
$$
A_\rho=M_\pi^2+\frac{m_a^2-s+y}{2}-M_\rho^2,\qquad
B_\pi=\frac{\sigma_\pi(s)\sqrt{\lambda(s,m_a^2,y)}}{2}.
$$
Then
$$
I(A_\rho,B_\pi)=\frac34\int_{-1}^{1}
\frac{1-z^2}{A_\rho-B_\pi z+i0}\,dz
=\frac{3\left[2A_\rho B_\pi+(A_\rho^2-B_\pi^2)
\{\log(A_\rho-B_\pi+i0)-\log(A_\rho+B_\pi+i0)\}\right]}
{4B_\pi^3}.
$$
The projected retained tree amplitude is therefore
$$
\boxed{f_1(s,y)=\frac{N_s(s,y)}{s-M_\rho^2+i0}
+2\frac{\partial N_t(r,y)}{\partial r}
+2N_t(M_\rho^2,y)I(A_\rho,B_\pi).}
$$
The continuous $B_\pi\to0$ limit is
$I(A_\rho,B_\pi)=1/(A_\rho+i0)$. Thus the retained tree
P wave is explicit once the source coefficients, vector residues, and pole
prescription are selected. The sign and common phase must still be matched
to Messerli:2025rnv's partial-wave definition before using its $192\pi$
spectral prefactor. The calculated $f_1$ includes the electromagnetic factor
$e$ and is the scalar in the crossed decay tensor above. Messerli:2025rnv's
current partial wave has an explicit current phase and its own electromagnetic
normalization; conversion between that quantity and $f_1$ is not yet
established here.

In a strictly elastic approximation, if the production kernel has no
additional right-hand-cut phase, Watson's theorem implies
$$
\mathcal F_{a\pi\pi\gamma^*}(s,y+i0)
=e^{i\delta_1^1(s)}\,\widehat{\mathcal F}_a(s,y),
$$
where $\widehat{\mathcal F}_a$ has no elastic phase. Left-hand cuts and
inelastic channels require a more general representation.

# Use in the spectral density

After matching the conversion from $f_1$ to Messerli:2025rnv's current
partial-wave convention, the density has the form
$$
\rho_a^{I=1}(s,y)=\frac{s\,\sigma_\pi^3(s)}{192\pi}
\operatorname{Im}\!\left[(F_\pi^V(s))^*
\mathcal F_{a\pi\pi\gamma^*}(s,y)\right],\qquad
\sigma_\pi(s)=\sqrt{1-\frac{4M_\pi^2}{s}}.
$$
If both amplitudes carry precisely the same elastic phase, those phases
cancel in their product. The remaining imaginary part depends on their
relative phase and on the virtuality dependence; it is not automatically
zero or nonzero. In particular, a real production kernel dressed by exactly
the same elastic phase as $F_\pi^V$ gives a real product in this strict
elastic limit. If a proposed density is nonzero in that setup, first check
the discontinuity convention, the partial-wave normalization, and any
additional phases before assigning spectral strength.

This construction supplies only the component supported by the stated
two-pion cut. Subtraction terms, single-spectral contributions, higher
energies, other open cuts, and matching that prevents double counting of rho
dynamics remain separate inputs [@Messerli:2025rnv].

# Assumptions and verification status

- Isospin symmetry and elastic two-pion unitarity are assumed in the phase
  statements; inelasticity and electromagnetic corrections are omitted.
- The Omnes expression fixes the phase representation, not a numerical form
  factor. A phase input, subtraction prescription, and high-energy
  continuation are needed for evaluation.
- The retained odd-resonance tree amplitude is source-resolved but the
  preliminary odd fit contains an unmeasured coefficient direction that
  affects the virtual-photon P wave. A numerical prediction must specify a
  representative or constrain that direction independently.
- The local anomalous current amplitude and rescattering completion have
  not been supplied. Also, the normalization and phase of $f_1$ must be
  matched directly to Messerli:2025rnv before inserting the quoted spectral
  prefactor.

# Reproducibility

Run `.venv/bin/python scripts/derive_alp_virtual_photon_pion_amplitude.py`
to print only the abstract $U(3)$ trace numerators, matrix amplitude,
crossed quotient and residue, and analytic $P$ wave. The component
calculation is retained internally for comparison, not displayed. Six
focused checks compare the trace contractions and amplitude to the component
result, compare the $y=0$ result to the existing real-photon amplitude,
reproduce the earlier fit-null projections, check the charged-pion
flavor/isospin structure, and compare the logarithmic angular integral to
direct quadrature. The checks were run with
`.venv/bin/python scripts/derive_alp_virtual_photon_pion_amplitude.py --test`
and passed. These checks validate the retained tree construction, not the
missing anomalous term, rescattering, or its normalization match to
Messerli:2025rnv.

# Inputs and provenance

- The scalar TFF and two-pion spectral normalization come from the
  [isovector spectral construction](alp-isovector-tff-spectral-construction.md)
  [@Messerli:2025rnv].
- The ALP source basis and virtual-current vertex conventions come from the
  [source-resolved doubly virtual TFF record](alp-doubly-virtual-tff.md)
  and the [odd photon plus two-pseudoscalar amplitude](../photon-two-pseudoscalar/alp-photon-pp-amplitude.md)
  [@Alda:2026nrd; @Ecker:1988te; @Ecker:1989yg; @Ruiz-Femenia:2003jdx].
- The pion phase and dispersive form-factor representation follow
  [@Watson:1952ji; @Omnes:1958hv; @GomezDumm:2013sib].
