---
title: Tree-level vector-exchange topologies for an ALP to four pseudoscalars
status: mixed
created: 2026-10-22
updated: 2026-10-22
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and inputs

This record catalogs tree-level diagrams for $a(p_a)\to P_1(p_1)P_2(p_2)
P_3(p_3)P_4(p_4)$ that contain at least one internal vector, at first order
in the ALP sources and at zero electromagnetic coupling. The no-mixing
matrices $c_q^A,c_q^P,c_q^{VS}$ are used as supplied; no further ALP--meson
rotation is made. The pseudoscalar generators obey
$\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$ and vector generators
have unit trace normalization. Natural units, metric $(+---)$, and
$\epsilon^{0123}=+1$ are used.

The retained interactions are the original seven $VJP$ and four $VVP$
odd-parity operators [@Ruiz-Femenia:2003jdx], the leading even-parity
$G_V$ interaction [@Ecker:1988te], and the fitted vector-mass spurion terms
proportional to $e_m^V$ and $\kappa_V$. Electromagnetic interactions, loops,
independent $VPPP$ operators, and operators outside the preliminary-fit
truncation are excluded. The notation and fit-null caveats are in
[the preliminary-fit record](../../rcht/odd-parity/odd-preliminary-fit.md).

The executable derivation is
[`scripts/derive_alp_four_pseudoscalar.py`](../../../scripts/derive_alp_four_pseudoscalar.py).
Run it with

```text
.venv/bin/python scripts/derive_alp_four_pseudoscalar.py
.venv/bin/python scripts/derive_alp_four_pseudoscalar.py --test
```

# Allowed vector-tree families

For five external legs including the incoming ALP, a connected tree made from
vertices of at least three fields has at most three vertices. Enumerating
source vertices, strong vertices, and internal line species gives these
families (the middle two source classes split into distinct parity/source
terms):

| Source and strong vertices | Internal lines | Content |
|---|---|---|
| $aVV+VPP+VPP$ | two vectors | odd $aVV$ source followed by two even $G_V$ decays |
| $aVP+VPP+VVP$ | two vectors | axial $G_V$ source with one odd $VVP$ vertex |
| $aVPP+VPP$ | one vector | vector-source contact and one even $G_V$ vertex |
| $aVV+VPP+VPP$ | two vectors | scalar-source mass variation of the vector quadratic term |
| $aPP+VPP+VPP$ | one pseudoscalar, one vector | leading chiral $aPP$ source and two $G_V$ vertices |

The three $aVV$ vertices are distinct: the odd $VVP$ source, the even
vector-source contribution in the unrotated variables, and the scalar
mass-spurion variation in the rotated variables. The unrotated vector-source
graphs and the rotated scalar-source graphs are alternative representations
of the same $c_q^{VS}$ physics. In particular, do not add the unrotated
$aVPP$ graph to the rotated scalar-source graphs; that would count the same
source twice. The five rows are field/source families, not five independent
fitted operators. All external labeled assignments and all flavor-allowed
vector generators must be summed coherently. Identical-particle factors
belong to phase space, not the amplitude.

# Odd $aVV$ numerator

For vector generators $T_R,T_S$ and virtualities $s_R=q_R^2$ and
$s_S=q_S^2$, define
$$
D_{RS}=\frac12\operatorname{tr}\!\left(\{T_R,T_S\}\mathcal U_a\right),
\qquad
P_{RS}=\frac12\operatorname{tr}\!\left(\{T_R,T_S\}\mathcal Y_a\right),
$$
where
$$
\mathcal U_a=\frac{2c_q^A}{f_a},\qquad
\mathcal Y_a=-\frac{2}{f_a}\{c_q^P,\chi\},\qquad \chi=2B_0m_q.
$$
Including both vector assignments at the odd vertex, the fitted
antisymmetric-tensor operators give
$$
\mathcal N_{RS}(a;s_R,s_S)
=4D_{RS}\left[d_1m_a^2-d_3(m_a^2-s_R-s_S)\right]-16d_2P_{RS}.
$$
The routed $d_4$ contribution vanishes. The executable coefficient
projection also verifies that this numerator annihilates all four original
fit-null directions. This statement is specific to the odd $aVV$ numerator;
it is not a fit to or uncertainty estimate for the complete four-body
amplitude.

# Tensor response and pairing signs

For a pair of outgoing pseudoscalars, the $VPP$ momentum current is
proportional to the wedge $p_i^\mu p_j^\nu-p_j^\mu p_i^\nu$. If its full
flavor and coupling coefficient is denoted by $J_{ij}$, inversion of the
antisymmetric-tensor Hessian gives
$$
V_R^{\mu\nu}(q)=\frac{2J_{ij}}{q^2-M_R^2}
\left(p_i^\mu p_j^\nu-p_j^\mu p_i^\nu\right),\qquad q=p_i+p_j.
$$
The full tensor propagator response is derived by the script from the
antisymmetric-tensor Hessian. It must be used on each $VPP$ current; replacing
it by a Proca pole alone omits the local tensor contribution. For the four
external labeled momenta, the three pair partitions have relative epsilon
signs
$$
(12)(34):+1,\qquad (13)(24):-1,\qquad (14)(23):+1
$$
relative to $\epsilon(p_1,p_2,p_3,p_4)$. These signs follow from determinant
antisymmetry and do not include flavor traces, coupling constants, propagator
denominators, or the Lorentz contraction coefficients.

# Scalar source vertices

The allowed flavor rotation that removes a pure-gradient $c_q^{VS}$ source
shifts $\chi_+$ by
$$
\delta\chi_+=-\frac{2ia}{f_a}[c_q^{VS},\chi].
$$
It yields the $aPP$ coefficient
$$
\frac{i}{4f_a}\operatorname{tr}\!\left(\{\lambda_i^\dagger,
\lambda_j^\dagger\}[c_q^{VS},\chi]\right),
$$
and the $aVV$ mass-term coefficient, with Lorentz indices contracted,
$$
\frac{2ie_m^V}{f_a}\operatorname{tr}\!\left(\{T_R,T_S\}
[c_q^{VS},\chi]\right)
-\frac{i\kappa_V}{2f_a}\left[
\operatorname{tr}(T_R)\operatorname{tr}\!\left(T_S[c_q^{VS},\chi]\right)
+\operatorname{tr}(T_S)\operatorname{tr}\!\left(T_R[c_q^{VS},\chi]\right)
\right].
$$
The invariant $M_V^2$ and $\Delta_V$ terms do not vary. The source rotation
is a change of variables: its induced vertices must be used consistently,
without also adding the removed pure-gradient-source vertices.

# Completion boundary

## Axial and pseudoscalar sources in the no-mixing basis

For this subsection $c_q^A$ and $c_q^P$ mean the input matrices already in
the no-ALP--meson-mixing basis of [the Appendix-A extract](../../extracts/Alda:2026nrd/appendix-a-expansion.md)
[@Alda:2026nrd]. They are not rotated again. The leading chiral source
orders show that neither source supplies an $aPP$ vertex. The axial part of
$u_\mu$ is even in $\Phi$, while the source-free part is odd; their cross
term in $\langle u_\mu u^\mu\rangle$ is therefore odd in $\Phi$. The
pseudoscalar source enters $\chi_+$ through the sine of the anticommutator
map and is also odd in $\Phi$. Consequently the $a\Phi^2$ coefficient is
zero. The possible $aP$ bilinear is removed by the prescribed no-mixing
basis; this cancellation is a basis input and is not re-derived for
arbitrary independent matrices here. The $aPP$--pseudoscalar--$VPP$--$VPP$
family in the topology inventory is thus not an axial or pseudoscalar
contribution: it belongs to the $c_q^{VS}$ scalar-source representation.

The nonzero axial even-parity source with one vector and one pseudoscalar
comes from the cross term of
$$
\mathcal L_{G_V}=\frac{iG_V}{\sqrt2}\langle V_{\mu\nu}u^\mu u^\nu\rangle.
$$
Using $u_\mu^{(a,A)}=2(\partial_\mu a)c_q^A/f_a$ and
$u_\mu^{(P)}=-(\partial_\mu\Phi)/f_\pi+O(\Phi^3)$, antisymmetry of
$V_{\mu\nu}$ combines the two cross terms. For a vector component
$V_{\mu\nu}=T_R V_{R,\mu\nu}$ and meson component
$\Phi=\lambda_i^\dagger P_i$, the result is
$$
\left.\mathcal L_{G_V}\right|_{aVP}=
-\frac{2iG_V}{\sqrt2 f_af_\pi}\,
V_R^{\mu\nu}\partial_\mu a\,\partial_\nu P_i\,
\operatorname{tr}\!\left(T_R[c_q^A,\lambda_i^\dagger]\right).
$$
Thus the coefficient of $V_R^{\mu\nu}\partial_\mu a\partial_\nu P_i$ is
$$
-\frac{2iG_V}{\sqrt2 f_af_\pi}
\operatorname{tr}\!\left(T_R[c_q^A,\lambda_i^\dagger]\right).
$$
The expression is a Lagrangian coefficient before Fourier and vertex
factors. It feeds the tree $aVP$--$VVP$--$VPP$ topology. The odd $aVV$
numerator already recorded above contains both source matrices through
$$
D_{RS}=\frac{1}{f_a}\operatorname{tr}\!\left(\{T_R,T_S\}c_q^A\right),
\qquad
P_{RS}=-\frac{1}{f_a}\operatorname{tr}\!\left(\{T_R,T_S\}\{c_q^P,\chi\}\right),
\quad \chi=2B_0m_q,
$$
in the definitions of that numerator. These supply the $aVV$--$VPP$--$VPP$
tree family. Its momentum routings, like those of the axial $aVP$ family,
must be summed over labeled external assignments and allowed vector flavors.

## Sewn scalar amplitudes

The full tensor responses on the two $VPP$ currents and on the axial $aVP$
current are sewn to the canonical odd tensor amplitude
$\mathcal N_{RS}\epsilon(\varepsilon_R,\varepsilon_S,q_R,q_S)/(M_RM_S)$.
Fourier signs are fixed by an incoming $a$ wave $e^{-ip_a x}$ and outgoing
meson waves $e^{+ip_i x}$. Define $s_{ij}=(p_i+p_j)^2$ and
$\epsilon(p_1,p_2,p_3,p_4)=\epsilon_{\mu\nu\rho\sigma}
p_1^\mu p_2^\nu p_3^\rho p_4^\sigma$. Each vector index $R,S$ is summed
over the real orthonormal vector nonet, with $\operatorname{tr}(T_RT_S)=
\delta_{RS}$. A denominator below means $q^2-M_R^2+i0$; no widths are
inserted.

For the odd $aVV$ source followed by two $VPP$ vertices, the result is
$$
\mathcal M_{aVV}=\epsilon(p_1,p_2,p_3,p_4)
\sum_{(ij|kl)=(12|34),(13|24),(14|23)}\sigma_{ij|kl}
\sum_{R,S}\frac{-G_V^2\,\mathcal N^a_{RS}(s_{ij},s_{kl})
\operatorname{tr}(T_R[\lambda_i^\dagger,\lambda_j^\dagger])
\operatorname{tr}(T_S[\lambda_k^\dagger,\lambda_l^\dagger])}
{2f_\pi^4(s_{ij}-M_R^2+i0)(s_{kl}-M_S^2+i0)},
\qquad (\sigma_{12|34},\sigma_{13|24},\sigma_{14|23})=(+1,-1,+1).
$$
Here $\mathcal N^a_{RS}$ is the previously defined odd numerator with
$D^a_{RS}=\operatorname{tr}(\{T_R,T_S\}c_q^A)/f_a$ and
$P^a_{RS}=-\operatorname{tr}(\{T_R,T_S\}\{c_q^P,\chi\})/f_a$.
The numerator already includes both vector assignments at its odd vertex.

For the axial $aVP$ source followed by an odd $VVP$ vertex and a $VPP$
vertex, let $i,j$ be ordered distinct labels and let $k<l$ be the remaining
two labels. Then
$$
\mathcal M_{aVP}=\epsilon(p_1,p_2,p_3,p_4)
\sum_{\substack{i\ne j\\\{k,l\}=\{1,2,3,4\}\setminus\{i,j\}\\k<l}}
\operatorname{sgn}(i,j,k,l)\sum_{R,S}
\frac{G_V^2\,\mathcal N^j_{RS}((p_a-p_i)^2,s_{kl})
\operatorname{tr}(T_R[c_q^A,\lambda_i^\dagger])
\operatorname{tr}(T_S[\lambda_k^\dagger,\lambda_l^\dagger])}
{f_\pi^3f_a((p_a-p_i)^2-M_R^2+i0)(s_{kl}-M_S^2+i0)}.
$$
The permutation sign is the sign of the ordered list $(i,j,k,l)$ relative
to $(1,2,3,4)$. The strong odd numerator $\mathcal N^j_{RS}$ is the same
VVP polynomial with its source replaced by the outgoing meson $P_j$ and its
source mass squared set to $m_j^2$:
$$
D^j_{RS}=-\frac{1}{2f_\pi}
\operatorname{tr}(\{T_R,T_S\}\lambda_j^\dagger),\qquad
P^j_{RS}=\frac{1}{2f_\pi}
\operatorname{tr}(\{T_R,T_S\}\{\lambda_j^\dagger,\chi\}),
\qquad \chi=2B_0m_q.
$$
The exchanged vector momentum at the odd vertex is crossed,
$q_R=-(p_a-p_i)$, as included in the routing sign above. The complete
axial/pseudoscalar-source amplitude in this retained vector-exchange sector
is $\mathcal M_{aVV}+\mathcal M_{aVP}$. The first term contains both
$c_q^A$ and $c_q^P$; the second contains $c_q^A$ only. The absent $aPP$
family is not included.

The executable calculation in
[`scripts/derive_alp_four_pseudoscalar.py`](../../../scripts/derive_alp_four_pseudoscalar.py)
derives the two scalar kernels, Lorentz contractions, source traces, and
routing signs; its focused checks pass. These are generic labeled-channel
amplitudes. A numerical prediction still requires specifying the four
external generators and masses, evaluating the vector-flavor sums, and
integrating the four-body phase space with the appropriate identical-particle
factor.
