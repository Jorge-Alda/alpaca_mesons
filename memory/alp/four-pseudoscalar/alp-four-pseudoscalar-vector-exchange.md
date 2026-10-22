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

The executable record establishes the complete field-counting inventory,
the odd off-shell numerator, the tensor response on a $VPP$ current, the
pairing signs, and the scalar-source vertex polynomials. It does **not** yet
contract all five families into a complete Lorentz-scalar channel amplitude
with every momentum routing, flavor sum, and propagator denominator. In
particular, no four-body width or numerical channel prediction follows from
this record. That assembly is required before describing these ingredients
as the full $a\to4P$ amplitude.
