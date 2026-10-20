---
title: Real-photon ALP resonance contribution from vector mixing
status: derived
created: 2026-10-20
updated: 2026-10-20
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Bauer:2020jbp, Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# Scope and inputs

This record derives the original-eleven-operator tree contribution from vector
photon conversion to the real-photon ALP coefficient. It uses the
[Bauer-normalized coefficient](alp-photon-invariant-from-mixing.md),
the [no-mixing two-vector kernels](../../rcht/odd-parity/alp-vv-amplitude.md),
and the [quadratic photon-vector source convention](../../rcht/foundations/rcht-photon-vector-diagonalization.md).
The selected even-sector and odd-sector inputs are the
[central negative-strange single-trace provisional fit](../../rcht/odd-parity/odd-preliminary-fit.md).
The fit's nuisance parameters and canonical odd-coefficient representative
are used as stored; no uncertainty is assigned along the four coefficient
null directions.

The user supplies dimensionless diagonal $c_q^A,c_q^P$ in the already
no-mixing basis, together with $m_a,f_a$ in GeV. No ALP field rotation or
additional meson-pole/WZW term is applied. Set
$$
U_a=\frac{2c_q^A}{f_a},\qquad
Y_a=-\frac{2}{f_a}\{c_q^P,\chi\},\qquad
\chi=\operatorname{diag}(2B_0m_u,2B_0m_d,2B_0m_s),
$$
and $Q=\operatorname{diag}(2/3,-1/3,-1/3)$. The fitted orthonormal neutral
generators $T_i$ obey $\operatorname{tr}(T_iT_j)=\delta_{ij}$, and the signed
one-trace residues $R_i=\operatorname{tr}(T_iS_\gamma)$ include the fitted
$F_V,z_V$ source. Define the real-photon response
$$
W_i(0)=\frac{R_i}{M_i^2}.
$$

# Analytic result

Use $q_1^2=q_2^2=0$ and $q_1\cdot q_2=m_a^2/2$. The ordered VJP scalar for a
vector leg $T_i$ and an explicit photon leg is
$$
V_{i\gamma}(0,0)=\frac{4}{M_V}\left[
\frac{m_a^2}{2}(c_1+c_2-c_5)\operatorname{tr}(T_i Q U_a)
-2c_3\operatorname{tr}(T_i QY_a)\right].
$$
The $c_6$ term vanishes because the photon virtuality is zero; neutral
commutator and on-shell wedge terms $c_4,c_7$ vanish. There are two photon
placements. With the signed mixing convention of the photon-vector source,
their scalar contribution, after factoring $i e^2 E$, is
$$
G_{\rm VJP}=-4\sqrt2\sum_i W_i(0)V_{i\gamma}(0,0).
$$

For VVP, both vector legs convert. The ordered numerator is
$$
P_{ij}(0,0)=4m_a^2(d_1-d_3)\operatorname{tr}(T_iT_jU_a)
-16d_2\operatorname{tr}(T_iT_jY_a),
$$
where the $d_4$ neutral contraction vanishes. The ordered sum includes both
assignments when $i\ne j$:
$$
G_{\rm VVP}=2\sum_{i,j}W_i(0)W_j(0)P_{ij}(0,0).
$$
The scalar convention is $i\mathcal M=i e^2 E G$, with
$E=\epsilon(q_1,\varepsilon_1,q_2,\varepsilon_2)$ and the canonical
$\epsilon^{0123}=+1$. In the Bauer convention the local coefficient has
$G_{\rm local}=-c_\gamma/(4\pi^2f_a)$. Therefore
$$
\boxed{\Delta c_\gamma^{\rm vector-mixing}
=-4\pi^2f_a\left(G_{\rm VJP}+G_{\rm VVP}\right)},\qquad
c_\gamma^{\rm eff}=c_\gamma+\Delta c_\gamma^{\rm vector-mixing}.
$$
The term in the Lagrangian is
$\alpha_{\rm em}\Delta c_\gamma^{\rm vector-mixing}aF\widetilde F/(4\pi f_a)$.
It is only the retained resonance contribution; any independent local photon
coefficient is supplied separately. In this no-mixing basis the meson-pole
piece already belongs to the transformed photon coefficient and is not added
again.

# Null-direction check and limits

The [finite real-photon null audit](alp-photon-rcht-blind-direction-audit.md)
finds zero projection of the retained real-photon ALP resonance rows onto the
current structural blind directions. The [measured-null audit](../../rcht/odd-parity/odd-alp-measured-null-audit.md)
also reports rank zero for the $a\gamma\gamma$ image of the original-eleven
fit nullspace at its fixed fitted nuisance points. Thus those odd-coefficient
nulls do not affect this retained projection. This is a fixed-input,
tree-level statement, not a claim about joint nuisance refits, higher
operators, loops, or all orders.

The calculation keeps leading electromagnetic order and stable tree vector
propagators at zero photon virtuality. It omits widths, loops, kinetic
breaking, all-orders coupled-propagator resummation and additional odd
operators. The executable is
[`scripts/implementation/alp_photon_coupling.py`](../../../scripts/implementation/alp_photon_coupling.py).
For example:

```text
.venv/bin/python scripts/implementation/alp_photon_coupling.py --axial 1 0 0 --pseudoscalar 0 0 0 --alp-mass 1 --alp-decay-constant 1000 --c-gamma 0
```

The script returns the VJP and VVP pieces separately and uses the fitted
generators, masses, signed residues, and one-trace $F_V,z_V$ source from the
provisional fit.
