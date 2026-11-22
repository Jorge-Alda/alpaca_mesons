---
title: One-loop power-counting estimate for an ALP and three mesons
status: mixed
created: 2026-11-22
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd, Manohar:1983md, Kaiser:1998ds, Herrera-Siklody:1996tqr]
provenance: analytical-derivation
---

# One-loop power-counting estimate for $a\to P_{n_1}P_{n_2}P_{n_3}$

Use the [dressed Lagrangian and vertex derivatives](alp-dressed-vertex-power-counting-invariants.md#result-1-vertices-in-the-invariants-new-form-of-l1-verified) and the [exact tree pole and amplitude formulas](alp-finite-mixing-three-meson-tree.md#result). The ALP-like pole has mass $m_a$, ALP component $1/\sqrt{N(m_a^2)}$, and meson component $2\epsilon\mathbb{M}_{m_a^2}^{-1}(\mathcal C)/\sqrt{N(m_a^2)}$. All vertex momenta are incoming. The leading explicit single-ALP-line term uses the unmixed meson projectors on the other lines. Corrections to those projectors start two powers of $\epsilon$ higher for that term. The exact tree amplitude stays in the denominator of the [crossing criterion](alp-loop-crossing-exact-tree-condition.md). The estimate uses the [natural-size assumption](alp-loop-naturalness-assumption.md) away from thresholds, internal tree poles and zeros of the tree amplitude.

## Direct six-leg tadpole

Use $V_{n_a;\lambda_1\cdots\lambda_{n_P}}$ from the general vertex formula, including the kinetic derivative for $n_a=0$. In the following projection each vertex includes the three daughter meson directions. The additional meson directions, where present, are all $\mathbb{M}_{m_a^2}^{-1}(\mathcal C)$. The three remaining legs of the six-leg vertex are the external ALP and both ends of the tadpole. Independent differentiation in the three pole directions gives

$$
\frac{8\epsilon^3 V_{0;6}+12\epsilon^2 V_{1;5}+6\epsilon V_{2;4}+V_{3;3}}
     {N(m_a^2)^{3/2}} .
$$

Here the semicolon counts meson legs and the notation suppresses their directions. Each $V_{n_a;6-n_a}$ already contains $\epsilon^{n_a}$. The four potential terms use $\widehat S_0$, $\widehat P_1$, $\widehat S_2$ and $\widehat P_3$, respectively, through the symmetrized flavour traces of the general vertex. There is no topological vertex with this many meson legs. The potential part of $V_{0;6}$ survives in general. The kinetic part follows by taking six directional derivatives of

$$
\frac{1}{1440f_\pi^4}
\left\langle
[\Phi,[\Phi,\partial_\mu\Phi]]
[\Phi,[\Phi,\partial^\mu\Phi]]
\right\rangle.
$$

The kinetic numerator is at most quadratic in the tadpole momentum $k$. Odd terms vanish, and a quadratic term reduces by Lorentz symmetry to its scalar $k^2$ moment. In dimensional regularization the polynomial part of $k^2/(k^2-m_a^2)$ is scaleless. The remaining part is $m_a^2/(k^2-m_a^2)$. Thus after this reduction the estimate for the magnitude of the direct contribution is

$$
\frac{m_a^2}{32\pi^2 N(m_a^2)^{3/2}}
\left|8\epsilon^3 V_{0;6}+12\epsilon^2 V_{1;5}
      +6\epsilon V_{2;4}+V_{3;3}\right|.
$$

The $1/2$ in the prefactor is the identical-line tadpole symmetry factor. The estimate uses the reduced kinetic numerator in $V_{0;6}$. The vertex has mass dimension $-2$, so this term is dimensionless, as required for the four-point amplitude. Relative to a nonzero tree term linear in $\epsilon$, it starts at order $\epsilon^2$.

## Two-quartic bubble

Split the four external legs into the three unordered pair partitions. For each partition sum over the second internal pole, a meson. At each of the two quartic vertices, insert one ALP-like direction and one meson direction for the internal legs, with momenta opposite at the two ends of each propagator. The external ALP direction is inserted at the vertex carrying that leg. The product of the two projected vertices gives the bubble numerator, with symmetry factor one for distinct internal species.

After factoring $\epsilon/\sqrt{N(m_a^2)}$ from each ALP-like direction, the potential parts of the quartic vertices containing one and two such directions, respectively, are

$$
V_{0;L_1L_2L_3L_4}+\epsilon^{-1}V_{1;L_2L_3L_4},
\qquad
V_{0;L_1L_2L_3L_4}
+2\epsilon^{-1}V_{1;L_2L_3L_4}
+\epsilon^{-2}V_{2;L_3L_4}.
$$

Here the $L_i$ are the meson arguments of $V$, with $L_1=2\mathbb{M}_{m_a^2}^{-1}(\mathcal C)$ in the first expression and $L_1=L_2=2\mathbb{M}_{m_a^2}^{-1}(\mathcal C)$ in the second. The other arguments are the relevant meson projectors. Let the indices on $\mathcal V$ specify physical pole projections, including the extracted factors when present. Counting the magnitudes of the momentum-independent bubble products separately gives
$$
\frac{1}{16\pi^2}
\sum_{(i,j,k)\,\mathrm{cyc}}\sum_n
\left|\mathcal V_{a n_i a n}\mathcal V_{n_j n_k a n}\right|.
$$
The $\mathcal V$ factors in this expression are the potential parts. The loop function is a bubble at $(p_a+p_{n_i})^2$ with masses $m_a,m_{n,0}$. The first term of each projected vertex also admits the [four-meson kinetic vertex](alp-finite-mixing-three-meson-tree.md#result), evaluated with the same meson directions and routed momenta. Its products use the [tensor reduction for the bubble](alp-loop-momentum-reduction.md#angular-average), including the scalar moments with fewer denominators produced by numerator cancellation. The overall projected numerator has the same $\epsilon^3/N(m_a^2)^{3/2}$ factor as the direct tadpole.

## Assembly at one loop

The direct four-point terms also include cubic–quintic bubbles and cubic–cubic–quartic triangles. The one-loop four-point one-particle-irreducible valences are
$$
(6),\quad (3,5),\quad (4,4),\quad (3,3,4),\quad (3,3,3,3).
$$
A pure-meson cubic vanishes in the dressed Lagrangian. For the explicit one-ALP-line box, one external ALP and one internal ALP provide only three ALP-like ends, insufficient for four nonzero leading cubic vertices. A box can contribute at higher $\epsilon$ through further overlaps or ALP-like lines. For the other valences, the projected vertices are the corresponding derivatives of the same Lagrangian. All contributing products have to be retained whenever they share the leading power.

For the potential part of a cubic–quintic bubble, let $i$ run over the four external pole labels and let $j,k,l$ denote its complement. The corresponding magnitude estimate with one internal ALP-like line is
$$
\frac{1}{16\pi^2}
\sum_i \sum_n
\left|\mathcal V_{i a n}\mathcal V_{jkl a n}\right|.
$$
The two propagating species are distinct, so there is no identical-line factor. For the potential part of the triangle, let $\{k,l\}$ run over the six external pairs, and let $i,j$ denote the complementary single legs in a fixed order. Sum over internal pole labels $r,s,t$ with exactly one equal to $a$. The analogous estimate is
$$
\frac{1}{16\pi^2}
\sum_{\{k,l\}}\,
\sum_{\substack{r,s,t\\ \delta_{ra}+\delta_{sa}+\delta_{ta}=1}}
\frac{\left|\mathcal V_{i r s}\mathcal V_{j s t}\mathcal V_{kl t r}\right|}{m^2}.
$$
The scale $m^2$ is chosen separately for each loop as in the [naturalness record](alp-loop-naturalness-assumption.md). Each external block is labeled, so each internal species assignment around this fixed cycle is counted once. The momentum-independent products are evaluated by the invariant vertex derivatives with the same pole projections as the two-quartic bubble. The kinetic quartic in the triangle has a single quadratic numerator, which uses the [triangle case of the tensor reduction](alp-loop-momentum-reduction.md#angular-average).

The [one-loop Hessian identity](../../nlo/loop-corrections-to-alp-meson-lsz-mixing.md#mixed-loop-self-energy) generates the projected one-, two-, three- and four-point corrections from the same vertex derivatives. Writing $H$ for that Euclidean Hessian, the part of its expansion needed for up to four external legs is

$$
\frac12\sum_{r=1}^{4}\frac{(-1)^{r+1}}{r}
\operatorname{Tr}\!\left[\big(H[0]^{-1}(H-H[0])\big)^r\right].
$$

In each trace, assign the external legs to nonempty ordered blocks around the loop. Each block differentiates the Hessian and supplies a vertex with two propagating legs. The trace includes the momentum routing and the internal pole sums. The projected Euclidean vertices use the negative of the Minkowski potential derivative and contract the Fourier derivatives with the Euclidean scalar product. Continue the resulting invariant momenta and convert the effective vertices to the inverse-propagator and vertex conventions of the tree record before assembling the on-shell expression. For the explicit ALP-line contribution retain the cycles with exactly one ALP-like propagator. This generates the vertices needed both for the direct terms and for the corrections in an exchange channel. With $D$ and the cubic and quartic vertices defined in the tree record, the first-order correction in a fixed channel before varying the external legs is

$$
\delta\mathcal V_4
-\delta\mathcal V_3 D^{-1}\mathcal V_3
-\mathcal V_3 D^{-1}\delta\mathcal V_3
+\mathcal V_3 D^{-1}\delta D D^{-1}\mathcal V_3 .
$$

The direct $\delta\mathcal V_4$ is included once, and the last three terms are summed over the three tree exchange channels. For the external projection, let $v^{(r)}$ be a real unit tree pole vector in the dressed basis and $\Sigma$ the loop change of $D$. Linearizing the [pole equation and residue condition](../../nlo/loop-corrections-to-alp-meson-lsz-mixing.md#pole-residue-and-amplitudes) at a nondegenerate pole gives

$$
\delta m_r^2=-v^{(r)T}\Sigma(m_r^2)v^{(r)},
$$
$$
\delta v^{(r)}
=-\sum_{n\ne r}
  \frac{v^{(n)}v^{(n)T}\Sigma(m_r^2)v^{(r)}}{m_r^2-m_n^2}
-\frac12 v^{(r)}
  v^{(r)T}\Sigma'(m_r^2)v^{(r)} .
$$

These signs use $D+\Sigma$ for the full inverse propagator and $D'(s)=\mathbf 1$. Vary each external leg once, and include the shift of the on-shell kinematics where required. For a fixed physical ALP pole mass, use the same pole-fixing convention in the two-point correction as in the tree comparison. With a nonvanishing class-6 one-point source, one can evaluate about the shifted stationary point. The displacement is

$$
\delta\Phi_I=-\sum_J [D(0)^{-1}]_{IJ}\delta\mathcal V_J .
$$

Here $\delta\mathcal V_J$ is the one-point derivative in the conventions of $D(0)$, and $\Phi_I$ is the full quadratic multiplet. This derivation uses the displacement to incorporate vacuum stationarity. In each vertex correction add the derivative of the tree vertex along $\delta\Phi_I$. In the two-point correction add the derivative of $D$ along the same displacement. The derivative of each tree quantity is another derivative of the dressed Lagrangian with one zero-momentum leg. This includes the same first-order term as explicitly attaching the one-point insertion through $D(0)^{-1}$, so include it once. If the stationary point has exact flat directions, the inverse and the choice of stationary point require a separate analysis.

Fix the physical ALP mass at one loop by using the two-point correction after the vacuum displacement in
$$
\delta m_{a,0}^2=\frac{v^{(a)T}\Sigma(m_a^2)v^{(a)}}{(v_a^{(a)})^2}.
$$
Subtract this change from the bare ALP diagonal entry of the inverse propagator before applying the pole and residue formulas. Then $\delta m_a^2=0$ while the external ALP projection and the internal exchange propagators can still change.

With the preceding changes included in $\delta D$ and $\delta\mathcal V$, the on-shell one-loop expression, prior to replacing the scalar moments by sizes, is
$$
\begin{aligned}
\delta\mathcal M
={}&\delta\mathcal V_{a n_1n_2n_3}
-\sum_{i=1}^3\Big[
  \delta\mathcal V_{a n_i}D^{-1}\mathcal V_{n_j n_k}
 +\mathcal V_{a n_i}D^{-1}\delta\mathcal V_{n_j n_k}
 -\mathcal V_{a n_i}D^{-1}\delta D D^{-1}\mathcal V_{n_j n_k}
\Big]\\
&+\sum_{r\in\{a,n_1,n_2,n_3\}}
\left[
\delta v^{(r)}\cdot\nabla_{v^{(r)}}+
\delta m_r^2\,\partial_{m_r^2}
\right]\mathcal M^{(0)} .
\end{aligned}
$$
The cubic vertices with two pole labels have an open internal index, and the associated propagators are evaluated at $(p_{n_j}+p_{n_k})^2$. In the external-leg sum use one term per external leg, including repeated species. The mass derivative acts on the on-shell kinematics with two independent Dalitz invariants held fixed. The last invariant then follows from the sum of external squared masses. In that derivative hold the pole vectors fixed, since their change has its own term.

## Extracting and estimating the leading coefficient

For the class with an explicit internal ALP-like pole, select exactly one such propagator in each cycle and project the corresponding external and internal legs. Meson-only cycles give an ordinary chiral correction at order $\epsilon$ whose ratio to tree stays finite as $\epsilon\to0$. For their additional mixing-dependent term, expand the meson-like part of the spectral sum using the [meson pole and projector expansion](alp-three-meson-loop-mixing-expansion.md). The external ALP-like normalization, daughter mass and overlap variations, double denominators and momentum derivatives all enter before the coefficient is extracted. The two procedures distinguish the spectral origins of the terms and avoid counting the ALP pole in both. With both contributions included, take
$$
\left.\frac{\epsilon^3}{3!}\frac{\partial^3\delta\mathcal M}{\partial\epsilon^3}\right|_{\epsilon=0}
$$
before assigning loop sizes. Take the derivative at fixed physical $m_a$ and fixed independent kinematic invariants on the smooth perturbative branch.

For each retained term build the product of the projected vertex numerators and perform the [scalar-moment reduction](alp-loop-momentum-reduction.md#scalar-moments-and-derivatives). For equal scalar moments at the same masses and invariants, combine their weights from the one-particle-irreducible, exchange, vacuum and pole terms before taking magnitudes. Assign a scalar moment with $r$ remaining denominators the size $(m^2)^{2-r}/(16\pi^2)$, including the raised powers from momentum or mass derivatives. The reference scale is the one specified by the [naturalness prescription](alp-loop-naturalness-assumption.md) for that class. Use the sum of the resulting magnitudes in the crossing, under the assumption that contributions with different scalar moments do not substantially cancel. This is a magnitude prescription. The signed physical one-loop amplitude depends on the actual loop functions.

## Symbolic checks

The [six-leg derivation](../../../scripts/derive_alp_three_meson_loop_estimate.py) derives the three ALP pole factors, expands the kinetic trace through six fields and checks the tadpole numerator division. The [topology derivation](../../../scripts/derive_alp_three_meson_loop_topologies.py) checks the direct valences, bubble partitions and symmetry factors, the quartic potential projections, and the first variations of the exchange, pole vectors, pole residues and vacuum. The [cycle assembly](../../../scripts/derive_alp_three_meson_loop_assembly.py) implements the Hessian expansion, projected kinetic vertices, momentum routing, and on-shell linear variation. The [projected cycle reduction](../../../scripts/derive_alp_three_meson_loop_cycle_reduction.py) builds the scalar-moment coefficients from these vertices and checks this interface on a symbolic cubic bubble and a six-leg tadpole. The scripts check the ingredients and their interfaces. They leave the full summed expression as the projection formula above, with all physical inputs symbolic.

```sh
.venv/bin/python scripts/derive_alp_three_meson_loop_estimate.py
.venv/bin/python scripts/derive_alp_three_meson_loop_topologies.py
PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_three_meson_loop_estimate derive_alp_three_meson_loop_topologies
```
