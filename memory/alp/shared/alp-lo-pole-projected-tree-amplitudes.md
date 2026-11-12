---
title: Pole-projected leading-order tree amplitudes
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Pole-projected leading-order tree amplitudes

## Inputs and assumptions

Use natural units, the fields $(P_i,a)$ and conventions of [@Alda:2026nrd]. The real symmetric, constant matrices $K$ and $M^2$ give $D(s)=sK-M^2$, with positive-definite $K$, as in the [pole solver](alp-lo-pole-solver.md). The source functions and their Taylor derivatives, aligned-vacuum assumption, and CCWZ construction of each amputated vertex are those of the [joint-vertex derivation](alp-lo-joint-vertex-expansion.md); no expansion in $f_\pi/f_a$ is made. Include the general mass-source second derivative in $M^2_{aa}$ as stated there. The explicit flavour-space inverse in the [operator-inverse record](alp-lo-finite-mixing-operator-inverse.md) assumes a linear mass spurion; when its other assumptions hold, its overlap can be used below, with the additional $M^2_{aa}$ term included in the pole equation for a nonlinear mass source. For poles with zero ALP component its overlap representation does not apply and the full spectrum solver supplies the pole vector.

Momenta are incoming; $s$ and entries of $M^2$ are squared masses in a common unit, for example GeV$^2$. A vertex $\mathcal V_{I_1\cdots I_n}$ is the derivative of the LO Lagrangian and omits the overall Feynman-rule factor $i$. External legs are stable, tree-level, simple poles. Momentum conservation and the pole conditions are imposed on the input scalar products. These results exclude loops, widths and complex poles. Near degeneracies and large admixtures require the qualifications in the [truncation record](alp-lo-truncation-limits.md); a degenerate subspace has a well-defined summed residue, while individual pole vectors depend on the chosen basis.

## Method and analytical expressions

1. The solver returns $v^{(n)}$ with $v^{(m)\mathsf T}Kv^{(n)}=\delta_{mn}$. Project every external leg, including the ALP component of each daughter. For three legs the reduced amplitude is
$$
\mathcal M_3=\mathcal V_{IJK}v_I^{(n_1)}v_J^{(n_2)}v_K^{(n_3)}.
$$
The pole vectors already contain the square roots of the residues; no additional external normalization is applied. A complex linear combination of vectors in a degenerate subspace may be supplied for a charged state. The directions use all-incoming fields: cross each outgoing momentum and choose the conjugate field direction appropriate to the outgoing charged particle.

2. For four legs, contract the four-point vertex with all four pole vectors. Add the three distinct partitions $\{12|34,13|24,14|23\}$ of external legs. For the partition $12|34$, put $q=-(p_1+p_2)=p_3+p_4$ and add
$$
-\bigl(\mathcal V_{IJK}(p_1,p_2,q)v_I^{(n_1)}v_J^{(n_2)}\bigr)
[D(q^2)^{-1}]_{KL}
\bigl(\mathcal V_{LMN}(-q,p_3,p_4)v_M^{(n_3)}v_N^{(n_4)}\bigr).
$$
The minus sign follows from the two $i\mathcal V_3$ vertices and the propagator $iD^{-1}$ relative to the four-point factor $i\mathcal V_4$. There is no extra factorial: each vertex is already symmetrized. The internal momenta are assigned with all momenta incoming at each vertex. The full inverse is used, without expanding the mixing.

3. Partition $D$ into its meson block $D_{PP}$ and ALP row and column. At values where $D_{PP}$ is invertible, its Schur complement and unnormalized pole direction are
$$
D_{aa}-D_{aP}D_{PP}^{-1}D_{Pa},\qquad
\begin{pmatrix}-D_{PP}^{-1}D_{Pa}\\1\end{pmatrix}.
$$
The full analytical inverse is the sum of $D_{PP}^{-1}$ embedded in the meson block and the outer product of this direction divided by the Schur complement. This expresses the same construction as the [flavour operator inverse](alp-lo-finite-mixing-operator-inverse.md), but is applicable to the general quadratic blocks of the joint-vertex derivation. At a simple zero, the derivative of the Schur complement is the quadratic form of $K$ on that direction. Apparent singularities at zeros of $\det D_{PP}$ can be removable; the direct inverse is used for exchange evaluation. A genuine exchanged pole requires a prescription beyond this stable tree-level calculation.

4. For this finite-dimensional LO theory the exact expression is
$$
D(s)^{-1}=\sum_n\frac{v^{(n)}v^{(n)\mathsf T}}{s-m_n^2}.
$$
It holds at **every non-pole** $s$, not merely on shell. At a pole both sides diverge, so their residues are compared instead. This equality does not mean that retaining only one pole is sufficient away from that pole.

5. For a constant invertible field change $\Phi=T\Phi'$, every vertex index transforms with $T$, $D'=T^{\mathsf T}DT$, and $v^{(n)\prime}=T^{-1}v^{(n)}$. Thus the contact and exchange contractions are individually unchanged, with external state phases transported consistently; see the [congruence record](alp-inverse-propagator-congruence.md).

## Reproducibility and checks

[`scripts/derive_alp_pole_projected_amplitudes.py`](../../../scripts/derive_alp_pole_projected_amplitudes.py) uses the spectrum solver and the CCWZ vertex generator. It takes a Gram matrix of all-incoming momenta, including their conservation and on-shell relations, and evaluates the three- or four-leg expression with exact SymPy arithmetic. Its example derives the Schur inverse and compares it with the complete pole sum. The test suite checks the symbolic Schur inverse and derivative normalization; the overlap and Schur complement against the explicit flavour operator inverse for exact generic Hermitian sources; the added mass-curvature term; full inverse and residues against the spectrum solver; a projected CCWZ three-leg vertex and four-leg contact plus exchanges; and an exact on-shell random constant field change with all tensor indices transformed, nonzero daughter ALP components and off-diagonal propagation. One additional test computes the quadratic matrices and the contact and exchange vertices from the same CCWZ sources in a two-field restriction, solves their poles, and evaluates an on-shell amplitude with nonzero contact and exchange contributions. This restricted algebraic example does not assert that the omitted meson directions decouple from the full theory. The illustrative matrices and source coefficients are test inputs, not a phenomenological prediction; the independently constructed matrices in the other amplitude tests are not claimed to be the Hessian of their illustrative source functions.

```sh
.venv/bin/python scripts/derive_alp_pole_projected_amplitudes.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_pole_projected_amplitudes.py
```
