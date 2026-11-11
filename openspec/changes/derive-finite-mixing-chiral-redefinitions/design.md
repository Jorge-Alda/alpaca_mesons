# Design

## Context

See proposal.md for motivation. The partonic transformation laws, transport matrices and dressed $\widehat M(x)$ exist in [alp-general-chiral-redefinition-invariants.md](../../../memory/alp/shared/alp-general-chiral-redefinition-invariants.md) and `scripts/derive_alp_invariant_coefficient_functions.py`. The finite-mixing record aligns one pole with an axial group element $L_a=\exp[-i a\,w_i^{(a)}\lambda_i/(2f_\pi)]$, $R_a=L_a^\dagger$, which is a constant-generator case of $g_{R,L}(x)$ with $G_R=-G_L=-(f_a/2f_\pi)\,w_i\lambda_i$. It treats the quadratic Hessian and one-ALP tangent-space vertices; the nonlinear map and the anomaly are only described in words. All conventions follow Alda:2026nrd. The results hold for the restricted Lagrangian of the partonic record: no dipoles, no higher-derivative bilinears, at most one $\partial_\mu a$ per interaction. Loops and higher chiral orders are outside the scope (assumption).

## Goals / Non-Goals

**Goals:**
- One derivation chain: partonic $g_{R,L}(x)$ → CCWZ sources and anomaly → meson coordinate map → quadratic form $D(s)$ → pole alignment → transformed vertices.
- Statements that hold at all orders in $1/f_a$ and in the meson fields, with explicit conditions for when they fail.
- Symbolic checks that exercise noncommuting matrices.

**Non-Goals:**
- Solving the spectrum numerically, decay widths, or fits (alp-sjp.2, .3, .5, .6).
- Any resummation of omitted chiral operators, loops or resonance effects.
- Operators outside the restricted truncation.

## Decisions

1. **Use $\widehat U=W_R^\dagger UW_L$ as the source-level bridge.** The transport matrices remove $r_\mu,\ell_\mu$ exactly, because the pullback connections depend on the single coordinate $x$ and are flat. The remaining ALP dependence sits in $\widehat M(x)$ and in the anomaly terms. Alternative: expand the CCWZ blocks order by order in $a$; rejected because it does not give all-orders statements (the existing expansion record covers that use).

2. **Treat anomaly through the Jacobian of the singlet part, using a general $\theta(x)$.** The trace of the generator $\operatorname{tr}(G_R-G_L)$ multiplies $\Theta$, and the matching meson-level term is the $X_{TC}$ shift, in Alda:2026nrd's normalization. The partonic-to-Alda conversion $\Theta_{\rm Alda}=(4\pi/\alpha_s)\Theta_{\rm raw}$ is taken from the partonic record. Alternative: ignore the singlet part; rejected since the issue lists the anomaly explicitly.

3. **Describe the quadratic problem as a congruence, $D'(s)=T^{\mathsf T}D(s)T$.** At the vacuum the Hessian depends only on the linear part $T$ of the coordinate map; invariance of poles and residues then follows from $\det D'=\det T^2\det D$ and $D'^{-1}=T^{-1}D^{-1}T^{-\mathsf T}$. The shift $P=P'+wa$ is the special case with $T$ triangular. Alternative: work only with shifts; rejected because it cannot describe alignment of several poles.

4. **Pole alignment is a statement about one column of $T$.** For a selected pole $s_k$ with null vector $(w^{(k)},1)$, choose $T$ whose last column is that vector. Several poles: the shift form has $N$ free parameters per column and the alignment conditions $D_{i\hat a}(s_k)+D_{ij}(s_k)w_j=0$ for different $k$ are consistent only when the $w^{(k)}(s_k)$ coincide; a general $T$ built from all null vectors diagonalizes $D(s)$ at every $s$, since $D^{-1}$ then has one simple pole per coordinate only if the residue structure allows it. The derivation must check that last claim rather than assume it; $K$ and $M^2$ simultaneous diagonalization is a statement about the pencil $(K,M^2)$ and holds for $K>0$.

5. **Nonlinear terms come from the group action, not from a tangent-space truncation.** The meson map is $U'=g_RUg_L^\dagger$ with the same $g_{R,L}(x)$ and, for a linear-in-$a$ pole shift, with $g$ built from the pole vector. Its compensator $h$ and the transformed $u_\mu$, $\chi_\pm$ give the full vertex; the script expands in meson fields to the order needed for the $a\to PP$ and $a\to PPP$ amplitudes. The meson coordinates $P_i'$ in the pole-aligned system are defined through this map, so the linear relation $P_i=P_i'+w_ia$ is only its leading term.

6. **Verification by exact algebra with random rational matrices, with `sympy`.** Noncommuting $U(2)$ and $U(3)$ generators, Taylor order up to 3 for partonic laws, and exact $\det$ and residue comparisons for the quadratic form. Scripts are written by the `coder` agent under `scripts/`.

## Risks / Trade-offs

- [Nonlinear map is large for $U(3)$ to the order needed] → verify symbolically for $U(2)$ first, $U(3)$ at fixed numeric rational data, and limit to the vertex orders needed by alp-sjp.3.
- [Transport matrices for noncommuting $C_{R,L}(x)$ have no closed form] → work with Taylor coefficients and the recursion; closed forms only for commuting or constant-generator cases.
- [The shifted pole vector diverges near degeneracy] → state the regularity condition and use the pencil formulation; do not hide it in a small-mixing expansion.
- [A general $T$ may mix mesons among themselves, leaving the Alda:2026nrd coordinates] → keep $T$ as a coordinate change only and map back with $T^{-1}$ before comparing amplitudes.
- [Anomaly normalization mismatch between partonic and Alda:2026nrd records] → one explicit conversion check in the script, as in the partonic record.

## Open Questions

- None that change the specs or the task breakdown. Notation for any new named quantity (for example the pole null vector) will be proposed to the user for approval as required by AGENTS.md; the existing $w_i^{(a)}$ is reused.
