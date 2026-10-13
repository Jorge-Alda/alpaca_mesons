---
title: One-ALP two-pseudoscalar vector- and scalar-source terms at LO in chiPT
status: derived
created: 2026-09-30
updated: 2026-09-30
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# One-ALP Two-Pseudoscalar Terms from $c_q^V$ and $c_q^S$

## Scope, conventions, and units

This record derives the terms linear in the ALP $a$ and quadratic in the light pseudoscalar matrix $\Phi$ from

$$
\mathcal L_{\mathrm{LO}}=\frac{f_\pi^2}{4}
\left\langle u_\mu u^\mu+\chi_+\right\rangle.
$$

The executable calculation is [`scripts/derive_alp_phi2_lo_chipt.py`](../../../scripts/derive_alp_phi2_lo_chipt.py). Natural units $\hbar=c=1$ are used. Thus $a$, $\Phi$, $f_a$, $f_\pi$, $B_0$, and $m_q$ carry their standard mass dimensions, and every displayed Lagrangian term has mass dimension four. Flavour brackets denote a trace. The commutator-map convention is

$$
\phi=\frac{\Phi}{2f_\pi},\qquad
\mathcal C_\phi[X]=[X,\phi],\qquad
\mathcal A_\phi[X]=\{X,\phi\}.
$$

These definitions, the LO Lagrangian, and the all-orders $u_\mu$ and $\chi_+$ inputs are those recorded in the [Appendix-A extract](../../extracts/Alda:2026nrd/appendix-a-expansion.md) [@Alda:2026nrd]. The physical light-field generator normalization used for the channel audit is the repository convention $\operatorname{tr}(\lambda_a\lambda_b)=2\delta_{ab}$ for the unmixed basis, consistent with the existing three-pseudoscalar implementation.

## Assumptions and exclusions

- The expansion is at leading chiral order and keeps exactly $O(a\Phi^2/f_a)$.
- Flavour matrices are not assumed to commute during the derivation.
- The physical-channel audit assumes an electrically neutral ALP source: its flavour matrix is block diagonal between $u$ and the equal-charge $(d,s)$ sector. Both diagonal and $d$--$s$ sources are retained.
- The symbolic $d$--$s$ source and its Hermitian conjugate are classified separately; Hermiticity relates their coefficients but is not used to erase either direction.
- The $\eta$--$\eta'$ rotation is kept at a generic real angle $\theta_P$; no numerical mixing approximation is chosen.
- Isospin symmetry and equal pseudoscalar masses are not assumed in the derivation. Exact particle--antiparticle mass equality is invoked only when identifying the on-shell zero of flavour-diagonal vector-source channels.
- This is the even-intrinsic-parity LO chiral Lagrangian. Wess-Zumino-Witten terms, vector mesons, and all heavier mesons are excluded.
- The equations are Lagrangian coefficients. No extra factorial from functional differentiation into a vertex rule is included.

## Stepwise derivation

### 0. Lowest allowed pseudoscalar multiplicity

The vector-source series in $u_\mu$ is odd in $\Phi$, while the ordinary mesonic vielbein is also odd, so their cross term contains only even powers $a\Phi^{2n}$. The scalar source enters $\chi_+$ through $\cos(\mathcal A_\phi)$ and likewise produces only $a\Phi^{2n}$. At zero pseudoscalar order its flavour trace is

$$
\left\langle[c_q^S,m_q]\right\rangle=0.
$$

Neither source generates an $a\Phi$ term. Therefore the lowest nontrivial mesonic interactions are the quadratic terms derived below, and the lowest-multiplicity physical candidates are two-body decays $a\to PP$.

### 1. Vector source in the kinetic term

The required terms from the [Appendix-A formulas](../../extracts/Alda:2026nrd/appendix-a-expansion.md) are [@Alda:2026nrd]

$$
u_\mu^{(0)}=-2\partial_\mu\phi+O(\phi^3),
\qquad
u_\mu^{(a,V)}=-2i\frac{\partial_\mu a}{f_a}
\mathcal C_\phi[c_q^V]+O(a\phi^3).
$$

SymPy forms both matrix orderings in $u_\mu u^\mu$, extracts one ALP and two powers of $\Phi$, and applies trace cyclicity. The result is

$$
\boxed{
\left.\mathcal L_{\mathrm{kin}}\right|_{a\Phi^2,c_q^V}
=\frac{i\,\partial_\mu a}{2f_a}
\left\langle
\partial^\mu\Phi\,[c_q^V,\Phi]
\right\rangle .
}
$$

The cancellation of $f_\pi$ is exact and is not an approximation.

### 2. Scalar source in the symmetry-breaking term

The ALP-dependent scalar source is

$$
\left.\chi_S\right|_a
=2iB_0\frac{a}{f_a}[c_q^S,m_q].
$$

Since

$$
2\cos(\mathcal A_\phi)[\chi_S]
=2\chi_S-\mathcal A_\phi^2[\chi_S]+O(\phi^4),
$$

coefficient extraction gives

$$
\left.\mathcal L_\chi\right|_{a\Phi^2,c_q^S}
=-\frac{iB_0a}{8f_a}
\left\langle
\{\Phi,\{\Phi,[c_q^S,m_q]\}\}
\right\rangle.
$$

The trace identity

$$
\left\langle\{\Phi,\{\Phi,X\}\}\right\rangle
=4\left\langle\Phi^2X\right\rangle
$$

then yields

$$
\boxed{
\left.\mathcal L_\chi\right|_{a\Phi^2,c_q^S}
=-\frac{iB_0a}{2f_a}
\left\langle
\Phi^2[c_q^S,m_q]
\right\rangle .
}
$$

In particular, this term vanishes when $[c_q^S,m_q]=0$.

### 3. Generic component form

For two distinct commuting fields in $\Phi=\phi_i\lambda_i+\phi_j\lambda_j+\cdots$, define

$$
F_{ij}^V=\left\langle c_q^V[\lambda_i,\lambda_j]\right\rangle,
\qquad
F_{ij}^S=\left\langle
\{\lambda_i,\lambda_j\}[c_q^S,m_q]
\right\rangle.
$$

Their interaction is

$$
\boxed{
\mathcal L_{a\phi_i\phi_j}
=\frac{i}{2f_a}F_{ij}^V\,\partial_\mu a
\left(\phi_i\partial^\mu\phi_j-\phi_j\partial^\mu\phi_i\right)
-\frac{iB_0}{2f_a}F_{ij}^S\,a\phi_i\phi_j.
}
$$

The vector coefficient is antisymmetric, so it vanishes for two identical fields. The scalar coefficient is symmetric. For an identical field $i$, its scalar Lagrangian coefficient follows directly from the matrix result as

$$
-\frac{iB_0a}{2f_a}\phi_i^2
\left\langle\lambda_i^2[c_q^S,m_q]\right\rangle.
$$

## Physical channel classification

The script audits all total-charge-zero pairs drawn from

$$
\{\pi^+,\pi^-,\pi^0,K^+,K^-,K^0,\bar K^0,\eta,\eta'\}.
$$

The classification concerns nonzero off-shell flavour coefficients:

| source structure | channels |
|---|---|
| diagonal $c_q^V$ | $\pi^+\pi^-$, $K^+K^-$, $K^0\bar K^0$ |
| $d$--$s$ $c_q^V$ selecting $S=+1$ | $K^+\pi^-$, $K^0\pi^0$, $K^0\eta$, $K^0\eta'$ |
| diagonal $c_q^S$ with diagonal $m_q$ | none, because $[c_q^S,m_q]=0$ |
| $d$--$s$ $c_q^S$ selecting $S=+1$ | $K^+\pi^-$, $K^0\pi^0$, $K^0\eta$, $K^0\eta'$ |

The Hermitian-conjugate $s$--$d$ sources generate the charge-conjugate $S=-1$ channels $K^-\pi^+$, $\bar K^0\pi^0$, $\bar K^0\eta$, and $\bar K^0\eta'$. For the vector interaction, on-shell two-body kinematics gives

$$
p_a\cdot(p_j-p_i)=m_j^2-m_i^2.
$$

Consequently the three flavour-diagonal particle--antiparticle vertices listed above give zero on-shell amplitudes when each particle and antiparticle has its exact common mass. They remain nonzero off-shell interaction terms. The flavour-changing $K\pi$, $K\eta$, and $K\eta'$ families are not removed by this kinematic identity because their masses differ.

Under the no-ALP--meson-mixing prescription recorded in the [Appendix-A extract](../../extracts/Alda:2026nrd/appendix-a-expansion.md),

$$
\widehat c_q^V=c_q^{VS},\qquad \widehat c_q^S=0.
$$

Thus the scalar-source term vanishes in that parametrization, while the vector formula remains with $c_q^V$ replaced by $c_q^{VS}$ [@Alda:2026nrd]. This is a field-basis statement, not an additional dynamical approximation.

## Verification

Run

```text
.venv/bin/python scripts/derive_alp_phi2_lo_chipt.py
```

The script:

1. derives both raw nested expressions from the CCWZ series by symbolic coefficient extraction;
2. verifies both trace reductions using fixed noncommuting matrices;
3. checks the complete total-charge-zero light-channel sets for diagonal and $d$--$s$ sources;
4. verifies the exact on-shell reduction to $m_j^2-m_i^2$.

Successful execution ends with `Verification: all symbolic trace, channel-classification, and on-shell checks passed.`
