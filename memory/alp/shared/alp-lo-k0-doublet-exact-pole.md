---
title: Exact LO poles of the $K^0$ doublet for diagonal $\chi_0$, with CP-conserving, aligned-phase and CP-violating sources
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Exact LO poles of the $K^0$ doublet for diagonal $\chi_0$, with CP-conserving, aligned-phase and CP-violating sources

## Assumptions, inputs and notation

- Exact pole equations and canonical frame of the [exact pole-equation record](alp-lo-exact-pole-equations.md) (Steps 2, 3, 5, 7); $\mathbb M_s$, $\mathcal C$, $\lambda_n$, $M^2_{\hat a\hat a,\rm can}$, $\epsilon=f_\pi/f_a$ as there; conventions of Alda:2026nrd [@Alda:2026nrd]. No expansion in $\epsilon$.
- $\chi_0=\operatorname{diag}(x_u,x_d,x_s)$, $m_d\neq m_s$, sources Hermitian and commuting with $Q$. $m_{K^0,0}^2=\tfrac12(x_d+x_s)$ and the CP convention are those of the [neutral-block record](alp-lo-neutral-block-reduction.md).
- The $d$-$s$ part of a $3\times3$ matrix $X$ is the matrix whose only nonzero entries are $X_{ds}$ and $X_{sd}=X_{ds}^*$. The $K^0$ doublet is the space of $d$-$s$ parts, spanned by $\lambda_6$, $\lambda_7$.

## Step 1. Doublet block and its contribution to the secular sum

$\mathbb M_s$ acts on $d$-$s$ parts as multiplication by $m_{K^0,0}^2-s$, and the $d$-$s$ part of $\mathcal C$ is $(\mathcal C)_{ds}=-m_{K^0,0}^2(c_q^{AP})_{ds}$. For any orthonormal pair of generators $\lambda_n$ of the doublet (for example $\lambda_6$, $\lambda_7$),
$$
\sum_{n\in K^0\ \rm doublet}\operatorname{tr}(\lambda_n\mathcal C)^2=\operatorname{tr}(\lambda_6\mathcal C)^2+\operatorname{tr}(\lambda_7\mathcal C)^2=4|(\mathcal C)_{ds}|^2=4m_{K^0,0}^4|(c_q^{AP})_{ds}|^2,
$$
and the doublet contribution to the secular sum $-2\operatorname{tr}(\mathcal C\,\mathbb M_s^{-1}(\mathcal C))$ of the exact record is
$$
\boxed{\,\sum_{n\in K^0\ \rm doublet}\frac{\operatorname{tr}(\lambda_n\mathcal C)^2}{s-m_{K^0,0}^2}=\frac{4m_{K^0,0}^4|(c_q^{AP})_{ds}|^2}{s-m_{K^0,0}^2}\,}
$$
equal to $-2\operatorname{tr}$ of the $d$-$s$ part of $\mathcal C$ times $\mathbb M_s^{-1}$ of the $d$-$s$ part of $\mathcal C$.

## Step 2. Coupled combination and exact pole

For $(c_q^{AP})_{ds}\neq0$, choose the two doublet generators as the $d$-$s$ matrices with $(d,s)$ entry $(\mathcal C)_{ds}/|(\mathcal C)_{ds}|$ (coupled combination) and $i(\mathcal C)_{ds}/|(\mathcal C)_{ds}|$ (decoupled combination). Both are eigen-generators of $\mathbb M_0$ with eigenvalue $m_{K^0,0}^2$, with $\operatorname{tr}(\lambda\lambda')=2\delta$. They give
$$
\operatorname{tr}(\lambda_{\rm coupled}\,\mathcal C)^2=4m_{K^0,0}^4|(c_q^{AP})_{ds}|^2,\qquad \operatorname{tr}(\lambda_{\rm decoupled}\,\mathcal C)=0 .
$$
Here "coupled" and "decoupled" are labels of the two combinations, not new symbols. In the canonical frame the decoupled combination has no entry to the ALP. Its pole is therefore exactly
$$
\boxed{\,s=m_{K^0,0}^2\,}
$$
for any $\epsilon$ and any phases of the $d$-$s$ sources, with zero ALP component (Step 7 of the exact record). The coupled combination satisfies Step 5 of the exact record:
$$
\boxed{\,(s-m_{K^0,0}^2)\Big[s-M^2_{\hat a\hat a,\rm can}-\epsilon^2\sum_{n\in\{\pi^0,\eta,\eta'\}}\frac{\operatorname{tr}(\lambda_n\mathcal C)^2}{s-m_{n,0}^2}\Big]=4\epsilon^2m_{K^0,0}^4|(c_q^{AP})_{ds}|^2\,}
$$
with the residues $\operatorname{tr}(\lambda_n\mathcal C)^2$ of the [diagonal-$\chi_0$ record](alp-lo-diagonal-chi0-meson-poles.md). For $(c_q^{AP})_{ds}=0$ both combinations decouple and $m_{K^0,0}^2$ is a double pole.

Consistency with the exact pole of the neutral-block record, (c). There, the null vector of $D(m_{K^0,0}^2)$ in the bare frame has components $(-D_{7a},D_{6a})$ on $(P_6,P_7)$. Since $D_{ia}(m_{K^0,0}^2)=\epsilon\operatorname{tr}(\lambda_i\mathcal C)$ on the doublet (Step 3 of the second-order shift record, at $s=m_{K^0,0}^2$), the corresponding generator $-\operatorname{tr}(\lambda_7\mathcal C)\lambda_6+\operatorname{tr}(\lambda_6\mathcal C)\lambda_7$ is proportional to the decoupled combination. The canonical field change shifts the meson fields only by multiples of $a$, so a vector with zero ALP component is the same in both frames. The two descriptions give the same pole and the same state.

## Step 3. Reductions

- **CP-conserving sources** (real $d$-$s$ entries of $c_q^A$, $c_q^P$). Then $(\mathcal C)_{ds}$ is real, the coupled combination is $\lambda_6$ (up to sign) and the decoupled one is $\lambda_7$. $K_{7a}=M^2_{7a}=0$ and $\operatorname{tr}(\lambda_7\mathcal C)=0$, so $P_7$ decouples in both frames: one pole exactly at $m_{K^0,0}^2$ ($P_7$), and the $P_6$ pole from the boxed coupled equation.
- **Aligned phases** ($(c_q^A)_{ds}$ and $(c_q^P)_{ds}$ with the same phase $w$, so $\operatorname{Im}((c_q^A)_{ds}(c_q^P)_{ds}^*)=0$). The $s$-quark rephasing $V=\operatorname{diag}(1,1,w)$ of the neutral-block record is a constant vector rotation. It keeps $\chi_0$ diagonal, makes both $d$-$s$ entries real, and leaves $\det D(s)$ unchanged. The case reduces to the CP-conserving one, with $(P_6,P_7)$ rotated by the phase of $w$.
- **CP-violating sources** ($\operatorname{Im}((c_q^A)_{ds}(c_q^P)_{ds}^*)\neq0$). In the bare frame $K_{6a}M^2_{7a}-K_{7a}M^2_{6a}\neq0$, so no combination of $(P_6,P_7)$ has zero $K_{ia}$ and $M^2_{ia}$ (neutral-block record, (b)). In the canonical frame the decoupled combination of Step 2 still exists, and the pole at $m_{K^0,0}^2$ with zero ALP component remains. At fixed $M^2_{\hat a\hat a,\rm can}$ the poles depend on the $d$-$s$ sources only through $|(c_q^{AP})_{ds}|^2$: a rephasing of $\mathcal C$ that makes $(\mathcal C)_{ds}$ real leaves $\det(s\mathbb 1-M^2_{\rm can})$ unchanged. At fixed bare $\operatorname{tr}S_2$ instead, the relative phase of $(c_q^A)_{ds}$ and $(c_q^P)_{ds}$ also enters through $\operatorname{tr}\widehat S_2$ in $M^2_{\hat a\hat a,\rm can}$.
- **No $d$-$s$ sources** ($(c_q^{AP})_{ds}=0$). Double pole at $m_{K^0,0}^2$, as in the minimal-block table of the neutral-block record.

## Verification

[`scripts/derive_alp_k0_doublet_exact_pole.py`](../../../scripts/derive_alp_k0_doublet_exact_pole.py):

- Symbolic (diagonal $\chi_0$, complex $(c_q^{AP})_{ds}$, symbolic $c_{GA}$, $m_0^2$, $s$, $\epsilon$):
  - $(\mathcal C)_{ds}=-m_{K^0,0}^2(c_q^{AP})_{ds}$, and $\mathbb M_s$ on $d$-$s$ parts;
  - the sum of squares and the boxed contribution of Step 1, against the trace definition with $\mathbb M_s^{-1}$ solved as a matrix equation;
  - orthonormality, eigen-equation and traces of the coupled and decoupled generators;
  - proportionality of the bare null vector $(-D_{7a},D_{6a})$ to the decoupled generator, for symbolic complex $(c_q^A)_{ds}$.
- Exact rationals, $\epsilon=1/2$ and $2$, for CP-conserving, aligned-phase ($w=(3+4i)/5$), CP-violating and $d$-$s$-free test sources ($m_{a,0}^2=3/5$ GeV$^2$, $B_0=27/10$ GeV, $\operatorname{tr}S_2=1/20$ GeV):
  - $s-m_{K^0,0}^2$ divides $\det(sK-M^2)$ once (twice without $d$-$s$ sources);
  - the exact null space of $D(m_{K^0,0}^2)$ (`inverse_propagator` of [`solve_alp_meson_spectrum.py`](../../../scripts/solve_alp_meson_spectrum.py)) has that dimension and zero ALP components;
  - `solve_spectrum_numerical` finds the pole at $m_{K^0,0}^2$ with ALP component below $10^{-10}$.
- Reductions: $K_{7a}=M^2_{7a}=\operatorname{tr}(\lambda_7\mathcal C)=0$ for CP-conserving sources. For aligned phases, the rephasing makes the entries real, leaves $\det(sK-M^2)$ unchanged and gives $K_{6a}M^2_{7a}-K_{7a}M^2_{6a}=0$. For CP-violating sources that combination is nonzero, and $\det(s\mathbb 1-M^2_{\rm can})$ is unchanged when $\mathcal C$ is rephased to a real $d$-$s$ entry.

Units: $x_q$, $m_{K^0,0}^2$, $s$ in GeV$^2$; $(\mathcal C)_{ds}$ in GeV$^2$; couplings and $\epsilon$ dimensionless. Inputs are test values, not fits.

```bash
.venv/bin/python scripts/derive_alp_k0_doublet_exact_pole.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_k0_doublet_exact_pole.py
```

All checks pass.
