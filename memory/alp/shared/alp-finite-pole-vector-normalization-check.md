---
title: Check of the finite-mixing pole-vector normalization in trace form
status: derived
created: 2026-11-20
updated: 2026-11-20
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Pole-vector normalization $v^{(r)T}Kv^{(r)}=1$ in trace form

## Conventions and units

Natural units. $\langle\bm A\rangle=\operatorname{tr}\bm A$. The $U(3)$ basis is $\bm{\lambda_0}=\sqrt{2/3}\,\bm 1$ and the Gell-Mann matrices $\bm{\lambda_1},\dots,\bm{\lambda_8}$, with $\langle\bm{\lambda_i}\bm{\lambda_j}\rangle=2\delta_{ij}$, as in `Alda:2026nrd` [@Alda:2026nrd]. $f_\pi$ and $f_a$ have mass dimension one; only $f_\pi/f_a$ enters. $\bm{c_q^A}$ and $\mathbb M_a^{-1}(\bm{\mathcal C})$ are dimensionless. No expansion in $f_\pi/f_a$ is made.

## Inputs

- Kinetic blocks $K_{ij}=\delta_{ij}$, $K_{i\hat a}=-(f_\pi/f_a)\langle\bm{c_q^A}\bm{\lambda_i}\rangle$, $K_{\hat a\hat a}=1+(f_\pi^2/f_a^2)\sum_i\langle\bm{c_q^A}\bm{\lambda_i}\rangle^2$: [finite ALP–meson mixing](finite-alp-meson-mixing.md), eq:finite-kinetic-blocks of `latex/finite_mixing.tex`.
- Pole vector $v_i^{(r)}=f_\pi\langle\bm{\lambda_i}[\bm{c_q^A}+\mathbb M_a^{-1}(\bm{\mathcal C})]\rangle/(f_a\sqrt N)$, $v_{\hat a}^{(r)}=1/\sqrt N$, and $N(s)=1+(2f_\pi^2/f_a^2)\langle[\mathbb M_a^{-1}(\bm{\mathcal C})]^2\rangle$: [operator inverse and normalization](alp-lo-finite-mixing-operator-inverse.md), eqs. finite-invariant-daughter-pole-vector and finite-pole-norm.
- $\bm{\lambda_r}=\sum_iv_i^{(r)}\bm{\lambda_i}$: eq:finite-field-expansion.
- **Assumption.** $\bm{c_q^A}$ and $\bm X\equiv\mathbb M_a^{-1}(\bm{\mathcal C})|_{s=m_r^2}$ are Hermitian and lie in the span of the generators kept as meson fields. The check treats $\bm X$ as a generic Hermitian matrix; nothing about its explicit form is used.

## Method

The script [check_finite_pole_vector_normalization.py](../../../scripts/check_finite_pole_vector_normalization.py) writes $\bm{c_q^A}=\tfrac12\sum_ic_i\bm{\lambda_i}$ and $\bm X=\tfrac12\sum_iX_i\bm{\lambda_i}$ with real symbols, builds $K$ as an explicit matrix and simplifies with `sympy`. This is done for two generator sets: the full $U(3)$ basis, and the generators commuting with the quark charge, $\{0,3,6,7,8\}$ ($\pi^0$, $\eta_8$, $\eta_0$ and the $K^0$, $\bar K^0$ directions).

1. For arbitrary real $v_i$, $v_{\hat a}$, the completeness relation $\sum_i\langle\bm A\bm{\lambda_i}\rangle^2=2\langle\bm A^2\rangle$ for Hermitian $\bm A$ in the span gives
$$
v^TKv=\tfrac12\langle\bm{\lambda_r}^2\rangle-\frac{2f_\pi}{f_a}v_{\hat a}\langle\bm{c_q^A}\bm{\lambda_r}\rangle+v_{\hat a}^2\left[1+\frac{2f_\pi^2}{f_a^2}\langle(\bm{c_q^A})^2\rangle\right].
$$
With $v_{\hat a}=1/\sqrt N$ this is the left side of eq:finite-pole-vector-normalization.
2. Inserting $v_i$ gives $\bm{\lambda_r}=2f_\pi[\bm{c_q^A}+\bm X]/(f_a\sqrt N)$, as stated in eq:finite-field-expansion.
3. Keeping $N$ as a free symbol and substituting step 2,
$$
N\times[\text{left side}]=\frac{2f_\pi^2}{f_a^2}\left(\langle(\bm c+\bm X)^2\rangle-2\langle\bm c(\bm c+\bm X)\rangle+\langle\bm c^2\rangle\right)+1=1+\frac{2f_\pi^2}{f_a^2}\langle\bm X^2\rangle,
$$
with $\bm c=\bm{c_q^A}$. All dependence on $\bm{c_q^A}$ cancels. The normalization therefore holds exactly when $N$ is given by eq:finite-pole-norm.
4. For $\bm{c_q^A}=0$ the identity reduces to $\tfrac12\langle\bm{\lambda_r}^2\rangle+1/N=1$.
5. As an independent check, $v^TKv-1$ is computed directly from the matrix $K$ and the vector $v$ with $N$ from eq:finite-pole-norm.

## Result

Eq:finite-pole-vector-normalization holds identically and is consistent with eqs. finite-pole-norm and finite-field-expansion, for generic Hermitian $\bm{c_q^A}$ (including its trace part) and generic Hermitian $\bm X$. The residuals of steps 1–5 are zero in both generator sets.

**Caveat (derived).** Step 1 needs the completeness relation in the set of meson fields that are kept. If $\bm{c_q^A}$ had components along the charged generators $\bm{\lambda_{1,2,4,5}}$ but only the neutral fields were kept, the residual would be
$$
v^TK v-[\text{left side}]=-\frac{f_\pi^2}{f_a^2}v_{\hat a}^2\left(c_1^2+c_2^2+c_4^2+c_5^2\right).
$$
This case is excluded by the statement in `latex/finite_mixing.tex` that the sources commute with the charge matrix.

## Reproduction and verification

```bash
.venv/bin/python scripts/check_finite_pole_vector_normalization.py
cd scripts && ../.venv/bin/python -m unittest test_check_finite_pole_vector_normalization
```

The unit test [test_check_finite_pole_vector_normalization.py](../../../scripts/test_check_finite_pole_vector_normalization.py) checks the generator normalization, steps 1–5 in both generator sets, that $N$ drops out of $N\times$[left side], and the caveat residual.
