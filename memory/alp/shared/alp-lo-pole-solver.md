---
title: Poles, pole vectors and residues of the LO ALP-meson system
status: mixed
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Poles, pole vectors and residues of the LO ALP-meson system

## Assumptions and inputs

- Fields $\Phi_I=(P_i,a)$, $\mathcal L_2=\tfrac12\partial_\mu\Phi^{\mathsf T}K\partial^\mu\Phi-\tfrac12\Phi^{\mathsf T}M^2\Phi$, $D(s)=sK-M^2$, as in the [finite-mixing record](finite-alp-meson-mixing.md) and [@Alda:2026nrd]. This record writes the kinetic mixing as $K_{ia}$; the finite-mixing record writes the same quantity as $-v_i$.
- $K$ and $M^2$ are real symmetric constant matrices and $K$ is positive definite (at LO see the [stability record](alp-lo-stability-conditions.md)). The solver needs no other property, so blocks with NLO local terms can be supplied unchanged. A $D(s)$ that is not linear in $s$ is not covered.
- LO inputs: the neutral $6\times6$ blocks of the [neutral-reduction record](alp-lo-neutral-block-reduction.md) (charged mesons decouple by the [decoupling record](alp-lo-charged-meson-decoupling.md)).
- The numerical inputs in the scripts are illustrative exact rationals, not fits.

## Method

1. Poles. The poles $m_n^2$ are the zeros of $\det D(s)$, equivalently the solutions of $M^2v=m^2Kv$.
2. Cholesky reduction. Write $K=LL^{\mathsf T}$. Then $\det D(s)=\det K\,\det(s-L^{-1}M^2L^{-\mathsf T})$, so the $m_n^2$ are the eigenvalues of the real symmetric matrix $L^{-1}M^2L^{-\mathsf T}$, which are real. If $y^{(n)}$ are its orthonormal eigenvectors, $v^{(n)}=L^{-\mathsf T}y^{(n)}$ satisfy $M^2v^{(n)}=m_n^2Kv^{(n)}$ and $v^{(m)\mathsf T}Kv^{(n)}=\delta_{mn}$.
3. Residues. Since $\sum_nv^{(n)}v^{(n)\mathsf T}=L^{-\mathsf T}\big(\sum_ny^{(n)}y^{(n)\mathsf T}\big)L^{-1}=K^{-1}$ (completeness) and $M^2=\sum_nm_n^2\,Kv^{(n)}v^{(n)\mathsf T}K$,
$$
D(s)^{-1}=\sum_n\frac{v^{(n)}v^{(n)\mathsf T}}{s-m_n^2},
\qquad \operatorname{Res}_{s=m_n^2}D^{-1}=v^{(n)}v^{(n)\mathsf T}\ \text{(simple pole)}.
$$
The composition of pole $n$ in $(P_i,a)$ is the list of components $v^{(n)}_I$.
4. Link to the $(w,1)$ normalization. The vector $(w_i^{(a)}(m_n^2),1)$ of the finite-mixing record is a null vector of $D(m_n^2)$ with ALP component $1$, so
$$
v^{(n)}=\pm\frac{(w^{(n)},1)}{\sqrt{N_a(m_n^2)}},\qquad
N_a(s)=K_{aa}+2K_{ia}w_i+K_{ij}w_iw_j,
$$
hence $w_i^{(n)}=v^{(n)}_i/v^{(n)}_a$, $(v^{(n)}_a)^2N_a(m_n^2)=1$, and $\big(\operatorname{Res}D^{-1}\big)_{aa}=1/N_a(m_n^2)$. The LSZ factor $N_a^{-1/2}$ of the finite-mixing record is $|v^{(n)}_a|$. The relation requires $v^{(n)}_a\neq0$ (for a pole with zero ALP component, such as the exact $K^0$ pole, $w^{(n)}$ is undefined).
5. Sign convention (a choice, not physics). Each $v^{(n)}$ is multiplied by $-1$ if needed so that its ALP component (the last field) is positive; if that component vanishes, the component of largest modulus is made positive.
6. Degenerate poles. For a degenerate pole the vectors are an arbitrary $K$-orthonormal basis of the null space of $D(m_n^2)$. Only the sum of their residues, $\sum_{n\in G}v^{(n)}v^{(n)\mathsf T}$, and the projector times $K$ are unique; individual compositions are not.
7. Invariance under a constant field change. For $\Phi=T\Phi'$ the [congruence record](alp-inverse-propagator-congruence.md) gives $D'=T^{\mathsf T}DT$. The poles are unchanged, $v^{(n)\prime}=\pm T^{-1}v^{(n)}$ for a simple pole (the sign is fixed by the convention of step 5 in the new coordinates), and $\operatorname{Res}D'^{-1}=T^{-1}(\operatorname{Res}D^{-1})T^{-\mathsf T}$, summed over a degenerate pole.

## Exact and numerical solvers

The numerical solver follows steps 2-3 in floating point. The exact solver uses SymPy: it finds the roots of $\det(sK-M^2)$, computes the null space of $D(m_n^2)$ exactly and $K$-orthonormalizes it by Gram-Schmidt. It does not form $L$ explicitly, because $L$ introduces square roots for rational $K$. It requires that all poles are found in closed form (for example rational poles) and raises an error otherwise. The positivity of $K$ is checked exactly with the leading principal minors.

## Verification

[`scripts/solve_alp_meson_spectrum.py`](../../../scripts/solve_alp_meson_spectrum.py) holds both solvers. [`scripts/verify_alp_lo_spectrum.py`](../../../scripts/verify_alp_lo_spectrum.py) checks, with exact rational data, random seeded symmetric $(K,M^2)$ and the LO neutral blocks (flavour-diagonal, CP-conserving and generic neutral sources): the null-vector equation, $K$-orthonormality, completeness, $D^{-1}(s)$ against a direct inverse, agreement of the ALP-like pole with $w^{(a)}(m_a^2)$ and $N_a$ (including the symbolic expressions of `derive_alp_finite_mixing.py`), invariance under random rational $T$ (three seeds, reusing the data of the congruence record), and agreement of the $6\times6$ system with the $10\times10$ system (whose spectrum is the union of the neutral and charged spectra). The ten-field exact $K^0$ pole is found with zero ALP and zero $(P_3,P_8,P_0)$ components.

```bash
.venv/bin/python scripts/solve_alp_meson_spectrum.py
.venv/bin/python scripts/verify_alp_lo_spectrum.py
```

The LO neutral spectrum for the illustrative inputs is obtained with (run from the repository root; the poles are in GeV$^2$ and `composition(4)` is the ALP-like pole):

```bash
PYTHONPATH=scripts .venv/bin/python -c "
import numpy as np
from derive_alp_lo_neutral_blocks import closed_form_blocks, illustrative_lo_inputs, NEUTRAL_FIELD_LABELS
from solve_alp_meson_spectrum import solve_spectrum_numerical
K, M2 = closed_form_blocks(**illustrative_lo_inputs('neutral'))
spectrum = solve_spectrum_numerical(np.array(K.evalf(30), dtype=float), np.array(M2.evalf(30), dtype=float), field_labels=NEUTRAL_FIELD_LABELS)
print([round(x, 6) for x in spectrum.pole_squared_masses]); print(spectrum.composition(4))"
```

For the illustrative neutral inputs the ALP-like pole is at $m^2\simeq0.769\ \mathrm{GeV}^2$ with $v_a\simeq0.643$ and a pole vector with $O(1)$ meson components; this point has large mixing and is not a physical prediction (see the [truncation record](alp-lo-truncation-limits.md)).
