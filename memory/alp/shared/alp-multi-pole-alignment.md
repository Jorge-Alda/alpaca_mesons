---
title: Alignment of several ALP-meson poles
status: derived
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Alignment of several ALP-meson poles

## Inputs and assumptions

- $D(s)=sK-M^2$, $D'=T^{\mathsf T}DT$, pole vectors $v^{(k)}$ and vertex transformation $\mathcal V'=T^{\mathsf T}\mathcal V$: [congruence record](alp-inverse-propagator-congruence.md).
- Single-pole shift and its regularity: [single-pole record](alp-single-pole-alignment.md); $w_i^{(k)}$ as in the [finite-mixing record](finite-alp-meson-mixing.md) and `Alda:2026nrd` [@Alda:2026nrd].
- $K$ real symmetric positive definite, $M^2$ real symmetric, poles $m_k^2$ distinct unless stated.

**Definition.** Pole $k$ is aligned with coordinate $J$ of $\Phi'$ when column $J$ of $D'(m_k^2)$ vanishes. Since $D'(m_k^2)e_J=T^{\mathsf T}D(m_k^2)Te_J$, this holds exactly when column $J$ of $T$ is proportional to $v^{(k)}$. For the ALP coordinate and a shift this is $D'_{i\hat a}(m_k^2)=0$.

## Shift-type maps $P_i=P_i'+w_ia'$

1. For regular poles, alignment of pole $k$ fixes $w=w^{(k)}(m_k^2)$. Two poles are aligned by one shift if and only if $w^{(k)}(m_k^2)=w^{(l)}(m_l^2)$.
2. With $K$ positive definite and $m_k^2\neq m_l^2$ this never happens: a common $(w,1)$ would satisfy $D(m_k^2)(w,1)^{\mathsf T}-D(m_l^2)(w,1)^{\mathsf T}=(m_k^2-m_l^2)K(w,1)^{\mathsf T}=0$, contradicting $K>0$. The argument does not use regularity, so it also covers degenerate poles. A shift aligns at most one pole with $a'$.
3. A shift cannot make $D(s)$ diagonal at all $s$. It leaves $K_{PP}$ and $M^2_{PP}$ unchanged, and removing both mixed blocks requires $w=-K_{PP}^{-1}K_{P\hat a}$ and $w=-M_{PP}^{2\,-1}M^2_{P\hat a}$ at once. At LO the first is the canonical shift $w_i=(f_\pi/f_a)\operatorname{tr}(c_q^A\lambda_i)$ of the finite-mixing record; the second differs in general.

## General $T$

1. Pole vectors of distinct poles are $K$-orthogonal: $v^{(k)\mathsf T}M^2v^{(l)}=m_l^2\,v^{(k)\mathsf T}Kv^{(l)}=m_k^2\,v^{(k)\mathsf T}Kv^{(l)}$.
2. Take $T$ with columns $v^{(k)}/\sqrt{v^{(k)\mathsf T}Kv^{(k)}}$. Then $T^{\mathsf T}KT=\mathbf 1$, $T^{\mathsf T}M^2T=\operatorname{diag}(m_k^2)$ and
$$
D'(s)=\operatorname{diag}(s-m_k^2)\quad\text{for every }s .
$$
This is the simultaneous diagonalization of the pair $(K,M^2)$. For $K>0$ it always exists (write $K=LL^{\mathsf T}$ and diagonalize the symmetric matrix $L^{-1}M^2L^{-\mathsf T}$ by an orthogonal matrix); for a degenerate pole choose a $K$-orthonormal basis of its null space. The primed fields are then the mass eigenstates, not the meson coordinates of `Alda:2026nrd`; vertices must be transformed with $\mathcal V'=T^{\mathsf T}\mathcal V$ before comparing amplitudes.
3. Aligning only some poles: the aligned coordinates are decoupled from each other at every $s$, because $D'(s)e_{J_k}=(s-m_k^2)T^{\mathsf T}Kv^{(k)}$ for $Te_{J_k}=v^{(k)}$ and the $v^{(k)}$ are $K$-orthogonal. They still mix off the poles with the remaining coordinates, unless the remaining columns of $T$ are chosen $K$-orthogonal to the aligned pole vectors.

## Check of design decision 4

The change design states that a general $T$ built from all pole vectors diagonalizes $D(s)$ at every $s$ "only if the residue structure allows it". The check gives:

- For $K$ positive definite, the statement holds without further conditions (item 2 above).
- The extra condition matters only if $K$ is not positive definite, which is outside the physical case. Example: $K=\begin{pmatrix}0&1\\1&0\end{pmatrix}$, $M^2=\operatorname{diag}(1,0)$ gives $\det D=-s^2$, a pole vector $(0,1)$ with $v^{\mathsf T}Kv=0$, and a $1/s^2$ term in $D^{-1}$. A diagonal $D'(s)$ with entries linear in $s$ has only simple poles in its inverse, and $D'^{-1}=T^{-1}D^{-1}T^{-\mathsf T}$ keeps the double pole, so no $T$ diagonalizes $D(s)$.
- A singular $K$ (a coordinate with no kinetic term) is not covered.

## Cases where alignment fails

- Two or more distinct poles with one shift-type map, always (for $K>0$).
- A single pole with a shift, when its pole vectors have no ALP component (see the single-pole record).
- Diagonalization of $D(s)$ at all $s$ by any shift, in general.
- Full diagonalization by a general $T$, when $K$ is not positive definite and $D^{-1}$ has a higher-order pole.

A general $T$ is a linear field redefinition. It is not the linear part of a chiral redefinition with $g_{R,L}(0)=\mathbf 1$: such a redefinition has a shift-type linear part only, see the [chiral-generator record](alp-pole-alignment-chiral-generator.md). Amplitudes are still unchanged under a general $T$ (congruence record).

## Verification

[`derive_alp_multi_pole_alignment.py`](../../../scripts/derive_alp_multi_pole_alignment.py) uses the exact rational $N_f=3$ data of the congruence record. For each of the six pairs of poles it checks that the $w^{(k)}$ differ, the stacked alignment equations are inconsistent (rank 3, augmented rank 4), $D(m_k^2)-D(m_l^2)=(m_k^2-m_l^2)K$, and that the shift aligning one pole leaves the other misaligned. It checks that the shifts removing $K_{P\hat a}$ and $M^2_{P\hat a}$ differ, $K$-orthogonality, $T^{\mathsf T}KT=\mathbf 1$, $T^{\mathsf T}M^2T=\operatorname{diag}(m_k^2)$, a diagonal $D'(s)$, the partial-alignment statements, and the indefinite-$K$ counterexample. A floating-point check with `scipy.linalg.eigh` for generic $K>0$ and symmetric $M^2$ with irrational poles gives deviations below $2\times10^{-15}$. All checks pass.

```sh
.venv/bin/python scripts/derive_alp_multi_pole_alignment.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_multi_pole_alignment.py
```
