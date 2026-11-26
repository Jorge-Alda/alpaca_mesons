---
title: Lowest-order weak meson-baryon Lagrangian
status: mixed
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Borasoy:2000ki, Alda:2026nrd]
provenance: analytical-derivation
---

# Lowest-order weak meson-baryon Lagrangian

The lowest-order non-leptonic weak meson-baryon Lagrangian in Borasoy:2000ki, Eq. (8), is written here with the notation requested for this project as

$$
\mathcal{L}_W^{(0)} = d_w \langle \bar{\bm{\mathcal{B}}} \{\bm{w}_+, \bm{\mathcal{B}}\}\rangle + f_w \langle \bar{\bm{\mathcal{B}}} [\bm{w}_+, \bm{\mathcal{B}}]\rangle .
$$

The matrices are

$$
\bm{w}_\pm = \bm{u} \bm{w}\bm{u}^\dagger \pm \bm{u} \bm{w}^\dagger \bm{u}^\dagger,
\qquad
(\bm{w})^a{}_b = \delta^a_2\delta^3_b.
$$

The weak couplings $d_w$ and $f_w$ correspond to $d$ and $f$ in Borasoy:2000ki. The trace $\langle\cdots\rangle$, baryon matrix $\bm{\mathcal{B}}$, and chiral field $\bm{u}$ follow Alda:2026nrd. The meson dependence is contained in $\bm{u}$; this expression is the meson-baryon contribution and does not specify a purely mesonic weak Lagrangian.

## Chiral convention and spurion transformation

Borasoy:2000ki uses $U\to g_L U g_R^\dagger$ and writes the corresponding weak matrix in Eq. (9) with $u^\dagger$ on the left. Alda:2026nrd instead uses $\bm{U}\to\bm{g_R}\bm{U}\bm{g_L}^\dagger$, with $\bm{U}=\bm{u}^2$ and, in Eq. (III.3),

$$
\bm{u}\to\bm{h}\bm{u}\bm{g_L}^\dagger
=\bm{g_R}\bm{u}\bm{h}^\dagger.
$$

Here $\bm{h}$ is the compensator of Alda:2026nrd. The weak matrix is treated as a left-handed spurion, with the transformation law

$$
\bm{w}\to\bm{g_L}\bm{w}\bm{g_L}^\dagger,
\qquad
\bm{w}^\dagger\to\bm{g_L}\bm{w}^\dagger\bm{g_L}^\dagger.
$$

Using the left-handed form of the transformation of $\bm{u}$ and unitarity of $\bm{g_L}$ gives, term by term,

$$
\bm{u}\bm{w}\bm{u}^\dagger
\to (\bm{h}\bm{u}\bm{g_L}^\dagger)
(\bm{g_L}\bm{w}\bm{g_L}^\dagger)
(\bm{g_L}\bm{u}^\dagger\bm{h}^\dagger)
=\bm{h}(\bm{u}\bm{w}\bm{u}^\dagger)\bm{h}^\dagger,
\qquad
\bm{w}_\pm\to\bm{h}\bm{w}_\pm\bm{h}^\dagger.
$$

Thus $\bm{w}_\pm$ transform like the baryon matter matrix. The original expression with $\bm{u}^\dagger$ on the left belongs to the opposite chiral convention; using it with the left-handed spurion and the conventions of Alda:2026nrd does not give this covariance. The transformation law applies to the spurion before it is set to the fixed physical matrix shown above. Once fixed, the weak interaction breaks the corresponding chiral symmetry explicitly.

The symbolic calculation in [weak_spurion_chiral_transformation.py](../../scripts/weak_spurion_chiral_transformation.py) implements these steps for both signs, treating $\bm{w}$ and $\bm{w}^\dagger$ as matrices with the same transformation law and inverses as adjoints of unitary matrices. Run it with `python scripts/weak_spurion_chiral_transformation.py`.
