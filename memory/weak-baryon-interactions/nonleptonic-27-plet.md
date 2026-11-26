---
title: Non-leptonic 27-plet and baryon interactions
status: mixed
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, Salone:2026byq, Tandean:1999nv]
provenance: analytical-derivation
---

# Non-leptonic 27-plet and baryon interactions

This note explores the 27-plet without deciding whether to include it. The octet calculation remains in [non-leptonic operator matching](nonleptonic-operator-matching.md).

## Inputs and conventions

Work with three light flavors ordered $(u,d,s)$, natural units, and the left-handed rotation of Alda:2026nrd, $\bm g_L=\exp[ia(\bm\delta-\bm\kappa)/f_a]$. The flavor tensor $t^{ik}{}_{jl}$ and its symmetric, traceless 27-plet component $T^{ik}{}_{jl}$ follow the index conventions of Salone:2026byq, App. A; the symbol $T$ for the 27-plet tensor is used in Tandean:1999nv; repeated indices are summed. Symmetrization parentheses carry unit weight. The leading octet-baryon coupling $\beta_{27}$ is the convention of Tandean:1999nv. The example is one color-singlet current-current operator, not the full weak Hamiltonian.

## Partonic operator and 27-plet projection

The representative weak Lagrangian of Salone:2026byq, Eq. (12), is

$$
\mathcal L^{\Delta S=1}=\frac{4G_F}{\sqrt2}V_{ud}V_{us}^*(\bar u\gamma^\mu P_Ls)(\bar d\gamma_\mu P_Lu)+\mathrm{h.c.}
$$

With $t^{ik}{}_{jl}=\delta^i_1\delta^3_j\delta^k_2\delta^1_l$, its symmetric traceless 27-plet component is calculated as

$$
\begin{aligned}
T^{ik}{}_{jl}={}&t^{(ik)}{}_{(jl)}-\frac15\bigl(
\delta^i_j\,t^{(mk)}{}_{(ml)}+
\delta^i_l\,t^{(mk)}{}_{(mj)}+
\delta^k_j\,t^{(mi)}{}_{(ml)}+
\delta^k_l\,t^{(mi)}{}_{(mj)}\bigr).
\end{aligned}
$$

The first term is symmetric separately in $i,k$ and $j,l$. The subtraction removes all upper–lower traces. The nonzero components obtained from this normalized projection are

$$
\begin{aligned}
T^{12}{}_{13}=T^{12}{}_{31}=T^{21}{}_{13}=T^{21}{}_{31}&=\frac15,\\
T^{22}{}_{23}=T^{22}{}_{32}=T^{23}{}_{33}=T^{32}{}_{33}&=-\frac1{10}.
\end{aligned}
$$

All other components vanish. Inserting the projected tensor into the partonic Lagrangian gives, for this particular operator,

$$
\begin{aligned}
\left.\mathcal L^{\Delta S=1}\right|_{27}={}&
\frac{4G_F}{\sqrt2}V_{ud}V_{us}^*\,\frac15\bigg[
2(\bar u_L\gamma^\mu u_L)(\bar d_L\gamma_\mu s_L)
+2(\bar u_L\gamma^\mu s_L)(\bar d_L\gamma_\mu u_L)\\
&\hspace{22mm}-(\bar d_L\gamma^\mu d_L)(\bar d_L\gamma_\mu s_L)
-(\bar d_L\gamma^\mu s_L)(\bar s_L\gamma_\mu s_L)\bigg]+\mathrm{h.c.}
\end{aligned}
$$

Here the currents are color singlets and the products use the same flavor and color contraction as the starting operator. This is the projection of that operator, before QCD running and matching; it is not the complete low-energy 27-plet Lagrangian. The overall normalization is tied to this particular projection and is not a Wilson coefficient or baryon matrix element. This tensor includes both $\Delta I=1/2$ and $\Delta I=3/2$ pieces. For comparison, Salone:2026byq, App. A, gives a pure $\Delta I=3/2$ tensor with first four entries $1/2$, next two entries $-1/2$, and the final two zero; the $\Delta I=1/2$ tensor is given separately there. Isospin labels presume the isospin limit.

## Chiral field redefinition

Regard the 27-plet as a $9\times9$ matrix acting on the tensor product of two flavor spaces, with matrix elements $(\bm T)^{ik}{}_{jl}=T^{ik}{}_{jl}$. In the spurion convention of the [octet matching record](nonleptonic-operator-matching.md), its transformation is

$$
\bm T\to(\bm g_L\otimes\bm g_L)\,\bm T\,
(\bm g_L^\dagger\otimes\bm g_L^\dagger).
$$

To first order in the ALP field, this becomes

$$
\bm T\to\bm T+\frac{ia}{f_a}
\left[\bigl(\bm\delta-\bm\kappa\bigr)\otimes\bm 1
+\bm 1\otimes\bigl(\bm\delta-\bm\kappa\bigr),\bm T\right]
+\mathcal O\!\left(\frac{a^2}{f_a^2}\right).
$$

The transformation preserves the symmetry and trace conditions. For diagonal $\bm\delta-\bm\kappa$, every nonzero component above obtains the same factor $1+ia[(\delta-\kappa)_{dd}-(\delta-\kappa)_{ss}]/f_a$. The up-quark phases cancel in the full four-quark tensor. For non-diagonal rotations the tensor can acquire other components, and the full tensor-product formula is needed.

## Leading octet-baryon interaction

Tandean:1999nv gives a leading heavy-baryon operator for $|\Delta I|=3/2$ with an independent coupling $\beta_{27}$. In the chiral convention of Alda:2026nrd, its octet-baryon structure can be written using the 27-plet tensor above as

$$
\mathcal L_{W,27}^{(0)}=\beta_{27}\,T^{ik}{}_{jl}\,
(\bm u^\dagger\bar{\bm{\mathcal B}}\bm u)^j{}_i
(\bm u^\dagger\bm{\mathcal B}\bm u)^l{}_k+\mathrm{h.c.}
$$

Equivalently, in tensor-product notation,

$$
\mathcal L_{W,27}^{(0)}=\beta_{27}\,\operatorname{Tr}_{\!\mathrm{flavor}\otimes\mathrm{flavor}}\!\left[\bm T\left((\bm u^\dagger\bar{\bm{\mathcal B}}\bm u)\otimes(\bm u^\dagger\bm{\mathcal B}\bm u)\right)\right]+\mathrm{h.c.}
$$

The baryon spinor indices are contracted between the two factors. Since $\bm u\to\bm h\bm u\bm g_L^\dagger$ and $\bm{\mathcal B}\to\bm h\bm{\mathcal B}\bm h^\dagger$, each dressed baryon transforms by $\bm g_L(\cdots)\bm g_L^\dagger$; the contraction is invariant when the tensor is transformed as above. The tensor may represent either isospin component; Tandean:1999nv specifically studies $|\Delta I|=3/2$. The displayed partonic tensor and the baryon source tensor must be matched with their flavor-index orientations, CKM convention and normalization specified; one must not identify their numerical coefficients from this symmetry argument. In Tandean:1999nv, $\beta_{27}$ already includes the normalization of the weak interaction. A decuplet-baryon operator with a separate coupling is also given in Tandean:1999nv; it is not included in the octet-only expression above.

The symmetry-based form does not determine $\beta_{27}$ from $d_w$ and $f_w$, nor does it establish a general ALP-dependent weak coefficient. Independent ALP-dependent 27-plet sources can require additional matching coefficients and baryon low-energy constants. ALP emission may arise from the transformed spurion and meson mixing, and baryon-pole graphs must also be considered for a physical amplitude. The expression above is an available operator, not a decision to include the 27-plet.

## Reproduction and checks

[nonleptonic_27_flavour.py](../../scripts/nonleptonic_27_flavour.py) constructs the current-current tensor, symmetrizes it, subtracts traces, and calculates its first-order rotation. Run `.venv/bin/python scripts/nonleptonic_27_flavour.py`; `--test` checks the 27-plet symmetries, tracelessness, and rotation. The calculation does not perform QCD evolution or hadronic matching.
