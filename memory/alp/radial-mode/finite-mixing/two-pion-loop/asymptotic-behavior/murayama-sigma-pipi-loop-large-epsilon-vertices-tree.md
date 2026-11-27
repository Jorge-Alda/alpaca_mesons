---
title: Large epsilon at fixed sources - powers of the vertices by field-end pattern and of the tree in sigma to pi+ pi-
status: derived
created: 2026-11-26
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Large $\epsilon$: vertices and tree

OpenSpec change `estimate-sigma-pipi-one-loop-murayama` (task 2.3). Assumptions (fixed sources, $\lambda$ fixed, $m_\sigma^2=2\lambda f_\pi^2/\epsilon^2+\epsilon^2\kappa_\sigma$): [large-$\epsilon$ input record](murayama-sigma-pipi-loop-large-epsilon-inputs.md). Exact vertex rules: [vertex-rule record](../one-loop-estimate/murayama-sigma-pipi-loop-vertex-rules.md). No loop is evaluated.

## Vertices by field-end pattern

A field end is one of three kinds: a radial end ($\sigma$), an ALP-field end (the field $a$, with no meson component), or a meson end (a bare $U(3)$ generator). No pole vector, no normalization $N(m_a^2)$ and no residue enters here; the leg factors that attach to the ends are in the [large-$\epsilon$ input record](murayama-sigma-pipi-loop-large-epsilon-inputs.md). The script lists the powers of $\epsilon$ present in each of the 46 patterns with $n=3,4,5$ legs, with $\lambda$ fixed. Two variants are given: generic fixed sources, and "$r\le1$ sources" ($\widehat S_{r,s}=\widehat P_{r,s}=0$ for $r\ge2$, a letter substitution that contains the up-only model).

**Rule.** Every pattern carries one power of $\epsilon$ per radial end and per ALP-field end, and none per meson end, because the source terms depend on $\sigma/f_a=\epsilon\sigma/f_\pi$ and $a/f_a=\epsilon a/f_\pi$. The exceptions are:

| field ends | generic, fixed momenta | generic, loop momenta $k\cdot k\sim\epsilon^2$ | $r\le1$ sources |
|---|---|---|---|
| $\sigma\sigma\sigma$ | $\epsilon^{-1}$, $\epsilon^3$ | same | $\epsilon^{-1}$ |
| $\sigma\sigma\sigma\sigma$ | $\epsilon^0$, $\epsilon^4$ | same | $\epsilon^0$ |
| $\sigma aa$ | $\epsilon$, $\epsilon^3$ | $\epsilon^3$ | as generic |
| $\sigma\sigma aa$ | $\epsilon^2$, $\epsilon^4$ | $\epsilon^4$ | $\epsilon^2$ (fixed momenta), $\epsilon^4$ (loop momenta) |
| four mesons | $\epsilon^0$ | $\epsilon^0$, $\epsilon^2$ | as generic |

- $V(\sigma,\sigma,\sigma)=-6\lambda f_\pi/\epsilon+6B_0\epsilon^3\langle\widehat S_{3,0}\rangle/f_\pi$ and $V(\sigma^4)=-6\lambda+24B_0\epsilon^4\langle\widehat S_{4,0}\rangle/f_\pi^2$. With $\lambda$ fixed, the term $-6\lambda f_a$ is of order $\epsilon^{-1}$. The $\epsilon^3$ ($\epsilon^4$) term comes only from $\langle\widehat S_{3,0}\rangle$ ($\langle\widehat S_{4,0}\rangle$) and vanishes for $r\le1$ sources.
- In the kinetic $\sigma aa$ and $\sigma\sigma aa$ vertices the momentum products are external products ($\sim m_\sigma^2\sim\epsilon^2$, or $m_{\pi^\pm}^2$) or loop products ($m_{\rm term}^2$, the per-term scale of D3: $\sim\epsilon^2$ when the term contains a radial line, an ALP-like line or pieces (ii)–(iv) of a meson line, and $m_P^2\sim\epsilon^0$ only for a term in which every line is a bare generator sum; the kinetic vertices have ALP-field ends, so in the terms where they act $m_{\rm term}\ge m_a\sim\epsilon$). For vertices with external pion legs the "loop momenta" column is an upper bound.
- Three patterns vanish identically: three meson ends ($\widehat P_{0,0}=0$), $\sigma$ with three meson ends ($\widehat P_{1,0}=0$), five meson ends.
- For $r\le1$ sources every pattern with two or more radial ends vanishes, except $\sigma\sigma\sigma$, $\sigma\sigma\sigma\sigma$ (through $\lambda$) and $\sigma\sigma aa$ (kinetic).

In the kinetic terms $\frac12(1+\sigma/f_a)^2(\partial a)^2$ only the radial ends carry $\epsilon$, since $a$ is canonically normalized: $C(\sigma)=2\epsilon/f_\pi$, $C(\sigma,\sigma)=2\epsilon^2/f_\pi^2$.

## Tree

$$
M_0=-\frac{B_0\epsilon}{2f_\pi}\big\langle\widehat S_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\big\rangle\qquad(\text{leading power }\epsilon^1,\ \text{exact}) .
$$
It is the exact tree of the [tree and legs record](../one-loop-estimate/murayama-sigma-pipi-loop-tree-and-legs.md): one radial end, fixed sources. For the up-only sources $M_0=-B_0m_u\epsilon/f_\pi$.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_large_epsilon.py --markdown
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_large_epsilon
```

[Script](../../../../../../scripts/derive_sigma_pipi_loop_large_epsilon.py), [tests](../../../../../../scripts/test_derive_sigma_pipi_loop_large_epsilon.py) (powers per pattern, the $\lambda$ handling: $V(\sigma,\sigma,\sigma)\sim\epsilon^{-1}$ for $r\le1$ sources and $\epsilon^3$ from $\langle\widehat S_{3,0}\rangle$ only for generic sources).
