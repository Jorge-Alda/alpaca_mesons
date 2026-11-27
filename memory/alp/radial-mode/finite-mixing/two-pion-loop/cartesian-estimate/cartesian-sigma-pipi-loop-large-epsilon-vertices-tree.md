---
title: Cartesian large epsilon at fixed sources - powers of the vertices by field-end pattern and of the tree in sigma to pi+ pi-
status: derived
created: 2026-11-27
updated: 2026-11-27
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Cartesian large $\epsilon$: vertices and tree

Task 2.3 of the OpenSpec change `estimate-sigma-pipi-one-loop-murayama-cartesian`. The Cartesian estimate is an independent order-of-magnitude estimate; neither it nor the polar estimate is authoritative. No loop integral, scalar loop function, loop coefficient or loop amplitude is evaluated. No polar vertex is used as input.

## Inputs and assumptions

- Exact vertex rules: [Cartesian vertex record](cartesian-sigma-pipi-loop-vertex-rules.md), with $\lambda$ eliminated through $m_\sigma^2$ (decision D-A3).
- Regime and fixed quantities (D10, assumptions): [Cartesian large-$\epsilon$ input record](cartesian-sigma-pipi-loop-large-epsilon-inputs.md). In particular, the sources, $\lambda$, $f_\pi$, $B_0$, $m_1$ and $m_P$ are fixed.
- Tree: [Cartesian frame and tree record](cartesian-sigma-pipi-loop-frame-and-tree.md), Step 6.

## Method

A field end is a radial end ($s$), an ALP-field end (the field $p$, with no meson component) or a meson end (a bare $U(3)$ generator). No pole vector, normalization or residue enters a field-end vertex. For each of the 46 patterns with $n=3,4,5$ ends, the script builds the vertex with $\lambda$ eliminated and applies the assumptions. It then inserts exactly
$$
m_\sigma^2=\frac{2\lambda f_\pi^2}{\epsilon^2}+\epsilon^2\kappa_\sigma,\qquad \langle\widehat S^C_{1,0}\rangle=\frac{\kappa_\sigma}{B_0}+2\langle\widehat S^C_{2,0}\rangle,\qquad \langle\widehat S^C_{0,2}\rangle=\frac{\kappa_\sigma-A}{2B_0}+\langle\widehat S^C_{2,0}\rangle ,
$$
so that $\lambda$ is fixed, and lists the powers of $\epsilon$ present over all independent products of traces. Momentum products are counted as $\epsilon^0$ for fixed momenta. For loop momenta they are counted as $\epsilon^2$, where $k\cdot k\sim m_{\rm term}^2$; this is an upper bound, since $m_{\rm term}=m_P$ in a term with only bare generator sums. The up-only variant applies decision D-A4 with $m_\sigma^2=2\lambda f_\pi^2/\epsilon^2+\epsilon^2B_0m_u$.

## Result: vertices by field-end pattern

**Rule.** The source part gives $\epsilon^{r+s}$: one power per radial end and per ALP-field end, none per meson end, because the sources depend on $s/f_a=\epsilon s/f_\pi$ and $p/f_a=\epsilon p/f_\pi$. The exceptions come from $V_\varphi$ and from the four-meson kinetic vertex:

| field ends | generic, fixed momenta | generic, loop momenta | up-only |
|---|---|---|---|
| $sss$ | $\epsilon^{-1}$, $\epsilon^3$ | same | $\epsilon^{-1}$ |
| $spp$ | $\epsilon^{-1}$, $\epsilon^3$ | same | $\epsilon^{-1}$ |
| $ssss$ | $\epsilon^0$, $\epsilon^4$ | same | $\epsilon^0$ |
| $sspp$ | $\epsilon^0$, $\epsilon^4$ | same | $\epsilon^0$ |
| $pppp$ | $\epsilon^0$, $\epsilon^4$ | same | $\epsilon^0$ |
| four mesons | $\epsilon^0$ | $\epsilon^0$, $\epsilon^2$ | $\epsilon^0$; loop momenta $\epsilon^0$, $\epsilon^2$ |

- $V(s,s,s)=-6\lambda f_\pi/\epsilon+6B_0\epsilon^3\langle\widehat S^C_{3,0}\rangle/f_\pi$ and $V(s;p,p)=-2\lambda f_\pi/\epsilon+2B_0\epsilon^3\langle\widehat S^C_{1,2}\rangle/f_\pi$. With $\lambda$ fixed, $-6\lambda f_a$ and $-2\lambda f_a$ are of order $\epsilon^{-1}$. The quartic $V_\varphi$ vertices $-6\lambda$, $-2\lambda$, $-6\lambda$ are of order $\epsilon^0$. The $\epsilon^3$ and $\epsilon^4$ terms come only from $\langle\widehat S^C_{r,s}\rangle$ with $r+s=n$.
- Three patterns vanish identically: three meson ends ($\widehat P^C_{0,0}=0$), $s$ with three meson ends ($\widehat P^C_{1,0}=0$), and five meson ends ($\widehat P^C_{0,0}=0$).
- Only the four-meson vertex depends on the momenta. No vertex has a kinetic $p$ end.
- **Up-only.** Only nine patterns survive: $sss$, $spp$ ($\epsilon^{-1}$, from $\lambda$); $ssss$, $sspp$, $pppp$ ($\epsilon^0$, from $\lambda$); $sMM$, $pMMM$, $sMMMM$ ($\epsilon^1$, from $\widehat S^C_{1,0}$, $\widehat P^C_{0,1}$, $\widehat S^C_{1,0}$); and $MMMM$ ($\epsilon^0$, from $\widehat S^C_{0,0}$ and the kinetic vertex). Here $M$ is a meson end. Every other pattern with two or more ends that are radial or $p$ vanishes, because $\widehat M=m_uP_u\varphi^*/f_a$ is linear in the fields.

## Result: tree

$$
M_0=-\frac{B_0\epsilon}{2f_\pi}\big\langle\widehat S^C_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\big\rangle\qquad(\text{power }\epsilon^1,\ \text{exact}),\qquad M_0^{\text{up-only}}=-\frac{B_0m_u\epsilon}{f_\pi}.
$$

## Verification

Exact symbolic checks: the pattern count (46); the exceptions to $\epsilon^{r+s}$ are exactly $sss$, $spp$, $ssss$, $sspp$, $pppp$; the $\lambda$ handling ($V(s,s,s)$ and $V(s;p,p)$ above, exact after the rewriting); the list of up-only patterns that do not vanish; only the four-meson vertex changes when loop momenta are counted; the tree power and the up-only tree.

## Reproduction

```sh
PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_large_epsilon_vertices.py
PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_large_epsilon_vertices
```

[Script](../../../../../../scripts/derive_sigma_pipi_cartesian_loop_large_epsilon_vertices.py) (sections 4 and 5), [tests](../../../../../../scripts/test_derive_sigma_pipi_cartesian_loop_large_epsilon_vertices.py).
