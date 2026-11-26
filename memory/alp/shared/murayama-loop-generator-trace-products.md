---
title: Compact generator-sum expressions for loop diagram products
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Compact generator-sum expressions for loop diagram products

## Inputs and assumptions

- Conventions, legs $(\lambda_x,\xi_x)$, vertex rules, retained diagrams 1.1–6.7, loop scale and sizing prescription: [estimate record](murayama-radial-up-only-estimate-result.md), following Alda:2026nrd and DiLuzio:2026poh. Units are GeV. Momenta are incoming, and every product of two momenta that contains the loop momentum is set to $m_a^2$ (diagrams with an ALP line) or $m_\sigma^2$ (all others).
- The sources $\widehat S_{rs}$, $\widehat P_{rs}$ and the pole matrices are general noncommuting matrices. Nothing is assumed about their commutators.
- The generators are normalized by $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$, with the singlet $\lambda_0=\sqrt{2/3}\,\mathbf 1$ included. The sum over the nine generators is defined by the two completeness identities
$$
\sum_i\langle A\lambda_iB\lambda_i\rangle=2\langle A\rangle\langle B\rangle,\qquad \sum_i\langle A\lambda_i\rangle\langle\lambda_iB\rangle=2\langle AB\rangle ,
$$
  for any matrices $A$, $B$. They are applied before any expansion, so no generator and no matrix component appears in the results.
- Internal lines. ALP line: the ALP pole leg $v_a$ at both ends. Radial line: $\sigma$ at both ends. Meson line: $\mathbf 1_{10}-v_av_a^T$ in the ten-field space of the nine meson fields and the bare ALP field, that is, the completeness term (nine generators), plus the bare angular term (ALP component 1 at both ends), minus the ALP pole term ($v_a$ at both ends). The subtraction removes the ALP pole, which the ALP line already carries.
- The radial tree exchanges $\sigma_{\pi^+\pi^-}$ and $\sigma_{a\pi^0}$ carry $V(\sigma;x,y)/(m_\sigma^2-s_{+-})$; both pairs have the same invariant $s_{+-}$ by momentum conservation.
- No loop integral is evaluated. The diagram list is not shown to be complete for general sources.

## Method

The vertex rules are written in compact trace form in the [vertex script](../../../scripts/derive_radial_loop_trace_contractions.py):

- potential part: the formula of the estimate record, with the sum over orderings of the meson matrices kept as one object, which for two matrices is the anticommutator;
- four-meson kinetic part, from $\mathcal L_4=\langle[\Phi,\partial_\mu\Phi][\Phi,\partial^\mu\Phi]\rangle/(48f_\pi^2)$: a sum over the three ways to split the four legs into two pairs, $-\frac{1}{24f_\pi^2}(k_1-k_2)\cdot(k_3-k_4)\,\langle[\lambda_1,\lambda_2][\lambda_3,\lambda_4]\rangle$;
- six-meson kinetic part, from $\mathcal L_6=\langle[\Phi,[\Phi,\partial_\mu\Phi]][\Phi,[\Phi,\partial^\mu\Phi]]\rangle/(1440f_\pi^4)$: for each pair of legs $c,f$ carrying the derivatives, $-k_c\cdot k_f/(720f_\pi^4)$ times the sum over the orderings of the other four meson matrices in the four open positions of $\langle[\cdot,[\cdot,\lambda_c]][\cdot,[\cdot,\lambda_f]]\rangle$. This form follows from writing $U^\dagger\partial_\mu U$ as a series of nested commutators; it is checked against the direct expansion below;
- ALP kinetic part and $V_{\sigma^3}=-6\lambda f_a+6B_0f_\pi^2\langle\widehat S_{30}\rangle/f_a^3$, the second term being the $r=3$ potential part.

Completeness is applied directly to these compact traces. A trace linear in a generator is written as $\langle Y\lambda_i\rangle$, using cyclicity, $\langle W[A,B]\rangle=\langle[W,A]B\rangle$, and the fact that the sum over orderings of $M_1,\dots,M_m,\lambda_i$ placed after $W$ gives $\langle Y\lambda_i\rangle$ with $Y$ the sum over orderings of $M_1,\dots,M_m,W$. Then $\sum_i\langle Y_1\lambda_i\rangle\langle\lambda_iY_2\rangle=2\langle Y_1Y_2\rangle$. Commutators and sums over orderings therefore survive the contraction.

The [assembly script](../../../scripts/derive_radial_loop_generator_products.py) holds the diagram list (label, symmetry factor, external legs at each vertex, internal lines) and multiplies the vertices around each loop, for every choice of the three meson-line terms. The [notation module](murayama-loop-trace-permutation-grouping.md) collects equal traces and prints them.

## Result

For each diagram 1.1–6.7 the output gives the factors shared by all meson-line terms once, followed by the remaining factors of the completeness, bare angular and ALP pole terms, each with its sign. Every factor is a sum of single or double traces of products of sources, external matrices, commutators and sums over orderings, with exact coefficients. Symmetry and loop factors are excluded. The per-class scripts print the same products for one class each; the external-leg corrections 7 and 8 are in the [residue record](murayama-radial-loop-external-residues.md).

## Reproduction and verification

Run `PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_generator_products.py`. All checks are exact symbolic identities with generic noncommuting symbols and exact rational coefficients.

- `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_trace_contractions`:
  - builds the nine generators, checks $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$, and checks both completeness identities exactly for two general $3\times3$ matrices of nine independent symbols each; the sum over the eight octet generators alone fails the first identity, which shows that the singlet is needed. This is the only check with explicit matrices, and it prints nothing;
  - compares the compact contraction with the word rule $\langle X\lambda_i\rangle\langle\lambda_iY\rangle\to2\langle XY\rangle$ for traces with sums over orderings, commutators and nested commutators;
  - compares the vertex rules for four general pole legs, six meson legs, and one or two radial legs with the direct expansion of the Lagrangian described below.
- Direct expansion, [expand_radial_loop_lagrangian.py](../../../scripts/expand_radial_loop_lagrangian.py). It does not use the vertex rules. It expands $\frac{f_\pi^2}{4}\langle\chi_+\rangle+\frac{f_\pi^2}{4}\langle\partial_\mu U\partial^\mu U^\dagger\rangle+\frac12(1+\sigma/f_a)^2(\partial a)^2-\lambda f_a\sigma^3$, with $\chi_+=u^\dagger\chi u^\dagger+u\chi^\dagger u$, $U=u^2$ and $u=\exp(i\Phi/(2f_\pi))$, as power series in noncommuting symbols, and keeps the coefficient linear in each leg. Each meson line is the pair of placeholder matrices $L\ldots L$ removed by the rewrites $\langle ALBL\rangle\to2\langle A\rangle\langle B\rangle$ and $\langle AL\rangle\langle LB\rangle\to2\langle AB\rangle$, plus the bare angular term, minus the ALP pole term. Each diagram is then a sum of products of traces of plain words.
- Each class test (`test_derive_radial_loop_class_*`, `test_derive_radial_loop_three_one`, `test_derive_radial_loop_six_seven`) expands the printed compact products back into words and requires exact equality with the direct expansion. Words are identified only up to cyclic rotation, with no Cayley–Hamilton relation and no reversal of words. All 44 diagrams agree, and a changed kinetic coefficient is detected.
- `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_generator_products` checks the number of diagrams per class and the symmetry factors, that no generator placeholder remains, that no coefficient is a floating-point number, and that the printed output contains no matrix component or generator label.
