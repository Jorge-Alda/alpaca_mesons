---
title: Permutation grouping of loop source traces
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Permutation grouping of loop source traces

## Inputs and assumptions

The traces are the flavour traces of the vertex rules and diagram products of the [generator-sum record](murayama-loop-generator-trace-products.md), with the conventions of Alda:2026nrd and DiLuzio:2026poh. The matrices are noncommuting symbols; no commutator is assumed to vanish unless its two entries are equal. The notation changes no input, loop-size prescription or unit.

## Method and result

The [notation module](../../../scripts/group_radial_loop_trace_permutations.py) represents the diagram products with four exact objects:

- the trace $\langle X\rangle$, linear in $X$, invariant under cyclic rotation of the product of factors, and with $\langle\mathbf 1\rangle=3$. Each trace is stored in one fixed rotation, so equal traces combine;
- the sum over all $m!$ orderings of a list of matrices, with repeated entries counted once per position, as in $\sum_{\pi\in S_m}\lambda_{\pi_1}\cdots\lambda_{\pi_m}$ of the estimate record. For $m=2$ it is the anticommutator $\{X,Y\}$. A scalar entry multiplies the sum over the orderings of the other entries by $m$;
- the commutator $[X,Y]$, with $[X,X]=0$ and its two entries kept in a fixed order up to a sign;
- the sum over the orderings of a list of matrices placed at marked positions of a fixed trace word, used for the six-meson kinetic part.

Sources stay at fixed positions unless the completeness contraction moves them inside a sum over orderings, where they are permuted together with the meson matrices. The printed output uses $\left(\sum_{\rm perm}X_1\cdots X_m\right)$ for the sum over orderings and a bullet for the marked positions; each script prints this legend first. The module also expands any expression back into a sum of products of traces of plain words, each word in its smallest cyclic rotation. This expansion is the reference used by all checks.

## Reproduction and verification

Run `PYTHONPATH=scripts .venv/bin/python -m unittest test_group_radial_loop_trace_permutations`. The exact checks cover cyclic rotation and $\langle\mathbf 1\rangle=3$, the words of sums over orderings with distinct, repeated and scalar entries, the antisymmetry and words of commutators, the cyclic words of a product of two traces, and the expansion of a sum over marked positions. The diagram-level checks are in the [generator-sum record](murayama-loop-generator-trace-products.md).
