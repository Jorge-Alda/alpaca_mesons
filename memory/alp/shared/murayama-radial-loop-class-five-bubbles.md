---
title: Radial and ALP bubble diagram source-trace products
status: derived
created: 2026-11-26
updated: 2026-11-26
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: analytical-derivation
---

# Radial and ALP bubble diagram products

## Inputs and assumptions

Use the polar connection-free frame, incoming momenta and dimensional loop prescription of the [estimate record](murayama-radial-up-only-estimate-result.md), following Alda:2026nrd and DiLuzio:2026poh, and the compact vertex rules of the [generator-sum record](murayama-loop-generator-trace-products.md). The ALP-like pole legs contain both their meson matrices and angular components. Masses are in GeV. Internal ALP and radial lines use their pole legs; the radial exchange keeps its tree denominator $m_\sigma^2-s_{+-}$. The source coefficients are general noncommuting matrices within the retained graph selection. No loop integral is evaluated.

## Method and result

The [script](../../../scripts/derive_radial_loop_class_five_bubbles.py) takes diagrams 5.1–5.4 from the diagram list of the [assembly script](../../../scripts/derive_radial_loop_generator_products.py). Each vertex is the sum, over the assignments of every pole leg to its meson matrix or angular component, of the potential part, written with the sum over orderings of the meson matrices, plus the angular kinetic term when the vertex has two angular legs and one or two radial legs. There is no meson line, so each diagram is the product of its two vertex factors and, for 5.3 and 5.4, the charged-pion radial exchange. Products with the loop momentum are $m_a^2$. The output gives these factors in compact traces, without loop measure or symmetry factor.

## Reproduction and verification

Run `PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_class_five_bubbles.py`. The test `PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_class_five_bubbles` expands every printed product back into trace words and requires exact equality with the independent direct expansion of the Lagrangian described in the [generator-sum record](murayama-loop-generator-trace-products.md), in which each meson line is a pair of placeholder matrices removed by the completeness rewrites, plus the bare angular term, minus the ALP pole term. Words are compared up to cyclic rotation only. All coefficients are exact rationals times symbols. The check does not test diagrams outside the retained list, nor the sizing prescription as a loop integral.
