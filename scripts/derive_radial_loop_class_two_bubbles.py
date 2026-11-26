#!/usr/bin/env python3
"""Compact flavour-trace products for the ALP--meson bubbles 2.1--2.5 (one ALP line and one meson line); GeV units.

Vertices, internal-line rules (meson line = 1_10 - v_a v_a^T), completeness and loop scale are
those of derive_radial_loop_generator_products.  Symmetry and loop factors are excluded; no loop
integral is evaluated.

Run:  PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_class_two_bubbles.py
"""
import derive_radial_loop_generator_products as products

LABELS = ('2.1', '2.2', '2.3', '2.4', '2.5')


def diagram_products():
    """[(label, unexpanded compact product)] for the diagrams of this script."""
    return [(label, products.diagram_total(label)) for label in LABELS]


if __name__ == '__main__':
    products.print_diagrams(LABELS)
