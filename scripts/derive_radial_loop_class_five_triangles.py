#!/usr/bin/env python3
"""Compact flavour-trace products for the radial--ALP triangles 5.5--5.12 (one radial, one ALP and one meson line); GeV units.

Vertices, internal-line rules (meson line = 1_10 - v_a v_a^T), completeness and loop scale are
those of derive_radial_loop_generator_products.  Symmetry and loop factors are excluded; no loop
integral is evaluated.

Run:  PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_class_five_triangles.py
"""
import derive_radial_loop_generator_products as products

LABELS = ('5.5', '5.6', '5.7', '5.8', '5.9', '5.10', '5.11', '5.12')


def diagram_products():
    """[(label, unexpanded compact product)] for the diagrams of this script."""
    return [(label, products.diagram_total(label)) for label in LABELS]


if __name__ == '__main__':
    products.print_diagrams(LABELS)
