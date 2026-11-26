#!/usr/bin/env python3
"""Radial cubic vertices and the two-radial-line bubbles 6.1 and 6.7; GeV units.

Prints V(sigma; x, y) and V(sigma, sigma; x, y) for general pole legs x, y, the sigma^3 vertex
V_sigma3 = -6 lambda f_a + 6 B_0 f_pi^2 <S_30>/f_a^3, and the products of diagrams 6.1 and 6.7
with the rules of derive_radial_loop_generator_products (radial tree exchanges into (pi+, pi-)
and (a, pi0) with denominator m_sigma^2 - s_{+-}).  Symmetry and loop factors are excluded; no
loop integral is evaluated.

Run:  PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_six_seven.py
"""
import derive_radial_loop_generator_products as products
import derive_radial_loop_trace_contractions as common
from derive_radial_loop_three_one import GENERAL_DOTS, GENERAL_NAMES, leg_x, leg_y

LABELS = ('6.1', '6.7')
radial_vertex = common.vertex([None, leg_x, leg_y], 0, GENERAL_DOTS)
radial_pair_vertex = common.vertex([None, None, leg_x, leg_y], 0, GENERAL_DOTS)
radial_cubic = common.vertex([None, None, None])


def diagram_products():
    """[(label, unexpanded compact product)] for diagrams 6.1 and 6.7."""
    return [(label, products.diagram_total(label)) for label in LABELS]


if __name__ == '__main__':
    print(common.NOTATION + '\n')
    print('V(sigma; x, y):\n' + common.to_latex(radial_vertex, GENERAL_NAMES) + '\n')
    print('V(sigma, sigma; x, y):\n' + common.to_latex(radial_pair_vertex, GENERAL_NAMES) + '\n')
    print('V_sigma3:\n' + common.to_latex(radial_cubic) + '\n')
    products.print_diagrams(LABELS)
