#!/usr/bin/env python3
"""Quartic vertex with two internal ALP pole legs and the class-3 products 3.1--3.3; GeV units.

The vertex V(a_loop, a_loop, x, y) is printed for two general pole legs x = (lambda_x, xi_x) and
y = (lambda_y, xi_y) with incoming momenta k_x, k_y; products with a loop momentum are m_a^2.
It contains the potential part and the four-meson kinetic part in commutator form.  The class-3
products use the rules of derive_radial_loop_generator_products.  Symmetry and loop factors are
excluded; no loop integral is evaluated.

Run:  PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_three_one.py
"""
import sympy as sp

import derive_radial_loop_generator_products as products
import derive_radial_loop_trace_contractions as common

lambda_x, lambda_y = sp.symbols('lambda_x lambda_y', commutative=False)
xi_x, xi_y = sp.symbols('xi_x xi_y', real=True)
k_xy = sp.Symbol('k_x.k_y', real=True)
leg_x = (lambda_x, xi_x, 'x')
leg_y = (lambda_y, xi_y, 'y')
GENERAL_DOTS = {frozenset(('x', 'y')): k_xy}
GENERAL_NAMES = {lambda_x: r'\lambda_x', lambda_y: r'\lambda_y', xi_x: r'\xi_x', xi_y: r'\xi_y', k_xy: r'k_x\cdot k_y'}
LABELS = ('3.1', '3.2', '3.3')

quartic_vertex = common.vertex([common.internal, common.internal, leg_x, leg_y], common.ma**2, GENERAL_DOTS)


def diagram_products():
    """[(label, unexpanded compact product)] for the class-3 diagrams."""
    return [(label, products.diagram_total(label)) for label in LABELS]


if __name__ == '__main__':
    print(common.NOTATION + '\n')
    print('Quartic vertex V(a_loop, a_loop, x, y):\n' + common.to_latex(quartic_vertex, GENERAL_NAMES) + '\n')
    products.print_diagrams(LABELS)
