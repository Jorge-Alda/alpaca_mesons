#!/usr/bin/env python3
"""Derive polar field expansion and four-point exchange cancellation.

Natural units; all momenta incoming. Example:
    .venv/bin/python scripts/derive_polar_cartesian_four_point.py
The four-point calculation retains the full matrix propagator and ordered
row-vector, matrix, column-vector contractions for one channel.
"""
import sympy as sp


def derive():
    radial, pseudoscalar, decay_constant, order = sp.symbols(
        's p f_varphi t', real=True, positive=True)
    # A common auxiliary variable counts the total number of fields.
    polar_radial = sp.sqrt((decay_constant + order*radial)**2 +
                           (order*pseudoscalar)**2) - decay_constant
    polar_phase = decay_constant*sp.atan(order*pseudoscalar /
                                       (decay_constant + order*radial))
    polar_radial = sp.series(polar_radial, order, 0, 4).removeO().subs(order, 1).expand()
    polar_phase = sp.series(polar_phase, order, 0, 4).removeO().subs(order, 1).expand()

    # MatrixExpr preserves the order of the full mixed propagator factors.
    dimension = sp.Symbol('n', integer=True, positive=True)
    inverse_propagator = sp.MatrixSymbol('D', dimension, dimension)
    left_vertex = sp.MatrixSymbol('V_12', 1, dimension)
    right_vertex = sp.MatrixSymbol('V_34', dimension, 1)
    left_change = sp.MatrixSymbol('C_12', 1, dimension)
    right_change = sp.MatrixSymbol('C_34', dimension, 1)
    transformed_left = left_vertex + left_change*inverse_propagator
    transformed_right = right_vertex + inverse_propagator*right_change
    exchange = sp.expand(-transformed_left*inverse_propagator**(-1)*transformed_right)
    original_exchange = -left_vertex*inverse_propagator**(-1)*right_vertex
    contact_change = (left_change*right_vertex + left_vertex*right_change +
                      left_change*inverse_propagator*right_change)
    difference = sp.simplify(exchange - original_exchange + contact_change)
    return polar_radial, polar_phase, exchange, contact_change, difference


def main():
    labels = ('Polar radial field', 'Polar phase field', 'Transformed exchange',
              'Contact change', 'Difference of full channel')
    for label, expression in zip(labels, derive()):
        print(f'{label}:\n{sp.latex(expression)}\n')


if __name__ == '__main__':
    main()
