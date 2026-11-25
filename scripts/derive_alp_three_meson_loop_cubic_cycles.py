"""Extract the cubic mixing coefficient of an on-shell routed loop cycle.

Run with PYTHONPATH=scripts .venv/bin/python -m unittest
derive_alp_three_meson_loop_cubic_cycles. Euclidean momenta are used and
all masses are symbolic. No loop integrals are evaluated.
"""
import unittest
from unittest.mock import patch

import sympy as sp

import derive_alp_three_meson_loop_assembly as assembly
from derive_alp_loop_tensor_estimate import reduce_routed_numerator


def cubic_cycle(cycle, external_directions, external_momenta, pole_directions,
                pole_masses_squared, invariants, epsilon, external_gram, dimension,
                order=3):
    """Return a coefficient through cubic order grouped by raised denominators.

    Pole directions, masses, external momenta and the Gram matrix may all
    depend on epsilon. The returned tuple lists, for each distinct denominator
    exponent tuple, its Feynman-parameter combined mass and reduced scalar
    moment coefficients. Symmetry weights are supplied by the caller.
    """
    if order not in (0, 1, 2, 3):
        raise ValueError('The mixing order must be between zero and three.')
    number_of_momenta = external_gram.rows
    loop_squared = sp.Symbol('k^2')
    loop_products = sp.symbols(f'k_dot_p_0:{number_of_momenta}')
    full_gram = sp.Matrix([[loop_squared, *loop_products],
                           *[[loop_products[i], *list(external_gram.row(i))]
                             for i in range(number_of_momenta)]])

    def dot(first, second):
        return sp.expand((sp.Matrix(first).T * full_gram * sp.Matrix(second))[0])

    loop = (sp.S.One,) + (sp.S.Zero,) * number_of_momenta
    numerator, shifts = assembly.routed_cycle_numerator(
        cycle, external_directions, external_momenta, pole_directions, loop,
        dot, invariants, epsilon, euclidean=True)
    denominators = [
        dot(tuple(loop[index] + shift[index] for index in range(len(loop))),
            tuple(loop[index] + shift[index] for index in range(len(loop))))
        + sp.sympify(pole_masses_squared[pole])
        for shift, (_, pole) in zip(shifts, cycle)
    ]
    reciprocal = sp.symbols(f'g_0:{len(cycle)}')
    expression = sp.series(numerator, epsilon, 0, order + 1).removeO()
    for denominator, inverse in zip(denominators, reciprocal):
        difference = denominator - denominator.subs(epsilon, 0)
        expression *= sum((-difference * inverse)**power for power in range(order + 1)) * inverse
        expression = sp.series(expression, epsilon, 0, order + 1).removeO().expand()
    coefficient_at_order = sp.expand(expression).coeff(epsilon, order)
    polynomial = sp.Poly(coefficient_at_order, *reciprocal)
    gram_zero = external_gram.subs(epsilon, 0)
    parameters = sp.symbols(f'x_1:{len(cycle)}')
    weights = tuple(parameters) + (1 - sum(parameters),)
    zero_shifts = [sp.Matrix(shift[1:]).subs(epsilon, 0) for shift in shifts]
    average_shift = sum((weight * shift for weight, shift in zip(weights, zero_shifts)),
                        sp.zeros(number_of_momenta, 1))
    combined_mass = sum(
        weight * (sp.sympify(pole_masses_squared[pole]).subs(epsilon, 0)
                  + (shift.T * gram_zero * shift)[0])
        for weight, shift, (_, pole) in zip(weights, zero_shifts, cycle)
    ) - (average_shift.T * gram_zero * average_shift)[0]
    output = []
    for exponents, coefficient in polynomial.terms():
        if any(exponent < 1 for exponent in exponents) or coefficient == 0:
            continue
        moments = reduce_routed_numerator(
            coefficient, loop_squared, loop_products, tuple(average_shift),
            gram_zero, sp.expand(combined_mass), sum(exponents), dimension)
        if moments:
            output.append((exponents, sp.factor(combined_mass), parameters, moments))
    return output


class CubicCycleTest(unittest.TestCase):
    def test_all_orders(self):
        epsilon, mass, mass_shift, vertex_zero, vertex_one = sp.symbols(
            'epsilon m2 delta_m2 v0 v1')
        with patch.object(assembly, 'projected_vertex',
                          lambda legs, momenta, invariants, epsilon, dot, **kwargs:
                          vertex_zero + epsilon*vertex_one):
            values = [cubic_cycle((((1,), 'n'),), {1: None}, {1: (0, 0)}, {'n': None},
                                  {'n': mass + epsilon**2*mass_shift}, None, epsilon,
                                  sp.zeros(1), sp.Symbol('d'), order=order)
                      for order in range(4)]
        indexed = [{entry[0]: entry[3] for entry in result} for result in values]
        self.assertEqual(indexed[0], {(1,): {1: vertex_zero}})
        self.assertEqual(indexed[1], {(1,): {1: vertex_one}})
        self.assertEqual(indexed[2], {(2,): {2: -vertex_zero*mass_shift}})
        self.assertEqual(indexed[3], {(2,): {2: -vertex_one*mass_shift}})

    def test_projected_three_meson_tadpole(self):
        epsilon, mass, b0, decay_constant, pseudoscalar = sp.symbols(
            'epsilon m_a2 B_0 f_pi P_3')
        zero = sp.zeros(3)
        source = sp.diag(pseudoscalar, 0, 0)
        meson = sp.diag(1, 0, 0)
        invariants = {'S': [zero] * 7, 'P': [zero, zero, zero, source, zero, zero, zero],
                      'theta': [0] * 7, 'f_pi': decay_constant, 'B_0': b0,
                      'anomaly_mass_squared': sp.Symbol('m_0^2')}
        directions = {1: (zero, 1), 2: (meson, 0), 3: (meson, 0), 4: (meson, 0)}
        momenta = {1: (0, 1, 0, 0), 2: (0, 0, 1, 0), 3: (0, 0, 0, 1),
                   4: (0, -1, -1, -1)}
        result = cubic_cycle((((1, 2, 3, 4), 'a'),), directions, momenta,
                             {'a': (zero, 1)}, {'a': mass}, invariants, epsilon,
                             sp.eye(3), sp.Symbol('d'))
        self.assertEqual(len(result), 1)
        self.assertTrue(result[0][3][1].has(pseudoscalar, b0))
        self.assertFalse(result[0][3][1].has(epsilon))

    def test_external_vertex_and_mass_shift(self):
        epsilon, mass, mass_shift, vertex_one, vertex_three = sp.symbols(
            'epsilon m2 delta_m2 v1 v3')
        external = {1: None}
        momenta = {1: (0, 0)}
        poles = {'n': None}
        with patch.object(assembly, 'projected_vertex',
                          lambda legs, momenta, invariants, epsilon, dot, **kwargs:
                          epsilon*vertex_one + epsilon**3*vertex_three):
            result = cubic_cycle((((1,), 'n'),), external, momenta, poles,
                                 {'n': mass + epsilon**2*mass_shift}, None, epsilon,
                                 sp.zeros(1), sp.Symbol('d'))
        self.assertEqual({entry[0]: entry[3] for entry in result},
                         {(1,): {1: vertex_three}, (2,): {2: -vertex_one*mass_shift}})

    def test_routed_momentum_shift(self):
        epsilon, mass, momentum, correction = sp.symbols('epsilon m2 s delta_s')
        cycle = (((1,), 'n'), ((2,), 'n'))
        with patch.object(assembly, 'projected_vertex',
                          lambda legs, momenta, invariants, epsilon, dot, **kwargs: epsilon):
            result = cubic_cycle(cycle, {1: None, 2: None}, {1: (0, 1), 2: (0, -1)},
                                 {'n': None}, {'n': mass}, None, epsilon,
                                 sp.Matrix([[momentum + epsilon*correction]]), sp.Symbol('d'))
        self.assertTrue(any(sum(exponents) == 3 for exponents, _, _, _ in result))


if __name__ == '__main__':
    epsilon, mass, mass_shift, vertex_one, vertex_three = sp.symbols(
        'epsilon m2 delta_m2 v1 v3')
    with patch.object(assembly, 'projected_vertex',
                      lambda legs, momenta, invariants, epsilon, dot, **kwargs:
                      epsilon * vertex_one + epsilon**3 * vertex_three):
        for exponents, combined_mass, _, moments in cubic_cycle(
                (((1,), 'n'),), {1: None}, {1: (0, 0)}, {'n': None},
                {'n': mass + epsilon**2 * mass_shift}, None, epsilon,
                sp.zeros(1), sp.Symbol('d')):
            print('Denominator exponents:')
            print(sp.latex(sp.Tuple(*exponents)))
            print('Combined squared mass:')
            print(sp.latex(combined_mass))
            for power, coefficient in moments.items():
                print(f'Scalar denominator power {power}, coefficient:')
                print(sp.latex(coefficient))
