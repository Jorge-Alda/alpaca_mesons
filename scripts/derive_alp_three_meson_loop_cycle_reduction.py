"""Route and reduce a projected loop cycle without evaluating scalar integrals.

Run: .venv/bin/python -m unittest discover -s scripts -p derive_alp_three_meson_loop_cycle_reduction.py
Momenta and pole masses are Euclidean, with masses in GeV.
"""
import unittest
from unittest.mock import patch
from math import factorial

import sympy as sp

import derive_alp_three_meson_loop_assembly as assembly
from derive_alp_loop_tensor_estimate import reduce_routed_numerator


def reduce_cycle(cycle, external_directions, external_momenta, pole_directions,
                 pole_masses_squared, invariants, epsilon, external_gram,
                 feynman_parameters, dimension, denominator_exponents=None):
    """Construct the routed vertices and reduce their loop momentum numerator.

    Momentum tuples use the loop momentum as their first basis vector. The
    remaining components refer to the Gram matrix of independent momenta.
    The supplied Feynman parameters sum to one. The result is a dictionary
    of scalar denominator powers and their parameter-dependent coefficients.
    The caller supplies the simplex measure and cycle symmetry weight. Optional
    positive integer exponents support repeated propagators.
    """
    number_of_momenta = external_gram.rows
    if external_gram.cols != number_of_momenta or len(feynman_parameters) != len(cycle):
        raise ValueError('Gram matrix and Feynman parameters must match the cycle.')
    if sp.simplify(sum(feynman_parameters) - 1) != 0:
        raise ValueError('Feynman parameters must sum to one.')
    if denominator_exponents is None:
        denominator_exponents = (1,) * len(cycle)
    if len(denominator_exponents) != len(cycle) or any(
            not isinstance(value, int) or value < 1 for value in denominator_exponents):
        raise ValueError('Each propagator requires a positive integer exponent.')
    loop_squared = sp.Symbol('k^2')
    loop_products = sp.symbols(f'k_dot_p_0:{number_of_momenta}')
    full_gram = sp.Matrix([[loop_squared, *loop_products],
                           *[[loop_products[i], *list(external_gram.row(i))]
                             for i in range(number_of_momenta)]])

    def dot(first, second):
        return sp.expand((sp.Matrix(first).T * full_gram * sp.Matrix(second))[0])

    loop_momentum = (sp.S.One,) + (sp.S.Zero,) * number_of_momenta
    numerator, shifts = assembly.routed_cycle_numerator(
        cycle, external_directions, external_momenta, pole_directions,
        loop_momentum, dot, invariants, epsilon, euclidean=True)
    if any(sp.simplify(shift[0]) != 0 for shift in shifts):
        raise ValueError('External momenta must have no loop component.')
    spatial_shifts = [sp.Matrix(shift[1:]) for shift in shifts]
    average_shift = sum((parameter * shift for parameter, shift in
                         zip(feynman_parameters, spatial_shifts)), sp.zeros(number_of_momenta, 1))
    combined_mass = sum(
        parameter * (pole_masses_squared[pole] + (shift.T * external_gram * shift)[0])
        for parameter, shift, (_, pole) in zip(feynman_parameters, spatial_shifts, cycle)
    ) - (average_shift.T * external_gram * average_shift)[0]
    moments = reduce_routed_numerator(
        numerator, loop_squared, loop_products, tuple(average_shift), external_gram,
        sp.expand(combined_mass), sum(denominator_exponents), dimension)
    return numerator, sp.factor(combined_mass), moments


def simplex_moment_coefficients(moments, independent_parameters, denominator_exponents):
    """Integrate polynomial Feynman weights with the normalized simplex measure.

    The final denominator parameter equals one minus the independent parameters.
    Positive integer exponents supply the standard Dirichlet weight, allowing
    repeated propagators from derivatives. Scalar integrals stay unevaluated.
    Coefficients for the same denominator power remain signed and can be summed.
    """
    exponents = tuple(sp.sympify(value) for value in denominator_exponents)
    if len(exponents) != len(independent_parameters) + 1:
        raise ValueError('Supply one exponent per denominator.')
    if any(value.is_Integer is not True or value <= 0 for value in exponents):
        raise ValueError('Denominator exponents must be positive integers.')
    total_power = sum(exponents)
    final_parameter = 1 - sum(independent_parameters)
    weight = sp.prod(parameter**(exponent - 1) for parameter, exponent in
                     zip(independent_parameters, exponents[:-1]))
    weight *= final_parameter**(exponents[-1] - 1)
    normalization = sp.Integer(factorial(int(total_power - 1)))
    normalization /= sp.prod(factorial(int(exponent - 1)) for exponent in exponents)
    results = {}
    for denominator_power, coefficient in moments.items():
        if not independent_parameters:
            results[denominator_power] = sp.factor(coefficient)
            continue
        polynomial = sp.Poly(sp.expand(weight * coefficient), *independent_parameters)
        integral = sp.S.Zero
        for powers, term in polynomial.terms():
            integral += term * sp.prod(factorial(power) for power in powers) / factorial(
                len(independent_parameters) + sum(powers))
        results[denominator_power] = sp.factor(normalization * integral)
    return results


def scalar_moment_sizes(coefficients, scale_squared):
    """Return common-scale sizes after equal-power signed coefficients are combined.

    This treats scalar moments as approximately constant on the simplex. It is
    an order-of-magnitude assumption away from exceptional kinematics.
    """
    return {power: sp.Abs(coefficient) * scale_squared**(2 - power) / (16 * sp.pi**2)
            for power, coefficient in coefficients.items()}


def projected_cycle_moments(external_directions, external_momenta, pole_directions,
                            pole_masses_squared, invariants, epsilon, external_gram,
                            dimension, alp_pole=None, number_of_propagators=None):
    """Yield weighted signed scalar-moment coefficients for projected cycles.

    Parameters are indexed by the propagator position, with the last parameter
    eliminated using the simplex constraint. Each result contains the canonical
    cycle weight, cycle, ordered pole-mass and shift signature, combined mass,
    parameter symbols, and signed moment coefficients. Group compatible scalar
    moments only after their loop routing and parameter conventions match.
    """
    cycles = assembly.projected_vertices(tuple(external_directions), tuple(pole_directions),
                                       alp_pole)
    for weight, cycle in cycles:
        if number_of_propagators is not None and len(cycle) != number_of_propagators:
            continue
        parameters = sp.symbols(f'x_1:{len(cycle)}')
        feynman_parameters = tuple(parameters) + (1 - sum(parameters),)
        numerator, combined_mass, moments = reduce_cycle(
            cycle, external_directions, external_momenta, pole_directions,
            pole_masses_squared, invariants, epsilon, external_gram,
            feynman_parameters, dimension)
        if not moments:
            continue
        shift = sp.zeros(external_gram.rows, 1)
        signature = []
        for block, pole in cycle:
            signature.append((pole, tuple(shift)))
            for label in block:
                shift += sp.Matrix(external_momenta[label][1:])
        signed = {power: sp.factor(weight * value) for power, value in moments.items()}
        yield weight, cycle, tuple(signature), combined_mass, parameters, signed


class CycleReductionTest(unittest.TestCase):
    def test_projected_six_leg_tadpole(self):
        mass, b0, decay_constant, epsilon, pseudoscalar = sp.symbols(
            'm_a2 B_0 f_pi epsilon P_3')
        zero = sp.zeros(3)
        source = sp.diag(pseudoscalar, 0, 0)
        meson = sp.diag(1, 0, 0)
        invariants = {'S': [zero] * 7, 'P': [zero, zero, zero, source, zero, zero, zero],
                      'theta': [0] * 7, 'f_pi': decay_constant, 'B_0': b0,
                      'anomaly_mass_squared': sp.Symbol('m_0^2')}
        external_directions = {1: (zero, 1), 2: (meson, 0),
                               3: (meson, 0), 4: (meson, 0)}
        external_momenta = {1: (0, 1, 0, 0), 2: (0, 0, 1, 0),
                            3: (0, 0, 0, 1), 4: (0, -1, -1, -1)}
        cycles = list(projected_cycle_moments(
            external_directions, external_momenta, {'a': (zero, 1)}, {'a': mass},
            invariants, epsilon, sp.eye(3), sp.Symbol('d'), 'a', 1))
        self.assertEqual(len(cycles), 1)
        weight, _, signature, combined_mass, parameters, moments = cycles[0]
        self.assertEqual(weight, sp.Rational(1, 2))
        self.assertEqual(signature, (('a', (0, 0, 0)),))
        self.assertEqual(combined_mass, mass)
        self.assertEqual(parameters, ())
        self.assertNotEqual(moments[1], 0)
        self.assertTrue(moments[1].has(pseudoscalar))

    def test_normalized_simplex_and_raised_denominators(self):
        x, y, scale = sp.symbols('x y m^2', positive=True)
        self.assertEqual(simplex_moment_coefficients({2: x}, (x,), (1, 1)), {2: sp.Rational(1, 2)})
        self.assertEqual(simplex_moment_coefficients({3: x*y}, (x, y), (1, 1, 1)),
                         {3: sp.Rational(1, 12)})
        self.assertEqual(simplex_moment_coefficients({3: x}, (x,), (2, 1)), {3: sp.Rational(2, 3)})
        self.assertEqual(simplex_moment_coefficients({1: scale}, (), (1,)), {1: scale})
        self.assertEqual(scalar_moment_sizes({2: -sp.Symbol('c', positive=True)}, scale),
                         {2: sp.Symbol('c', positive=True) / (16 * sp.pi**2)})

    def test_projected_cubic_bubble(self):
        parameter, mass_a, mass_n, momentum_squared, dimension = sp.symbols('x m_a2 m_n2 s d')
        scalar_entry, decay_constant, b0, epsilon = sp.symbols('S_1 f_pi B_0 epsilon')
        zero = sp.zeros(3)
        scalar_first = sp.Matrix([[0, 0, 0], [0, 0, scalar_entry], [0, scalar_entry, 0]])
        meson = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
        external = sp.diag(0, 1, 0)
        invariants = {'S': [zero, scalar_first, zero], 'P': [zero, zero, zero],
                      'theta': [0, 0, 0], 'f_pi': decay_constant, 'B_0': b0,
                      'anomaly_mass_squared': sp.Symbol('m_0^2')}
        numerator, _, moments = reduce_cycle(
            (((1,), 'a'), ((2,), 'n')),
            {1: (external, 0), 2: (external, 0)}, {1: (0, 1), 2: (0, -1)},
            {'a': (zero, 1), 'n': (meson, 0)}, {'a': mass_a, 'n': mass_n},
            invariants, epsilon, sp.Matrix([[momentum_squared]]),
            (1 - parameter, parameter), dimension)
        self.assertNotEqual(numerator, 0)
        self.assertEqual(moments, {2: numerator})

    def test_bubble_with_routed_momentum_numerators(self):
        parameter, mass_a, mass_n, momentum_squared, dimension = sp.symbols('x m_a2 m_n2 s d')
        momenta = {1: (0, 1), 2: (0, -1)}
        cycle = (((1,), 'a'), ((2,), 'n'))

        def vertex(legs, momenta, invariants, epsilon, dot, **kwargs):
            return dot(momenta[-2], momenta[-2])

        with patch.object(assembly, 'projected_vertex', vertex):
            numerator, combined_mass, moments = reduce_cycle(
                cycle, {1: None, 2: None}, momenta, {'a': None, 'n': None},
                {'a': mass_a, 'n': mass_n}, None, None, sp.Matrix([[momentum_squared]]),
                (1 - parameter, parameter), dimension)
        self.assertEqual(sp.expand(combined_mass - ((1-parameter)*mass_a + parameter*mass_n
                                                   + parameter*(1-parameter)*momentum_squared)), 0)
        self.assertTrue(moments)
        self.assertEqual(sp.Poly(numerator, sp.Symbol('k^2')).degree(), 2)

    def test_tadpole_moment(self):
        mass, dimension = sp.symbols('m_a2 d')
        with patch.object(assembly, 'projected_vertex', lambda legs, momenta, invariants, epsilon, dot, **kwargs:
                          dot(momenta[-2], momenta[-2])):
            _, combined_mass, moments = reduce_cycle(
                (((1,), 'a'),), {1: None}, {1: (0, 0)}, {'a': None}, {'a': mass},
                None, None, sp.zeros(1), (sp.S.One,), dimension)
        self.assertEqual(combined_mass, mass)
        self.assertEqual(moments, {1: -mass})


if __name__ == '__main__':
    unittest.main()
