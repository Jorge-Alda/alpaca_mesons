"""Connect projected Euclidean loop cycles to an on-shell 1PI callback.

The scalar-moment callable is supplied by the caller and must preserve the
combined-mass and Feynman-parameter dependence when differentiation is needed.
All momenta are incoming coefficient vectors in one fixed Euclidean basis.
Run: PYTHONPATH=scripts .venv/bin/python -m unittest derive_alp_loop_1pi_callback
"""
import unittest
from itertools import permutations
from unittest.mock import patch
import sympy as sp

from derive_alp_three_meson_loop_assembly import projected_vertices
from derive_alp_three_meson_loop_cubic_cycles import cubic_cycle


def canonical_bubble_moment(combined_mass, parameter, coefficient, exponents):
    """Canonicalize a bubble under exchange of its two propagators.

    The coefficient includes all numerator and normalized simplex weights.
    The returned exponents preserve the corresponding parameter order.
    """
    reflected_mass = sp.expand(combined_mass.subs(parameter, 1 - parameter))
    reflected_coefficient = sp.expand(coefficient.subs(parameter, 1 - parameter))
    if sp.srepr(reflected_mass) < sp.srepr(sp.expand(combined_mass)):
        return reflected_mass, reflected_coefficient, tuple(reversed(exponents))
    return sp.expand(combined_mass), sp.expand(coefficient), exponents


def canonical_simplex_moment(combined_mass, parameters, coefficient, exponents):
    """Canonicalize bubble and triangle moments under propagator relabeling.

    Raised-denominator simplex weights are absorbed into the numerator, with
    their normalization retained. Resulting exponents are all one.
    """
    weights = parameters + (1 - sum(parameters),)
    normalization = sp.factorial(sum(exponents) - 1) / (
        sp.prod(sp.factorial(exponent - 1) for exponent in exponents)
        * sp.factorial(len(exponents) - 1))
    weighted = sp.expand(normalization * coefficient * sp.prod(
        weight**(exponent - 1) for weight, exponent in zip(weights, exponents)))
    candidates = []
    for permutation in permutations(weights):
        substitution = dict(zip(parameters, permutation[:-1]))
        mass = sp.expand(combined_mass.subs(substitution, simultaneous=True))
        numerator = sp.expand(weighted.subs(substitution, simultaneous=True))
        candidates.append((sp.srepr(mass), mass, numerator))
    _, mass, numerator = min(candidates, key=lambda candidate: candidate[0])
    return mass, numerator, (1,) * len(exponents)


def minkowski_one_loop(directions, momenta, basis_directions, pole_directions,
                       pole_masses_squared, invariants, epsilon, gram, dimension,
                       scalar_moment, order, alp_pole=None, gram_variation=None):
    """Return one epsilon coefficient of a Minkowski 1PI derivative.

    scalar_moment(power, combined_mass, parameters, exponents, coefficient)
    returns the corresponding formal Euclidean scalar integral with the
    coefficient inside its Feynman-parameter integral and with the normalized
    Feynman measure for exponents included.
    A non-None gram_variation is the matrix dG/ds for a two-point derivative,
    at fixed dimensionless momentum coefficients. The Euclidean result is
    continued and multiplied by minus one by this function.
    """
    if gram_variation is not None and len(directions) != 2:
        raise ValueError('A momentum derivative requires two external legs.')
    parameter = sp.Dummy('t', real=True)
    varied_gram = gram if gram_variation is None else gram + parameter * gram_variation
    labels = tuple(range(len(directions)))
    external_directions = {}
    external_momenta = {}
    for label, (direction, momentum) in enumerate(zip(directions, momenta)):
        meson = sum((direction[index] * basis_directions[index][0]
                     for index in range(len(basis_directions))), sp.zeros(3))
        alp = sum(direction[index] * basis_directions[index][1]
                  for index in range(len(basis_directions)))
        external_directions[label] = (meson, alp)
        external_momenta[label] = (sp.S.Zero,) + tuple(momentum)
    output = sp.S.Zero
    for weight, cycle in projected_vertices(labels, tuple(pole_directions), alp_pole):
        for exponents, combined_mass, parameters, coefficients in cubic_cycle(
                cycle, external_directions, external_momenta, pole_directions,
                pole_masses_squared, invariants, epsilon, varied_gram, dimension, order=order):
            for power, coefficient in coefficients.items():
                moment_mass = combined_mass
                moment_exponents = exponents
                if len(parameters) in (1, 2):
                    moment_mass, coefficient, moment_exponents = canonical_simplex_moment(
                        combined_mass, parameters, coefficient, exponents)
                output -= weight * scalar_moment(
                    power, moment_mass, parameters, moment_exponents, coefficient)
    if gram_variation is not None:
        output = sp.diff(output, parameter).subs(parameter, 0)
    return sp.expand(output)


class CallbackTest(unittest.TestCase):
    def test_triangle_parameter_permutation(self):
        x, y, mass_a, mass_b, mass_c = sp.symbols('x y m_a m_b m_c')
        mass = x*mass_a + y*mass_b + (1-x-y)*mass_c
        numerator = 1 + x - y
        first = canonical_simplex_moment(mass, (x, y), numerator, (1, 2, 1))
        substitution = {x: y, y: x}
        second = canonical_simplex_moment(
            mass.subs(substitution, simultaneous=True), (x, y),
            numerator.subs(substitution, simultaneous=True), (2, 1, 1))
        self.assertEqual(first, second)

    def test_bubble_parameter_reversal(self):
        parameter, mass_a, mass_b, momentum = sp.symbols('x m_a2 m_b2 s')
        combined_mass = parameter * mass_a + (1 - parameter) * mass_b + parameter * (1 - parameter) * momentum
        coefficient = 1 + parameter
        direct = canonical_bubble_moment(combined_mass, parameter, coefficient, (1, 2))
        reflected = canonical_bubble_moment(
            combined_mass.subs(parameter, 1 - parameter), parameter,
            coefficient.subs(parameter, 1 - parameter), (2, 1))
        self.assertEqual(direct, reflected)

    def test_symbolic_alp_tadpole(self):
        source, b0, decay_constant, epsilon, mass = sp.symbols('c B f e m_a2')
        scalar_sources = [sp.zeros(3)] * 7
        scalar_sources[3] = sp.diag(source, 0, 0)
        invariants = {'S': scalar_sources, 'P': [sp.zeros(3)] * 7,
                      'theta': [0] * 7, 'f_pi': decay_constant, 'B_0': b0,
                      'anomaly_mass_squared': sp.Symbol('m0')}
        moment = sp.Symbol('I')
        result = minkowski_one_loop(
            (sp.Matrix([1]),), (sp.Matrix([0]),), ((sp.zeros(3), 1),),
            {'a': (sp.zeros(3), 1)}, {'a': mass}, invariants, epsilon,
            sp.zeros(1), sp.Symbol('d'), lambda power, mass, params, exponents, coefficient:
            coefficient * moment, 3)
        self.assertEqual(result, b0 * source * moment / (2 * decay_constant))

    def test_weights_sign_and_moment_derivative(self):
        mass, gram_entry, derivative = sp.symbols('m_a^2 s j')
        formal_moment = sp.Function('I')
        basis = ((sp.eye(3), 0),)
        poles = {'a': (sp.eye(3), 0)}
        direction = sp.Matrix([1])
        momentum = sp.Matrix([1])

        def fake_cycle(cycle, *args, **kwargs):
            current_gram = args[-2]
            return [((1,) * len(cycle), mass + current_gram[0, 0], (), {1: derivative})]

        with patch('derive_alp_loop_1pi_callback.cubic_cycle', fake_cycle):
            value = minkowski_one_loop((direction,), (sp.Matrix([0]),), basis, poles,
                                       {'a': mass}, {}, sp.Symbol('epsilon'),
                                       sp.Matrix([[gram_entry]]), 4,
                                       lambda power, combined_mass, params, exponents, coefficient:
                                       coefficient * formal_moment(combined_mass), 0)
            self.assertEqual(value, -derivative * formal_moment(mass + gram_entry) / 2)
            value = minkowski_one_loop((direction, direction), (momentum, -momentum),
                                       basis, poles, {'a': mass}, {}, sp.Symbol('epsilon'),
                                       sp.Matrix([[gram_entry]]), 4,
                                       lambda power, combined_mass, params, exponents, coefficient:
                                       coefficient * formal_moment(combined_mass), 0,
                                       gram_variation=sp.ones(1))
            self.assertEqual(value, 0)


if __name__ == '__main__':
    unittest.main()
