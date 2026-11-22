"""Symbolic one-loop momentum reduction; run with .venv/bin/python.

Euclidean momenta and symbolic squared masses are used throughout.
No momentum integrals are evaluated.
"""
import unittest
import sympy as sp
from functools import lru_cache


def derive_parameter_shift():
    x, y = sp.symbols('x y', real=True)
    loop_squared, loop_dot_first, loop_dot_second = sp.symbols('k2 kq1 kq2', real=True)
    first_squared, second_squared, first_dot_second = sp.symbols('q1sq q2sq q1q2', real=True)
    mass_zero, mass_one, mass_two = sp.symbols('m0sq m1sq m2sq', real=True)
    first_denominator = loop_squared + mass_zero
    second_denominator = loop_squared + 2 * loop_dot_first + first_squared + mass_one
    third_denominator = loop_squared + 2 * loop_dot_second + second_squared + mass_two
    combined = sp.expand((1-x-y)*first_denominator + x*second_denominator + y*third_denominator)
    shift_squared = x**2 * first_squared + y**2 * second_squared + 2*x*y*first_dot_second
    shifted_squared = loop_squared + 2*x*loop_dot_first + 2*y*loop_dot_second + shift_squared
    remainder = sp.expand(combined - shifted_squared)
    return combined, shifted_squared, remainder


def derive_angular_average():
    first_constant, second_constant = sp.symbols('a0 b0')
    first_vector, second_vector = sp.symbols('ap bp')
    first_quadratic, second_quadratic = sp.symbols('a2 b2')
    momentum_squared, dimension, vector_product = sp.symbols('k2 d ab', nonzero=True)
    numerator = sp.expand((first_constant + first_vector + first_quadratic*momentum_squared)
                          * (second_constant + second_vector + second_quadratic*momentum_squared))
    angular_average = sp.expand(numerator.subs({first_vector: 0, second_vector: 0}, simultaneous=True)
                                + momentum_squared*vector_product/dimension)
    return numerator, angular_average


def derive_moment_relations():
    momentum_squared, mass_squared = sp.symbols('k2 Delta')
    denominator_power = sp.Symbol('n', integer=True, positive=True)
    feynman_parameter, momentum_invariant = sp.symbols('x s', real=True)
    first_mass_squared, second_mass_squared = sp.symbols('m0sq m1sq')
    combined_mass = ((1-feynman_parameter)*first_mass_squared
                     + feynman_parameter*second_mass_squared
                     + feynman_parameter*(1-feynman_parameter)*momentum_invariant)
    recurrence = sp.simplify(
        momentum_squared*(momentum_squared+mass_squared)**(-denominator_power)
        - (momentum_squared+mass_squared)**(1-denominator_power)
        + mass_squared*(momentum_squared+mass_squared)**(-denominator_power)
    )
    derivative = sp.diff((momentum_squared+combined_mass)**(-denominator_power),
                         momentum_invariant)
    expected = (-denominator_power*feynman_parameter*(1-feynman_parameter)
                * (momentum_squared+combined_mass)**(-denominator_power-1))
    return recurrence, sp.simplify(derivative-expected)


def reduce_routed_numerator(numerator, loop_squared, loop_dot_products,
                            shift_coefficients, external_gram, combined_mass,
                            denominator_power, dimension):
    """Return scalar-moment coefficients after shift and dimensional reduction.

    The incoming denominator shift is k -> k + sum_i shift_coefficients[i] p_i.
    Terms with no remaining denominator are scaleless in dimensional regularization.
    """
    if len(loop_dot_products) != len(shift_coefficients):
        raise ValueError('Each external momentum needs one shift coefficient.')
    if external_gram.shape != (len(loop_dot_products), len(loop_dot_products)):
        raise ValueError('The Gram matrix must match the external momentum basis.')
    centered_squared = sp.Dummy('centered_squared')
    centered_products = sp.symbols(f'centered_dot_0:{len(loop_dot_products)}')
    denominator = sp.Dummy('denominator')
    shift = sp.Matrix(shift_coefficients)
    gram_shift = external_gram * shift
    substitutions = {
        loop_squared: centered_squared
        - 2 * sum(shift[i] * centered_products[i] for i in range(len(shift)))
        + (shift.T * external_gram * shift)[0]
    }
    substitutions.update({
        loop_dot_products[i]: centered_products[i] - gram_shift[i]
        for i in range(len(shift))
    })
    shifted = sp.Poly(sp.expand(numerator.subs(substitutions, simultaneous=True)),
                      *centered_products)

    @lru_cache(None)
    def contractions(indices):
        if not indices:
            return sp.Integer(1)
        first = indices[0]
        return sum(external_gram[first, indices[j]]
                   * contractions(indices[1:j] + indices[j+1:])
                   for j in range(1, len(indices)))

    averaged = 0
    for powers, coefficient in shifted.terms():
        indices = tuple(i for i, power in enumerate(powers) for _ in range(power))
        if len(indices) % 2:
            continue
        rank = len(indices) // 2
        if rank > 2:
            raise ValueError('The angular average supports loop rank at most four.')
        angular_denominator = sp.prod(dimension + 2*j for j in range(rank))
        averaged += coefficient * centered_squared**rank * contractions(indices) / angular_denominator

    reduced = sp.Poly(sp.expand(averaged.subs(centered_squared, denominator-combined_mass)),
                      denominator)
    return {
        denominator_power - power[0]: sp.factor(coefficient)
        for power, coefficient in reduced.terms()
        if denominator_power - power[0] > 0 and coefficient != 0
    }


def main():
    combined, shifted, remainder = derive_parameter_shift()
    numerator, averaged = derive_angular_average()
    for label, expression in [('Combined denominators', combined),
                              ('Shifted momentum squared', shifted),
                              ('Mass and external momentum terms', remainder),
                              ('Vertex product', numerator), ('Angular average', averaged)]:
        print(f'{label}:\n{sp.latex(expression)}\n')
    for label, expression in zip(('Scalar recurrence residual', 'Momentum derivative residual'),
                                 derive_moment_relations()):
        print(f'{label}:\n{sp.latex(expression)}\n')


class MomentumReductionTest(unittest.TestCase):
    def test_shift_for_three_denominators(self):
        combined, shifted, remainder = derive_parameter_shift()
        self.assertEqual(sp.expand(combined-shifted-remainder), 0)
        symbols = {str(item): item for item in remainder.free_symbols}
        self.assertFalse({'k2', 'kq1', 'kq2'} & set(symbols))

    def test_angular_average_and_triangle(self):
        _, averaged = derive_angular_average()
        a0, b0, a2, b2, k2, d, ab = sp.symbols('a0 b0 a2 b2 k2 d ab', nonzero=True)
        # Compare coefficients by name because nonzero assumptions differ.
        substitutions = {symbol: sp.Symbol(str(symbol)) for symbol in averaged.free_symbols}
        result = sp.expand(averaged.xreplace(substitutions))
        a0, b0, a2, b2, k2, d, ab = sp.symbols('a0 b0 a2 b2 k2 d ab')
        self.assertEqual(sp.expand(result - (a0*b0 + k2*(a0*b2+b0*a2+ab/d) + k2**2*a2*b2)), 0)
        self.assertEqual(sp.expand(result.subs({b0: 1, b2: 0, ab: 0}) - a0-a2*k2), 0)

    def test_moment_recurrence_and_momentum_derivative(self):
        self.assertEqual(derive_moment_relations(), (0, 0))

    def test_simplex_monomial_weights(self):
        x, y = sp.symbols('x y', nonnegative=True)
        self.assertEqual(sp.integrate(2*x*y, (y, 0, 1-x), (x, 0, 1)), sp.Rational(1, 12))

    def test_moment_mass_dimension(self):
        dimension, moment_power, denominator_power = sp.symbols('d r n', integer=True)
        exponent = dimension + 2*moment_power - 2*denominator_power
        self.assertEqual(sp.expand(exponent.subs(denominator_power, denominator_power+1)
                                   - exponent), -2)

    def test_routed_numerator_with_linear_terms(self):
        k2, kp, shift, p2, mass, d = sp.symbols('k2 kp x p2 Delta d')
        a, b, c, e = sp.symbols('a b c e')
        numerator = (a + b*kp) * (c + e*kp)
        result = reduce_routed_numerator(numerator, k2, (kp,), (shift,),
                                         sp.Matrix([[p2]]), mass, 2, d)
        self.assertEqual(sp.simplify(result[1] - b*e*p2/d), 0)
        self.assertEqual(sp.simplify(result[2] - (
            (a-b*shift*p2)*(c-e*shift*p2)-b*e*p2*mass/d)), 0)

    def test_rank_four_and_scaleless_terms(self):
        k2, kp, p2, mass, d = sp.symbols('k2 kp p2 Delta d')
        result = reduce_routed_numerator(kp**4, k2, (kp,), (0,),
                                         sp.Matrix([[p2]]), mass, 2, d)
        self.assertEqual(sp.simplify(result[1] + 6*p2**2*mass/(d*(d+2))), 0)
        self.assertEqual(sp.simplify(result[2] - 3*p2**2*mass**2/(d*(d+2))), 0)


if __name__ == '__main__':
    main()
