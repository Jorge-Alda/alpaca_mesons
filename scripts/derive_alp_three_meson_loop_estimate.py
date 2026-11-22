"""Symbolic external and internal ALP pole contraction of a six-leg vertex.

Run: .venv/bin/python scripts/derive_alp_three_meson_loop_estimate.py
The derivatives refer to the general invariant vertices in the associated record.
"""
import unittest
import sympy as sp


def derive_contact():
    epsilon, normalization = sp.symbols(r'\epsilon N', positive=True)
    external, internal_one, internal_two = sp.symbols('t_1 t_2 t_3')
    alp_source, meson_source = sp.symbols('a p')
    vertex_zero, vertex_one, vertex_two, vertex_three = sp.symbols(
        r'V_{0;6} V_{1;5} V_{2;4} V_{3;3}'
    )
    total = external + internal_one + internal_two
    potential = sum(
        vertex * alp_source**number / sp.factorial(number)
        * meson_source**(3-number) / sp.factorial(3-number)
        for number, vertex in enumerate((vertex_zero, vertex_one, vertex_two, vertex_three))
    )
    projected = potential.subs({
        alp_source: total / sp.sqrt(normalization),
        meson_source: 2 * epsilon * total / sp.sqrt(normalization),
    })
    result = sp.diff(projected, external, internal_one, internal_two).subs({
        external: 0, internal_one: 0, internal_two: 0,
    }).expand()
    return result, (vertex_zero, vertex_one, vertex_two, vertex_three), epsilon, normalization


def main():
    result, _, _, _ = derive_contact()
    print('Six-leg vertex contracted with three ALP pole vectors:')
    print(sp.latex(result))
    print('The tadpole symmetry factor is one half. The vertex symbols include the three meson daughters.')


class SixLegContractionTest(unittest.TestCase):
    def test_three_pole_contractions(self):
        result, vertices, epsilon, normalization = derive_contact()
        expected = (8 * epsilon**3 * vertices[0] + 12 * epsilon**2 * vertices[1]
                    + 6 * epsilon * vertices[2] + vertices[3]) / normalization**sp.Rational(3, 2)
        self.assertEqual(sp.simplify(result - expected), 0)


def sixth_order_kinetic():
    """Derive the six-field kinetic trace using cyclicity of the flavour trace."""
    field = sp.Symbol(r'\Phi', commutative=False)
    derivative = sp.Symbol(r'\partial_\mu\Phi', commutative=False)
    decay_constant = sp.Symbol(r'f_\pi', positive=True)
    left = []
    right = []
    for degree in range(1, 6):
        derivative_power = sum(field**position * derivative * field**(degree - 1 - position)
                               for position in range(degree))
        left.append(sp.I**degree * derivative_power / sp.factorial(degree))
        right.append((-sp.I)**degree * derivative_power / sp.factorial(degree))
    expression = sp.expand(sum(left[degree - 1] * right[6 - degree - 1]
                               for degree in range(1, 6)) / (4 * decay_constant**4))
    double_commutator = field * (field * derivative - derivative * field) - (field * derivative - derivative * field) * field
    expected = sp.expand(double_commutator**2 / (1440 * decay_constant**4))

    def cyclic_trace(value):
        # SymPy does not canonicalize cyclic traces of noncommutative symbols.
        terms = {}
        for term in sp.Add.make_args(sp.expand(value)):
            coefficient, factors = term.as_coeff_mul(field, derivative)
            word = []
            for factor in factors:
                base, exponent = factor.as_base_exp()
                word.extend([str(base)] * int(exponent))
            rotations = [tuple(word[index:] + word[:index]) for index in range(len(word))]
            representative = min(rotations)
            terms[representative] = terms.get(representative, 0) + coefficient
        return {word: sp.simplify(value) for word, value in terms.items() if sp.simplify(value) != 0}

    return cyclic_trace(expression), cyclic_trace(expected)


class SixMesonKineticTest(unittest.TestCase):
    def test_double_commutator_coefficient(self):
        obtained, expected = sixth_order_kinetic()
        self.assertEqual(obtained, expected)

    def test_tadpole_momentum_reduction(self):
        momentum_squared, pole_mass_squared = sp.symbols('k^2 m_a^2')
        quotient, remainder = sp.div(momentum_squared, momentum_squared - pole_mass_squared,
                                     momentum_squared)
        self.assertEqual(quotient, 1)
        self.assertEqual(remainder, pole_mass_squared)


if __name__ == '__main__':
    main()
    obtained, expected = sixth_order_kinetic()
    print('Six-field kinetic trace cyclic-word residual:')
    print(sp.latex(len(set(obtained.items()) ^ set(expected.items()))))
    momentum_squared, pole_mass_squared = sp.symbols('k^2 m_a^2')
    quotient, remainder = sp.div(momentum_squared, momentum_squared - pole_mass_squared,
                                 momentum_squared)
    print('Polynomial quotient and remainder of the tadpole momentum numerator:')
    print(sp.latex(quotient))
    print(sp.latex(remainder))
