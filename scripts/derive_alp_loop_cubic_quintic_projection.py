"""Symbolic projected cubic and quintic vertices, with GeV mass units.

Run: PYTHONPATH=scripts .venv/bin/python scripts/derive_alp_loop_cubic_quintic_projection.py
"""
import unittest
import sympy as sp

from derive_alp_dressed_vertex_power_counting import contract_potential_vertex, symmetrized_product


def projected_vertices():
    epsilon, f_pi, B_0, m_0_squared, theta_2 = sp.symbols('epsilon f_pi B_0 m_0_squared theta_2', nonzero=True, real=True)
    q_u, q_d, q_s = sp.symbols('q_u q_d q_s', real=True)
    scalar_ds = sp.Symbol('scalar_ds', real=True)
    p_u, p_d, p_s = sp.symbols('p_u p_d p_s', real=True)
    overlap = sp.diag(q_u, q_d, q_s)
    scalar_first = sp.Matrix([[0, 0, 0], [0, 0, scalar_ds], [0, scalar_ds, 0]])
    pseudoscalar_second = sp.diag(p_u, p_d, p_s)
    generators = (sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]), sp.diag(0, 1, 0),
                  sp.eye(3), sp.eye(3))
    zero = sp.zeros(3)
    invariants = {'f_pi': f_pi, 'B_0': B_0, 'anomaly_mass_squared': m_0_squared,
                  'S': [zero, scalar_first], 'P': [zero, zero, pseudoscalar_second],
                  'theta': [0, 0, theta_2]}
    alp_leg = (2 * epsilon * overlap, sp.S.One)
    meson_legs = tuple((generator, 0) for generator in generators)
    cubic_one = contract_potential_vertex((alp_leg, *meson_legs[:2]), invariants, epsilon)
    cubic_two = contract_potential_vertex((alp_leg, alp_leg, meson_legs[2]), invariants, epsilon)
    quintic_one = contract_potential_vertex((alp_leg, *meson_legs), invariants, epsilon)
    quintic_two = contract_potential_vertex((alp_leg, alp_leg, *meson_legs[:3]), invariants, epsilon)
    expected_cubic_one = -epsilon * B_0 / (2 * f_pi) * sp.trace(scalar_first * (generators[0] * generators[1] + generators[1] * generators[0]))
    expected_cubic_two = epsilon**2 / f_pi * sp.trace(generators[2] * (B_0 * (pseudoscalar_second - 2 * (overlap * scalar_first + scalar_first * overlap)) - m_0_squared * theta_2 * sp.eye(3) / 3))
    expected_quintic_one = epsilon * B_0 / f_pi**3 * sp.trace(scalar_first * symmetrized_product(generators)) / sp.factorial(4)
    expected_quintic_two = epsilon**2 * B_0 / f_pi**3 * (4 * sp.trace(scalar_first * symmetrized_product((overlap, *generators[:3]))) / sp.factorial(4) - sp.trace(pseudoscalar_second * symmetrized_product(generators[:3])) / sp.factorial(3))
    return ((cubic_one, cubic_two, quintic_one, quintic_two),
            (expected_cubic_one, expected_cubic_two, expected_quintic_one, expected_quintic_two))


def main():
    calculated, expected = projected_vertices()
    for label, value in zip(('one-ALP cubic', 'two-ALP cubic', 'one-ALP quintic', 'two-ALP quintic'), calculated):
        print(label + ':\n' + sp.latex(sp.factor(value)) + '\n')
    print('First cubic-quintic bubble numerator:\n' + sp.latex(sp.factor(calculated[0] * calculated[3])))
    print('Second cubic-quintic bubble numerator:\n' + sp.latex(sp.factor(calculated[1] * calculated[2])))
    assert all(sp.simplify(result - expression) == 0 for result, expression in zip(calculated, expected))


class ProjectedCubicQuinticTest(unittest.TestCase):
    def test_invariant_contractions(self):
        calculated, expected = projected_vertices()
        for result, expression in zip(calculated, expected):
            self.assertEqual(sp.simplify(result - expression), 0)
        self.assertNotEqual(sp.simplify(calculated[0] * calculated[3]), 0)
        self.assertNotEqual(sp.simplify(calculated[1] * calculated[2]), 0)
        self.assertEqual(sp.expand(calculated[0] * calculated[3] - expected[0] * expected[3]), 0)
        self.assertEqual(sp.expand(calculated[1] * calculated[2] - expected[1] * expected[2]), 0)


if __name__ == '__main__':
    main()
