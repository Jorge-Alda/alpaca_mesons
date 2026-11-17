"""Compare the derived singlet vertices with the independent joint expansion."""
import unittest

import sympy as sp

from derive_alp_lo_joint_vertices import vertex
from derive_alp_three_pseudoscalar_topological_vertices import topological_vertex


class TopologicalVertices(unittest.TestCase):
    def test_against_joint_ccwz_expansion(self):
        matrices = [sp.diag(1, 2), sp.diag(-1, 3), sp.diag(2, -2), sp.diag(1, -3)]
        components = [sp.Rational(2, 3), sp.Rational(-1, 2), 1, sp.Rational(3, 5)]
        zero = sp.zeros(2)
        theta = [0, sp.Rational(2, 5), sp.Rational(3, 7), sp.Rational(-1, 3), 0]
        for count in (3, 4):
            with self.subTest(count=count):
                common = (matrices[:count], components[:count], [0] * count,
                          [zero] * 5, [zero] * 5, [zero] * 5)
                result = vertex(*common, theta, 2, 3, 5, 0)
                without_topological_source = vertex(*common, [0] * 5, 2, 3, 5, 0)
                expression = topological_vertex(count)
                values = {'f_pi': 2, 'f_a': 3, 'm_0_squared': 5, 'N_f': 2,
                          **{f'Theta_{n}': theta[n] for n in range(1, 4)},
                          **{f'xi_{n}': components[n] for n in range(count)},
                          **{f'trace_lambda_{n}': sp.trace(matrices[n]) for n in range(count)}}
                expected = expression.subs({symbol: values[symbol.name] for symbol in expression.free_symbols})
                self.assertEqual(sp.simplify(result - without_topological_source - expected), 0)


if __name__ == '__main__':
    unittest.main()
