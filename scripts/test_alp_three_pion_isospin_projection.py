"""Unit tests for the isospin Born projection and its flavour traces.

Run: .venv/bin/python -m unittest discover -s scripts -p test_alp_three_pion_isospin_projection.py
"""
import unittest
import sympy as sp
from derive_alp_three_pion_isospin_projection import derive


class ThreePionProjectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = derive()

    def test_trace_vertex_matches_components(self):
        self.assertEqual(self.results['trace_vertex_difference'], 0)

    def test_trace_projection_matches_components(self):
        self.assertEqual(self.results['trace_projection_difference'], 0)

    def test_mass_source_cancels(self):
        self.assertEqual(self.results['mass_source_coefficients'], (0, 0))

    def test_isoscalar_axial_source_vanishes(self):
        self.assertEqual(sp.factor(self.results['projected'].subs(sp.Symbol('c_Ad'), sp.Symbol('c_Au'))), 0)

    def test_exchange_antisymmetry(self):
        substitutions = {}
        for positive, negative in [('p_dot_pplus', 'p_dot_pminus'), ('q_dot_pplus', 'q_dot_pminus'), ('pplus_mu', 'pminus_mu')]:
            substitutions[sp.Symbol(positive)] = sp.Symbol(negative)
            substitutions[sp.Symbol(negative)] = sp.Symbol(positive)
        expression = self.results['projected']
        self.assertEqual(sp.factor(expression + expression.xreplace(substitutions)), 0)

    def test_ward_identity_from_vector_components(self):
        expression = self.results['projected'].subs({sp.Symbol('pplus_mu'): sp.Symbol('q_dot_pplus'), sp.Symbol('pminus_mu'): sp.Symbol('q_dot_pminus'), sp.Symbol('pzero_mu'): sp.Symbol('q_dot_p0'), sp.Symbol('q_mu'): sp.Symbol('q2')})
        self.assertEqual(sp.factor(expression), 0)



if __name__ == '__main__':
    unittest.main()
