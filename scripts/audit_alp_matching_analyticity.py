#!/usr/bin/env python3
"""Audit one-variable analytic matching examples; natural units. Run with --test."""
import sys
import unittest
import sympy as sp

virtuality, pole_mass_squared, second_pole_mass_squared = sp.symbols('Q2 M2 N2', positive=True)
asymptotic_coefficient, normalization = sp.symbols('C F0', nonzero=True)
residue_one, residue_two = sp.symbols('r_1 r_2')
correction = residue_one / (virtuality + pole_mass_squared) + residue_two / (virtuality + second_pole_mass_squared)
equations = [correction.subs(virtuality, 0), sp.limit(virtuality * correction, virtuality, sp.oo) - asymptotic_coefficient]
residues = sp.solve(equations, (residue_one, residue_two), dict=True)[0]
matched_correction = sp.factor(correction.subs(residues))
low_form = normalization * pole_mass_squared / (virtuality + pole_mass_squared)
high_form = asymptotic_coefficient / virtuality
jump = sp.factor(high_form - low_form)
matching_scale = sp.symbols('Q0', positive=True)

class TestMatching(unittest.TestCase):
    def test_normalization(self):
        self.assertEqual(sp.simplify(matched_correction.subs(virtuality, 0)), 0)
    def test_asymptotic_coefficient(self):
        self.assertEqual(sp.limit(virtuality * matched_correction, virtuality, sp.oo), asymptotic_coefficient)
    def test_distinct_pole_strengths(self):
        self.assertEqual(sp.simplify(residues[residue_one] + residues[residue_two] - asymptotic_coefficient), 0)
        self.assertNotEqual(residues[residue_one], 0)
        self.assertNotEqual(residues[residue_two], 0)
    def test_sharp_join_not_generically_continuous(self):
        self.assertNotEqual(sp.simplify(jump.subs(virtuality, matching_scale)), 0)

if __name__ == '__main__':
    if '--test' in sys.argv:
        sys.argv.remove('--test')
        unittest.main()
    else:
        for label, expression in [('Sharp join difference', jump), ('Zero-normalization condition', equations[0]), ('Asymptotic condition', equations[1]), ('Pole residues', sp.Tuple(residues[residue_one], residues[residue_two])), ('Correction', matched_correction)]:
            print(f'{label}:\n{sp.latex(expression)}\n')
