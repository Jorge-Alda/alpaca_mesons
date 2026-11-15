#!/usr/bin/env python3
"""Symbolic conditional matching variation, natural units. Run with --test."""
import sys
import unittest
import sympy as sp
import derive_alp_spectral_matching as matching
import derive_alp_low_energy_benchmarks as benchmarks

virtuality = matching.virtuality
matching_scale_squared = matching.matching_scale_squared
other_matching_scale_squared = sp.symbols('s_m_prime', positive=True)
first_mass_squared, second_mass_squared = sp.symbols('M1_squared M2_squared', positive=True)
first_residue, second_residue = sp.symbols('g1 g2', real=True)
pole_residue, pole_mass_squared = matching.pole_residue, matching.pole_mass_squared
continuum = matching.spectral_contribution
scale_difference = sp.factor(continuum.subs(matching_scale_squared, other_matching_scale_squared) - continuum)
scale_derivative = sp.diff(continuum, matching_scale_squared)
single_pole = pole_residue * pole_mass_squared**2 / (pole_mass_squared + virtuality)**2
two_poles = first_residue * first_mass_squared**2 / (first_mass_squared + virtuality)**2 + second_residue * second_mass_squared**2 / (second_mass_squared + virtuality)**2
residues = sp.solve([first_residue + second_residue - pole_residue, first_residue*first_mass_squared + second_residue*second_mass_squared - pole_residue*pole_mass_squared], (first_residue, second_residue), dict=True)[0]
pole_difference = sp.factor(two_poles.subs(residues) - single_pole)

class VariationTests(unittest.TestCase):
    def test_scale_derivative(self):
        self.assertEqual(sp.simplify(scale_derivative + matching.spectral_integrand.subs(matching.spectral_mass, matching_scale_squared)), 0)
    def test_scale_difference(self):
        self.assertEqual(sp.simplify(scale_difference.subs(other_matching_scale_squared, matching_scale_squared)), 0)
        self.assertEqual(sp.limit(virtuality * scale_difference, virtuality, sp.oo), 0)
    def test_pole_constraints(self):
        self.assertEqual(sp.simplify((first_residue + second_residue - pole_residue).subs(residues)), 0)
        self.assertEqual(sp.simplify((first_residue*first_mass_squared + second_residue*second_mass_squared - pole_residue*pole_mass_squared).subs(residues)), 0)
        self.assertEqual(sp.limit(pole_difference, virtuality, 0), 0)
        self.assertEqual(sp.limit(virtuality * pole_difference, virtuality, sp.oo), 0)
    def test_both_low_energy_benchmarks(self):
        for low_form in (benchmarks.tree, benchmarks.dispersive):
            low_equal = low_form.subs({benchmarks.s1: -virtuality, benchmarks.s2: -virtuality})
            first_completion = low_equal + continuum + single_pole
            second_completion = low_equal + continuum + two_poles.subs(residues)
            self.assertEqual(sp.simplify(sp.expand(second_completion - first_completion) - pole_difference), 0)
            scale_changed = first_completion.xreplace({matching_scale_squared: other_matching_scale_squared})
            self.assertEqual(sp.simplify(sp.expand(scale_changed - first_completion) - scale_difference), 0)

if __name__ == '__main__':
    if '--test' in sys.argv:
        sys.argv.remove('--test')
        unittest.main()
    else:
        for label, expression in [('Continuum', continuum), ('Scale derivative', scale_derivative), ('Finite scale difference', scale_difference), ('Two-pole residues', sp.Tuple(residues[first_residue], residues[second_residue])), ('Pole-model difference', pole_difference)]:
            print(f'{label}:\n{sp.latex(expression)}\n')
