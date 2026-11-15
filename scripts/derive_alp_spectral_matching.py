#!/usr/bin/env python3
"""Conditional spectral matching checks in natural units. Invoke with --test."""
import sys
import unittest
import sympy as sp

spectral_mass, matching_scale_squared, virtuality = sp.symbols('s s_m Q2', positive=True)
coefficient, low_normalization, normalization, low_single_coefficient, single_coefficient = sp.symbols('C L_0 F_0 L_infty C_infty', real=True)
pole_mass_squared, pole_residue = sp.symbols('M_eff2 R_eff', nonzero=True)
spectral_integrand = 3 * coefficient * virtuality**2 / (spectral_mass + virtuality)**4
spectral_contribution = sp.integrate(spectral_integrand, (spectral_mass, matching_scale_squared, sp.oo))
effective_pole = pole_residue * pole_mass_squared**2 / (pole_mass_squared + virtuality)**2
single_effective_pole = pole_residue * pole_mass_squared / (pole_mass_squared + virtuality)
matching_equations = [pole_residue + low_normalization - normalization,
                      pole_residue * pole_mass_squared + low_single_coefficient - single_coefficient]
matching_solution = sp.solve(matching_equations, (pole_residue, pole_mass_squared), dict=True)[0]

class MatchingTests(unittest.TestCase):
    def test_zero_normalization(self):
        self.assertEqual(sp.limit(spectral_contribution, virtuality, 0), 0)
        self.assertEqual(sp.simplify((pole_residue + low_normalization - normalization).subs(matching_solution)), 0)
    def test_symmetric_limit(self):
        self.assertEqual(sp.limit(virtuality * spectral_contribution, virtuality, sp.oo), coefficient)
        self.assertEqual(sp.limit(virtuality * effective_pole, virtuality, sp.oo), 0)
    def test_single_limit(self):
        self.assertEqual(sp.simplify((sp.limit(virtuality * single_effective_pole, virtuality, sp.oo) + low_single_coefficient - single_coefficient).subs(matching_solution)), 0)
    def test_cutoff_density(self):
        self.assertEqual(sp.simplify(sp.diff(spectral_contribution, matching_scale_squared) + spectral_integrand.subs(spectral_mass, matching_scale_squared)), 0)
    def test_unequal_hard_ratio(self):
        fraction = sp.symbols('u', real=True)
        spectral = sp.integrate(sp.apart(3 * 2 * 5 / ((spectral_mass + 2)**2 * (spectral_mass + 5)**2), spectral_mass), (spectral_mass, 0, sp.oo))
        light_cone = sp.integrate(6 * fraction * (1 - fraction) / (2 * fraction + 5 * (1 - fraction)), (fraction, 0, 1))
        self.assertEqual(sp.simplify(spectral - light_cone), 0)
    def test_disjoint_spectral_benchmark(self):
        low_density = sp.Piecewise((1, (spectral_mass >= 0) & (spectral_mass <= 1)), (0, True))
        high_density = sp.Piecewise((1, (spectral_mass >= 2) & (spectral_mass <= 3)), (0, True))
        self.assertEqual(sp.integrate(low_density * high_density, (spectral_mass, 0, 3)), 0)
        self.assertEqual(sp.integrate(low_density + high_density, (spectral_mass, 0, 3)), 2)
    def test_spacelike_continuity(self):
        self.assertEqual(sp.simplify(sp.limit(spectral_contribution, virtuality, 1, dir='+') - sp.limit(spectral_contribution, virtuality, 1, dir='-')), 0)

if __name__ == '__main__':
    if '--test' in sys.argv:
        sys.argv.remove('--test')
        unittest.main()
    else:
        for label, expression in [('Spectral integrand', spectral_integrand), ('Equal virtualities', spectral_contribution), ('Normalization and singly virtual constraints', sp.Tuple(*matching_equations)), ('Effective pole parameters', sp.Tuple(matching_solution[pole_residue], matching_solution[pole_mass_squared]))]:
            print(f'{label}:\n{sp.latex(expression)}\n')
