import unittest
import sympy as sp
from derive_aligned_mass_anticommutators import derive, nested_anticommutator, shifted_derivatives, scalar_pseudoscalar_derivatives


class TestAlignedMassAnticommutators(unittest.TestCase):
    def test_scalar_pseudoscalar_derivatives(self):
        anticommutator_map, scalar, pseudoscalar, derivatives = scalar_pseudoscalar_derivatives()
        for n, (scalar_result, pseudoscalar_result) in enumerate(derivatives):
            scalar_expected = sum(sp.binomial(n, k)*anticommutator_map**k*(sp.cos(sp.pi*k/2)*scalar[n-k]-sp.sin(sp.pi*k/2)*pseudoscalar[n-k]) for k in range(n+1))
            pseudoscalar_expected = sum(sp.binomial(n, k)*anticommutator_map**k*(sp.sin(sp.pi*k/2)*scalar[n-k]+sp.cos(sp.pi*k/2)*pseudoscalar[n-k]) for k in range(n+1))
            self.assertEqual(sp.expand(scalar_result-scalar_expected), 0)
            self.assertEqual(sp.expand(pseudoscalar_result-pseudoscalar_expected), 0)

    def test_shifted_derivatives(self):
        self.assertEqual(shifted_derivatives(), derive()[2])

    def test_formal_exponential_coefficients(self):
        x, anticommutator_map = sp.symbols('x A')
        mass = sp.symbols('M0:5')
        source = sum(x**n * mass[n] / sp.factorial(n) for n in range(5))
        series = sp.series(sp.exp(sp.I*x*anticommutator_map)*source, x, 0, 5).removeO()
        for n in range(5):
            coefficient = sp.diff(series, x, n).subs(x, 0)
            expected = sum(sp.binomial(n,k)*sp.I**k*anticommutator_map**k*mass[n-k] for k in range(n+1))
            self.assertEqual(sp.expand(coefficient - expected), 0)

    def test_derivatives_through_fourth_order(self):
        axial, mass, derivatives = derive()
        for n, derivative in enumerate(derivatives):
            expected = sum(sp.binomial(n, k)*sp.I**k*nested_anticommutator(axial, mass[n-k], k) for k in range(n+1))
            self.assertEqual(sp.expand(derivative - expected), 0)


if __name__ == '__main__':
    unittest.main()
