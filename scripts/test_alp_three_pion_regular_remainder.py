"""Unit test for the conditional retained-numerator pole subtraction."""
import unittest
import sympy as sp
import derive_radiative_tff_amplitude as tensors
from derive_alp_three_pion_regular_remainder import derive


class RegularRemainderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = derive()

    def test_regular_subtraction_reconstructs_each_numerator(self):
        mass = sp.Symbol('M_i', positive=True)
        for name in ('direct', 'pair'):
            numerator = self.results[name]
            reconstructed = numerator.subs(tensors.s1, mass**2) + (mass**2 - tensors.s1) * self.results[name + '_regular']
            self.assertEqual(sp.expand(numerator - reconstructed), 0)


if __name__ == '__main__':
    unittest.main()
