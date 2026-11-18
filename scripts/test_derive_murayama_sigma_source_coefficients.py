import unittest
import sympy as sp
from scripts.derive_murayama_sigma_source_coefficients import derive


class SourceCoefficientTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_source_matching(self):
        result = self.result
        self.assertEqual(result['scalar2'], result['pseudoscalar1'])
        self.assertEqual(result['scalar0_radial'], -result['pseudoscalar1']/result['fa'])

    def test_polynomial_and_coupling(self):
        result = self.result
        self.assertEqual(sp.simplify(result['meson_polynomial']-result['source_polynomial']), 0)
        self.assertEqual(sp.simplify(result['mass_bracket']-result['original_bracket']), 0)


if __name__ == '__main__':
    unittest.main()
