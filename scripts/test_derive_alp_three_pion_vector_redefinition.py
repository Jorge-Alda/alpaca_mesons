"""Checks for the finite-mixing three-pion coordinate comparison."""
import unittest
import sympy as sp
from derive_alp_three_pion_vector_redefinition import derive


class ThreePionVectorRedefinitionTest(unittest.TestCase):
    def test_nonzero_contact_and_exchanges_cancel(self):
        original, contact, exchange_t, exchange_u, transformed = derive()
        rotation = sp.Symbol('v', real=True, nonzero=True)
        self.assertEqual(original, sp.Rational(188832,28561))
        self.assertEqual(contact, -sp.Rational(732,169)*rotation**2)
        self.assertEqual(exchange_t, sp.Rational(366,169)*rotation**2)
        self.assertEqual(exchange_u, sp.Rational(366,169)*rotation**2)
        self.assertEqual(transformed, original)
        self.assertNotEqual(original+contact,original)


if __name__ == '__main__':
    unittest.main()
