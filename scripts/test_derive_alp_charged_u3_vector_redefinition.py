"""Checks for the charge-conserving U(3) vertex comparison."""
import unittest
import sympy as sp
from derive_alp_charged_u3_vector_redefinition import derive


class ChargedU3VectorRedefinitionTest(unittest.TestCase):
    def test_contact_and_nonzero_physical_vertex(self):
        contact, amplitude, on_shell_contact = derive()
        self.assertNotEqual(sp.simplify(contact), 0)
        self.assertEqual(on_shell_contact, 0)
        self.assertNotEqual(sp.simplify(amplitude), 0)
        self.assertTrue(amplitude.has(sp.Symbol('w_a', real=True, nonzero=True)))


if __name__ == '__main__':
    unittest.main()
