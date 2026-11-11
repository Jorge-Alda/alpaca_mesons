"""Checks for the vector-redefinition three-point example."""

import unittest
import sympy as sp
from derive_alp_two_meson_vector_redefinition import derive


class VectorRedefinitionTest(unittest.TestCase):
    def test_contact_and_physical_vertex(self):
        _, contact, contact_vertex, physical_on_shell = derive()
        self.assertNotEqual(contact, 0)
        self.assertNotEqual(contact_vertex, 0)
        self.assertNotEqual(physical_on_shell, 0)
        self.assertTrue(sp.simplify(contact_vertex).has(sp.Symbol("p_1_sq", real=True)))


if __name__ == "__main__":
    unittest.main()
