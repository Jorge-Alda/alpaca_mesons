"""Checks for the vector-redefinition three-point example."""

import unittest
import sympy as sp
from derive_alp_two_meson_vector_redefinition import derive, derive_mixed_pole_relation


class VectorRedefinitionTest(unittest.TestCase):
    def test_contact_and_physical_vertex(self):
        _, contact, contact_vertex, physical_on_shell = derive()
        self.assertNotEqual(contact, 0)
        self.assertNotEqual(contact_vertex, 0)
        self.assertNotEqual(physical_on_shell, 0)
        self.assertTrue(sp.simplify(contact_vertex).has(sp.Symbol("p_1_sq", real=True)))

    def test_mixed_pole_contact(self):
        off_diagonal, pole_piece, contact, transformed_vertex = derive_mixed_pole_relation()
        self.assertNotEqual(off_diagonal, 0)
        self.assertNotEqual(contact, 0)
        self.assertEqual(sp.simplify(transformed_vertex - pole_piece - contact), 0)
        self.assertNotEqual(sp.simplify(transformed_vertex - pole_piece), 0)


if __name__ == "__main__":
    unittest.main()
