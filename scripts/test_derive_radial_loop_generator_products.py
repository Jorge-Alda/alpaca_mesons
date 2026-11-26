"""Exact checks of the assembled diagram products: completeness removes every generator
placeholder, the coefficients are exact, and the printed output has no component notation.

Run:  PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_generator_products
"""
import contextlib
import io
import re
import unittest

import sympy as sp

import derive_radial_loop_generator_products as products


class AssembledProductChecks(unittest.TestCase):
    def test_diagram_list(self):
        classes = [label.split('.')[0] for label, *_ in products.DIAGRAMS]
        self.assertEqual([classes.count(str(k)) for k in range(1, 7)], [3, 5, 3, 14, 12, 7])
        for label, symmetry, blocks, lines in products.DIAGRAMS:
            self.assertEqual(len(blocks), len(lines), label)
            self.assertEqual(symmetry, sp.Rational(1, 2) if len(lines) == 1 or lines in (('A', 'A'), ('S', 'S')) else 1)

    def test_no_generator_and_exact_coefficients(self):
        for label, *_ in products.DIAGRAMS:
            with self.subTest(label=label):
                total = products.diagram_total(label)
                self.assertFalse(total.has(*products.GENERATORS))
                self.assertFalse(total.atoms(sp.Float))
                self.assertNotEqual(total, 0)

    def test_printed_output_has_no_component_notation(self):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            products.print_diagrams([label for label, *_ in products.DIAGRAMS])
        output = buffer.getvalue()
        for pattern in (r'\\lambda_\{?[ij0-8]\}?(?![a-z])', r'E_\{?[0-9ij]', r'G_[0-9]', r'_slot', r'lambda_[ij]\b'):
            self.assertIsNone(re.search(pattern, output), pattern)


if __name__ == '__main__':
    unittest.main()
