"""Tests for scripts/derive_alp_meson_mass_shifts_exact.py (exact LO pole equations).

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_meson_mass_shifts_exact.py
"""

import unittest

import sympy as sp

from derive_alp_lo_neutral_blocks import illustrative_lo_inputs
from derive_alp_meson_mass_shifts_exact import (
    SQUARED_MOMENTUM,
    canonical_mass_matrix_from_sources,
    derive,
    explicit_secular_polynomial,
    flavour_diagonal_meson_polynomial,
    rational_checks,
    inverse_meson_operator,
)
from derive_alp_meson_mass_shifts_second_order import rational_test_sources


class ExactPoleEquationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_all_checks_pass(self):
        for section, content in self.result.items():
            for label, passed in content["checks"].items():
                self.assertTrue(passed, f"{section}: {label}")

    def test_flavour_diagonal_sources(self):
        content = rational_checks(epsilon=sp.Rational(1, 3), base_sources=illustrative_lo_inputs("flavor_diagonal"))
        for label, passed in content["checks"].items():
            self.assertTrue(passed, label)
        kaon = sp.Rational(1, 2) * (sp.Rational(1, 40) + sp.Rational(47, 100))
        # c^AP_ds = 0: the K^0 doublet decouples, double root at m_K0,0^2
        self.assertTrue(content["polynomial"].rem(sp.Poly((SQUARED_MOMENTUM - kaon)**2, SQUARED_MOMENTUM)).is_zero)

    def test_zero_mixing_limit(self):
        s = SQUARED_MOMENTUM
        diagonal = (sp.Rational(3, 200), sp.Rational(1, 40), sp.Rational(47, 100))
        polynomial = explicit_secular_polynomial(diagonal, sp.diag(1, 2, 3), sp.Rational(4, 5), sp.Rational(1, 2),
                                                 0, sp.Rational(3, 5))
        kaon = (diagonal[1] + diagonal[2]) / 2
        expected = (s - kaon)**2 * flavour_diagonal_meson_polynomial(diagonal, sp.Rational(4, 5)) * (s - sp.Rational(3, 5))
        self.assertEqual(sp.expand(polynomial - expected), 0)

    def test_canonical_matrix_symmetric_and_trace(self):
        sources = rational_test_sources(sp.Rational(1, 2))
        condensate, trace = sp.Rational(27, 10), sp.Rational(1, 20)
        canonical = canonical_mass_matrix_from_sources(sources, condensate, trace)
        self.assertEqual((canonical - canonical.T).applyfunc(sp.simplify), sp.zeros(6))
        roots_sum = -sp.Poly(explicit_secular_polynomial(
            [sources["quark_spurion"][index, index] for index in range(3)],
            sources["axial_source"] - sources["pseudoscalar_source"], sources["anomaly_mass_squared"],
            sources["gluon_coupling"] + sp.trace(sources["axial_source"]), sp.Rational(1, 2), canonical[5, 5]),
            SQUARED_MOMENTUM).all_coeffs()[1]
        self.assertEqual(sp.simplify(roots_sum - canonical.trace()), 0)

    def test_inverse_meson_operator_without_anomaly_term(self):
        # m_0^2 = 0, chi_0 = diag(x_q): M_s^{-1} acts entrywise as 1/((x_q + x_r)/2 - s)
        s = SQUARED_MOMENTUM
        quark_spurion = sp.diag(1, 2, 4)
        argument = sp.Matrix([[5, 0, 0], [0, 7, 2 - sp.I], [0, 2 + sp.I, 11]])
        expected = sp.Matrix(3, 3, lambda row, column: argument[row, column]
                             / ((quark_spurion[row, row] + quark_spurion[column, column]) / 2 - s))
        self.assertEqual((inverse_meson_operator(argument, quark_spurion, 0) - expected).applyfunc(sp.simplify), sp.zeros(3))


if __name__ == "__main__":
    unittest.main()
