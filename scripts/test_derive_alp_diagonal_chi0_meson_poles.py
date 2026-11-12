"""Tests for scripts/derive_alp_diagonal_chi0_meson_poles.py (explicit pi^0, eta, eta' pole equations).

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_diagonal_chi0_meson_poles.py
"""

import unittest

import mpmath
import sympy as sp

from derive_alp_diagonal_chi0_meson_poles import (
    derive,
    explicit_inputs,
    explicit_residues,
    unmixed_flavour_diagonal_masses,
)
from derive_alp_lo_neutral_blocks import illustrative_lo_inputs
from derive_alp_meson_mass_shifts_exact import SQUARED_MOMENTUM, mixing_numerator_polynomial
from derive_alp_meson_mass_shifts_second_order import rational_test_sources


class DiagonalQuarkSpurionPoleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_all_checks_pass(self):
        for label, passed in self.result["checks"].items():
            self.assertTrue(passed, label)

    def test_unmixed_masses_without_anomaly_term(self):
        # m_0^2 = 0: the roots of P(s) are x_u, x_d, x_s
        sources = dict(rational_test_sources(sp.Rational(1, 2), illustrative_lo_inputs("flavor_diagonal")),
                       anomaly_mass_squared=0)
        inputs = explicit_inputs(sources)
        roots = unmixed_flavour_diagonal_masses(inputs)
        for root, value in zip(roots, sorted(inputs["quark_spurion_diagonal"])):
            self.assertLess(abs(root - mpmath.mpf(str(sp.N(value, 60)))), mpmath.mpf("1e-45"))

    def test_residues_sum_to_twice_leading_coefficient(self):
        # sum_n R_n = 2 x (coefficient of s^2 in N), from the partial fractions at large s
        sources = rational_test_sources(sp.Rational(1, 2), illustrative_lo_inputs("neutral"))
        inputs = explicit_inputs(sources)
        residues = explicit_residues(inputs, unmixed_flavour_diagonal_masses(inputs))
        numerator = mixing_numerator_polynomial(
            inputs["quark_spurion_diagonal"], [inputs["difference_source"][index, index] for index in range(3)],
            inputs["anomaly_mass_squared"], inputs["gluon_axial_combination"])
        leading = sp.Poly(numerator, SQUARED_MOMENTUM).coeff_monomial(SQUARED_MOMENTUM**2)
        self.assertLess(abs(mpmath.fsum(residues) - 2 * mpmath.mpf(str(sp.N(leading, 60)))), mpmath.mpf("1e-45"))


if __name__ == "__main__":
    unittest.main()
