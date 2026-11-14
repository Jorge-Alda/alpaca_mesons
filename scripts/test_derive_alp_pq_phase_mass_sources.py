"""Unit tests for the ALP mass sources from a Peccei-Quinn phase.

Run: .venv/bin/python -m unittest scripts.test_derive_alp_pq_phase_mass_sources
"""

import unittest

import sympy as sp

from scripts.derive_alp_pq_phase_mass_sources import (
    ALP_FIELD_RATIO,
    PQ_CHARGES,
    QUARK_MASSES,
    check_closed_forms,
    check_first_order_matching,
    check_numerical_examples,
    check_redefinition_invariance,
    check_sign_mapping,
    first_order_couplings,
    pq_charge_matrix,
    quark_mass_matrix,
    source_taylor_coefficients,
)

MAX_ORDER = 6


class PQPhaseMassSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scalar_coefficients, cls.pseudoscalar_coefficients = \
            source_taylor_coefficients(MAX_ORDER)

    def test_closed_forms(self):
        self.assertTrue(check_closed_forms(
            MAX_ORDER, scalar_coefficients=self.scalar_coefficients,
            pseudoscalar_coefficients=self.pseudoscalar_coefficients))

    def test_explicit_low_orders(self):
        mass_matrix = quark_mass_matrix()
        charge_matrix = pq_charge_matrix()
        self.assertEqual(self.scalar_coefficients[0], mass_matrix)
        self.assertEqual(self.pseudoscalar_coefficients[0], sp.zeros(3, 3))
        self.assertEqual(self.pseudoscalar_coefficients[1],
                         -mass_matrix * charge_matrix)
        self.assertEqual(self.scalar_coefficients[2],
                         -mass_matrix * charge_matrix ** 2)
        self.assertEqual(self.pseudoscalar_coefficients[3],
                         mass_matrix * charge_matrix ** 3)

    def test_sign_mapping(self):
        result = check_sign_mapping()
        mass_up, charge_up = QUARK_MASSES[0], PQ_CHARGES[0]
        self.assertEqual(sp.simplify(
            result["scalar"][0, 0]
            - mass_up * sp.cos(charge_up * ALP_FIELD_RATIO)), 0)
        self.assertEqual(sp.simplify(
            result["pseudoscalar"][0, 0]
            + mass_up * sp.sin(charge_up * ALP_FIELD_RATIO)), 0)

    def test_first_order_matching(self):
        pseudoscalar_coupling = check_first_order_matching(
            self.scalar_coefficients[1], self.pseudoscalar_coefficients[1])
        self.assertEqual(pseudoscalar_coupling, -pq_charge_matrix() / 2)
        self.assertEqual(first_order_couplings()[0], sp.zeros(3, 3))

    def test_wrong_sign_first_order_rejected(self):
        with self.assertRaises(AssertionError):
            check_first_order_matching(
                self.scalar_coefficients[1],
                -self.pseudoscalar_coefficients[1])

    def test_wrong_sign_closed_form_rejected(self):
        wrong_pseudoscalar = list(self.pseudoscalar_coefficients)
        wrong_pseudoscalar[1] = -wrong_pseudoscalar[1]
        with self.assertRaises(AssertionError):
            check_closed_forms(
                MAX_ORDER, scalar_coefficients=self.scalar_coefficients,
                pseudoscalar_coefficients=wrong_pseudoscalar)

    def test_redefinition_invariance_diagonal(self):
        connection = check_redefinition_invariance(MAX_ORDER)
        self.assertEqual(len(connection), MAX_ORDER + 1)

    def test_redefinition_invariance_non_diagonal(self):
        generator = sp.Matrix([[1, 2 - sp.I, 0],
                               [2 + sp.I, sp.Rational(-1, 3), sp.I],
                               [0, -sp.I, 2]])
        check_redefinition_invariance(4, generator)

    def test_numerical_examples(self):
        self.assertTrue(check_numerical_examples(MAX_ORDER))


if __name__ == "__main__":
    unittest.main()
