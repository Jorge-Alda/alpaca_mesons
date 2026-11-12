"""Tests for scripts/derive_alp_lo_mass_sum_rules.py (exact LO pole relations).

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_mass_sum_rules.py
"""

import random
import unittest

import numpy as np
import scipy.linalg
import sympy as sp

from derive_alp_lo_mass_sum_rules import (
    ALP_POSITION_IN_NEUTRAL_BASIS,
    RANDOM_SEEDS,
    RATIONAL_CONDENSATE,
    RATIONAL_SECOND_ORDER_SCALAR_TRACE,
    _random_rule_sources,
    _random_series,
    add_second_order_scalar_term,
    blocks_with_second_order_scalar,
    constant_rotation,
    derive,
    dressed_second_order_scalar_from_definition,
    dressed_second_order_scalar_trace,
    product_rule_invariant_form,
    solver_test_series,
    trace_sum_rule_invariant_form,
    kinetic_inverse_closed_form,
    neutral_kaon_mass_squared,
    neutral_kaon_null_vector,
    product_rule_closed_form,
    rational_sources,
    symbolic_second_order_scalar,
    symbolic_sources,
    trace_sum_rule_closed_form,
)
from derive_alp_lo_neutral_blocks import ALL_GENERATOR_INDICES, CHARGED_GENERATOR_INDICES, closed_form_blocks


class MassSumRulesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_all_checks_pass(self):
        for section, content in self.result.items():
            for label, passed in content["checks"].items():
                self.assertTrue(passed, f"{section}: {label}")

    def test_kinetic_inverse_for_rational_vector(self):
        kinetic_vector = sp.Matrix([sp.Rational(1, 3), -2, sp.Rational(5, 7), 0, sp.Rational(-1, 4)])
        kinetic = sp.eye(6)
        kinetic[:5, 5] = -kinetic_vector
        kinetic[5, :5] = -kinetic_vector.T
        kinetic[5, 5] = 1 + (kinetic_vector.T * kinetic_vector)[0]
        self.assertEqual(kinetic * kinetic_inverse_closed_form(kinetic_vector), sp.eye(6))
        self.assertEqual(kinetic.det(), 1)

    def test_sum_and_product_against_generalized_eigenvalues(self):
        # Independent of the project solver: scipy generalized symmetric eigenproblem.
        for ds_entry in (0, sp.Rational(1, 50) + sp.I * sp.Rational(1, 70)):
            sources = rational_sources(ds_entry)
            kinetic, mass_squared = closed_form_blocks(**sources)
            poles = scipy.linalg.eigh(np.array(mass_squared.evalf(30), dtype=float),
                                      np.array(kinetic.evalf(30), dtype=float), eigvals_only=True)
            self.assertAlmostEqual(poles.sum() / float(sp.expand(trace_sum_rule_closed_form(**sources))), 1.0, places=11)
            self.assertAlmostEqual(np.prod(poles) / float(sp.expand(product_rule_closed_form(**sources))), 1.0, places=9)

    def test_second_order_scalar_against_generalized_eigenvalues(self):
        # Same as above with the tr S_2 term in M^2_aa (scipy, independent of the project solver).
        second_order = {"condensate": RATIONAL_CONDENSATE,
                        "second_order_scalar_trace": RATIONAL_SECOND_ORDER_SCALAR_TRACE}
        for ds_entry in (0, sp.Rational(1, 50) + sp.I * sp.Rational(1, 70)):
            sources = rational_sources(ds_entry)
            kinetic, mass_squared = blocks_with_second_order_scalar(sources, **second_order)
            poles = scipy.linalg.eigh(np.array(mass_squared.evalf(30), dtype=float),
                                      np.array(kinetic.evalf(30), dtype=float), eigvals_only=True)
            self.assertAlmostEqual(poles.sum() / float(sp.expand(
                trace_sum_rule_closed_form(**sources, **second_order))), 1.0, places=11)
            self.assertAlmostEqual(np.prod(poles) / float(sp.expand(
                product_rule_closed_form(**sources, **second_order))), 1.0, places=9)

    def test_second_order_scalar_term_only_changes_alp_entry(self):
        sources = symbolic_sources()
        condensate, second_order_scalar_trace = symbolic_second_order_scalar()
        _, mass_squared = closed_form_blocks(**sources)
        shifted = add_second_order_scalar_term(mass_squared, sources["pion_decay_constant"],
                                               sources["alp_decay_constant"], condensate, second_order_scalar_trace)
        expected = sp.zeros(6)
        expected[5, 5] = -(sources["pion_decay_constant"] / sources["alp_decay_constant"])**2 * condensate * second_order_scalar_trace
        self.assertEqual(shifted - mass_squared, expected)

    def test_second_order_scalar_equals_bare_mass_shift(self):
        # In both closed forms tr S_2 acts as m_a0^2 -> m_a0^2 - r^2 B_0 tr S_2.
        sources = symbolic_sources(diagonal_quark_spurion=False)
        condensate, second_order_scalar_trace = symbolic_second_order_scalar()
        ratio = sources["pion_decay_constant"] / sources["alp_decay_constant"]
        shifted_sources = dict(sources)
        shifted_sources["bare_alp_mass_squared"] = (sources["bare_alp_mass_squared"]
                                                    - ratio**2 * condensate * second_order_scalar_trace)
        for closed_form in (trace_sum_rule_closed_form, product_rule_closed_form):
            difference = (closed_form(**sources, condensate=condensate,
                                      second_order_scalar_trace=second_order_scalar_trace)
                          - closed_form(**shifted_sources))
            self.assertEqual(sp.expand(difference), 0)

    def test_kaon_null_vector_depends_only_on_axial_minus_pseudoscalar(self):
        # Shifting the d-s entries of c_q^A and c_q^P by the same amount leaves (-D_7a, D_6a) unchanged.
        sources = symbolic_sources(diagonal_quark_spurion=True)
        shift = sp.Symbol("shift_ds_re", real=True) + sp.I * sp.Symbol("shift_ds_im", real=True)
        shift_matrix = sp.zeros(3)
        shift_matrix[1, 2], shift_matrix[2, 1] = shift, sp.conjugate(shift)
        shifted = dict(sources)
        shifted["axial_source"] = sources["axial_source"] + shift_matrix
        shifted["pseudoscalar_source"] = sources["pseudoscalar_source"] + shift_matrix
        _, original = neutral_kaon_null_vector(*closed_form_blocks(**sources), sources["quark_spurion"])
        _, transformed = neutral_kaon_null_vector(*closed_form_blocks(**shifted), shifted["quark_spurion"])
        self.assertEqual((original - transformed).applyfunc(lambda entry: sp.simplify(sp.expand(entry))), sp.zeros(6, 1))
        self.assertEqual(original[ALP_POSITION_IN_NEUTRAL_BASIS], 0)

    def test_invariant_forms_against_generalized_eigenvalues(self):
        # scipy generalized eigenproblem (independent of the project solver) vs the forms with tr hat S_2.
        second_order = {"condensate": RATIONAL_CONDENSATE,
                        "second_order_scalar_trace": RATIONAL_SECOND_ORDER_SCALAR_TRACE}
        for ds_entry in (0, sp.Rational(1, 50) + sp.I * sp.Rational(1, 70)):
            sources = rational_sources(ds_entry)
            kinetic, mass_squared = blocks_with_second_order_scalar(sources, **second_order)
            series = solver_test_series(sources, **second_order)
            new = {"condensate": RATIONAL_CONDENSATE,
                   "dressed_scalar_trace": dressed_second_order_scalar_from_definition(**series).trace()}
            poles = scipy.linalg.eigh(np.array(mass_squared.evalf(30), dtype=float),
                                      np.array(kinetic.evalf(30), dtype=float), eigvals_only=True)
            self.assertAlmostEqual(poles.sum() / float(sp.expand(
                trace_sum_rule_invariant_form(**sources, **new))), 1.0, places=11)
            self.assertAlmostEqual(np.prod(poles) / float(sp.expand(
                product_rule_invariant_form(**sources, **new))), 1.0, places=9)

    def test_invariant_forms_reject_raw_trace(self):
        # Inserting tr S_2 instead of tr hat S_2 into the invariant forms gives a different value (c^A, c^P != 0).
        second_order = {"condensate": RATIONAL_CONDENSATE,
                        "second_order_scalar_trace": RATIONAL_SECOND_ORDER_SCALAR_TRACE}
        sources = rational_sources()
        wrong = {"condensate": RATIONAL_CONDENSATE, "dressed_scalar_trace": RATIONAL_SECOND_ORDER_SCALAR_TRACE}
        self.assertNotEqual(sp.expand(trace_sum_rule_invariant_form(**sources, **wrong)
                                      - trace_sum_rule_closed_form(**sources, **second_order)), 0)
        self.assertNotEqual(sp.expand(product_rule_invariant_form(**sources, **wrong)
                                      - product_rule_closed_form(**sources, **second_order)), 0)

    def test_trace_identity_rejects_wrong_coefficient(self):
        series, _ = _random_series(random.Random(5))
        trace_from_definition = dressed_second_order_scalar_from_definition(**series).trace()
        axial = series["axial_couplings"][0]
        pseudoscalar, quark_mass = series["first_order_pseudoscalar"], series["quark_mass"]
        self.assertEqual(sp.expand(trace_from_definition - dressed_second_order_scalar_trace(
            series["second_order_scalar"].trace(), axial, pseudoscalar, quark_mass)), 0)
        wrong = (series["second_order_scalar"].trace() + 2 * sp.trace(axial * pseudoscalar)
                 - 4 * sp.trace(axial * axial * quark_mass))
        self.assertNotEqual(sp.expand(trace_from_definition - wrong), 0)

    def test_invariance_needs_gluon_coupling_shift(self):
        # Without c_G -> c_G - tr kappa the invariant trace form changes (tr kappa != 0).
        generator = random.Random(RANDOM_SEEDS[0])
        series, couplings = _random_series(generator)
        condensate, sources = _random_rule_sources(generator, couplings, series["quark_mass"])
        axial_shift = sp.diag(sp.Rational(1, 2), sp.Rational(1, 3), 0)
        rotated = constant_rotation(series, sp.zeros(3), axial_shift)
        shifted_sources = dict(sources, axial_source=sources["axial_source"] + axial_shift,
                               pseudoscalar_source=sources["pseudoscalar_source"] + axial_shift)
        original = trace_sum_rule_invariant_form(
            **sources, condensate=condensate,
            dressed_scalar_trace=dressed_second_order_scalar_from_definition(**series).trace())
        shifted_without_gluon = trace_sum_rule_invariant_form(
            **shifted_sources, condensate=condensate,
            dressed_scalar_trace=dressed_second_order_scalar_from_definition(**rotated).trace())
        self.assertNotEqual(sp.expand(shifted_without_gluon - original), 0)
        shifted_sources["gluon_coupling"] = sources["gluon_coupling"] - sp.trace(axial_shift)
        self.assertEqual(sp.expand(trace_sum_rule_invariant_form(
            **shifted_sources, condensate=condensate,
            dressed_scalar_trace=dressed_second_order_scalar_from_definition(**rotated).trace()) - original), 0)

    def test_kaon_pole_value(self):
        sources = rational_sources()
        self.assertEqual(neutral_kaon_mass_squared(sources["quark_spurion"]),
                         (sources["quark_spurion"][1, 1] + sources["quark_spurion"][2, 2]) / 2)

    def test_charged_block_with_down_strange_spurion_entry(self):
        # A d-s entry of chi_0 mixes pi^+ with K^+ inside the charged block only;
        # the charged poles are (chi_uu + eigenvalues of the d-s block of chi_0)/2.
        sources = rational_sources(sp.Rational(1, 50) + sp.I * sp.Rational(1, 70))
        _, mass_squared = closed_form_blocks(**sources, generator_indices=ALL_GENERATOR_INDICES)
        charged = list(CHARGED_GENERATOR_INDICES)
        charged_poles = sorted(np.linalg.eigvalsh(np.array(mass_squared.extract(charged, charged).evalf(30), dtype=float)))
        quark_spurion = np.array(sources["quark_spurion"].evalf(30), dtype=complex)
        down_strange = np.linalg.eigvalsh(quark_spurion[1:, 1:])
        expected = sorted([(quark_spurion[0, 0].real + value) / 2 for value in down_strange for _ in range(2)])
        np.testing.assert_allclose(charged_poles, expected, rtol=1e-12)


if __name__ == "__main__":
    unittest.main()
