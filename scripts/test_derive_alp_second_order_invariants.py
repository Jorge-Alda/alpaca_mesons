"""Unit tests for the second-order dressed sources (derivative convention F_n = d^nF/dx^n at 0).

Run: .venv/bin/python -m unittest scripts.test_derive_alp_second_order_invariants
"""

import unittest

import sympy as sp

from scripts.derive_alp_second_order_invariants import (
    LEFT_CURRENT,
    MASS,
    PSEUDOSCALAR,
    QUARK_MASS,
    RIGHT_CURRENT,
    SCALAR,
    SCALAR_COUPLING,
    VECTOR_COUPLING,
    alda_compact_forms,
    anticommutator,
    commutator,
    check_chiral_trace_forms,
    check_constant_connections,
    check_existing_transport_agreement,
    check_first_order_regression,
    check_groupings,
    check_groupings_on_matrices,
    check_hermiticity,
    check_invariance,
    check_numerical_transport,
    check_real_mass_on_matrices,
    check_vacuum_angle,
    derive_real_mass_specialization,
    real_mass_vector_axial_forms,
    closed_form_dressed_mass_two,
    derive_alda_specialization,
    derive_theta_two,
    derive_transport,
    grouping_chiral_nested,
    grouping_vector_axial_nested,
    is_zero,
    transport_scalar_pseudoscalar_two,
)


def perturbed_chiral_nested(s, p, right, left, bracket, anti_bracket):
    """Chosen grouping with the coefficient of {C_R1 - C_L1, P_0} changed from 1/2 to 1/4."""
    scalar, pseudoscalar = grouping_chiral_nested(s, p, right, left, bracket, anti_bracket)
    return scalar - anti_bracket(right[1] - left[1], p[0]) / 4, pseudoscalar


class TransportAndClosedFormTests(unittest.TestCase):
    def test_transport_derivatives(self):
        for chirality, values in derive_transport().items():
            for residual in (values["ode_residuals"] + values["unitarity_residuals"]
                             + values["closed_form_residuals"]):
                self.assertEqual(residual, 0, chirality)

    def test_agreement_with_existing_transport_recursion(self):
        for key, residual in check_existing_transport_agreement().items():
            self.assertEqual(residual, 0, key)

    def test_wrong_closed_form_is_detected(self):
        """Dropping the 2 C_{R,0} M_0 C_{L,0} term must give a nonzero residual."""
        dressed, _ = transport_scalar_pseudoscalar_two()
        wrong = closed_form_dressed_mass_two(MASS, RIGHT_CURRENT, LEFT_CURRENT) \
            - 2 * RIGHT_CURRENT[0] * MASS[0] * LEFT_CURRENT[0]
        self.assertNotEqual(sp.expand(dressed[2] - wrong), 0)

    def test_theta_coefficients(self):
        values = derive_theta_two()
        self.assertEqual(values["theta_residual"], 0)
        self.assertEqual(values["theta_photon_residual"], 0)


class GroupingTests(unittest.TestCase):
    def test_all_groupings_symbolic(self):
        values = check_groupings()
        for name, row in values["table"].items():
            self.assertEqual(row["scalar_residual"], 0, name)
            self.assertEqual(row["pseudoscalar_residual"], 0, name)

    def test_shortest_groupings_and_counts(self):
        values = check_groupings()
        self.assertEqual(values["shortest_chiral"], "chiral, nested")
        self.assertEqual(values["shortest_vector_axial"], "vector/axial, nested")
        for name in ("chiral, nested", "vector/axial, nested"):
            self.assertEqual(values["table"][name]["scalar_count"][0], 5)
            self.assertEqual(values["table"][name]["pseudoscalar_count"][0], 5)

    def test_chosen_forms_use_only_scalar_pseudoscalar_and_chiral_currents(self):
        row = check_groupings()["table"]["chiral, nested"]
        allowed = set(SCALAR + PSEUDOSCALAR + RIGHT_CURRENT + LEFT_CURRENT)
        for expression in row["display"]:
            self.assertTrue(expression.free_symbols <= allowed)

    def test_all_groupings_on_matrices(self):
        for name, passed in check_groupings_on_matrices().items():
            self.assertTrue(passed, name)

    def test_perturbed_grouping_is_detected(self):
        """A wrong coefficient must fail both the symbolic and the matrix check,
        and must not be selected as the shortest chiral grouping."""
        values = check_groupings(chiral_groupings={"perturbed": perturbed_chiral_nested},
                                 vector_axial_groupings={"vector/axial, nested":
                                                         grouping_vector_axial_nested})
        row = values["table"]["perturbed"]
        self.assertNotEqual(row["scalar_residual"], 0)
        self.assertIsNone(values["shortest_chiral"])
        matrices = check_groupings_on_matrices(groupings={"perturbed": (perturbed_chiral_nested,
                                                                        True)})
        self.assertFalse(matrices["perturbed"])


class AldaSpecializationTests(unittest.TestCase):
    def test_compact_forms_and_raw_dependence(self):
        values = derive_alda_specialization()
        for key in ("scalar_residual", "pseudoscalar_residual",
                    "scalar_invariant_residual", "pseudoscalar_invariant_residual"):
            self.assertEqual(values[key], 0, key)
        self.assertTrue(values["scalar_depends_on_raw"])
        self.assertTrue(values["pseudoscalar_depends_on_raw"])

    def test_perturbed_compact_form_is_detected(self):
        """Replacing 2c^S - c^V by c^S - c^V in hat S_2 must leave a residual."""
        def perturbed(bracket=commutator, anti_bracket=anticommutator):
            scalar, pseudoscalar = alda_compact_forms(bracket, anti_bracket)
            return (scalar - bracket(VECTOR_COUPLING, bracket(SCALAR_COUPLING, QUARK_MASS)),
                    pseudoscalar)

        values = derive_alda_specialization(compact_forms=perturbed)
        self.assertNotEqual(values["scalar_residual"], 0)
        self.assertEqual(values["pseudoscalar_residual"], 0)


def perturbed_real_mass_forms(mass, scalar, pseudoscalar, vector, axial,
                              bracket=commutator, anti_bracket=anticommutator):
    """Real-mass candidates with -{C_A1, m_q} replaced by -(1/2){C_A1, m_q} in hat P_2."""
    dressed_scalar, dressed_pseudoscalar = real_mass_vector_axial_forms(
        mass, scalar, pseudoscalar, vector, axial, bracket, anti_bracket)
    return dressed_scalar, dressed_pseudoscalar + anti_bracket(axial[1], mass) / 2


def wrong_cross_term(vector, axial, mass, bracket, anti_bracket):
    """Leibniz rewriting with the wrong factor 1 instead of 2 in front of {C_A0, [C_V0, m_q]}."""
    return (anti_bracket(bracket(vector, axial), mass)
            + anti_bracket(axial, bracket(vector, mass)))


class RealMassSpecializationTests(unittest.TestCase):
    def test_symbolic(self):
        values = derive_real_mass_specialization()
        for key in ("scalar_residual", "pseudoscalar_residual", "chiral_scalar_residual",
                    "chiral_pseudoscalar_residual"):
            self.assertEqual(values[key], 0, key)
        for name, residual in values["cross_residuals"].items():
            self.assertEqual(residual, 0, name)
        self.assertGreater(values["cross_monomial_number"],
                           values["largest_single_bracket_monomial_number"])

    def test_matrices(self):
        values = check_real_mass_on_matrices()
        self.assertTrue(values["noncommuting"])
        for key in ("scalar_residual", "pseudoscalar_residual", "chiral_scalar_residual",
                    "chiral_pseudoscalar_residual"):
            self.assertTrue(is_zero(values[key]), key)
        for name, residual in values["cross_residuals"].items():
            self.assertTrue(is_zero(residual), name)

    def test_perturbed_specialized_form_is_rejected(self):
        values = derive_real_mass_specialization(forms=perturbed_real_mass_forms,
                                                 cross_forms={"wrong": wrong_cross_term})
        self.assertEqual(values["scalar_residual"], 0)
        self.assertNotEqual(values["pseudoscalar_residual"], 0)
        self.assertNotEqual(values["cross_residuals"]["wrong"], 0)
        matrices = check_real_mass_on_matrices(forms=perturbed_real_mass_forms,
                                               cross_forms={"wrong": wrong_cross_term})
        self.assertFalse(is_zero(matrices["pseudoscalar_residual"]))
        self.assertFalse(is_zero(matrices["cross_residuals"]["wrong"]))


class IndependentCheckTests(unittest.TestCase):
    def assert_invariant(self, values, include_photon):
        self.assertTrue(values["group_is_unitary"])
        self.assertTrue(values["generators_noncommuting"])
        self.assertTrue(values["generator_traces_nonzero"])
        self.assertTrue(values["raw_scalar_changes"])
        self.assertTrue(values["raw_theta_changes"])
        for key, residual in values["residuals"].items():
            if key == "theta_photon_two" and not include_photon:
                continue
            self.assertTrue(is_zero(residual), key)

    def test_invariance_general_redefinition(self):
        self.assert_invariant(check_invariance(), include_photon=False)

    def test_invariance_general_redefinition_other_seeds(self):
        self.assert_invariant(check_invariance(seed=5, redefinition_seed=11),
                              include_photon=False)

    def test_invariance_charge_conserving_redefinition(self):
        values = check_invariance(charge_conserving=True)
        self.assert_invariant(values, include_photon=True)
        self.assertTrue(values["raw_theta_photon_changes"])
        self.assertTrue(all(r == 0 for r in values["photon_trace_form_residual"]))

    def test_constant_connections(self):
        for key, residual in check_constant_connections().items():
            self.assertTrue(is_zero(residual), key)

    def test_numerical_transport(self):
        values = check_numerical_transport()
        self.assertTrue(values["passed"], values["relative_deviation"])

    def test_first_order_regression(self):
        for key, residual in check_first_order_regression().items():
            self.assertEqual(residual, 0, key)

    def test_vacuum_angle(self):
        values = check_vacuum_angle()
        self.assertTrue(all(r == 0 for r in values["residuals"]))
        self.assertTrue(values["second_order_phase_changes"])

    def test_chiral_trace_forms(self):
        for key, residual in check_chiral_trace_forms().items():
            self.assertEqual(residual, 0, key)

    def test_wrong_trace_normalizations_are_detected(self):
        """Power-series-like normalizations (1/4, N_c/2, 4/3, 1/3) must leave residuals."""
        values = check_chiral_trace_forms(theta_coefficient=sp.Rational(1, 4),
                                          photon_coefficient=sp.Rational(3, 2),
                                          component_coefficients=(sp.Rational(4, 3),
                                                                  sp.Rational(1, 3)))
        for key, residual in values.items():
            self.assertNotEqual(residual, 0, key)

    def test_hermiticity(self):
        values = check_hermiticity()
        self.assertEqual(values["symbolic_scalar"], 0)
        self.assertEqual(values["symbolic_pseudoscalar"], 0)
        self.assertTrue(values["numeric_scalar"])
        self.assertTrue(values["numeric_pseudoscalar"])


if __name__ == "__main__":
    unittest.main()
