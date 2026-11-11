"""Unit tests for the dressed topological coefficient hat(Theta).

Run: .venv/bin/python -m unittest scripts.test_derive_alp_dressed_topological_invariant
"""

import unittest

from scripts.derive_alp_dressed_topological_invariant import (
    check_first_order_vacuum_angle,
    check_liouville,
    check_meson_topological_combination,
    check_redefinition_invariance,
    check_tadpole_source,
    check_taylor_form,
    check_vacuum_angle,
)


class DressedTopologicalInvariantTests(unittest.TestCase):
    def test_liouville(self):
        for chirality, values in check_liouville().items():
            self.assertTrue(values["ode_is_zero"], chirality)
            self.assertEqual(values["determinant_residual"], 0, chirality)

    def test_redefinition_invariance(self):
        for traceless in (False, True):
            values = check_redefinition_invariance(traceless=traceless)
            self.assertEqual(values["anomaly_shift_residual"], 0)
            self.assertEqual(values["theta_changes"], not traceless)
            self.assertTrue(values["transport_covariance_is_zero"])
            self.assertEqual(values["dressed_theta_residual"], 0)

    def test_redefinition_invariance_other_seeds(self):
        values = check_redefinition_invariance(seed=3, redefinition_seed=5)
        self.assertEqual(values["dressed_theta_residual"], 0)

    def test_taylor_form(self):
        values = check_taylor_form()
        self.assertTrue(all(residual == 0 for residual in values["taylor_residuals"]))
        self.assertEqual(values["first_order_residual"], 0)
        self.assertEqual(values["shift_residual"], 0)

    def test_vacuum_angle(self):
        for key, value in check_vacuum_angle().items():
            self.assertEqual(value, 0, key)

    def test_first_order_vacuum_angle(self):
        for case, values in check_first_order_vacuum_angle().items():
            for key, value in values.items():
                self.assertEqual(value, 0, f"{case}: {key}")

    def test_meson_topological_combination(self):
        values = check_meson_topological_combination()
        self.assertTrue(values["goldstone_is_unitary"])
        for key in ("unitary_determinant_residual", "factorization_residual",
                    "trace_log_shift_residual", "topological_combination_residual"):
            self.assertEqual(values[key], 0, key)

    def test_tadpole_source(self):
        for key, value in check_tadpole_source().items():
            self.assertEqual(value, 0, key)


if __name__ == "__main__":
    unittest.main()
