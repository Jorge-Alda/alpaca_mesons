"""Unit tests for the dressed photon coefficient hat(Theta_gamma).

Run: .venv/bin/python -m unittest scripts.test_derive_alp_dressed_photon_invariant
"""

import unittest

from scripts.derive_alp_dressed_photon_invariant import (
    check_application_identity,
    check_charge_mixing_failure,
    check_charge_preservation,
    check_independence_from_gluons,
    check_redefinition_invariance,
    check_taylor_form,
)


class DressedPhotonInvariantTests(unittest.TestCase):
    def test_charge_preservation(self):
        for key, value in check_charge_preservation().items():
            self.assertTrue(value, key)

    def test_redefinition_invariance(self):
        for mode, theta_gamma_expected in (("nonlinear", True), ("traceless", True),
                                     ("block_traceless", False)):
            values = check_redefinition_invariance(mode=mode)
            self.assertTrue(values["generators_commute_with_charge"], mode)
            self.assertTrue(values["ds_block_noncommuting"], mode)
            self.assertTrue(all(residual == 0 for residual in values["block_trace_residuals"]))
            self.assertEqual(values["anomaly_law_residual"], 0, mode)
            self.assertEqual(values["theta_gamma_changes"], theta_gamma_expected, mode)
            self.assertTrue(values["transport_covariance_is_zero"], mode)
            self.assertEqual(values["dressed_theta_gamma_residual"], 0, mode)
            self.assertEqual(values["dressed_theta_residual"], 0, mode)

    def test_redefinition_invariance_other_seeds(self):
        values = check_redefinition_invariance(seed=3, redefinition_seed=5)
        self.assertEqual(values["dressed_theta_gamma_residual"], 0)

    def test_taylor_form(self):
        values = check_taylor_form()
        self.assertTrue(all(residual == 0 for residual in values["taylor_residuals"]))
        for key in ("first_order_residual", "shift_residual", "explicit_form_residual"):
            self.assertEqual(values[key], 0, key)
        self.assertTrue(values["only_diagonal_entries"])

    def test_charge_mixing_failure(self):
        values = check_charge_mixing_failure()
        for key in ("right_charge_changed", "left_charge_changed", "charges_chiral",
                    "charge_x_dependent", "couplings_violate_charge"):
            self.assertTrue(values[key], key)
        for key in ("dressed_theta_gamma_residual", "homogeneous_residual",
                    "inhomogeneous_residual"):
            self.assertNotEqual(values[key], 0, key)
        self.assertEqual(values["lowest_residual_power"], 3)
        self.assertEqual(values["dressed_theta_residual"], 0)

    def test_independence_from_gluons(self):
        values = check_independence_from_gluons()
        self.assertEqual(values["weight_rank"], 2)
        self.assertEqual(values["shift_rank"], 2)
        self.assertEqual(values["ds_weight_rank"], 1)
        self.assertEqual(values["ds_proportionality_residual"], 0)
        self.assertTrue(values["full_not_proportional"])

    def test_application_identity(self):
        for key, value in check_application_identity().items():
            self.assertEqual(value, 0, key)


if __name__ == "__main__":
    unittest.main()
