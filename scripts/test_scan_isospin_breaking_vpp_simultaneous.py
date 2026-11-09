"""Tests of the simultaneous r_1, r_2, r_3 scan of the isospin-breaking VPP model.

Run: .venv/bin/python -m unittest scripts.test_scan_isospin_breaking_vpp_simultaneous
"""
import sys
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import derive_isospin_breaking_vpp_operators as derivation  # noqa: E402
import scan_isospin_breaking_vpp_simultaneous as scan  # noqa: E402


class LinearReduction(unittest.TestCase):
    def test_isospin_limit_splits(self):
        splits = scan.symbolic_linear_reduction()
        self.assertEqual(sp.simplify(splits['kstar'] - (scan.R1 - scan.R2/3)*(scan.XU - scan.XD)), 0)
        self.assertEqual(sp.simplify(splits['phi'] - (scan.R2 - scan.R3)*(scan.XU - scan.XD)), 0)

    def test_code_inputs_gradient_kstar_r1_and_r3(self):
        gradients = scan.numerical_linear_gradients(derivation.vacuum_spurion_values())
        self.assertAlmostEqual(gradients['kstar'][0], 1.0, places=4)
        self.assertAlmostEqual(gradients['kstar'][2], 0.0, places=8)


class ExactFamily(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spurions = derivation.vacuum_spurion_values()
        cls.target = scan.targets()
        cls.coefficients = scan.trace_coefficients(cls.spurions)

    def test_family_points_satisfy_both_splits(self):
        rows = scan.family_table(self.spurions, self.target, np.linspace(-2, 2, 21))
        self.assertEqual(len(rows), 42)
        for row in rows:
            self.assertAlmostEqual(row['kstar_ratio'], self.target['kstar'], places=10)
            self.assertAlmostEqual(row['phi_change'], self.target['phi'], places=10)
            if max(abs(row['r2']), abs(row['r3'])) <= 5:  # published pi0 only; grows with |r2|
                self.assertLess(abs(row['rho_split'] - 1), 2e-3)

    def test_family_contains_record_two_coefficient_solutions(self):
        for r1, r2_expected, r3_expected in ((0.0, 3.7362, 1.2482), (-0.8415, -0.7871, 0.0)):
            rows = scan.family_table(self.spurions, self.target, [r1])
            closest = min(rows, key=lambda row: abs(row['r2'] - r2_expected))
            self.assertAlmostEqual(closest['r2'], r2_expected, places=3)
            self.assertAlmostEqual(closest['r3'], r3_expected, places=3)

    def test_record_single_coefficient_values(self):
        values = scan.observables(self.coefficients, 0, 3.736157, 0)
        self.assertAlmostEqual(values['kstar_ratio'], self.target['kstar'], places=5)
        self.assertAlmostEqual(values['phi_change'], -0.1625, places=3)
        self.assertAlmostEqual(values['record_phi_over_kstar'], 0.8232, places=3)

    def test_kstar_ratio_independent_of_r3(self):
        first = scan.observables(self.coefficients, -0.3, 1.2, 0.0)['kstar_ratio']
        second = scan.observables(self.coefficients, -0.3, 1.2, 4.0)['kstar_ratio']
        self.assertAlmostEqual(first, second, places=13)

    def test_rho_charge_ratio_isospin_limit(self):
        traces = derivation.effective_traces((1, -1, 0), 0)
        self.assertEqual(sp.simplify(traces['rho+_pi+_pi0'] + traces['rho0_pi+_pi-']), 0)


if __name__ == '__main__':
    unittest.main()
