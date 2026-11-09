"""Tests of the isospin-breaking VPP operator derivation.

Run: .venv/bin/python -m unittest scripts.test_derive_isospin_breaking_vpp_operators
"""
import sys
import unittest
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import derive_isospin_breaking_vpp_operators as derivation  # noqa: E402

R1, R2, R3 = (derivation.SYMBOLS[k] for k in ('r1', 'r2', 'r3'))


class StructureCounting(unittest.TestCase):
    def test_mass_spurion_ranks(self):
        for dimension in (3, 5):
            rank, _ = derivation.structure_rank(('V', 'P1', 'P2', 'S'), dimension)
            single, _ = derivation.structure_rank(('V', 'P1', 'P2', 'S'), dimension, single_trace_only=True)
            self.assertEqual((rank, single), (5, 2))

    def test_chi_minus_operator_is_redundant(self):
        self.assertTrue(derivation.check_chi_minus_reduction())


class ZeroCoefficientReduction(unittest.TestCase):
    def test_flavor_factors_recovered(self):
        differences = derivation.check_zero_coefficient_reduction()
        self.assertLess(max(differences.values()), 1e-13)

    def test_rho_charge_states_equal_in_isospin_limit(self):
        traces = derivation.effective_traces((1, -1, 0), 0)
        ratio = sp.simplify(traces['rho+_pi+_pi0']/traces['rho0_pi+_pi-'])
        self.assertEqual(ratio, -1)

    def test_leading_phi_ratio_matches_endpoint(self):
        self.assertAlmostEqual(derivation.check_phi_leading_ratio(), 0.5084953/0.3355123, places=5)

    def test_linear_splits(self):
        splits = derivation.linear_isospin_limit_splits()
        self.assertEqual(splits['rho'], 0)


class CoefficientSolve(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.target = derivation.required_kstar_ratio()
        cls.ratio, _ = derivation.kstar_width_ratio_function()

    def test_vector_singlet_term_leaves_kstar_ratio(self):
        self.assertEqual(sp.simplify(sp.diff(self.ratio, R3)), 0)

    def test_single_coefficient_roots(self):
        expected = {'r1': -0.695084, 'r2': 3.736157}
        for name, value in expected.items():
            root, _ = derivation.solve_single_coefficient(name, self.target)
            self.assertAlmostEqual(root, value, places=5)
            others = {s: 0 for s in (R1, R2, R3) if s != derivation.SYMBOLS[name]}
            achieved = float(self.ratio.subs(others).subs(derivation.SYMBOLS[name], root))
            self.assertAlmostEqual(achieved, self.target, places=12)


if __name__ == '__main__':
    unittest.main()
