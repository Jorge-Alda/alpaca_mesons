"""Tests of the K pi saturation freedom derivation for the K*(892) couplings.

Run: .venv/bin/python -m unittest scripts.test_derive_kstar_kpi_saturation_freedom
"""
import sys
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import derive_kstar_kpi_saturation_freedom as derivation
import fit_babar_kkpi_global_vector_sector as joint_fit


class SymbolicScaling(unittest.TestCase):
    def test_width_matches_code_formula_and_vertex_roots(self):
        symbolic = derivation.symbolic_width_scaling()
        kappa = sp.Symbol('kappa', positive=True)
        self.assertEqual(set(symbolic['vertex_roots']), {sp.sqrt(kappa), -sp.sqrt(kappa)})

    def test_isospin_limit_traces(self):
        traces = derivation.symbolic_flavor_traces()['isospin_limit']
        self.assertEqual(traces['kstar_charged_K+_pi0'], sp.Rational(1, 2))
        self.assertEqual(traces['kstar_neutral_K+_pi-'], 1)


class AmplitudeSubstitution(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.values = derivation.accepted_endpoint()
        cls.model = joint_fit.GlobalLikelihood(1, 3, 4)
        cls.state = cls.model.unpack(cls.values)
        cls.baseline_residual, _ = cls.model._babar_blocks_full(cls.values, cls.state)

    def test_unit_kappa_is_baseline(self):
        _, residual = derivation.babar_amplitudes(self.model, self.values, self.state)
        self.assertTrue(np.array_equal(residual, self.baseline_residual))
        self.assertIs(joint_fit.running_propagator, derivation.BASELINE_RUNNING_PROPAGATOR)

    def test_neutral_kappa_absent_from_charged_kaon_channel(self):
        totals = self.state[4]
        fixed = (totals['Kstarplus'], totals['Kstar0'])
        reference, _ = derivation.babar_amplitudes(self.model, self.values, self.state, kstar_widths=fixed)
        varied, _ = derivation.babar_amplitudes(self.model, self.values, self.state, 1., 2.,
                                                kstar_widths=fixed)
        np.testing.assert_array_equal(reference['kkpi0'], varied['kkpi0'])
        self.assertGreater(np.max(abs(reference['kskpi']-varied['kskpi'])), 0.)

    def test_saturated_widths_reproduce_totals(self):
        _, _, mk, gv, totals = self.state
        widths = derivation.saturated_widths(self.model, mk, gv, derivation.symbolic_width_scaling())
        for name in derivation.KSTAR_NAMES:
            self.assertAlmostEqual(widths[name]/totals[name], 1., delta=1e-13)


if __name__ == '__main__':
    unittest.main()
