"""Bitwise equality of IsospinBreakingVppLikelihoodWithReuse with IsospinBreakingVppLikelihood.

Run: .venv/bin/python -m unittest scripts.test_joint_vector_isospin_breaking_vpp_reuse
Fit quadrature (2,16,16) with four-bin blocks, default option switches.
Points: restart endpoint (run 2 of 'scan_point' in the fit ledger), the same
point with r = (0, 0, 0), and the same point with r at the scan point.
"""
import json
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import fit_joint_vector_isospin_breaking_vpp as driver  # noqa: E402
import joint_vector_isospin_breaking_vpp_option as vpp_option  # noqa: E402
import joint_vector_isospin_breaking_vpp_reuse as reuse_module  # noqa: E402

LEDGER = 'data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json'
FIT_ORDERS = (2, 16, 16)
ORIGINAL_RUNNING = joint_fit.running_propagator


def bits(array):
    """IEEE-754 bit patterns of a float or complex array."""
    return np.ascontiguousarray(array).view(np.uint64)


class BitwiseEqualityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        options = dict(bin_block_size=4, **driver.MODEL_OPTIONS)
        cls.plain = vpp_option.IsospinBreakingVppLikelihood(*FIT_ORDERS, **options)
        cls.reuse = reuse_module.IsospinBreakingVppLikelihoodWithReuse(*FIT_ORDERS, **options)
        _, run = driver.restart_run(json.loads((driver.ROOT/LEDGER).read_text()))
        endpoint = np.asarray(run['endpoint'])
        cls.points = dict(restart=endpoint, r_zero=np.r_[endpoint[:32], 0., 0., 0.],
                          scan_point=np.r_[endpoint[:32], vpp_option.SCAN_POINT_GEV_MINUS2])

    def test_residual_and_jacobian_bitwise(self):
        for name, point in self.points.items():
            for method in ('residual', 'jacobian'):
                with self.subTest(point=name, method=method):
                    expected = getattr(self.plain, method)(point)
                    result = getattr(self.reuse, method)(point)
                    self.assertEqual(result.shape, expected.shape)
                    self.assertTrue(np.array_equal(bits(result), bits(expected)))

    def test_driver_model_uses_reuse(self):
        self.assertTrue(issubclass(driver.FitLikelihood, reuse_module.IsospinBreakingVppLikelihoodWithReuse))

    def test_globals_restored(self):
        self.reuse.jacobian(self.points['restart'])
        self.assertIs(vpp_option.RUNNING_PROPAGATOR, ORIGINAL_RUNNING)
        self.assertIs(joint_fit.running_propagator, ORIGINAL_RUNNING)
        self.assertIs(vpp_option.s_channel_factors, reuse_module.S_CHANNEL_FACTORS)
        self.assertIsNone(self.reuse.jacobian_point)


class GroundPropagatorMemoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        model = joint_fit.GlobalLikelihood(1, 3, 4)
        cls.reference = dict(model.preliminary, m_rho_GeV=.775, m_Kstar_GeV=.8955,
                             Gamma_total_omega_GeV=.00868, Gamma_total_phi_GeV=.004249)
        cls.invariant = np.linspace(.5, 3., 41).reshape(41, 1)

    def test_full_and_zero_excitation_calls_bitwise(self):
        for name, mass, width in (('kstar', .8955, .0473), ('rho', .775, .149)):
            memory = reuse_module.GroundPropagatorMemory()
            arguments = (mass, width, self.reference)
            full = memory(name, self.invariant, (.3, -.2), *arguments)
            ground = memory(name, self.invariant, (0., 0.), *arguments)
            self.assertEqual((memory.evaluated_calls, memory.reused_calls), (1, 1))
            self.assertTrue(np.array_equal(bits(full), bits(ORIGINAL_RUNNING(
                name, self.invariant, (.3, -.2), *arguments))))
            self.assertTrue(np.array_equal(bits(ground), bits(ORIGINAL_RUNNING(
                name, self.invariant, (0., 0.), *arguments))))

    def test_no_reuse_when_call_differs(self):
        memory = reuse_module.GroundPropagatorMemory()
        memory('kstar', self.invariant, (.3, -.2), .8955, .0473, self.reference)
        other = self.invariant.copy()
        for arguments in (('kstar', other, (0., 0.), .8955, .0473, self.reference),
                          ('kstar', self.invariant, (0., 0.), .8955, .05, self.reference),
                          ('rho', self.invariant, (0., 0.), .8955, .0473, self.reference)):
            memory.stored = ('kstar', self.invariant, .8955, .0473, self.reference, None, None)
            result = memory(*arguments)
            self.assertTrue(np.array_equal(bits(result), bits(ORIGINAL_RUNNING(*arguments))))
        self.assertEqual(memory.reused_calls, 0)


class FactorMemoryTests(unittest.TestCase):
    def test_copies_and_values(self):
        model = reuse_module.IsospinBreakingVppLikelihoodWithReuse(1, 3, 4, isospin_breaking_vpp=True)
        plain = vpp_option.IsospinBreakingVppLikelihood(1, 3, 4, isospin_breaking_vpp=True)
        coefficients = np.array(vpp_option.SCAN_POINT_GEV_MINUS2)
        first = model.factors(coefficients)
        first['changed'] = 1.
        self.assertEqual(model.factors(coefficients), plain.factors(coefficients))
        generators = [np.diag([1., 1., 0.])/np.sqrt(2), np.diag([0., 0., 1.])]
        self.assertEqual(model.s_channel_factors(coefficients, model.spurion_diagonal, *generators),
                         vpp_option.s_channel_factors(coefficients, model.spurion_diagonal, *generators))


if __name__ == '__main__':
    unittest.main()
