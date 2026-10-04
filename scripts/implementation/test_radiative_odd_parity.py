"""Tests for the on-shell omega -> pi0 gamma RChT predictor.

Invocation:
    .venv/bin/python -m unittest scripts.implementation.test_radiative_odd_parity
"""

import math
import unittest
from dataclasses import FrozenInstanceError, replace

import numpy as np

from scripts.implementation.odd_parity_rcht import OddParityCouplings, VJPCouplings, VVPCouplings
from scripts.implementation.radiative_odd_parity import (
    OmegaPiGammaInputs,
    G_omega_pi_gamma,
    constrained_G_omega_pi_gamma,
    omega_pi_gamma_direct_coupling,
    omega_pi_gamma_rho_coupling,
    omega_pi_gamma_width,
    omega_pi_gamma_branching_fraction_predictor,
    omega_pi_gamma_width_predictor,
    vector_to_pseudoscalar_gamma_width,
)


class OmegaPiGammaTests(unittest.TestCase):
    def setUp(self):
        self.inputs = OmegaPiGammaInputs(0.3028, 0.092, 0.145, 0.770, 0.783, 0.135)
        self.couplings = OddParityCouplings(
            VJPCouplings(0.01, -0.02, 0.03, 91.0, 0.04, -0.05, -73.0),
            VVPCouplings(0.06, -0.01, 0.02, 111.0),
        )

    def test_direct_eq27_transcription_and_absent_couplings(self):
        x = self.inputs
        c, d = self.couplings.vjp, self.couplings.vvp
        expected_direct = 2 * math.sqrt(2) * x.e / (x.Momega * x.MV * x.F) * (
            (c.c2 - c.c1 + c.c5 - 2 * c.c6) * x.Momega**2
            + (c.c1 + c.c2 + 8 * c.c3 - c.c5) * x.mpi**2
        )
        expected_rho = -4 * x.e * x.FV / (x.MV**2 * x.Momega * x.F) * (
            d.d3 * x.Momega**2 + (d.d1 + 8 * d.d2 - d.d3) * x.mpi**2
        )
        self.assertAlmostEqual(omega_pi_gamma_direct_coupling(self.couplings, x), expected_direct)
        self.assertAlmostEqual(omega_pi_gamma_rho_coupling(self.couplings, x), expected_rho)
        self.assertAlmostEqual(G_omega_pi_gamma(self.couplings, x), expected_direct + expected_rho)
        changed = OddParityCouplings(
            replace(c, c4=-1e7, c7=2e7), replace(d, d4=-3e7)
        )
        self.assertEqual(G_omega_pi_gamma(changed, x), G_omega_pi_gamma(self.couplings, x))

    def test_all_five_ope_substitutions_equal_eq28(self):
        x, n_c = self.inputs, 3
        c1, c5, d2 = 0.071, -0.023, 0.014
        short_c = n_c * x.MV / (64 * math.pi**2 * math.sqrt(2) * x.FV)
        short_d = n_c * x.MV**2 / (64 * math.pi**2 * x.FV**2)
        couplings = OddParityCouplings(
            VJPCouplings(c1, c1 + c5, -c1 / 4, 9.0, c5, c5 - short_c, -8.0),
            VVPCouplings(
                -short_d + x.F**2 / (4 * x.FV**2) - 8 * d2,
                d2,
                -short_d + x.F**2 / (8 * x.FV**2),
                13.0,
            ),
        )
        self.assertAlmostEqual(
            G_omega_pi_gamma(couplings, x), constrained_G_omega_pi_gamma(x, n_c), places=14
        )

    def test_width_reproduces_eq29(self):
        x = replace(self.inputs, FV=math.sqrt(2) * self.inputs.F)
        width = vector_to_pseudoscalar_gamma_width(
            constrained_G_omega_pi_gamma(x), x.Momega, x.mpi
        )
        alpha = x.e**2 / (4 * math.pi)
        expected = alpha / 192 * x.Momega * (1 - x.mpi**2 / x.Momega**2) ** 3 * (
            3 / (4 * math.pi**2) * x.Momega**2 / x.F**2
            - x.Momega**2 / x.MV**2 * (1 + x.mpi**2 / x.Momega**2)
        ) ** 2
        self.assertAlmostEqual(width, expected, places=16)

    def test_scaling_and_signed_direct_rho_interference(self):
        direct = omega_pi_gamma_direct_coupling(self.couplings, self.inputs)
        rho = omega_pi_gamma_rho_coupling(self.couplings, self.inputs)
        self.assertGreater(direct, 0.0)
        self.assertLess(rho, 0.0)
        self.assertAlmostEqual(G_omega_pi_gamma(self.couplings, self.inputs), direct + rho)
        scale = -2.5
        scaled = OddParityCouplings(
            VJPCouplings(*(scale * value for value in self.couplings.vjp.as_array())),
            VVPCouplings(*(scale * value for value in self.couplings.vvp.as_array())),
        )
        self.assertAlmostEqual(G_omega_pi_gamma(scaled, self.inputs), scale * (direct + rho))
        self.assertAlmostEqual(
            omega_pi_gamma_width(scaled, self.inputs), scale**2 * omega_pi_gamma_width(self.couplings, self.inputs)
        )

    def test_predictor_shape_and_validation(self):
        parameters = {
            **{f"c{i}": value for i, value in enumerate(self.couplings.vjp.as_array(), 1)},
            **{f"d{i}": value for i, value in enumerate(self.couplings.vvp.as_array(), 1)},
        }
        config = {name: getattr(self.inputs, name) for name in ("e", "F", "FV", "MV", "Momega", "mpi")}
        prediction = omega_pi_gamma_width_predictor(parameters, config)
        self.assertEqual(prediction.shape, (1,))
        self.assertFalse(prediction.flags.writeable)
        self.assertEqual(prediction[0], omega_pi_gamma_width(self.couplings, self.inputs))
        with self.assertRaisesRegex(ValueError, "missing required keys"):
            omega_pi_gamma_width_predictor({key: value for key, value in parameters.items() if key != "c7"}, config)
        with self.assertRaisesRegex(ValueError, "missing required keys"):
            omega_pi_gamma_width_predictor(parameters, {key: value for key, value in config.items() if key != "MV"})
        with self.assertRaisesRegex(ValueError, "unexpected keys"):
            omega_pi_gamma_width_predictor(parameters, {**config, "Zomega": 1.0})
        with self.assertRaisesRegex(ValueError, "finite"):
            omega_pi_gamma_width_predictor({**parameters, "d4": np.nan}, config)

        total_width = 0.00849
        branching_config = {**config, "GammaOmegaTotal": total_width}
        branching = omega_pi_gamma_branching_fraction_predictor(parameters, branching_config)
        self.assertEqual(branching.shape, (1,))
        np.testing.assert_allclose(branching, [omega_pi_gamma_width(self.couplings, self.inputs) / total_width])
        self.assertFalse(branching.flags.writeable)
        with self.assertRaisesRegex(ValueError, "GammaOmegaTotal"):
            omega_pi_gamma_branching_fraction_predictor(parameters, config)
        for invalid in (0.0, -1.0, np.inf):
            with self.assertRaisesRegex(ValueError, "positive|finite"):
                omega_pi_gamma_branching_fraction_predictor(
                    parameters, {**config, "GammaOmegaTotal": invalid}
                )

    def test_frozen_and_nonfinite_or_unphysical_inputs_rejected(self):
        with self.assertRaises(FrozenInstanceError):
            self.inputs.MV = 1.0
        with self.assertRaisesRegex(ValueError, "finite"):
            OmegaPiGammaInputs(np.inf, 0.092, 0.13, 0.77, 0.783, 0.135)
        with self.assertRaisesRegex(ValueError, "open decay"):
            replace(self.inputs, mpi=self.inputs.Momega)
        with self.assertRaisesRegex(ValueError, "finite scalar"):
            vector_to_pseudoscalar_gamma_width(complex(np.nan, 0), 1.0, 0.1)
        with self.assertRaisesRegex(ValueError, "width must be finite"):
            vector_to_pseudoscalar_gamma_width(1e308, 1.0, 0.1)


if __name__ == "__main__":
    unittest.main()
