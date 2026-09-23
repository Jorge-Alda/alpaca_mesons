"""Tests for literature K-pi elastic phase parametrizations."""

import unittest

import numpy as np

from scripts.implementation.kpi_phases import PelaezRodas2016ElasticPhase
from scripts.implementation.omnes import OmnesFunction


class KPiPhaseTests(unittest.TestCase):
    def test_threshold_signs_and_domains(self) -> None:
        half = PelaezRodas2016ElasticPhase("1/2")
        three_half = PelaezRodas2016ElasticPhase("3/2")
        self.assertEqual(half(half.threshold), 0.0)
        self.assertEqual(three_half(three_half.threshold), 0.0)
        self.assertGreater(half(0.8), 0.0)
        self.assertLess(three_half(0.8), 0.0)
        with self.assertRaisesRegex(ValueError, "valid only"):
            half(half.maximum + 1e-5)

    def test_cfd_and_ufd_are_distinct(self) -> None:
        cfd = PelaezRodas2016ElasticPhase("1/2", "CFD")
        ufd = PelaezRodas2016ElasticPhase("1/2", "UFD")
        self.assertNotAlmostEqual(cfd(0.9), ufd(0.9), places=8)

    def test_phase_builds_normalized_omnes_function(self) -> None:
        phase = PelaezRodas2016ElasticPhase("1/2")
        endpoint = phase.maximum
        omega = OmnesFunction(
            phase=phase,
            threshold=phase.threshold,
            phase_input_max=endpoint,
            subtraction_point=0.0,
            tail="constant",
            # Assumption: freeze the phase above the K-eta matching point.
            asymptotic_phase=phase(endpoint),
        )
        self.assertEqual(omega(0.0), 1.0 + 0.0j)
        value = omega(0.8)
        self.assertTrue(np.isfinite(value))
        self.assertAlmostEqual(np.angle(value), phase(0.8), places=9)


if __name__ == "__main__":
    unittest.main()
