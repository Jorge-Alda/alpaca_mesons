"""Tests for pole tracking along m_a0^2 and f_pi/f_a (task 2.4)."""

import unittest

import numpy as np
import sympy as sp

from derive_alp_lo_neutral_blocks import illustrative_lo_inputs
from diagnose_alp_meson_spectrum import (
    UnstableVacuumError,
    isolated_two_level_expectation,
    isolated_two_level_inputs,
    lo_bare_mass_stability_bound,
)
from solve_alp_meson_spectrum import solve_spectrum_numerical
from track_alp_meson_poles import (
    avoided_crossings,
    derive,
    fast_lo_blocks_builder,
    lo_blocks_builder,
    track_poles,
)


class TrackingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = derive()
        cls.isolated = isolated_two_level_inputs()
        cls.isolated_builder = staticmethod(fast_lo_blocks_builder(cls.isolated, "bare_alp_mass_squared"))
        cls.isolated_track = track_poles(np.linspace(0.002, 1.5, 31), cls.isolated_builder)

    def test_all_scan_checks(self):
        for name, summary in self.results.items():
            for label, passed in summary["checks"].items():
                self.assertTrue(passed, f"{name}: {label}")

    def test_fast_builder_matches_exact_blocks(self):
        inputs = illustrative_lo_inputs("neutral")
        for name, value in (("bare_alp_mass_squared", 0.41), ("ratio", 0.13)):
            fast = fast_lo_blocks_builder(inputs, name)(value)
            exact = lo_blocks_builder(inputs, name)(value)
            for left, right in zip(fast, exact):
                np.testing.assert_allclose(left, right, rtol=1e-13, atol=1e-15)

    def test_tracked_vectors_solve_the_eigenproblem(self):
        tracked = self.isolated_track
        for index, parameter in enumerate(tracked.parameter_values):
            kinetic, mass_squared = self.isolated_builder(parameter)
            vectors = tracked.pole_vectors[index]
            np.testing.assert_allclose(vectors.T @ kinetic @ vectors, np.eye(6), atol=1e-10)
            np.testing.assert_allclose(mass_squared @ vectors, kinetic @ vectors * tracked.poles[index], atol=1e-10)
            np.testing.assert_allclose(np.sort(tracked.poles[index]),
                                       solve_spectrum_numerical(kinetic, mass_squared).pole_squared_masses, atol=1e-12)

    def test_continuity_of_tracked_vectors(self):
        tracked = self.isolated_track
        for index in range(len(tracked.parameter_values) - 1):
            metric = (tracked.kinetic[index] + tracked.kinetic[index + 1]) / 2
            overlaps = np.diag(tracked.pole_vectors[index].T @ metric @ tracked.pole_vectors[index + 1])
            self.assertGreater(np.min(overlaps), 0.8)

    def test_alp_like_branch_and_two_level_gap(self):
        tracked = self.isolated_track
        alp_like = tracked.alp_like_branch()
        self.assertEqual(alp_like[0], 0)
        self.assertEqual(alp_like[-1], 1)
        self.assertEqual(tracked.alp_like_switches()[0][1:], (0, 1))
        # Branch 1 is the ALP-like state at the end: it crossed the decoupled K0 doublet, eta and eta'.
        self.assertAlmostEqual(tracked.poles[-1, 1], max(tracked.poles[-1]), places=12)
        self.assertLess(tracked.poles[0, 1], tracked.poles[0, 2])
        crossing = avoided_crossings(tracked, self.isolated_builder)
        self.assertEqual(len(crossing), 1)
        splitting = float(isolated_two_level_expectation(self.isolated)["splitting"])
        self.assertAlmostEqual(crossing[0]["minimal_gap"], splitting, places=10)
        # Both pole branches through the avoided crossing stay ordered (no label crossing).
        self.assertTrue(np.all(tracked.poles[:, 0] < tracked.poles[:, 1]))

    def test_decoupled_kaon_doublet_keeps_its_plane(self):
        tracked = self.isolated_track
        for branch in (2, 3):
            components = tracked.pole_vectors[:, :, branch]
            self.assertLess(np.max(np.abs(components[:, [0, 1, 2, 5]])), 1e-12)

    def test_coarse_grid_without_refinement_swaps_labels(self):
        coarse = track_poles(np.linspace(0.002, 0.05, 3), self.isolated_builder, minimal_step=1.0)
        self.assertTrue(coarse.unresolved_crossings())
        refined = track_poles(np.linspace(0.002, 0.05, 3), self.isolated_builder)
        self.assertFalse(refined.unresolved_crossings())
        self.assertGreater(refined.refinement_count, 0)
        self.assertEqual(refined.order_exchanges_between_alp_carrying_branches(), [])

    def test_unstable_point_stops_tracking(self):
        inputs = illustrative_lo_inputs("neutral")
        bound = float(lo_bare_mass_stability_bound(inputs)["bound"])
        builder = fast_lo_blocks_builder(inputs, "bare_alp_mass_squared")
        with self.assertRaises(UnstableVacuumError) as context:
            track_poles(np.linspace(0.6, bound - 0.01, 5), builder)
        self.assertLess(context.exception.parameter, bound)
        self.assertIn("negative eigenvalue", str(context.exception))

    def test_epsilon_scan_starts_at_unmixed_spectrum(self):
        inputs = illustrative_lo_inputs("flavor_diagonal")
        builder = fast_lo_blocks_builder(inputs, "ratio")
        tracked = track_poles(np.linspace(1e-4, 0.3, 20), builder)
        kinetic, mass_squared = builder(1e-4)
        unmixed = np.sort(np.append(np.linalg.eigvalsh(mass_squared[:5, :5]), float(inputs["bare_alp_mass_squared"])))
        np.testing.assert_allclose(np.sort(tracked.poles[0]), unmixed, atol=1e-6)
        self.assertFalse(tracked.unresolved_crossings())
        self.assertEqual(tracked.order_exchanges_between_alp_carrying_branches(), [])
        self.assertAlmostEqual(float(sp.Rational(3, 10)), tracked.poles[0, tracked.alp_like_branch()[0]], places=6)


if __name__ == "__main__":
    unittest.main()
