"""Tests for the stability and near-degeneracy diagnostics (tasks 2.2 and 2.3)."""

import unittest

import numpy as np
import sympy as sp

from derive_alp_lo_neutral_blocks import NEUTRAL_FIELD_LABELS, closed_form_blocks, illustrative_lo_inputs
from diagnose_alp_meson_spectrum import (
    KINETIC_CONDITION,
    MASS_CONDITION,
    POLE_CONDITION,
    UnstableVacuumError,
    derive,
    exact_degenerate_groups,
    exact_stability_diagnostics,
    find_near_degeneracies,
    isolated_two_level_expectation,
    isolated_two_level_inputs,
    lo_bare_mass_stability_bound,
    lo_blocks_float,
    lo_kinetic_matrix_checks,
    single_pole_formula_singular,
    solve_physical_spectrum,
    solve_physical_spectrum_exact,
    stability_diagnostics,
    two_level_estimates,
)
from solve_alp_meson_spectrum import solve_spectrum_exact, solve_spectrum_numerical
from verify_alp_lo_spectrum import random_symmetric_quadratic_form


class StabilityTest(unittest.TestCase):
    def test_lo_kinetic_matrix(self):
        for label, passed in lo_kinetic_matrix_checks()["checks"].items():
            self.assertTrue(passed, label)

    def test_lo_kinetic_determinant_numerically(self):
        generator = np.random.default_rng(5)
        for _ in range(20):
            mixing = generator.normal(scale=3.0, size=5)
            kinetic = np.eye(6)
            kinetic[:5, 5] = kinetic[5, :5] = mixing
            kinetic[5, 5] = 1 + mixing @ mixing
            diagnostics = stability_diagnostics(kinetic, np.eye(6))
            self.assertAlmostEqual(diagnostics.kinetic_determinant, 1.0, delta=1e-9 * (1 + mixing @ mixing)**2)
            self.assertGreater(diagnostics.kinetic_smallest_eigenvalue, 0)
            self.assertTrue(diagnostics.kinetic_positive_definite)

    def test_all_derive_checks(self):
        for group, checks in derive().items():
            for label, value in checks.items():
                if isinstance(value, bool):
                    self.assertTrue(value, f"{group}: {label}")

    def test_unstable_vacuum_names_condition(self):
        inputs = illustrative_lo_inputs("neutral")
        bound = lo_bare_mass_stability_bound(inputs)["bound"]
        unstable = dict(inputs, bare_alp_mass_squared=bound - sp.Rational(1, 100))
        with self.assertRaises(UnstableVacuumError) as context:
            solve_physical_spectrum(*lo_blocks_float(unstable))
        self.assertIn(MASS_CONDITION, context.exception.violated_conditions)
        self.assertNotIn(KINETIC_CONDITION, context.exception.violated_conditions)
        self.assertIn("negative eigenvalue", str(context.exception))
        self.assertLess(context.exception.diagnostics.mass_squared_smallest_eigenvalue, 0)

    def test_unstable_vacuum_exact(self):
        inputs = isolated_two_level_inputs()
        bound = lo_bare_mass_stability_bound(inputs)["bound"]
        kinetic, mass_squared = closed_form_blocks(**dict(inputs, bare_alp_mass_squared=bound - sp.Rational(1, 1000)))
        diagnostics = exact_stability_diagnostics(kinetic, mass_squared)
        self.assertEqual(diagnostics.violated_conditions, [MASS_CONDITION, POLE_CONDITION])
        with self.assertRaises(UnstableVacuumError):
            solve_physical_spectrum_exact(kinetic, mass_squared)
        kinetic, mass_squared = closed_form_blocks(**dict(inputs, bare_alp_mass_squared=bound))
        self.assertTrue(exact_stability_diagnostics(kinetic, mass_squared).stable)

    def test_nonpositive_kinetic_matrix(self):
        kinetic = np.diag([1.0, 1.0, -0.5])
        with self.assertRaises(UnstableVacuumError) as context:
            solve_physical_spectrum(kinetic, np.diag([1.0, 2.0, 3.0]))
        self.assertIn(KINETIC_CONDITION, context.exception.violated_conditions)
        singular = np.diag([1.0, 0.0])
        with self.assertRaises(UnstableVacuumError) as context:
            solve_physical_spectrum(singular, np.eye(2))
        self.assertIn(KINETIC_CONDITION, context.exception.violated_conditions)
        exact = exact_stability_diagnostics(sp.Matrix([[1, 2], [2, 1]]), sp.eye(2))
        self.assertIn(KINETIC_CONDITION, exact.violated_conditions)
        self.assertEqual(exact.kinetic_determinant, -3)

    def test_random_inputs_consistency(self):
        # For K > 0 the number of negative poles equals that of negative eigenvalues of M^2.
        for seed in range(10):
            kinetic, mass_squared = random_symmetric_quadratic_form(5, seed)
            diagnostics = stability_diagnostics(kinetic, mass_squared)
            self.assertTrue(diagnostics.kinetic_positive_definite)
            self.assertEqual(diagnostics.mass_squared_positive_semidefinite, diagnostics.poles_real_and_nonnegative)
            negative_poles = sum(1 for pole in diagnostics.poles if pole < 0)
            self.assertEqual(negative_poles, int(np.sum(np.linalg.eigvalsh(mass_squared) < 0)))

    def test_massless_state_is_stable(self):
        diagnostics = stability_diagnostics(np.eye(2), np.diag([0.0, 1.0]))
        self.assertTrue(diagnostics.stable)


class DegeneracyTest(unittest.TestCase):
    def test_two_level_formula_away_from_degeneracy(self):
        base = isolated_two_level_inputs()
        for offset in (sp.Rational(-1, 100), sp.Rational(-1, 1000), sp.Rational(1, 500), sp.Rational(3, 100)):
            inputs = dict(base, bare_alp_mass_squared=base["bare_alp_mass_squared"] + offset)
            kinetic, mass_squared = lo_blocks_float(inputs)
            spectrum = solve_spectrum_numerical(kinetic, mass_squared)
            estimates = two_level_estimates(kinetic, mass_squared)
            pion = min(estimates, key=lambda entry: abs(entry["meson_mass_squared"] - 0.02))
            for pole in pion["poles"]:
                self.assertLess(np.min(np.abs(np.array(spectrum.pole_squared_masses) - pole)), 1e-12)
            report = find_near_degeneracies(spectrum, relative_tolerance=0.05)
            pair = [group for group in report.groups if group.alp_carrying_state_count == 2]
            self.assertEqual(len(pair), 1)
            self.assertAlmostEqual(pair[0].mixing_angle, pion["mixing_angle"], places=10)
            self.assertAlmostEqual(pair[0].spread, pion["splitting"], places=12)

    def test_canonical_mixing_closed_form(self):
        inputs = isolated_two_level_inputs()
        expectation = isolated_two_level_expectation(inputs)
        estimates = two_level_estimates(*lo_blocks_float(inputs))
        coupled = [entry for entry in estimates if entry["canonical_mixing"] != 0.0]
        self.assertEqual(len(coupled), 1)
        self.assertAlmostEqual(abs(coupled[0]["canonical_mixing"]), float(abs(expectation["canonical_mixing"])), places=14)
        self.assertAlmostEqual(coupled[0]["canonical_alp_mass_squared"], float(expectation["meson_mass_squared"]), places=14)
        self.assertAlmostEqual(coupled[0]["splitting_at_degeneracy"], float(expectation["splitting"]), places=14)

    def test_alp_near_decoupled_kaons_is_not_alp_meson_degeneracy(self):
        # ALP-like pole placed exactly at m_K0^2 in the isolated system: a threefold degenerate pole.
        base = isolated_two_level_inputs()
        kinetic, mass_squared = lo_blocks_float(base)
        kaon = float((base["quark_spurion"][1, 1] + base["quark_spurion"][2, 2]) / 2)
        estimate = [entry for entry in two_level_estimates(kinetic, mass_squared) if entry["canonical_mixing"] != 0.0][0]
        mixing, pion = estimate["canonical_mixing"], estimate["meson_mass_squared"]
        canonical_target = kaon - mixing**2 / (kaon - pion)
        shift = estimate["canonical_alp_mass_squared"] - mass_squared[5, 5]
        mass_squared[5, 5] = canonical_target - shift
        spectrum = solve_spectrum_numerical(kinetic, mass_squared, field_labels=NEUTRAL_FIELD_LABELS)
        report = find_near_degeneracies(spectrum, relative_tolerance=1e-3)
        kaon_groups = [group for group in report.groups if abs(group.poles[0] - kaon) < 1e-6]
        self.assertEqual(len(kaon_groups), 1)
        group = kaon_groups[0]
        self.assertEqual(len(group.pole_indices), 3)
        self.assertTrue(group.exactly_degenerate)
        self.assertTrue(group.contains_alp_like_pole)
        self.assertEqual(group.alp_carrying_state_count, 1)
        self.assertFalse(group.alp_meson_near_degeneracy)
        np.testing.assert_allclose(group.basis.T @ kinetic @ group.basis, np.eye(3), atol=1e-12)
        np.testing.assert_allclose(group.projector @ group.projector, group.projector, atol=1e-12)
        self.assertAlmostEqual(np.trace(group.projector), 3.0, places=12)
        # ALP weight equals the ALP component squared of the two-level state at m_K0^2.
        alp_component_squared = (kaon - pion)**2 / ((kaon - pion)**2 + mixing**2)
        self.assertAlmostEqual(group.alp_weight, alp_component_squared, places=9)
        self.assertEqual(report.alp_meson_near_degeneracies, [])

    def test_alp_like_pole_near_exact_kaon_pole_neutral_sources(self):
        # Generic neutral sources: the exact m_K0^2 pole has zero ALP component; its
        # degeneracy with the ALP-like pole must be reported only through the coupled states.
        inputs = illustrative_lo_inputs("neutral")
        kinetic, mass_squared = lo_blocks_float(inputs)
        spectrum = solve_spectrum_numerical(kinetic, mass_squared)
        kaon = float((inputs["quark_spurion"][1, 1] + inputs["quark_spurion"][2, 2]) / 2)
        singular = single_pole_formula_singular(spectrum)
        for entry in singular:
            self.assertEqual(entry["singular"], abs(entry["pole"] - kaon) < 1e-9, entry)
        report = find_near_degeneracies(spectrum, relative_tolerance=1e-2)
        for group in report.groups:
            kaon_members = [position for position, pole in enumerate(group.poles) if abs(pole - kaon) < 1e-9]
            for position in kaon_members:
                self.assertLess(abs(group.alp_components[position]), 1e-10)

    def test_projector_is_basis_independent(self):
        spectrum = solve_spectrum_numerical(*lo_blocks_float(illustrative_lo_inputs("flavor_diagonal")))
        report = find_near_degeneracies(spectrum, relative_tolerance=1e-6)
        doublet = [group for group in report.groups if group.exactly_degenerate]
        self.assertEqual(len(doublet), 1)
        group = doublet[0]
        expected = np.diag([0, 0, 0, 1.0, 1.0, 0])
        np.testing.assert_allclose(group.projector, expected, atol=1e-12)
        angle = 0.7
        rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        rotated = group.basis @ rotation
        np.testing.assert_allclose(rotated @ rotated.T @ spectrum.kinetic, expected, atol=1e-12)
        self.assertEqual(group.alp_carrying_state_count, 0)
        self.assertFalse(group.alp_meson_near_degeneracy)

    def test_exact_degenerate_groups(self):
        kinetic, mass_squared = closed_form_blocks(**isolated_two_level_inputs())
        spectrum = solve_spectrum_exact(kinetic, mass_squared, field_labels=NEUTRAL_FIELD_LABELS)
        groups = exact_degenerate_groups(spectrum)
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0]["pole"], sp.Rational(49, 200))
        self.assertEqual(groups[0]["projector"], sp.diag(0, 0, 0, 1, 1, 0))
        self.assertEqual(groups[0]["alp_weight"], 0)

    def test_tolerance_controls_grouping(self):
        spectrum = solve_spectrum_numerical(*lo_blocks_float(illustrative_lo_inputs("cp_conserving")))
        loose = find_near_degeneracies(spectrum, relative_tolerance=1e-2)
        tight = find_near_degeneracies(spectrum, relative_tolerance=1e-8)
        self.assertGreaterEqual(len(loose.groups), len(tight.groups))
        self.assertEqual(tight.groups, [])
        self.assertAlmostEqual(loose.minimal_gap, np.min(np.diff(spectrum.pole_squared_masses)))


if __name__ == "__main__":
    unittest.main()
