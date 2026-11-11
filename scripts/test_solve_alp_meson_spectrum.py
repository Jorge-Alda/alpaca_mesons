"""Tests for the general pole solver of D(s) = s K - M^2."""

import unittest

import numpy as np
import sympy as sp

from solve_alp_meson_spectrum import (
    KineticMatrixNotPositiveDefiniteError,
    solve_spectrum_exact,
    solve_spectrum_numerical,
)
from verify_alp_lo_spectrum import (
    exact_and_numerical_agree,
    exact_spectrum_checks,
    numerical_spectrum_checks,
    random_rational_quadratic_form,
    random_symmetric_quadratic_form,
)


class ExactSolverTest(unittest.TestCase):
    def test_random_rational_data_with_known_poles(self):
        for seed in (1, 20261111):
            data = random_rational_quadratic_form(size=4, seed=seed)
            spectrum = solve_spectrum_exact(data["kinetic"], data["mass_squared"])
            for label, passed in exact_spectrum_checks(spectrum, data["pole_squared_masses"]).items():
                self.assertTrue(passed, f"seed {seed}: {label}")
            for label, passed in exact_and_numerical_agree(spectrum).items():
                self.assertTrue(passed, f"seed {seed}: {label}")

    def test_degenerate_pole_gives_kinetic_orthonormal_basis(self):
        conversion = sp.Matrix([[1, sp.Rational(1, 2), 0], [0, 1, sp.Rational(-1, 3)], [sp.Rational(1, 4), 0, 1]])
        kinetic = conversion.T * conversion
        mass_squared = conversion.T * sp.diag(sp.Rational(1, 5), sp.Rational(1, 5), 2) * conversion
        spectrum = solve_spectrum_exact(kinetic, mass_squared)
        self.assertEqual(spectrum.pole_squared_masses, [sp.Rational(1, 5), sp.Rational(1, 5), 2])
        for label, passed in exact_spectrum_checks(spectrum, symbolic_propagator=False).items():
            self.assertTrue(passed, label)

    def test_sign_convention(self):
        data = random_rational_quadratic_form(size=4, seed=3)
        spectrum = solve_spectrum_exact(data["kinetic"], data["mass_squared"])
        for index in range(4):
            vector = spectrum.pole_vector(index)
            reference = vector[3] if vector[3] != 0 else max(vector, key=lambda entry: abs(sp.N(entry)))
            self.assertGreater(sp.N(reference), 0)

    def test_composition_labels(self):
        data = random_rational_quadratic_form(size=3, seed=5)
        spectrum = solve_spectrum_exact(data["kinetic"], data["mass_squared"], field_labels=("pi", "eta", "a"))
        self.assertEqual(list(spectrum.composition(0)), ["pi", "eta", "a"])
        self.assertEqual(spectrum.composition(1)["a"], spectrum.pole_vector(1)[2])

    def test_kinetic_not_positive_definite(self):
        with self.assertRaises(KineticMatrixNotPositiveDefiniteError):
            solve_spectrum_exact(sp.Matrix([[1, 2], [2, 1]]), sp.eye(2))
        with self.assertRaises(KineticMatrixNotPositiveDefiniteError):
            solve_spectrum_exact(sp.diag(1, 0), sp.eye(2))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            solve_spectrum_exact(sp.Matrix([[1, 1], [0, 1]]), sp.eye(2))
        with self.assertRaises(ValueError):
            solve_spectrum_exact(sp.eye(2), sp.eye(3))

    def test_poles_without_closed_form_raise(self):
        mass_squared = sp.Matrix(5, 5, lambda row, column: 2 if row == column else (1 if abs(row - column) == 1 else 0))
        mass_squared[0, 0] = 3
        mass_squared[4, 4] = sp.Rational(1, 2)
        with self.assertRaises(ValueError):
            solve_spectrum_exact(sp.eye(5), mass_squared)


class NumericalSolverTest(unittest.TestCase):
    def test_random_symmetric_data(self):
        for size, seed in ((2, 0), (4, 1), (10, 2), (10, 3), (16, 4)):
            kinetic, mass_squared = random_symmetric_quadratic_form(size, seed)
            spectrum = solve_spectrum_numerical(kinetic, mass_squared)
            self.assertTrue(np.all(np.diff(spectrum.pole_squared_masses) >= 0))
            for label, passed in numerical_spectrum_checks(spectrum).items():
                self.assertTrue(passed, f"size {size}, seed {seed}: {label}")

    def test_kinetic_not_positive_definite(self):
        with self.assertRaises(KineticMatrixNotPositiveDefiniteError) as context:
            solve_spectrum_numerical(np.array([[1.0, 2.0], [2.0, 1.0]]), np.eye(2))
        self.assertIn("smallest eigenvalue -1", str(context.exception))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            solve_spectrum_numerical(np.array([[1.0, 0.5], [0.0, 1.0]]), np.eye(2))
        with self.assertRaises(ValueError):
            solve_spectrum_numerical(np.eye(2), np.eye(3))


if __name__ == "__main__":
    unittest.main()
