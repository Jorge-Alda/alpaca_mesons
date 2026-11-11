"""Tests for the weak-mixing limit f_pi/f_a -> 0 of the exact LO spectrum (task 2.1)."""

import unittest

import numpy as np
import sympy as sp

from derive_alp_lo_neutral_blocks import illustrative_lo_inputs
from solve_alp_meson_spectrum import solve_spectrum_numerical
from verify_alp_weak_mixing_limit import (
    EXPECTED_EXPONENTS,
    alda_first_order_components,
    generic_mass_basis_model,
    lo_blocks_at_ratio,
    numerical_weak_mixing_checks,
    symbolic_alda_identity_checks,
    symbolic_weak_mixing_checks,
)


class SymbolicSeriesTest(unittest.TestCase):
    def test_two_meson_series(self):
        result = symbolic_weak_mixing_checks(meson_count=2)
        for label, passed in result["checks"].items():
            self.assertTrue(passed, label)

    def test_model_structure(self):
        model = generic_mass_basis_model(2)
        epsilon = sp.Symbol("epsilon", positive=True)
        kinetic = model["kinetic"]
        self.assertEqual(sp.expand(kinetic.det()), 1)
        self.assertEqual(kinetic.subs(epsilon, 0), sp.eye(3))

    def test_alda_identity(self):
        for label, passed in symbolic_alda_identity_checks()["checks"].items():
            self.assertTrue(passed, label)


class NumericalLimitTest(unittest.TestCase):
    def test_residual_exponents(self):
        for source_structure in ("flavor_diagonal", "neutral"):
            result = numerical_weak_mixing_checks(source_structure)
            for label, passed in result["checks"].items():
                self.assertTrue(passed, f"{source_structure}: {label}; exponents {result['exponents']}")
            self.assertEqual(set(result["exponents"]), set(EXPECTED_EXPONENTS))

    def test_first_order_components_match_solver_at_tiny_ratio(self):
        ratio = 1e-4
        inputs = illustrative_lo_inputs("neutral")
        kinetic, mass_squared = lo_blocks_at_ratio(ratio, "neutral")
        spectrum = solve_spectrum_numerical(kinetic, mass_squared)
        bare = float(inputs["bare_alp_mass_squared"])
        index = int(np.argmin(np.abs(np.array(spectrum.pole_squared_masses) - bare)))
        components, _ = alda_first_order_components(inputs)
        np.testing.assert_allclose(spectrum.pole_vector(index)[:5] / ratio, components, rtol=1e-5, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
