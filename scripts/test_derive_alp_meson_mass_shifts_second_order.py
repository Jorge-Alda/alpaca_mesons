"""Tests for scripts/derive_alp_meson_mass_shifts_second_order.py (O(epsilon^2) meson shifts).

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_meson_mass_shifts_second_order.py
"""

import unittest

import numpy as np
import sympy as sp

from derive_alp_lo_neutral_blocks import closed_form_blocks, illustrative_lo_inputs
from derive_alp_meson_mass_shifts_second_order import (
    EPSILON_VALUES,
    MESON_COUNT,
    _block_inputs,
    derive,
    fitted_exponent,
    high_precision_trace,
    invariant_mixing_vector,
    mass_eigenstate_generators,
    meson_mass_operator,
    mixing_source_matrix,
    neutral_components,
    second_order_meson_poles,
    rational_test_sources,
    solve_neutral_matrix_equation,
)


class SecondOrderMesonShiftTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_all_checks_pass(self):
        for section, content in self.result.items():
            for label, passed in content["checks"].items():
                self.assertTrue(passed, f"{section}: {label}")

    def test_exponents_reported(self):
        self.assertAlmostEqual(self.result["solver (neutral)"]["second_order_exponent"], 4, delta=0.15)
        self.assertAlmostEqual(self.result["solver (flavour diagonal)"]["zeroth_order_exponent"], 2, delta=0.15)

    def test_predicted_shifts_invariant_under_rational_kappa(self):
        sources = rational_test_sources(sp.Rational(1, 10))
        kappa = sp.Matrix([[sp.Rational(2, 3), 0, 0], [0, sp.Rational(-1, 4), sp.Rational(1, 5) + sp.I / 3],
                           [0, sp.Rational(1, 5) - sp.I / 3, sp.Rational(1, 2)]])
        shifted = dict(sources, axial_source=sources["axial_source"] + kappa,
                       pseudoscalar_source=sources["pseudoscalar_source"] + kappa,
                       gluon_coupling=sources["gluon_coupling"] - sp.trace(kappa))
        original = np.array(second_order_meson_poles(sources, sp.Rational(1, 10)))
        transformed = np.array(second_order_meson_poles(shifted, sp.Rational(1, 10)))
        np.testing.assert_allclose(original, transformed, rtol=1e-12)
        self.assertEqual(invariant_mixing_vector(**_block_inputs(sources)).applyfunc(sp.simplify),
                         invariant_mixing_vector(**_block_inputs(shifted)).applyfunc(sp.simplify))

    def test_no_shift_without_alp(self):
        sources = rational_test_sources(sp.Rational(1, 10))
        for unmixed, predicted in second_order_meson_poles(sources, 0):
            self.assertEqual(unmixed, predicted)

    def test_level_repulsion_sign(self):
        sources = rational_test_sources(sp.Rational(1, 10))
        alp = float(sources["bare_alp_mass_squared"])
        for unmixed, predicted in second_order_meson_poles(sources, sp.Rational(1, 10)):
            self.assertGreaterEqual((predicted - unmixed) * (unmixed - alp), 0.0)

    def test_mixed_entry_from_blocks(self):
        sources = rational_test_sources(sp.Rational(1, 3), illustrative_lo_inputs("cp_conserving"))
        kinetic, mass_squared = closed_form_blocks(**sources)
        squared_momentum = sp.Rational(7, 10)
        entry = squared_momentum * kinetic[:MESON_COUNT, MESON_COUNT] - mass_squared[:MESON_COUNT, MESON_COUNT]
        reduced = neutral_components(mixing_source_matrix(squared_momentum, **_block_inputs(sources)))
        self.assertEqual((entry - sp.Rational(1, 3) * reduced).applyfunc(sp.simplify), sp.zeros(MESON_COUNT, 1))

    def test_fitted_exponent_of_power_law(self):
        values = np.array([float(value)**3 for value in EPSILON_VALUES])
        self.assertAlmostEqual(fitted_exponent(EPSILON_VALUES, values), 3.0, places=10)

    def test_neutral_matrix_equation_solver(self):
        quark_spurion = sp.diag(1, 2, 3)
        right_side = sp.Matrix([[1, 0, 0], [0, 2, 1 + sp.I], [0, 1 - sp.I, 3]])
        solution = solve_neutral_matrix_equation(
            lambda argument: meson_mass_operator(argument, quark_spurion, 0, sp.Rational(1, 2)), right_side)
        self.assertEqual((meson_mass_operator(solution, quark_spurion, 0, sp.Rational(1, 2)) - right_side)
                         .applyfunc(sp.simplify), sp.zeros(3))
        self.assertEqual(solution[0, 0], 2)

    def test_mass_eigenstate_generators_without_anomaly_term(self):
        # m_0^2 = 0 and chi_0 diagonal: the eigenvalues are x_u, x_d, x_s and (x_d + x_s)/2 twice
        sources = dict(rational_test_sources(sp.Rational(1, 10), illustrative_lo_inputs("flavor_diagonal")),
                       anomaly_mass_squared=0)
        diagonal = [sources["quark_spurion"][index, index] for index in range(3)]
        expected = sorted(diagonal + [(diagonal[1] + diagonal[2]) / 2] * 2)
        eigenbasis = mass_eigenstate_generators(sources)
        for (value, generator), reference in zip(sorted(eigenbasis, key=lambda pair: pair[0]), expected):
            self.assertAlmostEqual(float(value), float(reference), places=30)
            self.assertAlmostEqual(float(high_precision_trace(generator * generator).real), 2.0, places=30)


if __name__ == "__main__":
    unittest.main()
