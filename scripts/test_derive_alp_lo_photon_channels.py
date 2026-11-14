"""Tests for scripts/derive_alp_lo_photon_channels.py (LO finite-mixing photon channels).

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_photon_channels.py
"""

import itertools
import random
import unittest
from dataclasses import replace

import numpy as np
import sympy as sp

from derive_alp_lo_photon_channels import (
    CHARGE_MATRIX, NEUTRAL_GENERATORS, NUMBER_OF_COLOURS, WZW_ORIENTATION, alp_component_squared_inverse,
    alp_like_pole, alp_secular_function, canonical_frame_matrix, dressed_first_order_invariants,
    electric_charge, float_channel_comparison, illustrative_theory, inverse_mass_operator,
    inverse_on_tadpole_argument, leg_flavour_matrix, leg_matrix_from_pole_vector, lo_blocks, mass_operator,
    photon_two_pseudoscalar_factor, poles_with_alp_component, random_charge_commuting_hermitian,
    random_redefinition, redefinition_taylor_coefficients, schur_mixing_vector, six_term_flavour_trace,
    solve_float_spectrum, tadpole_argument, transform_theory, two_photon_bracket,
    two_photon_bracket_from_pole_vector, two_photon_factor_from_lagrangian, two_photon_weak_mixing_coefficients,
    unmixed_meson_squared_masses, wzw_amplitude_prefactor, wzw_orientation_from_neutral_pion_matching,
    wzw_vertex_sum)
from derive_wzw_alp_photon_two_pseudoscalar import (
    alp_decay_constant as record_alp_decay_constant, amplitude_coefficient as record_amplitude_coefficient,
    antisymmetric_trace, electric_charge as record_electric_charge, number_colors as record_colour_number,
    orientation as record_orientation, outgoing_generators, pion_decay_constant as record_pion_decay_constant)
from symbolic_utils import levi_civita_contract

FLOAT_TOLERANCE = 1e-10
squared_momentum = sp.Symbol("s")
ratio_symbol = sp.Symbol("epsilon", positive=True)


def generic_matrix(name):
    return sp.Matrix(3, 3, lambda row, column: sp.Symbol(f"{name}_{row}{column}"))


def is_zero(matrix_or_value):
    if isinstance(matrix_or_value, sp.MatrixBase):
        return all(sp.cancel(entry) == 0 for entry in matrix_or_value)
    return sp.cancel(matrix_or_value) == 0


def complex_matrix(matrix):
    return np.array(sp.Matrix(matrix).evalf(), dtype=complex)


class NormalizationAndOrientation(unittest.TestCase):
    def test_orientation_is_fixed_by_neutral_pion_matching(self):
        self.assertEqual(wzw_orientation_from_neutral_pion_matching(), -1)
        self.assertEqual(WZW_ORIENTATION, -1)

    def test_two_photon_factor(self):
        self.assertEqual(sp.simplify(two_photon_factor_from_lagrangian() + electric_charge**2 / (4 * sp.pi**2)), 0)

    def test_wzw_prefactor(self):
        """A/F = -orientation N_c e/(24 pi^2 f_pi^3), derived from the gamma Phi^3 coefficient of the record."""
        expected = (-record_orientation * record_colour_number * record_electric_charge
                    / (24 * sp.pi**2 * record_pion_decay_constant**3))
        self.assertEqual(sp.simplify(wzw_amplitude_prefactor() - expected), 0)

    def test_six_term_trace(self):
        first, second, third, charge = (generic_matrix(name) for name in "ABCQ")
        explicit = six_term_flavour_trace(charge, first, second, third)
        self.assertEqual(sp.expand(explicit - antisymmetric_trace(charge, first, second, third)), 0)
        self.assertEqual(sp.expand(explicit + six_term_flavour_trace(charge, first, third, second)), 0)
        self.assertEqual(sp.expand(six_term_flavour_trace(charge, first, second, second)), 0)


class IdentitiesFromAntisymmetry(unittest.TestCase):
    def test_ward_identities_and_two_alp_legs(self):
        vectors = [sp.Matrix(sp.symbols(f"{name}_0:4")) for name in ("q1", "q2", "e2", "k", "p1", "p2", "pa")]
        first, second, polarization, photon, meson_one, meson_two, alp = vectors
        self.assertEqual(sp.expand(levi_civita_contract(first, first, second, polarization)), 0)
        self.assertEqual(sp.expand(levi_civita_contract(photon, meson_one, meson_two, photon)), 0)
        legs = [generic_matrix("X"), generic_matrix("X"), generic_matrix("Y")]
        self.assertEqual(wzw_vertex_sum(legs, [alp, alp, meson_one], polarization), 0)


class ExactFlavourIdentities(unittest.TestCase):
    """Exact rational arithmetic, symbolic s and epsilon."""

    @classmethod
    def setUpClass(cls):
        cls.theory = illustrative_theory()
        cls.symbolic = replace(cls.theory, decay_constant_ratio=ratio_symbol)
        cls.argument = tadpole_argument(cls.theory)
        cls.inverse = inverse_on_tadpole_argument(squared_momentum, cls.theory)

    def test_tadpole_sign_against_exact_pole_record(self):
        theory = self.theory
        difference = theory.axial_source - theory.pseudoscalar_source
        record = (-anticommutator(theory.condensate_parameter * theory.quark_mass, difference)
                  - theory.anomaly_mass_squared / 3 * (theory.gluon_coefficient + theory.axial_source.trace())
                  * sp.eye(3))
        self.assertTrue(is_zero(self.argument - record))
        invariants = dressed_first_order_invariants(theory)
        self.assertTrue(is_zero(invariants["dressed_pseudoscalar_source"]
                                + anticommutator(theory.quark_mass, difference)))

    def test_inverse_mass_operator(self):
        self.assertTrue(is_zero(mass_operator(self.inverse, squared_momentum, self.theory) - self.argument))

    def test_canonical_frame(self):
        kinetic, mass_squared = lo_blocks(self.symbolic)
        transform = canonical_frame_matrix(self.symbolic)
        self.assertTrue(is_zero(transform.T * kinetic * transform - sp.eye(6)))
        mixed = (transform.T * mass_squared * transform)[:5, 5]
        expected = [-ratio_symbol * (generator * self.argument).trace() for generator in NEUTRAL_GENERATORS]
        self.assertTrue(is_zero(mixed - sp.Matrix(expected)))

    def test_pole_direction_normalization_and_vertices(self):
        mixing = schur_mixing_vector(squared_momentum, self.symbolic)
        expected = [ratio_symbol * (generator * (self.theory.axial_source + self.inverse)).trace()
                    for generator in NEUTRAL_GENERATORS]
        self.assertTrue(is_zero(mixing - sp.Matrix(expected)))
        direction = sp.Matrix(expected + [1])
        kinetic, _ = lo_blocks(self.symbolic)
        normalization = alp_component_squared_inverse(squared_momentum, self.symbolic)
        self.assertTrue(is_zero((direction.T * kinetic * direction)[0] - normalization))
        self.assertTrue(is_zero(sp.diff(alp_secular_function(squared_momentum, self.symbolic), squared_momentum)
                                - normalization))
        self.assertTrue(is_zero(two_photon_bracket_from_pole_vector(direction, self.symbolic)
                                - two_photon_bracket(squared_momentum, self.symbolic)))
        self.assertTrue(is_zero(leg_matrix_from_pole_vector(direction, self.symbolic)
                                - 2 * ratio_symbol * self.inverse))

    def test_secular_function_equals_schur_complement(self):
        kinetic, mass_squared = lo_blocks(self.symbolic)
        propagator = squared_momentum * kinetic - mass_squared
        schur = propagator[5, 5] - (propagator[5, :5] * propagator[:5, :5].LUsolve(propagator[:5, 5]))[0]
        self.assertTrue(is_zero(schur - alp_secular_function(squared_momentum, self.symbolic)))


def anticommutator(first, second):
    return first * second + second * first


class SolverComparison(unittest.TestCase):
    def test_alp_pole_mixing_vector_normalization_and_two_photon(self):
        for ratio in (sp.Rational(1, 2), sp.Rational(1, 10)):
            comparison = float_channel_comparison(illustrative_theory(ratio))
            self.assertLess(comparison["mixing_vector_difference"], FLOAT_TOLERANCE)
            self.assertLess(comparison["normalization_difference"], FLOAT_TOLERANCE)
            vector_value, invariant_value = comparison["two_photon"]
            self.assertLess(abs(vector_value - invariant_value), FLOAT_TOLERANCE)
            for _, _, difference in comparison["legs"].values():
                self.assertLess(difference, FLOAT_TOLERANCE)

    def test_aligned_frame_has_alp_pole_without_meson_component(self):
        theory = illustrative_theory(sp.Rational(1, 2))
        pole = alp_like_pole(theory)
        kinetic, mass_squared = lo_blocks(theory)
        inverse = inverse_on_tadpole_argument(pole, theory)
        transform = sp.eye(6)
        for index, generator in enumerate(NEUTRAL_GENERATORS):
            transform[index, 5] = theory.decay_constant_ratio * (generator * (theory.axial_source + inverse)).trace()
        aligned = transform.T * (pole * kinetic - mass_squared) * transform
        self.assertLess(max(abs(sp.N(entry)) for entry in aligned[:, 5]), 1e-35)
        self.assertLess(abs(sp.N(aligned[5, 5])), 1e-35)
        self.assertLess(abs(sp.N((transform.T * kinetic * transform)[5, 5]
                                 - alp_component_squared_inverse(pole, theory))), 1e-35)


class WeakMixingLimits(unittest.TestCase):
    def test_two_photon_limit_and_epsilon_squared_residual(self):
        theory = illustrative_theory(sp.Rational(1, 1000))
        zeroth, first, _ = two_photon_weak_mixing_coefficients(theory)
        invariants = dressed_first_order_invariants(theory)
        reference = invariants["dressed_photon_coefficient"] + 2 * NUMBER_OF_COLOURS * (
            CHARGE_MATRIX**2 * inverse_on_tadpole_argument(theory.bare_alp_mass_squared, theory)).trace()
        self.assertTrue(is_zero(zeroth - reference))
        deviations = []
        for ratio in (sp.Rational(1, 100), sp.Rational(1, 1000)):
            theory = illustrative_theory(ratio)
            pole = alp_like_pole(theory)
            value = two_photon_bracket(pole, theory) / sp.sqrt(alp_component_squared_inverse(pole, theory))
            residual = sp.re(sp.N((value - zeroth) / ratio**2, 30))
            deviations.append(abs(float(residual - first)))
        # (value - zeroth)/epsilon^2 - first is O(epsilon^2): it drops by 100 from 1/100 to 1/1000.
        self.assertLess(abs(deviations[0] / deviations[1] / 100 - 1), 0.01)

    def test_photon_two_pseudoscalar_limits(self):
        """Charged channel -> WZW table; neutral legs -> unmixed generators; residuals O(epsilon^2)."""
        charged_residuals, neutral_residuals = [], []
        for ratio in (sp.Rational(1, 100), sp.Rational(1, 1000)):
            theory = illustrative_theory(ratio)
            limit = inverse_on_tadpole_argument(theory.bare_alp_mass_squared, theory)
            poles = poles_with_alp_component(theory)
            alp_pole = min(poles, key=lambda item: item[1])[0]
            parent = leg_flavour_matrix(alp_pole, theory)
            charged = six_term_flavour_trace(CHARGE_MATRIX, parent, outgoing_generators["pi+"],
                                             outgoing_generators["pi-"]) / (2 * ratio)
            charged_residuals.append(abs(sp.N(charged + 2 * limit[0, 0])) / float(ratio)**2)
            unmixed = unmixed_meson_squared_masses(theory)
            worst = 0
            for pole, _ in poles:
                if pole == alp_pole:
                    continue
                leg = leg_flavour_matrix(pole, theory)
                nearest = min(unmixed, key=lambda value: abs(sp.N(value) - pole))
                eigen_residual = mass_operator(leg, nearest, theory)
                worst = max(worst, max(abs(sp.N(entry)) for entry in eigen_residual),
                            abs(sp.N((leg * leg).trace() - 2)))
            neutral_residuals.append(worst / float(ratio)**2)
        self.assertLess(abs(charged_residuals[1] / charged_residuals[0] - 1), 0.05)
        self.assertLess(neutral_residuals[1], 2 * neutral_residuals[0])

    def test_no_mixing_table_and_prefactor(self):
        up, down, strange, down_strange, strange_down = sp.symbols("S_u S_d S_s S_ds S_sd")
        solution = sp.Matrix([[up, 0, 0], [0, down, down_strange], [0, strange_down, strange]])
        expected = {("pi+", "pi-"): -2 * up, ("K+", "K-"): -2 * up, ("K0", "K0bar"): 2 * (down - strange),
                    ("K0", "pi0"): sp.sqrt(2) * down_strange, ("K0", "eta8"): -sp.sqrt(6) * down_strange,
                    ("K0bar", "pi0"): -sp.sqrt(2) * strange_down, ("K0bar", "eta8"): sp.sqrt(6) * strange_down}
        for (first, second), value in expected.items():
            trace = six_term_flavour_trace(CHARGE_MATRIX, solution, outgoing_generators[first],
                                           outgoing_generators[second])
            self.assertEqual(sp.simplify(trace - value), 0)
        for name, generator in outgoing_generators.items():
            self.assertEqual(sp.expand(six_term_flavour_trace(CHARGE_MATRIX, solution, generator, generator)), 0)
        theory = replace(illustrative_theory(), decay_constant_ratio=ratio_symbol)
        mine = photon_two_pseudoscalar_factor(2 * ratio_symbol * solution, outgoing_generators["pi+"],
                                              outgoing_generators["pi-"], theory)
        record = record_amplitude_coefficient.subs({
            record_orientation: WZW_ORIENTATION, record_colour_number: NUMBER_OF_COLOURS,
            record_electric_charge: electric_charge, record_pion_decay_constant: theory.pion_decay_constant,
            record_alp_decay_constant: theory.pion_decay_constant / ratio_symbol}) * (-2 * up)
        self.assertEqual(sp.simplify(mine - record), 0)


def amplitude_set(theory):
    """Pole-vector (original-field) amplitudes of the ALP-like pole: gamma gamma bracket and gamma P P traces."""
    spectrum = solve_float_spectrum(theory)
    vectors = spectrum.pole_vectors
    alp_index = int(np.argmax(np.abs(vectors[5, :])))
    legs = [leg_matrix_from_pole_vector(list(vectors[:, index]), theory) for index in range(6)]
    channels = {"gamma gamma": complex(two_photon_bracket_from_pole_vector(list(vectors[:, alp_index]), theory)),
                "pi+ pi-": complex(six_term_flavour_trace(CHARGE_MATRIX, legs[alp_index], outgoing_generators["pi+"],
                                                          outgoing_generators["pi-"]))}
    others = [index for index in range(6) if index != alp_index]
    for first, second in itertools.combinations(others, 2):
        channels[(first, second)] = complex(six_term_flavour_trace(CHARGE_MATRIX, legs[alp_index], legs[first],
                                                                   legs[second]))
    return np.array(spectrum.pole_squared_masses), channels


class ChiralRedefinitionInvariance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.theory = illustrative_theory(sp.Rational(1, 2))
        cls.reference_masses, cls.reference_channels = amplitude_set(cls.theory)

    def assert_invariant(self, transformed):
        before = dressed_first_order_invariants(self.theory)
        after = dressed_first_order_invariants(transformed)
        for key in ("dressed_pseudoscalar_source", "dressed_gluon_coefficient", "dressed_photon_coefficient"):
            self.assertTrue(is_zero(before[key] - after[key]), key)
        self.assertTrue(is_zero(tadpole_argument(self.theory) - tadpole_argument(transformed)))
        rational_s = sp.Rational(7, 11)
        self.assertTrue(is_zero(two_photon_bracket(rational_s, self.theory) - two_photon_bracket(rational_s, transformed)))
        self.assertTrue(is_zero(alp_component_squared_inverse(rational_s, self.theory)
                                - alp_component_squared_inverse(rational_s, transformed)))
        masses, channels = amplitude_set(transformed)
        self.assertLess(np.max(np.abs(masses - self.reference_masses)), FLOAT_TOLERANCE)
        for key, value in self.reference_channels.items():
            self.assertLess(abs(channels[key] - value), 1e-9, key)

    def test_constant_axial_redefinition(self):
        generator = random.Random(5)
        rotation = random_charge_commuting_hermitian(generator)
        transformed, _ = transform_theory(self.theory, redefinition_taylor_coefficients(rotation, sp.zeros(3)),
                                          redefinition_taylor_coefficients(-rotation, sp.zeros(3)))
        self.assertTrue(is_zero(transformed.axial_source - self.theory.axial_source - rotation))
        self.assertTrue(is_zero(transformed.pseudoscalar_source - self.theory.pseudoscalar_source - rotation))
        self.assertTrue(is_zero(transformed.photon_coefficient - self.theory.photon_coefficient
                                + 2 * NUMBER_OF_COLOURS * (CHARGE_MATRIX**2 * rotation).trace()))
        self.assertTrue(is_zero(transformed.gluon_coefficient - self.theory.gluon_coefficient + rotation.trace()))
        self.assert_invariant(transformed)

    def test_random_nonlinear_redefinitions(self):
        for seed in (11, 12, 13):
            transformed, _ = transform_theory(self.theory, *random_redefinition(seed))
            self.assert_invariant(transformed)

    def test_second_order_redefinition_terms_do_not_enter(self):
        generator = random.Random(21)
        right_first, left_first = (random_charge_commuting_hermitian(generator) for _ in range(2))
        outputs = []
        for _ in range(2):
            right_second, left_second = (random_charge_commuting_hermitian(generator) for _ in range(2))
            outputs.append(transform_theory(self.theory, redefinition_taylor_coefficients(right_first, right_second),
                                            redefinition_taylor_coefficients(left_first, left_second)))
        (first_theory, first_second_order), (second_theory, second_second_order) = outputs
        self.assertEqual(first_theory, second_theory)
        self.assertFalse(is_zero(first_second_order["second_mass_derivative"]
                                 - second_second_order["second_mass_derivative"]))


if __name__ == "__main__":
    unittest.main()
