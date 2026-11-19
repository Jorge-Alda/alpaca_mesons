"""Unit tests: a -> B_1 bar(B_2) in the only-up model (GeV units).

Run: .venv/bin/python -m unittest scripts.test_derive_murayama_baryon_decays
"""
import math
import unittest

import numpy as np
import sympy as sp
from scipy.linalg import expm

from scripts.derive_murayama_baryon_decays import (
    BARYON_MASSES, D_BARYON, F_BARYON, S_BARYON, TOPOLOGICAL_SINGLET_BARYON,
    channel_widths, decay_width, derive_building_blocks, derive_kinetic_normalization,
    derive_only_up_diagonal_inverse, derive_source_identification, derive_spin_sum, inverse_mass_operator,
    lagrangian_couplings, neutral_flavour_conserving_channels, numerical_spin_sum, octet_generators,
    only_up_alp_data, original_basis_neutral_pole, physical_couplings, section_five_couplings,
    spin_summed_squared_amplitude, threshold_kappa)

RANDOM = np.random.default_rng(2026)


def random_hermitian():
    matrix = RANDOM.normal(size=(3, 3)) + 1j * RANDOM.normal(size=(3, 3))
    return (matrix + matrix.conj().T) / 2


class SymbolicTest(unittest.TestCase):
    def test_invariant_source_is_b0_p01(self):
        invariant, expected = derive_source_identification()
        self.assertEqual(sp.simplify(invariant - expected), sp.zeros(3))

    def test_building_blocks_match_section_five(self):
        blocks = derive_building_blocks()
        self.assertEqual(sp.simplify(blocks["vielbein"] - blocks["expected_vielbein"]), sp.zeros(3))
        self.assertEqual(sp.simplify(blocks["connection"]), sp.zeros(3))
        self.assertEqual(sp.simplify(blocks["topological_gradient"] - blocks["expected_topological_gradient"]), 0)

    def test_kinetic_normalization(self):
        normalization, expected = derive_kinetic_normalization()
        self.assertEqual(sp.simplify(normalization - expected), 0)

    def test_diagonal_inverse_closed_form(self):
        solution, unknowns, closed = derive_only_up_diagonal_inverse()
        for unknown, value in zip(unknowns, closed):
            self.assertEqual(sp.simplify(solution[unknown] - value), 0)

    def test_spin_sum_trace_equals_closed_form(self):
        spin = derive_spin_sum()
        self.assertEqual(sp.simplify(sp.expand(spin["trace"] - spin["closed"])), 0)

    def test_diagonal_width_reduces_to_paper(self):
        spin = derive_spin_sum()
        m1, m2, alp_mass, f_a, vector, axial = spin["symbols"]
        mass = sp.Symbol("m_B", positive=True)
        width = spin["width"].subs({m1: mass, m2: mass, vector: 0})
        paper = alp_mass * mass**2 / (2 * sp.pi * f_a**2) * axial * sp.conjugate(axial) \
            * sp.sqrt(1 - 4 * mass**2 / alp_mass**2)
        ratio = sp.simplify((width / paper).subs(alp_mass, 3 * mass))
        self.assertEqual(ratio, 1)


class TraceAndOperatorTest(unittest.TestCase):
    def test_section_five_traces_against_lagrangian_noncommuting(self):
        axial_source, vector_source = random_hermitian(), random_hermitian()
        gluon_axial = 0.37
        vielbein = -2 * axial_source
        connection = -1j * vector_source  # nabla B|_a = -i (da/f_a)[c^VS, B] means Gamma|_a = -i c^VS
        topological = -2 * (gluon_axial + np.trace(axial_source))
        self.assertGreater(np.linalg.norm(axial_source @ vector_source - vector_source @ axial_source), 0.1)
        names = list(octet_generators())
        for first in names:
            for second in names:
                expected = section_five_couplings(first, second, axial_source, vector_source, gluon_axial)
                direct = lagrangian_couplings(first, second, vielbein, connection, topological)
                np.testing.assert_allclose(direct, expected, atol=1e-13)
                # Hermiticity c_{B2B1} = c_{B1B2}^*
                reverse = section_five_couplings(second, first, axial_source, vector_source, gluon_axial)
                np.testing.assert_allclose(np.conj(reverse), expected, atol=1e-13)

    def test_numerical_chiral_shift_noncommuting(self):
        shift = random_hermitian()
        f_a, step = 0.7, 1e-6

        def rotation(value):
            return expm(1j * value * shift / f_a)
        point = 0.3
        derivative = (rotation(point + step) - rotation(point - step)) / (2 * step)
        u, u_dagger = rotation(point), rotation(point).conj().T
        derivative_dagger = derivative.conj().T
        vielbein = 1j * (u_dagger @ derivative - u @ derivative_dagger)
        connection = (u_dagger @ derivative + u @ derivative_dagger) / 2
        np.testing.assert_allclose(vielbein, -2 * shift / f_a, atol=1e-8)
        np.testing.assert_allclose(connection, 0, atol=1e-8)

    def test_inverse_operator_light_limit(self):
        # Alda:2026nrd Eq. Ma_light_limit, exact at m_a = 0.
        data = only_up_alp_data(0.5)
        masses, B_0, singlet = data["quark_masses"], data["B_0"], data["m_1_squared"]
        source = B_0 * np.diag([-masses[0], 0., 0.])
        result = inverse_mass_operator(source, 0., B_0, masses, singlet)
        axial_minus_pseudoscalar = np.diag([0.5, 0., 0.])
        gluon_pseudoscalar = -0.5
        inverse_masses = np.diag(1 / masses)
        expected = (-axial_minus_pseudoscalar - gluon_pseudoscalar * inverse_masses
                    / (np.trace(inverse_masses) + 6 * B_0 / singlet))
        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_inverse_operator_noncommuting(self):
        masses = np.array([0.0022, 0.0047, 0.093])
        source = random_hermitian()
        result = inverse_mass_operator(source, 0.31, 2.8, masses, 0.73)
        image = (2.8 * (np.diag(masses) @ result + result @ np.diag(masses))
                 + 0.73 / 3 * np.trace(result) * np.eye(3) - 0.31 * result)
        np.testing.assert_allclose(image, source, atol=1e-12)


class BasisTest(unittest.TestCase):
    def test_unmixed_basis_equals_original_pole_projection(self):
        """Original polar fields: u_mu = -dPhi/f_pi, d X_TC = -<dPhi>/f_pi on the canonical pole."""
        for kappa in (0.003, 0.05, 0.3, 0.9, 1.0):
            data = only_up_alp_data(kappa)
            mass_squared, meson_matrix, _ = original_basis_neutral_pole(kappa)
            self.assertAlmostEqual(mass_squared, data["m_a"]**2, delta=1e-9 * data["m_a"]**2)
            vielbein = -meson_matrix / data["f_pi"] * data["f_a"]  # coefficient of d a_phys/f_a
            topological = -np.trace(meson_matrix) / data["f_pi"] * data["f_a"]
            for first, second in neutral_flavour_conserving_channels():
                direct = lagrangian_couplings(first, second, vielbein, np.zeros((3, 3)), topological)
                unmixed = physical_couplings(data, first, second)
                np.testing.assert_allclose(direct, unmixed, rtol=1e-8, atol=1e-14)

    def test_small_mixing_limit(self):
        """N - 1 vanishes like kappa_u^2 and the coupling approaches the Section V light-ALP form."""
        small, smaller = only_up_alp_data(2e-3), only_up_alp_data(1e-3)
        self.assertLess(small["normalization"] - 1, 1e-3)
        ratio = (small["normalization"] - 1) / (smaller["normalization"] - 1)
        self.assertAlmostEqual(ratio, 4., delta=0.01)
        for name in ("p", "n", "Sigma0"):
            vector, axial = physical_couplings(smaller, name, name)
            unnormalized = section_five_couplings(name, name, smaller["axial_source"], np.zeros((3, 3)), 0.)
            self.assertAlmostEqual(axial.real, unnormalized[1].real / math.sqrt(smaller["normalization"]), 14)
            self.assertAlmostEqual(axial.real, unnormalized[1].real, delta=1e-4 * abs(axial))
        # Section V couplings with Alda:2026nrd Eq. Ma_light_limit (exact at m_a = 0, finite m_0):
        # agreement up to O(m_a^2/m_pi^2) and O(f_pi^2/f_a^2) corrections.
        masses, B_0, singlet = smaller["quark_masses"], smaller["B_0"], smaller["m_1_squared"]
        inverse_masses = np.diag(1 / masses)
        light_source = -np.diag([0.5, 0., 0.]) + 0.5 * inverse_masses / (np.trace(inverse_masses) + 6 * B_0 / singlet)
        for name in octet_generators():
            light = section_five_couplings(name, name, light_source, np.zeros((3, 3)), 0.)[1].real
            self.assertAlmostEqual(physical_couplings(smaller, name, name)[1].real, light, delta=1e-3 * abs(light))
        # Alda:2026nrd Eq. effcoupl_fc (large m_0, m_{u,d}/m_s neglected), c_uu^AP = 1/2, c_GP = -1/2.
        # These approximations are of order 6 B_0/(m_1^2 <m_q^-1>) ~ 3% and m_u/m_s ~ 2%; allow 10%.
        delta_u, delta_d = D_BARYON + F_BARYON + S_BARYON, S_BARYON
        m_u, m_d = masses[:2]
        proton = 0.5 * delta_u - 0.5 * (m_u * delta_d + m_d * delta_u) / (m_u + m_d)
        neutron = 0.5 * delta_d - 0.5 * (m_u * delta_u + m_d * delta_d) / (m_u + m_d)
        self.assertAlmostEqual(physical_couplings(smaller, "p", "p")[1].real, proton, delta=0.1 * abs(proton))
        self.assertAlmostEqual(physical_couplings(smaller, "n", "n")[1].real, neutron, delta=0.1 * abs(neutron))


class WidthTest(unittest.TestCase):
    def test_spin_sum_explicit_spinors(self):
        cases = [("p", "p", 0., -0.0013), ("Lambda", "Sigma0", 0.0021 - 0.0004j, -0.0007 + 0.0002j),
                 ("Sigma0", "Lambda", 0.0011, 0.0009)]
        for first, second, vector, axial in cases:
            masses = BARYON_MASSES[first], BARYON_MASSES[second]
            for alp_mass in (sum(masses) + 0.05, 2.9):
                explicit = numerical_spin_sum(vector, axial, *masses, alp_mass, 0.0035)
                closed = spin_summed_squared_amplitude(vector, axial, *masses, alp_mass, 0.0035)
                self.assertAlmostEqual(explicit / closed, 1., places=12)

    def test_width_explicit_spinors_model_channels(self):
        data = only_up_alp_data(1.0)
        for first, second in (("p", "p"), ("Lambda", "Sigma0")):
            vector, axial = physical_couplings(data, first, second)
            masses = BARYON_MASSES[first], BARYON_MASSES[second]
            momentum = math.sqrt((data["m_a"]**2 - sum(masses)**2) * (data["m_a"]**2 - (masses[0] - masses[1])**2)) \
                / (2 * data["m_a"])
            explicit = momentum / (8 * math.pi * data["m_a"]**2) * numerical_spin_sum(
                vector, axial, *masses, data["m_a"], data["f_a"])
            self.assertAlmostEqual(decay_width(vector, axial, *masses, data["m_a"], data["f_a"]) / explicit, 1., places=12)

    def test_threshold_behaviour(self):
        masses = BARYON_MASSES["Lambda"], BARYON_MASSES["Sigma0"]
        threshold = sum(masses)
        self.assertEqual(decay_width(0.1, 0.1, *masses, threshold, 0.003), 0.)
        self.assertEqual(decay_width(0.1, 0.1, *masses, threshold - 1e-3, 0.003), 0.)
        ratios_axial, ratios_vector = [], []
        for excess in (1e-4, 1e-6, 1e-8):
            alp_mass = threshold + excess
            momentum = math.sqrt((alp_mass**2 - threshold**2) * (alp_mass**2 - (masses[0] - masses[1])**2)) / (2 * alp_mass)
            ratios_axial.append(decay_width(0., 1., *masses, alp_mass, 0.003) / momentum)
            ratios_vector.append(decay_width(1., 0., *masses, alp_mass, 0.003) / momentum**3)
        # S wave (axial) ~ |p|, P wave (vector) ~ |p|^3, both with finite nonzero limits.
        for ratios in (ratios_axial, ratios_vector):
            self.assertGreater(ratios[-1], 0)
            self.assertAlmostEqual(ratios[-1] / ratios[-2], 1., places=5)

    def test_model_channels_and_thresholds(self):
        channels = neutral_flavour_conserving_channels()
        self.assertEqual(len(channels), 10)
        self.assertIn(("Lambda", "Sigma0"), channels)
        self.assertIsNotNone(threshold_kappa("p", "p"))
        self.assertIsNone(threshold_kappa("Sigma+", "Sigma+"))
        kappa = threshold_kappa("p", "p")
        self.assertEqual(channel_widths(kappa * (1 - 1e-6))[("p", "p")], 0.)
        self.assertGreater(channel_widths(kappa * (1 + 1e-6))[("p", "p")], 0.)
        # Only-up couplings: vanishing vector part; Sigma0-Lambda coupling independent of tilde S_B.
        data = only_up_alp_data(1.0)
        for pair in channels:
            self.assertEqual(physical_couplings(data, *pair)[0], 0)
        np.testing.assert_allclose(physical_couplings(data, "Lambda", "Sigma0", 0.3),
                                   physical_couplings(data, "Lambda", "Sigma0", TOPOLOGICAL_SINGLET_BARYON))


if __name__ == "__main__":
    unittest.main()
