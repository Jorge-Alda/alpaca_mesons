"""Tests for the canonical vector-meson pole contribution.

Invocation:
    .venv/bin/python -m unittest scripts.implementation.test_vector_meson_exchange
"""

import unittest
from dataclasses import replace
from unittest.mock import patch

import numpy as np

from scripts.implementation import (
    build_vector_exchange_amplitude,
    build_vector_matching_amplitude,
)


class VectorMesonExchangeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.zero = np.zeros((3, 3))
        self.identity = np.eye(3)
        self.diagonal = np.diag([1.0, -0.3, 0.2])
        self.generic = np.array(
            [[0.4, 0.2j, 0.1], [-0.2j, -0.3, 0.05j], [0.1, -0.05j, 0.2]]
        )

    def build(self, c_a, c_p=None, channel=("pi+", "pi-", "pi0"), **kwargs):
        return build_vector_exchange_amplitude(
            1.5,
            c_a,
            self.zero if c_p is None else c_p,
            *channel,
            100.0,
            **kwargs,
        )

    def test_identity_spurion_has_no_contributions(self) -> None:
        model = self.build(2.7 * self.identity)
        self.assertEqual(model.contributions, ())
        self.assertEqual(model(0.3, 0.4), 0.0j)

    def test_diagonal_spurion_has_no_diagonal_neutral_vector(self) -> None:
        model = self.build(self.diagonal, channel=("K+", "K-", "pi0"))
        self.assertTrue(
            all(name not in {"rho0", "omega", "phi"} for name in model.exchange_names)
        )

    def test_pi_plus_pi_minus_eta_isospin_symmetric_diagonal_channel_is_zero(self) -> None:
        # This is diagonal but not proportional to identity.  Equal u,d ALP
        # couplings make both charged-pion spectator commutators vanish.
        isospin_symmetric = np.diag([1.0, 1.0, -0.3])
        model = self.build(isospin_symmetric, channel=("pi+", "pi-", "eta"))
        self.assertEqual(model.contributions, ())
        self.assertEqual(model(0.3, 0.4), 0.0j)

    def test_allowed_charged_rho_exchange_and_metadata(self) -> None:
        model = self.build(self.diagonal)
        self.assertEqual(model.exchange_names, ("rho-", "rho+"))
        self.assertEqual(
            [(item.spectator_index, item.pair_indices) for item in model.contributions],
            [(0, (1, 2)), (1, (0, 2))],
        )
        self.assertNotEqual(model(0.3, 0.4), 0.0j)

    def test_cqP_is_validated_but_does_not_enter(self) -> None:
        first = self.build(self.diagonal, self.zero)(0.3, 0.4)
        second = self.build(self.diagonal, np.diag([8.0, -2.0, 1.0]))(0.3, 0.4)
        self.assertEqual(first, second)
        nonhermitian = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        with self.assertRaisesRegex(ValueError, "Hermitian"):
            self.build(self.diagonal, nonhermitian)

    def test_amplitude_scales_with_cqA_and_inverse_fa(self) -> None:
        base = self.build(self.diagonal)(0.3, 0.4)
        doubled_coupling = self.build(2.0 * self.diagonal)(0.3, 0.4)
        doubled_fa = build_vector_exchange_amplitude(
            1.5, self.diagonal, self.zero, "pi+", "pi-", "pi0", 200.0
        )(0.3, 0.4)
        np.testing.assert_allclose(doubled_coupling, 2.0 * base, rtol=2e-14, atol=0.0)
        np.testing.assert_allclose(doubled_fa, 0.5 * base, rtol=2e-14, atol=0.0)

    def test_ordered_mandelstam_permutation_covariance(self) -> None:
        ordered = self.build(self.generic)
        cyclic = self.build(self.generic, channel=("pi0", "pi+", "pi-"))
        s12, s23 = 0.30, 0.40
        s13 = ordered.alp_mass**2 + sum(mass**2 for mass in ordered.masses) - s12 - s23
        np.testing.assert_allclose(ordered(s12, s23), cyclic(s13, s12), rtol=2e-14, atol=1e-15)

    def test_repeated_fields_keep_labeled_spectator_diagrams(self) -> None:
        channel = ("pi0", "pi0", "K0")
        ordered = self.build(self.generic, channel=channel)
        permuted = self.build(self.generic, channel=("pi0", "K0", "pi0"))
        self.assertEqual(len(ordered.contributions), 2)
        self.assertEqual(ordered.exchange_names, ("K*0", "K*0"))
        s12, s23 = 0.35, 0.55
        s13 = ordered.alp_mass**2 + sum(mass**2 for mass in ordered.masses) - s12 - s23
        np.testing.assert_allclose(ordered(s12, s23), permuted(s13, s23), rtol=2e-14, atol=1e-15)

    def test_pair_reversal_commutator_and_current_signs_cancel(self) -> None:
        model = self.build(self.diagonal)
        item = model.contributions[0]
        i, j = item.pair_indices
        k = item.spectator_index
        pair = {(0, 1): 0.30, (1, 2): 0.40}
        pair[(0, 2)] = model.alp_mass**2 + sum(m**2 for m in model.masses) - 0.70
        s, sik, sjk = pair[(i, j)], pair[tuple(sorted((i, k)))], pair[tuple(sorted((j, k)))]
        mi2, mj2, mk2 = model.masses[i] ** 2, model.masses[j] ** 2, model.masses[k] ** 2
        numerator = 0.25 * s * (
            (s + mj2 - mi2) * (sik - mi2 - mk2)
            - (s + mi2 - mj2) * (sjk - mj2 - mk2)
        )
        reversed_numerator = 0.25 * s * (
            (s + mi2 - mj2) * (sjk - mj2 - mk2)
            - (s + mj2 - mi2) * (sik - mi2 - mk2)
        )
        np.testing.assert_allclose(reversed_numerator, -numerator, rtol=0.0, atol=1e-15)
        np.testing.assert_allclose(
            (-item.vpp_trace) * item.alp_trace * reversed_numerator,
            item.flavor_factor * numerator,
            rtol=2e-14,
            atol=1e-15,
        )

    def test_default_width_keeps_amplitude_finite_at_a_pole(self) -> None:
        model = self.build(self.diagonal)
        rho_mass = model.contributions[0].pole.mass
        value = model(rho_mass**2, 0.40)
        self.assertTrue(np.isfinite(value.real))
        self.assertTrue(np.isfinite(value.imag))

    def test_pi0_kl_ks_cancels_before_pole_kinematics(self) -> None:
        model = self.build(
            self.diagonal,
            channel=("pi0", "KL", "KS"),
        )
        self.assertEqual(model.contributions, ())
        self.assertEqual(model.exchange_names, ())
        with patch(
            "scripts.implementation.vector_meson_exchange._ordered_pair_invariants",
            side_effect=AssertionError("kinematics evaluated"),
        ):
            self.assertEqual(model(0.30, 0.40), 0.0j)
        with self.assertRaisesRegex(ValueError, "finite"):
            model(np.nan, 0.40)

        flavor = self.build(
            self.diagonal,
            channel=("pi0", "K0", "K0bar"),
        )
        self.assertTrue(flavor.contributions)
        self.assertNotEqual(flavor(0.30, 0.40), 0.0j)


class VectorMatchingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.zero = np.zeros((3, 3))
        self.identity = np.eye(3)
        self.diagonal = np.diag([1.0, -0.3, 0.2])
        self.generic = np.array(
            [[0.4, 0.2j, 0.1], [-0.2j, -0.3, 0.05j], [0.1, -0.05j, 0.2]]
        )

    def build(self, c_a=None, c_p=None, channel=("pi+", "pi-", "pi0"), **kwargs):
        return build_vector_matching_amplitude(
            1.5,
            self.diagonal if c_a is None else c_a,
            self.zero if c_p is None else c_p,
            *channel,
            100.0,
            **kwargs,
        )

    @staticmethod
    def manual_value(model, s12, s23):
        invariant_sum = model.alp_mass**2 + sum(mass**2 for mass in model.masses)
        pair_invariants = {
            (0, 1): s12,
            (1, 2): s23,
            (0, 2): invariant_sum - s12 - s23,
        }
        result = 0.0j
        for item in model.contributions:
            i, j = item.pair_indices
            k = item.spectator_index
            s_ik = pair_invariants[tuple(sorted((i, k)))]
            s_jk = pair_invariants[tuple(sorted((j, k)))]
            mi2, mj2, mk2 = (
                model.masses[i] ** 2,
                model.masses[j] ** 2,
                model.masses[k] ** 2,
            )
            local_kinematic = 0.25 * (
                (s_ik - mi2 - mk2) * (model.alp_mass**2 + mj2 - s_ik)
                - (s_jk - mj2 - mk2) * (model.alp_mass**2 + mi2 - s_jk)
            )
            coefficient = model.g_v**2 / (
                model.f_pi**3 * model.f_a * item.pole.mass**2
            )
            result += coefficient * item.flavor_factor * local_kinematic
        return result

    def test_matches_manual_sum_of_labeled_contact_terms(self) -> None:
        model = self.build(self.generic)
        expected = self.manual_value(model, 0.30, 0.40)
        np.testing.assert_allclose(model(0.30, 0.40), expected, rtol=2e-14, atol=1e-15)

    def test_scales_linearly_with_cqA_and_inverse_fa(self) -> None:
        base = self.build()(0.30, 0.40)
        doubled_coupling = self.build(2.0 * self.diagonal)(0.30, 0.40)
        doubled_fa = build_vector_matching_amplitude(
            1.5, self.diagonal, self.zero, "pi+", "pi-", "pi0", 200.0
        )(0.30, 0.40)
        np.testing.assert_allclose(doubled_coupling, 2.0 * base, rtol=2e-14, atol=0.0)
        np.testing.assert_allclose(doubled_fa, 0.5 * base, rtol=2e-14, atol=0.0)

    def test_scales_quadratically_with_gv(self) -> None:
        base = self.build(g_v=0.06)(0.30, 0.40)
        doubled = self.build(g_v=0.12)(0.30, 0.40)
        np.testing.assert_allclose(doubled, 4.0 * base, rtol=2e-14, atol=0.0)

    def test_is_independent_of_width_at_fixed_vector_mass(self) -> None:
        masses = {
            item.exchange_name: item.pole.mass for item in self.build().contributions
        }
        zero_widths = {name: (mass, 0.0) for name, mass in masses.items()}
        broad_widths = {name: (mass, 2.5) for name, mass in masses.items()}
        first = self.build(pole_parameters=zero_widths)(0.30, 0.40)
        second = self.build(pole_parameters=broad_widths)(0.30, 0.40)
        self.assertEqual(first, second)

    def test_isolated_contribution_scales_with_inverse_mass_squared(self) -> None:
        model = self.build()
        target = next(item for item in model.contributions if item.exchange_name == "rho-")
        isolated = replace(model, contributions=(target,))
        changed_model = self.build(
            pole_parameters={
                "rho-": (2.0 * target.pole.mass, target.pole.width),
            }
        )
        changed_target = next(
            item for item in changed_model.contributions if item.exchange_name == "rho-"
        )
        changed_isolated = replace(changed_model, contributions=(changed_target,))
        np.testing.assert_allclose(
            changed_isolated(0.30, 0.40),
            0.25 * isolated(0.30, 0.40),
            rtol=2e-14,
            atol=0.0,
        )

    def test_ordered_mandelstam_permutation_covariance(self) -> None:
        ordered = self.build(self.generic)
        cyclic = self.build(self.generic, channel=("pi0", "pi+", "pi-"))
        s12, s23 = 0.30, 0.40
        s13 = ordered.alp_mass**2 + sum(mass**2 for mass in ordered.masses) - s12 - s23
        np.testing.assert_allclose(ordered(s12, s23), cyclic(s13, s12), rtol=2e-14, atol=1e-15)

    def test_repeated_fields_are_bose_symmetric_but_two_variable(self) -> None:
        model = self.build(self.generic, channel=("pi0", "pi0", "K0"))
        self.assertEqual(len(model.contributions), 2)
        s12, s23 = 0.35, 0.55
        invariant_sum = model.alp_mass**2 + sum(mass**2 for mass in model.masses)
        s13 = invariant_sum - s12 - s23
        np.testing.assert_allclose(model(s12, s23), model(s12, s13), rtol=2e-14, atol=1e-15)
        self.assertNotAlmostEqual(model(s12, s23), model(s12, 0.45), places=12)

    def test_three_identical_fields_and_singlet_spurion_vanish(self) -> None:
        identical = self.build(self.generic, channel=("pi0", "pi0", "pi0"))
        singlet = self.build(2.7 * self.identity)
        self.assertEqual(identical.contributions, ())
        self.assertEqual(identical(0.30, 0.40), 0.0j)
        self.assertEqual(singlet.contributions, ())
        self.assertEqual(singlet(0.30, 0.40), 0.0j)

    def test_cqP_is_validated_but_does_not_enter(self) -> None:
        base = self.build(c_p=self.zero)(0.30, 0.40)
        changed = self.build(c_p=np.diag([8.0, -2.0, 1.0]))(0.30, 0.40)
        self.assertEqual(base, changed)
        nonhermitian = np.array(
            [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        )
        with self.assertRaisesRegex(ValueError, "Hermitian"):
            self.build(c_p=nonhermitian)

    def test_dalitz_inputs_use_pole_callable_validation(self) -> None:
        model = self.build()
        with self.assertRaisesRegex(TypeError, "real scalars"):
            model("bad", 0.40)
        with self.assertRaisesRegex(ValueError, "finite"):
            model(np.inf, 0.40)

    def test_pi0_kl_ks_cancels_before_matching_kinematics(self) -> None:
        model = self.build(channel=("pi0", "KL", "KS"))
        self.assertEqual(model.contributions, ())
        self.assertEqual(model.exchange_names, ())
        with patch(
            "scripts.implementation.vector_meson_exchange._ordered_pair_invariants",
            side_effect=AssertionError("kinematics evaluated"),
        ):
            self.assertEqual(model(0.30, 0.40), 0.0j)
        with self.assertRaisesRegex(ValueError, "finite"):
            model(0.30, np.nan)

        flavor = self.build(channel=("pi0", "K0", "K0bar"))
        self.assertTrue(flavor.contributions)
        self.assertNotEqual(flavor(0.30, 0.40), 0.0j)


if __name__ == "__main__":
    unittest.main()
