"""Tests of the finite-mixing LO WZW coefficient of a -> gamma P P.

Inputs common to all tests unless stated (GeV units):
B_0 = 2.5 GeV; f_pi = 0.093 GeV and m_0_squared = 0.67124 GeV^2 (Alda:2026nrd
snapshot); e = sqrt(4 pi alpha_em) with alpha_em of alpaca_chiral vector
inputs; hat_S[0] = m_q = diag(B0mq snapshot)/B_0; hat_S[2] = 0;
hat_P[1] = -{m_q, c^AP}; hat_Theta[1] = c_GA = 0.4; m_a0 = 1.5 GeV;
epsilon_K = EPSILON_K of three_meson_width.py (2.228e-3 exp(i 43.52 deg)).
c^AP presets: DIAGONAL = diag(1, 0.5, -0.3);
FLAVOUR_CHANGING = DIAGONAL plus (d,s) entry 0.2 + 0.3 i.
The no-mixing scale is N_c e/(12 pi^2 f_a f_pi^2).

Run from the repository root::

    .venv/bin/python -m unittest scripts/implementation/test_alp_photon_pp_wzw_finite_mixing.py
"""
import json
import math
import unittest

import numpy as np

from alpaca_chiral.pseudoscalar_diag import ALPcouplings
from alpaca_chiral.three_pseudoscalar import F_PI
from alpaca_chiral.vector_inputs import _read_vector_inputs
from .alp_photon_pp_wzw import photon_two_pseudoscalar_wzw_amplitude
from .alp_photon_pp_wzw_finite_mixing import (
    DAUGHTER_LABELS, EPSILON_K, FLAVOUR_KAON_OUTGOING_GENERATORS, NEUTRAL_GELL_MANN_GENERATORS, SNAPSHOT,
    _all_legs, _kaon_mixing_coefficients, alternating_trace, external_legs, flavour_kaon_leg,
    tadpole_source, unmixed_neutral_generators, wzw_coefficient)

SNAPSHOT_VALUES = json.loads(SNAPSHOT.read_text())
B_0 = 2.5
QUARK_MASSES = np.diag(SNAPSHOT_VALUES['B0mq_gev2'])/B_0
# Isospin-symmetric light masses, used where the unmixed pi0 must be lambda_3.
ISOSPIN_QUARK_MASSES = np.diag([.5*(QUARK_MASSES[0, 0] + QUARK_MASSES[1, 1])]*2 + [QUARK_MASSES[2, 2]])
ELECTRIC_COUPLING = math.sqrt(4*math.pi*_read_vector_inputs()['alpha_em'])
F_PI_SNAPSHOT = SNAPSHOT_VALUES['fpi_gev']
M_0_SQUARED = SNAPSHOT_VALUES['m0_squared_gev2']
CONSTANTS = dict(f_pi=F_PI_SNAPSHOT, B_0=B_0, m_0_squared=M_0_SQUARED, e=ELECTRIC_COUPLING)
HAT_THETA_1 = .4
M_A0 = 1.5

DIAGONAL = np.diag([1., .5, -.3]).astype(complex)
FLAVOUR_CHANGING = np.array([[1., 0., 0.], [0., .5, .2+.3j], [0., .2-.3j, -.3]])
LAMBDA_3, LAMBDA_8, LAMBDA_0, LAMBDA_6, LAMBDA_7 = NEUTRAL_GELL_MANN_GENERATORS
K0_OUTGOING = math.sqrt(2)*np.array([[0, 0, 0], [0, 0, 0], [0, 1, 0]], dtype=complex)
K0BAR_OUTGOING = K0_OUTGOING.T.copy()
NEUTRAL_KAON_LABELS = ('K0', 'K0bar', 'K_L', 'K_S')


def hat_P_1(axial_minus_pseudoscalar, quark_masses=QUARK_MASSES):
    """First-order invariant hat_P_1 = -{m_q, c^AP}."""
    return -(quark_masses @ axial_minus_pseudoscalar + axial_minus_pseudoscalar @ quark_masses)


def source_of(axial_minus_pseudoscalar, quark_masses=QUARK_MASSES):
    """C = B_0 hat_P_1 - m_0^2/3 hat_Theta_1, computed here independently of the module."""
    return B_0*hat_P_1(axial_minus_pseudoscalar, quark_masses) - M_0_SQUARED/3*HAT_THETA_1*np.eye(3)


def coefficient(meson_1, meson_2, f_a, axial_minus_pseudoscalar, quark_masses=QUARK_MASSES,
                m_a0=M_A0, hat_Theta_1=HAT_THETA_1, **overrides):
    """Coefficient with the common inputs of the module docstring."""
    return wzw_coefficient(
        meson_1, meson_2, m_a0, f_a, [quark_masses, np.zeros((3, 3)), np.zeros((3, 3))],
        [np.zeros((3, 3)), hat_P_1(axial_minus_pseudoscalar, quark_masses)],
        [0., hat_Theta_1], **{**CONSTANTS, **overrides})


def all_legs(f_a, axial_minus_pseudoscalar, quark_masses=QUARK_MASSES, epsilon_K=EPSILON_K):
    """Legs of every label, internal kaon pole legs and pole masses, with the common inputs."""
    return _all_legs(M_A0, f_a, [quark_masses], [np.zeros((3, 3)), hat_P_1(axial_minus_pseudoscalar, quark_masses)],
                     [0., HAT_THETA_1], F_PI_SNAPSHOT, B_0, M_0_SQUARED, epsilon_K)


def no_mixing_scale(f_a, f_pi=F_PI_SNAPSHOT):
    """N_c e/(12 pi^2 f_a f_pi^2) in GeV^-3."""
    return 3*ELECTRIC_COUPLING/(12*math.pi**2*f_a*f_pi**2)


def mass_operator_inverse(quark_masses, squared_mass, source, m_0_squared=M_0_SQUARED):
    """M_s^{-1}(source), solved as a 9x9 linear system, independently of the module."""
    quark_mass_term = B_0*quark_masses
    operator = np.zeros((9, 9), dtype=complex)
    for column in range(9):
        unit = np.zeros(9, dtype=complex)
        unit[column] = 1
        unit = unit.reshape(3, 3)
        operator[:, column] = (quark_mass_term @ unit + unit @ quark_mass_term
                               + m_0_squared/3*np.trace(unit)*np.eye(3) - squared_mass*unit).reshape(9)
    return np.linalg.solve(operator, source.reshape(9)).reshape(3, 3)


def unmixed_mixing_solution(axial_minus_pseudoscalar, quark_masses=QUARK_MASSES, m_a0=M_A0):
    """S = M_{m_a0^2}^{-1}(C) at the unmixed ALP mass."""
    return mass_operator_inverse(quark_masses, m_a0**2, source_of(axial_minus_pseudoscalar, quark_masses))


def unmixed_diagonal_generators(quark_masses):
    """pi0, eta, eta' unmixed generators and squared masses from the 3x3 block, signs as in the module."""
    chi = 2*B_0*np.diag(quark_masses)
    basis = (LAMBDA_3, LAMBDA_8, LAMBDA_0)
    block = np.array([[np.real(np.trace((a @ b + b @ a) @ np.diag(chi))/4
                               + M_0_SQUARED*np.trace(a).real*np.trace(b).real/6)
                       for b in basis] for a in basis])
    squared_masses, vectors = np.linalg.eigh(block)
    generators = []
    for column in range(3):
        generator = sum(vectors[row, column]*basis[row] for row in range(3))
        generators.append(generator if np.trace(generator @ basis[column]).real > 0 else -generator)
    return generators, squared_masses


def unmixed_legs(axial_minus_pseudoscalar, quark_masses=QUARK_MASSES, epsilon_K=EPSILON_K):
    """epsilon -> 0 limits of the daughter legs, built here without the module.

    Charged mesons and K0, K0bar: lambda^dagger. pi0, eta, eta': the boxed leg
    2 eps X/sqrt(1 + 2 eps^2 tr X^2) tends to sign(tr(U_j C) (m_a0^2 - mu_j)) U_j,
    since its ALP component is positive and tr(U_j C)/(mu_j - s_j) has the sign
    of tr(U_j C) D_j, D_j -> m_a0^2 - mu_j. K_S, K_L: p* K0 +- q* K0bar.
    """
    source = source_of(axial_minus_pseudoscalar, quark_masses)
    generators, squared_masses = unmixed_diagonal_generators(quark_masses)
    legs = {'K0': K0_OUTGOING, 'K0bar': K0BAR_OUTGOING}
    for label, generator, squared_mass in zip(('pi0', 'eta', 'eta_prime'), generators, squared_masses):
        legs[label] = np.sign(np.trace(generator @ source).real*(M_A0**2 - squared_mass))*generator
    p_coefficient, q_coefficient = _kaon_mixing_coefficients(epsilon_K)
    legs['K_S'] = np.conj(p_coefficient)*K0_OUTGOING + np.conj(q_coefficient)*K0BAR_OUTGOING
    legs['K_L'] = np.conj(p_coefficient)*K0_OUTGOING - np.conj(q_coefficient)*K0BAR_OUTGOING
    return legs


class TableConventionTests(unittest.TestCase):
    """The WZW table F(Q; S, lambda^dagger, lambda^dagger) with the module's alternating_trace."""

    def test_table_entries(self):
        mixing = unmixed_mixing_solution(FLAVOUR_CHANGING)
        pi_plus = np.zeros((3, 3), dtype=complex)
        pi_plus[1, 0] = math.sqrt(2)
        cases = (
            ((mixing, pi_plus, pi_plus.T), -2*mixing[0, 0]),
            ((mixing, K0_OUTGOING, K0BAR_OUTGOING), 2*(mixing[1, 1] - mixing[2, 2])),
            ((mixing, K0_OUTGOING, LAMBDA_3), math.sqrt(2)*mixing[1, 2]),
            ((mixing, K0_OUTGOING, LAMBDA_8), -math.sqrt(6)*mixing[1, 2]),
            ((mixing, K0BAR_OUTGOING, LAMBDA_3), -math.sqrt(2)*mixing[2, 1]),
            ((mixing, K0BAR_OUTGOING, LAMBDA_8), math.sqrt(6)*mixing[2, 1]),
            ((mixing, K0_OUTGOING, LAMBDA_0), 0.),
            ((mixing, LAMBDA_3, LAMBDA_8), 0.),
        )
        for legs, expected in cases:
            self.assertAlmostEqual(alternating_trace(*legs), expected, places=14)
        np.testing.assert_array_equal(FLAVOUR_KAON_OUTGOING_GENERATORS['K0'], K0_OUTGOING)
        np.testing.assert_array_equal(FLAVOUR_KAON_OUTGOING_GENERATORS['K0bar'], K0BAR_OUTGOING)


class LabelTests(unittest.TestCase):
    """Public labels: K0, K0bar, K_L, K_S; the internal pole labels K0_1, K0_2 are rejected."""

    def test_internal_kaon_labels_rejected(self):
        self.assertNotIn('K0_1', DAUGHTER_LABELS)
        self.assertNotIn('K0_2', DAUGHTER_LABELS)
        for label in NEUTRAL_KAON_LABELS:
            self.assertIn(label, DAUGHTER_LABELS)
        for pair in (('K0_1', 'pi0'), ('pi0', 'K0_2'), ('KS', 'pi0')):
            with self.assertRaisesRegex(ValueError, 'unknown daughter'):
                coefficient(*pair, 10., FLAVOUR_CHANGING)

    def test_returned_masses_have_no_internal_kaon_labels(self):
        _, masses = external_legs('K0', 'K_S', M_A0, 10., [QUARK_MASSES], [np.zeros((3, 3)), hat_P_1(DIAGONAL)],
                                  [0., HAT_THETA_1], f_pi=F_PI_SNAPSHOT, B_0=B_0, m_0_squared=M_0_SQUARED)
        self.assertFalse({'K0_1', 'K0_2'} & set(masses))
        self.assertIn('a', masses)


class WeakMixingLimitTests(unittest.TestCase):
    """(a) epsilon = f_pi/f_a -> 0 against the no-mixing WZW table, residual of order epsilon^2.

    S = M_{m_a0^2}^{-1}(C) at the unmixed ALP mass m_a0 = 1.5 GeV. Residuals are
    relative to the table value; f_a = 1e2, 1e3, 1e4 GeV give ratios of
    successive residuals of 100 (epsilon^2 scaling).
    """

    f_a_values = (1e2, 1e3, 1e4)

    def assert_epsilon_squared_limit(self, evaluate, expected):
        residuals = [abs(evaluate(f_a)/no_mixing_scale(f_a) - expected)/abs(expected)
                     for f_a in self.f_a_values]
        self.assertLess(residuals[-1], 1e-6)
        for larger, smaller in zip(residuals, residuals[1:]):
            self.assertAlmostEqual(larger/smaller, 100., delta=2.)

    def test_charged_pairs(self):
        """pi+ pi-, K+ K-: -2 S_u (reversed order: +2 S_u); also against alp_photon_pp_wzw."""
        for axial_minus_pseudoscalar in (DIAGONAL, FLAVOUR_CHANGING):
            mixing = unmixed_mixing_solution(axial_minus_pseudoscalar)
            for pair, sign in ((('pi+', 'pi-'), 1), (('K+', 'K-'), 1), (('pi-', 'pi+'), -1)):
                self.assert_epsilon_squared_limit(
                    lambda f_a: coefficient(*pair, f_a, axial_minus_pseudoscalar), -2*sign*mixing[0, 0])
            # No-mixing script with hat c_q^A = -S (ALPcouplings holds twice the Alda matrix).
            couplings = ALPcouplings(0., np.diag([-2*mixing[0, 0], 0., 0.]), np.zeros((3, 3)),
                                     -2*mixing[1:, 1:], np.zeros((3, 3)))
            reference = photon_two_pseudoscalar_wzw_amplitude(
                M_A0, 1e3, couplings, 'pi+', 'pi-', orientation=-1)*(F_PI/F_PI_SNAPSHOT)**2
            self.assertAlmostEqual(reference/(no_mixing_scale(1e3)*-2*mixing[0, 0]), 1., places=12)

    def test_flavour_kaon_pairs_with_isospin_symmetric_masses(self):
        """K0 pi0 -> s_pi0 sqrt2 S_ds, K0bar pi0 -> -s_pi0 sqrt2 S_sd, K0 K0bar -> 2(S_d - S_s).

        m_u = m_d makes the unmixed pi0 exactly lambda_3. s_pi0 =
        sign(tr(lambda_3 C)(m_a0^2 - m_pi0,0^2)) = -1 for these inputs is the
        sign of the pi0 state fixed by the boxed leg (positive ALP component).
        K0 and K0bar carry no sign choice. K0 eta -> s_eta a_8 (-sqrt6 S_ds),
        with a_8 the lambda_8 part of the unmixed eta.
        """
        mixing = unmixed_mixing_solution(FLAVOUR_CHANGING, ISOSPIN_QUARK_MASSES)
        legs = unmixed_legs(FLAVOUR_CHANGING, ISOSPIN_QUARK_MASSES)
        pion_sign = np.trace(legs['pi0'] @ LAMBDA_3).real/2
        np.testing.assert_allclose(legs['pi0'], pion_sign*LAMBDA_3, atol=1e-15)
        self.assertEqual(pion_sign, -1.)
        eta_octet_part = np.trace(legs['eta'] @ LAMBDA_8).real/2

        def evaluate(pair):
            return lambda f_a: coefficient(*pair, f_a, FLAVOUR_CHANGING, ISOSPIN_QUARK_MASSES)

        self.assert_epsilon_squared_limit(evaluate(('K0', 'pi0')), pion_sign*math.sqrt(2)*mixing[1, 2])
        self.assert_epsilon_squared_limit(evaluate(('K0bar', 'pi0')), -pion_sign*math.sqrt(2)*mixing[2, 1])
        self.assert_epsilon_squared_limit(evaluate(('K0', 'eta')), eta_octet_part*-math.sqrt(6)*mixing[1, 2])
        self.assert_epsilon_squared_limit(evaluate(('K0', 'K0bar')), 2*(mixing[1, 1] - mixing[2, 2]))
        self.assert_epsilon_squared_limit(evaluate(('K0bar', 'K0')), -2*(mixing[1, 1] - mixing[2, 2]))

    def test_neutral_pairs_with_snapshot_masses(self):
        """Snapshot quark masses: every pair with a neutral kaon -> F(Q; S, unmixed leg 1, unmixed leg 2)."""
        mixing = unmixed_mixing_solution(FLAVOUR_CHANGING)
        legs = unmixed_legs(FLAVOUR_CHANGING)
        for pair in (('K0', 'pi0'), ('K0bar', 'eta'), ('K0', 'eta_prime'), ('K0', 'K0bar'),
                     ('K_L', 'K_S'), ('K_S', 'pi0'), ('K_L', 'pi0'), ('K_S', 'eta'), ('K_L', 'eta_prime')):
            expected = alternating_trace(mixing, legs[pair[0]], legs[pair[1]])
            self.assertGreater(abs(expected), 1e-3)
            self.assert_epsilon_squared_limit(lambda f_a: coefficient(*pair, f_a, FLAVOUR_CHANGING), expected)


class FlavourKaonLegTests(unittest.TestCase):
    """K0 and K0bar legs as combinations of the two real kaon pole legs."""

    def test_legs_tend_to_outgoing_generators(self):
        """f_a = 1e8 GeV: K0, K0bar legs equal lambda^dagger to 1e-12 (residual of order epsilon^2 ~ 1e-18)."""
        legs, _, _ = all_legs(1e8, FLAVOUR_CHANGING)
        np.testing.assert_allclose(legs['K0'], K0_OUTGOING, atol=1e-12)
        np.testing.assert_allclose(legs['K0bar'], K0BAR_OUTGOING, atol=1e-12)

    def test_diagonal_sources_give_exact_flavour_legs(self):
        """(C)_ds = 0: both kaon poles decouple, so K0, K0bar legs are lambda^dagger at any f_a."""
        for f_a in (1., 10., 1e4):
            legs, _, _ = all_legs(f_a, DIAGONAL)
            np.testing.assert_allclose(legs['K0'], K0_OUTGOING, atol=1e-15)
            np.testing.assert_allclose(legs['K0bar'], K0BAR_OUTGOING, atol=1e-15)

    def test_independent_of_pole_signs_and_matches_direct_pole_legs(self):
        """f_a = 1, 3, 10 GeV: the combination does not change under L_k -> -L_k, and agrees with
        the combination built from the 9x9 solution 2 eps X/sqrt(1 + 2 eps^2 tr X^2) for the coupled
        pole and from the decoupled generator -tr(lambda_7 C) lambda_6 + tr(lambda_6 C) lambda_7."""
        source = source_of(FLAVOUR_CHANGING)
        decoupled = -np.trace(LAMBDA_7 @ source).real*LAMBDA_6 + np.trace(LAMBDA_6 @ source).real*LAMBDA_7
        decoupled = decoupled/math.sqrt(np.trace(decoupled @ decoupled).real/2)
        for f_a in (1., 3., 10.):
            mixing_ratio = F_PI_SNAPSHOT/f_a
            legs, pole_legs, pole_masses = all_legs(f_a, FLAVOUR_CHANGING)
            solution = mass_operator_inverse(QUARK_MASSES, pole_masses['K0_1']**2, source)
            coupled = 2*mixing_ratio*solution/np.sqrt(1 + 2*mixing_ratio**2*np.trace(solution @ solution))
            for label, outgoing in (('K0', K0_OUTGOING), ('K0bar', K0BAR_OUTGOING)):
                for first_sign in (1, -1):
                    for second_sign in (1, -1):
                        np.testing.assert_allclose(
                            flavour_kaon_leg(outgoing, [first_sign*pole_legs['K0_1'],
                                                        second_sign*pole_legs['K0_2']]), legs[label], atol=1e-14)
                np.testing.assert_allclose(flavour_kaon_leg(outgoing, [coupled, decoupled]), legs[label],
                                           atol=1e-9)

    def test_flavour_states_orthonormal(self):
        """The coefficients of |K0>, |K0bar> on the two pole states form a unitary 2x2 matrix."""
        _, pole_legs, _ = all_legs(1., FLAVOUR_CHANGING)
        rows = []
        for outgoing in (K0_OUTGOING, K0BAR_OUTGOING):
            row = []
            for label in ('K0_1', 'K0_2'):
                doublet = sum(np.trace(generator @ pole_legs[label])/2*generator for generator in (LAMBDA_6, LAMBDA_7))
                row.append(np.trace(outgoing @ doublet)/2/math.sqrt(np.trace(doublet @ doublet).real/2))
            rows.append(row)
        rows = np.array(rows)
        np.testing.assert_allclose(rows @ rows.conj().T, np.eye(2), atol=1e-14)


class NeutralKaonSelectionRuleTests(unittest.TestCase):
    """Odd-photon selection-rule record: K_L K_S = 2 p* q* A(K0 K0bar); equal-label pairs vanish."""

    def test_long_short_pair(self):
        """f_a = 1, 10, 1e4 GeV; epsilon_K = EPSILON_K and 0.01 + 0.02 i; both source presets."""
        for axial_minus_pseudoscalar in (DIAGONAL, FLAVOUR_CHANGING):
            for f_a in (1., 10., 1e4):
                flavour = coefficient('K0', 'K0bar', f_a, axial_minus_pseudoscalar)
                for epsilon_K in (EPSILON_K, .01+.02j):
                    p_coefficient, q_coefficient = _kaon_mixing_coefficients(epsilon_K)
                    value = coefficient('K_L', 'K_S', f_a, axial_minus_pseudoscalar, epsilon_K=epsilon_K)
                    self.assertAlmostEqual(abs(value - 2*np.conj(p_coefficient*q_coefficient)*flavour)
                                           / no_mixing_scale(f_a), 0., places=14)
                    self.assertGreater(abs(value), .1*no_mixing_scale(f_a))

    def test_equal_labels_vanish(self):
        """K0 K0, K0bar K0bar, K_S K_S, K_L K_L are zero (identical legs), both source presets."""
        for axial_minus_pseudoscalar in (DIAGONAL, FLAVOUR_CHANGING):
            for f_a in (1., 10.):
                for label in NEUTRAL_KAON_LABELS:
                    self.assertEqual(abs(coefficient(label, label, f_a, axial_minus_pseudoscalar)), 0.)

    def test_diagonal_sources_conserve_strangeness(self):
        """DIAGONAL c^AP: a single neutral kaon with pi0, eta, eta' gives exactly zero."""
        for f_a in (1., 10.):
            for kaon in NEUTRAL_KAON_LABELS:
                for other in ('pi0', 'eta', 'eta_prime'):
                    self.assertEqual(abs(coefficient(kaon, other, f_a, DIAGONAL)), 0.)

    def test_zero_epsilon_K(self):
        """epsilon_K = 0: K_S = (K0 + K0bar)/sqrt2, K_L = (K0 - K0bar)/sqrt2 for legs and amplitudes."""
        for f_a in (1., 10.):
            legs, _, _ = all_legs(f_a, FLAVOUR_CHANGING, epsilon_K=0.)
            np.testing.assert_allclose(legs['K_S'], (legs['K0'] + legs['K0bar'])/math.sqrt(2), atol=1e-15)
            np.testing.assert_allclose(legs['K_L'], (legs['K0'] - legs['K0bar'])/math.sqrt(2), atol=1e-15)
            for other in ('pi0', 'eta', 'eta_prime', 'K0'):
                kaon_0 = coefficient('K0', other, f_a, FLAVOUR_CHANGING)
                kaon_0_bar = coefficient('K0bar', other, f_a, FLAVOUR_CHANGING)
                for label, sign in (('K_S', 1), ('K_L', -1)):
                    value = coefficient(label, other, f_a, FLAVOUR_CHANGING, epsilon_K=0.)
                    self.assertAlmostEqual(abs(value - (kaon_0 + sign*kaon_0_bar)/math.sqrt(2))
                                           / no_mixing_scale(f_a), 0., places=14)


class BoseSymmetryTests(unittest.TestCase):
    """(b) A(P_1, P_2) = -A(P_2, P_1); identical daughters give zero.

    f_a = 1, 10, 1e6 GeV with c^AP = FLAVOUR_CHANGING. Tolerance 1e-13 times
    the no-mixing scale at that f_a.
    """

    def test_exchange_and_identical_daughters(self):
        for f_a in (1., 10., 1e6):
            scale = no_mixing_scale(f_a)
            for index, first in enumerate(DAUGHTER_LABELS):
                self.assertLess(abs(coefficient(first, first, f_a, FLAVOUR_CHANGING)), 1e-13*scale)
                for second in DAUGHTER_LABELS[index + 1:]:
                    forward = coefficient(first, second, f_a, FLAVOUR_CHANGING)
                    backward = coefficient(second, first, f_a, FLAVOUR_CHANGING)
                    self.assertLess(abs(forward + backward), 1e-13*scale)


class ConstantRedefinitionTests(unittest.TestCase):
    """(c) c_q^{A,P} -> c_q^{A,P} + kappa, c_G -> c_G - tr kappa, tr S_2 -> tr S_2 - 4 tr(kappa P_1) - 4 tr(kappa^2 m_q).

    The invariant inputs hat_P_1 = -{m_q, c^A - c^P}, hat_Theta_1 = c_G + tr c^A
    and tr hat_S_2 = tr S_2 + 4 tr(c^A {c^P, m_q}) - 4 tr((c^A)^2 m_q)
    (sum-rule record, divided by B_0) are rebuilt from the transformed
    couplings and must equal the original ones; then every amplitude agrees.
    Inputs: c^P = diag(0.3, 0, 0.1), c^A = FLAVOUR_CHANGING + c^P, c_G = 0.4 - tr c^A,
    tr S_2 = 0.05 GeV (as hat_S[2] = diag(0.05, 0, 0)), kappa with (d,s) entry 0.1 - 0.4 i,
    f_a = 3 GeV.
    """

    @staticmethod
    def invariant_inputs(axial, pseudoscalar, gluonic, trace_scalar_second):
        quark_masses = QUARK_MASSES
        trace_hat_S_2 = (trace_scalar_second
                         + 4*np.trace(axial @ (pseudoscalar @ quark_masses + quark_masses @ pseudoscalar))
                         - 4*np.trace(axial @ axial @ quark_masses))
        return (hat_P_1(axial - pseudoscalar), gluonic + np.trace(axial), trace_hat_S_2)

    def test_constant_axial_redefinition(self):
        kappa = np.array([[.7, 0., 0.], [0., -.2, .1-.4j], [0., .1+.4j, .9]])
        pseudoscalar = np.diag([.3, 0., .1]).astype(complex)
        axial = FLAVOUR_CHANGING + pseudoscalar
        gluonic = HAT_THETA_1 - np.trace(axial)
        trace_scalar_second = .05
        quark_pseudoscalar = pseudoscalar @ QUARK_MASSES + QUARK_MASSES @ pseudoscalar
        transformed_trace = (trace_scalar_second - 4*np.trace(kappa @ quark_pseudoscalar)
                             - 4*np.trace(kappa @ kappa @ QUARK_MASSES))
        original = self.invariant_inputs(axial, pseudoscalar, gluonic, trace_scalar_second)
        rotated = self.invariant_inputs(axial + kappa, pseudoscalar + kappa,
                                        gluonic - np.trace(kappa), transformed_trace)
        np.testing.assert_allclose(rotated[0], original[0], atol=1e-15)
        self.assertAlmostEqual(rotated[1], original[1], places=14)
        self.assertAlmostEqual(rotated[2], original[2], places=14)

        def evaluate(pair, inputs):
            invariant_P_1, invariant_Theta_1, trace_hat_S_2 = inputs
            return wzw_coefficient(*pair, M_A0, 3., [QUARK_MASSES, np.zeros((3, 3)),
                                                     np.diag([trace_hat_S_2.real, 0, 0])],
                                   [np.zeros((3, 3)), invariant_P_1], [0., invariant_Theta_1.real],
                                   **CONSTANTS)
        for pair in (('pi+', 'pi-'), ('K0', 'pi0'), ('K0bar', 'eta'), ('K0', 'K0bar'), ('K_L', 'K_S'),
                     ('K_S', 'pi0'), ('K_L', 'eta_prime')):
            value = evaluate(pair, original)
            self.assertGreater(abs(value), .01*no_mixing_scale(3.))
            self.assertAlmostEqual(abs(evaluate(pair, rotated) - value)/abs(value), 0., places=12)


class LargeDecayConstantTests(unittest.TestCase):
    """(d) f_a f_pi^2 A -> N_c e/(12 pi^2) F(Q; S, unmixed legs), correction of order epsilon^2.

    Snapshot masses and c^AP = FLAVOUR_CHANGING. f_a from 1e2 to 1e10 GeV.
    The relative difference to the limit falls by 100 per decade until it
    reaches the rounding level (1e-12 bound).
    """

    def test_scaling(self):
        mixing = unmixed_mixing_solution(FLAVOUR_CHANGING)
        legs = unmixed_legs(FLAVOUR_CHANGING)
        limits = {('pi+', 'pi-'): -2*mixing[0, 0], ('K+', 'K-'): -2*mixing[0, 0]}
        for pair in (('K0', 'pi0'), ('K0bar', 'eta'), ('K0', 'K0bar'), ('K_L', 'K_S'), ('K_S', 'pi0')):
            limits[pair] = alternating_trace(mixing, legs[pair[0]], legs[pair[1]])
        f_a_values = [10.**power for power in range(2, 11)]
        for pair, limit in limits.items():
            residuals = [abs(coefficient(*pair, f_a, FLAVOUR_CHANGING)*f_a*F_PI_SNAPSHOT**2
                             / (3*ELECTRIC_COUPLING/(12*math.pi**2)) - limit)/abs(limit)
                         for f_a in f_a_values]
            self.assertAlmostEqual(residuals[0]/residuals[1], 100., delta=2.)
            self.assertLess(max(residuals[4:]), 1e-12)


class PoleLegFormulaTests(unittest.TestCase):
    """Real pole legs equal the record's 2 eps X/sqrt(1 + 2 eps^2 tr X^2), X = M_{m^2}^{-1}(C), with its sign.

    X is solved here as a 9x9 linear system at the physical mass from
    calculate_masses; f_a = 1, 3, 10 GeV, c^AP = FLAVOUR_CHANGING. No sign
    adjustment: the module uses the boxed sign (positive ALP component).
    """

    def test_legs_match_direct_solution(self):
        source = source_of(FLAVOUR_CHANGING)
        for f_a in (1., 3., 10.):
            mixing_ratio = F_PI_SNAPSHOT/f_a
            legs, pole_legs, masses = all_legs(f_a, FLAVOUR_CHANGING)
            for label, leg in (('a', legs['a']), ('pi0', legs['pi0']), ('eta', legs['eta']),
                               ('eta_prime', legs['eta_prime']), ('K0_1', pole_legs['K0_1'])):
                solution = mass_operator_inverse(QUARK_MASSES, masses[label]**2, source)
                direct = 2*mixing_ratio*solution/np.sqrt(1 + 2*mixing_ratio**2*np.trace(solution @ solution))
                np.testing.assert_allclose(leg, direct, atol=1e-9)
                self.assertAlmostEqual(np.trace(leg @ leg).real/2 + 1/(1 + 2*mixing_ratio**2*np.trace(
                    solution @ solution).real), 1., places=9)

    def test_decoupled_kaon_leg(self):
        """K0_2: M_{m_K^2}(U) = 0, tr(U C) = 0, tr U^2 = 2, (d,s) entry i C_ds/|C_ds|; pole exactly at mu_K."""
        source = tadpole_source(hat_P_1(FLAVOUR_CHANGING), HAT_THETA_1, B_0, M_0_SQUARED)
        squared_masses, generators = unmixed_neutral_generators(B_0*QUARK_MASSES, M_0_SQUARED, source)
        _, pole_legs, masses = all_legs(3., FLAVOUR_CHANGING)
        leg = pole_legs['K0_2']
        np.testing.assert_array_equal(leg, generators['K0_2'])
        image = (B_0*QUARK_MASSES @ leg + leg @ (B_0*QUARK_MASSES)
                 + M_0_SQUARED/3*np.trace(leg)*np.eye(3) - squared_masses['K0_2']*leg)
        np.testing.assert_allclose(image, 0, atol=1e-15)
        self.assertAlmostEqual(abs(np.trace(leg @ source)), 0., places=15)
        self.assertAlmostEqual(np.trace(leg @ leg).real, 2., places=15)
        self.assertAlmostEqual(leg[1, 2], 1j*source[1, 2]/abs(source[1, 2]), places=15)
        self.assertAlmostEqual(generators['K0_1'][1, 2], source[1, 2]/abs(source[1, 2]), places=15)
        self.assertAlmostEqual(masses['K0_2']**2/squared_masses['K0_2'], 1., places=14)

    def test_constants_are_used(self):
        """f_pi and m_0_squared reach calculate_masses: the eta' mass changes with m_0_squared."""
        arguments = ('pi0', 'pi+', M_A0, 10., [QUARK_MASSES], [np.zeros((3, 3)), hat_P_1(DIAGONAL)],
                     [0., HAT_THETA_1])
        default = external_legs(*arguments, f_pi=F_PI_SNAPSHOT, B_0=B_0, m_0_squared=M_0_SQUARED)[1]
        changed = external_legs(*arguments, f_pi=.0924, B_0=B_0, m_0_squared=.5)[1]
        self.assertGreater(abs(default['eta_prime'] - changed['eta_prime']), 1e-2)
        self.assertNotEqual(default['a'], changed['a'])


class DiagonalNeutralPairTests(unittest.TestCase):
    """Pairs from {pi0, eta, eta'} with the ALP give zero; pairs with a neutral kaon do not.

    Cause: the ALP leg is a real combination of lambda_3, lambda_8, lambda_0
    and the coupled kaon generator, whose (d,s) entry is C_ds/|C_ds|; the
    daughter legs are real combinations of the same matrices. After the
    s-quark rephasing diag(1, 1, C_ds*/|C_ds|) (commutes with Q and m_q) all
    are real symmetric, tr(Q A B C) = tr(Q C B A), and each permutation
    cancels its reverse. Inputs: c^AP = FLAVOUR_CHANGING, f_a = 1, 10, 1e4 GeV;
    bound 1e-13 times the no-mixing scale; for comparison K0 pi0 is above
    1e-2 times it.
    """

    def test_zero_and_nonzero_pairs(self):
        for f_a in (1., 10., 1e4):
            scale = no_mixing_scale(f_a)
            for pair in (('pi0', 'eta'), ('pi0', 'eta_prime'), ('eta', 'eta_prime')):
                self.assertLess(abs(coefficient(*pair, f_a, FLAVOUR_CHANGING)), 1e-13*scale)
            self.assertGreater(abs(coefficient('K0', 'pi0', f_a, FLAVOUR_CHANGING)), 1e-2*scale)


if __name__ == '__main__':
    unittest.main()
