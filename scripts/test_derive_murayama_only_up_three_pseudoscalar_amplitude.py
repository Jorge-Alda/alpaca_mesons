"""Unit tests: a_phys -> P_i P_j P_k in the model of DiLuzio:2026poh (GeV units).

Cubic-vertex parity rule, sigma exchange channel by channel, quartic vertex in invariant form
(original fields, original source coefficients) and its on-shell independence of the field basis,
vanishing contact for three distinct neutral daughters.

Run: .venv/bin/python -m unittest scripts.test_derive_murayama_only_up_three_pseudoscalar_amplitude
"""
import itertools
import unittest

import numpy as np
import sympy as sp
from scipy.linalg import expm

from scripts.derive_alp_three_pseudoscalar_radial_exchange import pole_normalization, random_parameters, vertex
from scripts.derive_murayama_only_up_three_pseudoscalar_amplitude import (
    LAMBDA, NUMBER_OF_FIELDS, PAIRINGS, assemble, axial_rotation_parameters, charge_model_parameters,
    column_from_direction, cosine_direction_masses_squared, dalitz_invariants, decay_gram,
    derive_charge_source_parity, derive_cubic_trace_parity, derive_diagonal_forms, derive_kinetic_channel_identity,
    derive_neutral_pole_relation, derive_neutral_secular_product, derive_only_up_radial_vertex,
    derive_pole_leg_form, derive_quartic_lagrangian, derive_quartic_vertex, derive_symmetric_chiral_map,
    derive_three_term_identity, diagonal_charged_pair_quartic_vertex, diagonal_neutral_quartic_vertex,
    direct_pieces, final_states, invariant_legs, maximum_cubic_vertex, multilinear_complex, neutral_coordinates,
    neutral_legs, neutral_sector_model, off_shell_gram, only_up_radial_vertex_closed_form,
    original_basis_direct_quartic, alp_shift_matrix, flavour_direction, inverse_mass_operator_of_p01,
    no_mixing_field_change_component_form, no_mixing_field_change_invariant_form, pole_leg,
    quartic_vertex_from_legs, radial_vertex_from_direction, single_flavour_neutral_quartic_vertex, spectrum,
    unmixed_basis_direct_quartic, unmixed_basis_vertex, with_extra_sources)

QUARK_MASSES = (0.0022, 0.0047, 0.095)
MODEL = dict(B0=2.6, pion_decay_constant=0.0924, decay_constant=0.4, singlet_mass_squared=0.75)
POINT = (0.3, 0.5)  # (s_12, s_13) in GeV^2; algebraic checks, not restricted to the physical region


def only_up():
    parameters = charge_model_parameters((1, 0, 0))
    return parameters, spectrum(parameters)


def invariants_tuple(invariants):
    """(s_ij, s_ik, s_jk) from the dict keyed by daughter pairs."""
    return invariants[frozenset({1, 2})], invariants[frozenset({1, 3})], invariants[frozenset({2, 3})]


class SymbolicTest(unittest.TestCase):
    def test_charge_sources_have_definite_parity(self):
        table, violations = derive_charge_source_parity()
        self.assertEqual(violations, [])
        n = sp.Symbol('n', real=True)
        self.assertEqual(sp.simplify(table[(1, 1)][1] + n ** 2), 0)  # P_11 = -n^2 m (only-up: -m_u)
        self.assertEqual(sp.simplify(table[(0, 2)][0] + n ** 2 / 2), 0)

    def test_cubic_trace_words(self):
        general, parity, kinetic = derive_cubic_trace_parity()
        x = sp.Symbol('x', real=True)
        self.assertEqual(general, {('S3',): 2 * x ** 3, ('P0', 'Phi', 'Phi', 'Phi'): sp.Rational(-1, 3),
                                   ('Phi', 'Phi', 'S1'): -x, ('P2', 'Phi'): 2 * x ** 2})
        self.assertEqual(parity, {})
        self.assertEqual(kinetic, {})

    def test_symmetric_chiral_map_has_no_even_terms(self):
        X, Y = sp.symbols('X Y', commutative=False)
        commutator = lambda first, second: first * second - second * first
        orders = derive_symmetric_chiral_map()
        self.assertEqual(orders[1], 2 * X + Y)
        self.assertEqual(orders[2], 0)
        self.assertEqual(orders[4], 0)
        third = -sp.Rational(1, 6) * (commutator(X, commutator(X, Y)) + commutator(Y, commutator(X, Y)))
        self.assertEqual(sp.expand(orders[3] - third), 0)

    def test_quartic_lagrangian_and_vertex(self):
        self.assertEqual(derive_quartic_lagrangian(), {})
        self.assertEqual(derive_quartic_vertex(), {})
        classes, residual = derive_kinetic_channel_identity()
        self.assertEqual(classes, {})
        self.assertEqual(residual, 0)
        self.assertEqual(derive_pole_leg_form(), {})

    def test_diagonal_forms(self):
        self.assertEqual(derive_diagonal_forms(), [0, 0, 0])

    def test_neutral_sector_relations(self):
        residuals, up_bracket = derive_neutral_pole_relation()
        self.assertEqual(residuals, [0, 0, 0])
        self.assertEqual(up_bracket, 0)
        self.assertEqual(derive_neutral_secular_product(), [0, 0, 0])
        self.assertEqual(derive_three_term_identity(), (0, 0))

    def test_only_up_radial_vertex_closed_form(self):
        self.assertEqual(derive_only_up_radial_vertex()[1], 0)


class CubicVertexTest(unittest.TestCase):
    """V3 = 0 iff M(-a, sigma) = M(a, sigma)^dagger; flavour conservation is not the criterion."""

    CONSERVING = staticmethod(lambda: [charge_model_parameters((1, 0, 0)), charge_model_parameters((0, 1, 0)),
                                       charge_model_parameters((2, 0, 1)),
                                       axial_rotation_parameters(0.5 * LAMBDA[6].real + 0.3 * LAMBDA[4].real)])

    @staticmethod
    def violating():
        base = charge_model_parameters((1, 0, 0))
        return [with_extra_sources(base, {('P', (0, 2)): 0.002 * np.diag([0.0, 1.0, -1.0])}),
                with_extra_sources(base, {('S', (0, 1)): 0.002 * np.diag([1.0, -1.0, 0.0])})]

    def test_parity_conserving_sources_have_no_cubic_vertex(self):
        for parameters in self.CONSERVING():
            largest, radial = maximum_cubic_vertex(parameters)
            self.assertGreater(radial, 1e-3)
            self.assertLess(largest, 1e-14 * radial)

    def test_parity_violating_sources_have_cubic_vertex(self):
        for parameters in self.violating():
            largest, radial = maximum_cubic_vertex(parameters)
            self.assertGreater(largest, 0.1 * radial)
            data = spectrum(parameters)
            parent, neutral = neutral_legs(data)
            legs = [parent] + neutral[:3]
            pieces = direct_pieces(parameters, data, [leg[0] for leg in legs], [leg[1] for leg in legs], *POINT)
            self.assertGreater(abs(pieces['non_radial_exchange']), 1e-3 * abs(pieces['total']))

    def test_basis_without_mixing_keeps_the_rule(self):
        """Cubic vertices of the Lagrangian after the chiral shift, random (a, P') directions."""
        generator = np.random.default_rng(5)
        momenta = generator.normal(size=(3, 3))
        gram = momenta @ np.diag([1.0, -1.0, -1.0]) @ momenta.T
        triples = [[np.concatenate([[0.0], generator.normal(size=NUMBER_OF_FIELDS - 1)]) for _ in range(3)]
                   for _ in range(10)]

        def largest(parameters):
            data = spectrum(parameters)
            parent_mass_squared = data['masses_squared'][data['parent']]
            return max(abs(unmixed_basis_vertex(parameters, parent_mass_squared, triple, gram)) for triple in triples)

        violating = [largest(parameters) for parameters in self.violating()]
        self.assertGreater(min(violating), 1e-3)
        for parameters in self.CONSERVING():
            self.assertLess(largest(parameters), 1e-11 * min(violating))  # round-off of O(1) terms


class PoleVectorTest(unittest.TestCase):
    """eq:murayama-physical-fields: Hessian pole vectors are ((2 B0 f_pi/f_a) M^{-1}(P_01), 1)/sqrt(N)."""

    def check(self, parameters):
        data = spectrum(parameters)
        for n in data['mixed']:
            mass_squared = data['masses_squared'][n]
            meson_matrix, alp_component = pole_leg(parameters, mass_squared)
            expected = flavour_direction(meson_matrix, alp_component)
            self.assertLess(np.max(np.abs(data['vectors'][:, n] - expected)), 1e-10)

    def test_only_up_and_random_sources(self):
        self.check(only_up()[0])
        self.check(random_parameters(2026))


class QuarticVertexBasisTest(unittest.TestCase):
    """Invariant V4 (original fields) against the direct vertex in the original fields and in the
    basis without ALP-meson mixing: equal on shell, different off shell."""

    def compare(self, parameters, legs):
        directions = [leg[0] for leg in legs]
        masses_squared = [leg[1] for leg in legs]
        on_shell = decay_gram(masses_squared[0], masses_squared[1:],
                              dalitz_invariants(masses_squared[0], masses_squared[1:], *POINT))
        formula_legs = invariant_legs(parameters, legs)
        results = {}
        for kind, gram in (('on', on_shell), ('off', off_shell_gram())):
            original = original_basis_direct_quartic(parameters, directions, gram)
            unmixed = unmixed_basis_direct_quartic(parameters, masses_squared[0], directions, gram)
            formula = quartic_vertex_from_legs(parameters, formula_legs, gram)
            self.assertGreater(abs(original), 1e-3)
            # precision of M^{-1}_{m^2} near meson-like poles
            self.assertLess(abs(formula - original), 1e-10 * abs(original))
            results[kind] = abs(unmixed - original) / abs(original)
        return results

    def test_random_sources(self):
        parameters = random_parameters(2026)
        parent, others = neutral_legs(spectrum(parameters))
        for selection in ((0, 1, 2), (3, 4, 5), (2, 2, 4), (1, 5, 5)):
            with self.subTest(selection=selection):
                results = self.compare(parameters, [parent] + [others[n] for n in selection])
                self.assertLess(results['on'], 1e-12)
                self.assertGreater(results['off'], 1e-2)

    def test_only_up_charged_legs(self):
        parameters, data = only_up()
        parent, _ = neutral_legs(data)
        for name in ('pi0 pi+ pi-', 'eta K+ K-', 'pi0 K0 K0bar', 'pi0 pi0 eta'):
            with self.subTest(final_state=name):
                results = self.compare(parameters, [parent] + final_states(data)[name])
                self.assertLess(results['on'], 1e-12)


class OnlyUpAmplitudeTest(unittest.TestCase):
    """No exchange of mesons or the ALP, sigma channels, assembly on all final states."""

    @classmethod
    def setUpClass(cls):
        cls.parameters, cls.data = only_up()
        cls.pieces = {name: assemble(cls.parameters, cls.data, daughters, *POINT)
                      for name, daughters in final_states(cls.data).items()}

    def test_no_meson_or_alp_exchange(self):
        for name, pieces in self.pieces.items():
            with self.subTest(final_state=name):
                self.assertLess(abs(pieces['non_radial_exchange']), 1e-14 * abs(pieces['total']))

    def test_radial_channels(self):
        up = self.parameters['sources']['S'][(1, 0)][0, 0]
        for name, pieces in self.pieces.items():
            with self.subTest(final_state=name):
                direct, trace_form = pieces['radial_channels'], pieces['radial_channels_trace_form']
                self.assertLess(np.max(np.abs(direct - trace_form)), 1e-12 * abs(pieces['total']))
                if name in ('pi0 eta eta_prime', 'pi0 pi0 eta', 'pi0 pi0 pi0'):
                    self.assertTrue(np.all(np.abs(direct) > 1e-6))
                elif name == 'pi0 K0 K0bar':
                    self.assertLess(np.max(np.abs(direct)), 1e-25)
                else:  # only (neutral | charged pair)
                    self.assertGreater(abs(direct[0]), 1e-6)
                    self.assertLess(np.max(np.abs(direct[1:])), 1e-25)
        states = final_states(self.data)
        for name in ('pi0 pi+ pi-', 'pi0 K+ K-'):
            _, plus, minus = states[name]
            coupling = radial_vertex_from_direction(self.parameters, plus, minus, 0.7)
            self.assertAlmostEqual(coupling.real, self.parameters['B_0'] * up / self.parameters['f_a'], places=14)

    def test_closed_form_radial_vertex(self):
        parent, neutral = neutral_legs(self.data)
        for x in [parent] + neutral:
            for y in [parent] + neutral:
                for momentum_squared in (0.0, 0.77):
                    general = radial_vertex_from_direction(self.parameters, x, y, momentum_squared)
                    closed = only_up_radial_vertex_closed_form(self.parameters, x[1], y[1], momentum_squared)
                    # Precision limited by M_{m^2}^{-1} near meson-like poles (N(m^2) ~ 10^3): ~7e-10.
                    self.assertLess(abs(general - closed), 1e-8 * abs(closed))

    def test_three_distinct_neutral_daughters_have_no_contact(self):
        pieces = self.pieces['pi0 eta eta_prime']
        without_radial = pieces['total'] - pieces['radial_exchange']
        self.assertLess(abs(without_radial), 1e-12 * abs(pieces['radial_exchange']))
        self.assertLess(abs(pieces['invariant_quartic']),
                        1e-12 * abs(self.pieces['pi0 pi0 eta']['invariant_quartic']))

    def test_assembled_amplitude_and_on_shell_basis_independence(self):
        for name, pieces in self.pieces.items():
            with self.subTest(final_state=name):
                scale = abs(pieces['total'])
                self.assertLess(abs(pieces['assembled'] - pieces['total']), 1e-11 * scale)
                self.assertLess(abs(pieces['invariant_quartic'] - pieces['contact']), 1e-11 * scale)
                self.assertLess(abs(pieces['unmixed_basis_contact'] - pieces['contact']), 1e-12 * scale)

    def test_explicit_only_up_quartic_forms(self):
        states = final_states(self.data)
        parent, _ = neutral_legs(self.data)
        for name, pair in (('pi0 pi+ pi-', (0, 1)), ('eta pi+ pi-', (0, 1)), ('pi0 K+ K-', (0, 2)),
                           ('eta K+ K-', (0, 2)), ('pi0 K0 K0bar', (1, 2))):
            invariants = invariants_tuple(self.pieces[name]['invariants'])
            explicit = diagonal_charged_pair_quartic_vertex(self.parameters, (1, 0, 0), pair, parent[1],
                                                            states[name][0][1], invariants)
            with self.subTest(final_state=name):
                self.assertLess(abs(explicit - self.pieces[name]['contact']), 1e-12 * abs(self.pieces[name]['total']))
        for name in ('pi0 pi0 eta', 'pi0 pi0 pi0', 'pi0 eta eta_prime'):
            masses = [parent[1]] + [leg[1] for leg in states[name]]
            with self.subTest(final_state=name):
                scale = abs(self.pieces[name]['total'])
                self.assertLess(abs(diagonal_neutral_quartic_vertex(self.parameters, (1, 0, 0), masses)
                                    - self.pieces[name]['contact']), 1e-10 * scale)
                self.assertLess(abs(single_flavour_neutral_quartic_vertex(self.parameters, (1, 0, 0), masses)
                                    - self.pieces[name]['contact']), 1e-10 * scale)


class ContactVanishingTest(unittest.TestCase):
    """Predictions of the separable cosine sector and the three-term identity."""

    def no_radial_amplitude(self, charges):
        parameters = charge_model_parameters(charges)
        data = spectrum(parameters)
        parent, neutral = neutral_legs(data)
        legs = [parent] + neutral[:3]
        pieces = direct_pieces(parameters, data, [leg[0] for leg in legs], [leg[1] for leg in legs], *POINT)
        model = neutral_sector_model(charges, QUARK_MASSES, MODEL['B0'], MODEL['pion_decay_constant'],
                                     MODEL['decay_constant'], MODEL['singlet_mass_squared'])
        return pieces, model[2]([neutral_coordinates(leg[0]) for leg in legs]), model, legs

    def test_single_charged_quark_gives_zero(self):
        for charges in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), (-1, 0, 0)):
            pieces, model_contact, _, _ = self.no_radial_amplitude(charges)
            with self.subTest(charges=charges):
                self.assertLess(abs(pieces['total'] - pieces['radial_exchange']), 1e-12 * abs(pieces['radial_exchange']))
                self.assertLess(abs(model_contact), 1e-12 * abs(pieces['radial_exchange']))

    def test_two_charged_quarks_give_nonzero(self):
        for charges in ((1, 1, 0), (1, 0, 1)):
            pieces, model_contact, model, legs = self.no_radial_amplitude(charges)
            without_radial = pieces['total'] - pieces['radial_exchange']
            with self.subTest(charges=charges):
                self.assertGreater(abs(without_radial), 1e-2 * abs(pieces['radial_exchange']))
                self.assertLess(abs(model_contact - without_radial), 1e-10 * abs(without_radial))
                np.testing.assert_allclose(np.sort([leg[1] for leg in legs]), np.sort(model[0]), rtol=1e-10)

    def test_bare_mass_and_identical_daughters_give_nonzero(self):
        arguments = (QUARK_MASSES, MODEL['B0'], MODEL['pion_decay_constant'], MODEL['decay_constant'],
                     MODEL['singlet_mass_squared'])
        _, vectors, contact = neutral_sector_model((1, 0, 0), *arguments)
        distinct = abs(contact([vectors[:, n] for n in range(4)]))
        identical = abs(contact([vectors[:, 0], vectors[:, 1], vectors[:, 1], vectors[:, 3]]))
        self.assertLess(distinct, 1e-12 * identical)
        for bare_mass_squared in (1e-3, 0.05):
            _, vectors, contact = neutral_sector_model((1, 0, 0), *arguments, bare_mass_squared)
            self.assertGreater(abs(contact([vectors[:, n] for n in range(4)])), 1e-3 * identical)


# ---------------------------------------------------------------------------
# Equations of latex/murayama_matching.tex (sec:murayama-vertex and sec:murayama-three-body),
# coded literally in the invariant notation of the draft: original fields, original source
# coefficients, M^{-1}_{m^2}(P_01) and N(m^2).  M_s^{-1} is solved here on the nine complex
# entries of a 3 x 3 matrix, independently of the generator-component solver of the scripts.
# The results are compared with direct vertices of the exact Lagrangian.
# ---------------------------------------------------------------------------

def commutator(first, second):
    return first @ second - second @ first


def anticommutator(first, second):
    return first @ second + second @ first


def ordered_product(matrices):
    result = np.eye(3, dtype=complex)
    for matrix in matrices:
        result = result @ matrix
    return result


def draft_inverse_mass_operator(parameters, momentum_squared, matrix):
    """eq:murayama-operator inverted: Y with B0 {S_00, Y} + (m_1^2/3) <Y> 1 - s Y = X (GeV^-1 for X in GeV)."""
    B0, singlet = parameters['B_0'], parameters['m_1_squared']
    S00 = parameters['sources']['S'][(0, 0)]
    operator = np.zeros((9, 9), dtype=complex)
    for column in range(9):
        unit = np.zeros(9, dtype=complex)
        unit[column] = 1.0
        Y = unit.reshape(3, 3)
        operator[:, column] = (B0 * (S00 @ Y + Y @ S00) + singlet / 3 * np.trace(Y) * np.eye(3)
                               - momentum_squared * Y).reshape(9)
    return np.linalg.solve(operator, np.asarray(matrix, dtype=complex).reshape(9)).reshape(3, 3)


def draft_inverse_of_p01(parameters, mass_squared):
    """M^{-1}_{m^2}(P_01)."""
    return draft_inverse_mass_operator(parameters, mass_squared, parameters['sources']['P'][(0, 1)])


def draft_normalization(parameters, mass_squared):
    """eq:murayama-pole-norm: N(m^2) = 1 + (2 B0^2 f_pi^2/f_a^2) <M^{-1}_{m^2}(P_01)^2>."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    inverse = draft_inverse_of_p01(parameters, mass_squared)
    return (1 + 2 * B0 ** 2 * fpi ** 2 / fa ** 2 * np.trace(inverse @ inverse)).real


def draft_physical_fields(parameters, mass_squared):
    """eq:murayama-physical-fields for one pole r: coefficients of P_r in Phi (matrix) and in a,
    Phi = (2 B0 f_pi/f_a) M^{-1}_{m_r^2}(P_01) P_r/sqrt(N(m_r^2)),  a = P_r/sqrt(N(m_r^2))."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    root = np.sqrt(draft_normalization(parameters, mass_squared))
    return 2 * B0 * fpi / fa * draft_inverse_of_p01(parameters, mass_squared) / root, 1 / root


def draft_radial_mass_squared(parameters):
    """eq:murayama-potential-derivatives: m_sigma^2 = 2 lambda f_a^2 + (B0 f_pi^2/f_a^2)(<S_10> - 2<S_20>)."""
    S = parameters['sources']['S']
    zero = np.zeros((3, 3))
    return (2 * parameters['lambda'] * parameters['f_a'] ** 2 + parameters['B_0'] * parameters['f_pi'] ** 2
            / parameters['f_a'] ** 2 * np.trace(S[(1, 0)] - 2 * S.get((2, 0), zero)).real)


def draft_polar_vertex(parameters, mass_i, mass_j):
    """eq:murayama-polar-vertex:
    g_ij = (1/f_a){(2 B0 f_pi^2/f_a^2)[B0^2 <S_10 {M_i, M_j}> - B0 <P_11 (M_i + M_j)> - <S_12>]
                    + m_sigma^2 - m_i^2 - m_j^2} prod_{r=i,j}{1 + (2 B0^2 f_pi^2/f_a^2) <M_r^2>}^{-1/2},
    M_r = M^{-1}_{m_r^2}(P_01) in this docstring only."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    S, P = parameters['sources']['S'], parameters['sources']['P']
    inverse_i, inverse_j = draft_inverse_of_p01(parameters, mass_i), draft_inverse_of_p01(parameters, mass_j)
    brackets = (B0 ** 2 * np.trace(S[(1, 0)] @ anticommutator(inverse_i, inverse_j))
                - B0 * np.trace(P[(1, 1)] @ (inverse_i + inverse_j)) - np.trace(S[(1, 2)]))
    braces = 2 * B0 * fpi ** 2 / fa ** 2 * brackets + draft_radial_mass_squared(parameters) - mass_i - mass_j
    product = np.prod([(1 + 2 * B0 ** 2 * fpi ** 2 / fa ** 2 * np.trace(inverse @ inverse)) ** -0.5
                       for inverse in (inverse_i, inverse_j)])
    return (braces * product / fa).real


def draft_cartesian_vertex(parameters, mass_i, mass_j):
    """eq:murayama-cartesian-vertex."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    S, P = parameters['sources']['S'], parameters['sources']['P']
    zero = np.zeros((3, 3))
    inverse_i, inverse_j = draft_inverse_of_p01(parameters, mass_i), draft_inverse_of_p01(parameters, mass_j)
    braces = (2 * parameters['lambda'] * fa ** 2
              + 2 * B0 ** 3 * fpi ** 2 / fa ** 2 * np.trace(S[(1, 0)] @ anticommutator(inverse_i, inverse_j))
              - 2 * B0 ** 2 * fpi ** 2 / fa ** 2 * np.trace((P[(1, 1)] - P[(0, 1)]) @ (inverse_i + inverse_j))
              + B0 * fpi ** 2 / fa ** 2 * np.trace(S[(1, 0)] - 2 * S.get((2, 0), zero) + 4 * S[(0, 2)] - 2 * S[(1, 2)]))
    return (braces / fa / np.sqrt(draft_normalization(parameters, mass_i) * draft_normalization(parameters, mass_j))).real


def draft_polar_cartesian_difference(parameters, mass_i, mass_j):
    """Right-hand side of eq:murayama-polar-cartesian (before '= 0')."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    S, P = parameters['sources']['S'], parameters['sources']['P']
    zero = np.zeros((3, 3))
    value = (draft_radial_mass_squared(parameters) - 2 * parameters['lambda'] * fa ** 2
             - B0 * fpi ** 2 / fa ** 2 * np.trace(S[(1, 0)] - 2 * S.get((2, 0), zero)).real)
    for mass_squared in (mass_i, mass_j):
        value -= (mass_squared + 2 * B0 * fpi ** 2 / fa ** 2 * (
            np.trace(S[(0, 2)]) + B0 * np.trace(P[(0, 1)] @ draft_inverse_of_p01(parameters, mass_squared)))).real
    return value


def draft_offshell_vertex(parameters, mass_x, mass_y, momentum_squared):
    """eq:murayama-offshell-vertex, two legs with an ALP component."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    S, P = parameters['sources']['S'], parameters['sources']['P']
    inverse_x, inverse_y = draft_inverse_of_p01(parameters, mass_x), draft_inverse_of_p01(parameters, mass_y)
    right = (2 * B0 ** 3 * fpi ** 2 / fa ** 2 * np.trace(S[(1, 0)] @ anticommutator(inverse_x, inverse_y))
             - 2 * B0 ** 2 * fpi ** 2 / fa ** 2 * np.trace(P[(1, 1)] @ (inverse_x + inverse_y))
             - 2 * B0 * fpi ** 2 / fa ** 2 * np.trace(S[(1, 2)]) + momentum_squared - mass_x - mass_y)
    return (right / (fa * np.sqrt(draft_normalization(parameters, mass_x) * draft_normalization(parameters, mass_y)))).real


def draft_offshell_vertex_with_legs(parameters, first, second, momentum_squared, literal_text=False):
    """eq:murayama-offshell-vertex with the charged-leg rule; first, second = (m^2, generator or None).

    Rule of eq:murayama-quartic: for a leg x without ALP component, (2 B0 f_pi/f_a) M^{-1}_{m_x^2}(P_01)
    -> G_x, N(m_x^2) -> 1, and every term in which x enters through its ALP component is dropped
    (the P_11 term of the other leg, the S_12 term and q^2 - m_x^2 - m_y^2).
    literal_text=True instead replaces M^{-1}_{m_x^2}(P_01) itself by G_x, as the sentence after
    eq:murayama-offshell-vertex reads."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    S, P = parameters['sources']['S'], parameters['sources']['P']
    matrices, norms, alp = [], [], []
    for mass_squared, generator in (first, second):
        if generator is None:
            matrices.append(draft_inverse_of_p01(parameters, mass_squared))
            norms.append(draft_normalization(parameters, mass_squared))
            alp.append(1.0)
        else:
            matrices.append(generator if literal_text else fa / (2 * B0 * fpi) * generator)
            norms.append(1.0)
            alp.append(0.0)
    (mass_x, _), (mass_y, _) = first, second
    right = (2 * B0 ** 3 * fpi ** 2 / fa ** 2 * np.trace(S[(1, 0)] @ anticommutator(*matrices))
             - 2 * B0 ** 2 * fpi ** 2 / fa ** 2 * np.trace(P[(1, 1)] @ (alp[1] * matrices[0] + alp[0] * matrices[1]))
             + alp[0] * alp[1] * (-2 * B0 * fpi ** 2 / fa ** 2 * np.trace(S[(1, 2)]) + momentum_squared - mass_x - mass_y))
    return (right / (fa * np.sqrt(norms[0] * norms[1]))).real


def direct_radial_vertex(parameters, coordinates, first, second, momentum_squared):
    """g_xy(q^2) from the exact Lagrangian, L > -(1/2) sigma sum g_xy P_x P_y: minus the cubic vertex
    on (sigma, x, y), legs x, y on shell with outgoing k_x, k_y, sigma momentum -(k_x + k_y), q^2 = (k_x + k_y)^2."""
    (direction_x, mass_x), (direction_y, mass_y) = first, second
    product = (momentum_squared - mass_x - mass_y) / 2
    gram = np.array([[momentum_squared, -(mass_x + product), -(mass_y + product)],
                     [-(mass_x + product), mass_x, product],
                     [-(mass_y + product), product, mass_y]])
    radial = np.eye(NUMBER_OF_FIELDS)[0].astype(complex)
    return -multilinear_complex(
        lambda real: {'g': vertex(parameters, coordinates, real, np.eye(3), gram)},
        [radial, direction_x, direction_y])['g']


def draft_quartic_vertex(parameters, legs, invariants):
    """eq:murayama-quartic, legs = [(m_x^2, generator or None)] for a, i, j, k:
    V4 = (2 B0 f_pi^2/(3 f_a^4)) prod_x N(m_x^2)^{-1/2} { B0^4 sum_{24 orderings} <S_00 M M M M>
         - 2 B0^3 sum_{x_1} sum_{6 orderings of the others} <P_01 M M M> - 6 B0^2 sum_{pairs {x_1,x_2}} <S_02 {M_x3, M_x4}>
         + 18 B0 sum_{x_1} <P_03 M_x1> + 36 <S_04>
         - B0^3 sum_{(r|uv)} s_uv (<[M_u,M_a][M_v,M_r]> + (u <-> v)) },
    M_x = M^{-1}_{m_x^2}(P_01) in this docstring only.  Charged-leg rule as stated in the text: for a
    leg x without ALP component (generator G_x), (2 B0 f_pi/f_a) M_x -> G_x, N(m_x^2) -> 1, and every
    term in which x enters through its ALP component is dropped (x_1 in the second, x_1 or x_2 in the
    third, x_2, x_3, x_4 in the fourth, all legs in the fifth term).
    invariants = {(u, v): s_uv} for daughter positions 1, 2, 3."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    zero = np.zeros((3, 3))
    S, P = parameters['sources']['S'], parameters['sources']['P']
    matrices, norms, has_alp = [], [], []
    for mass_squared, generator in legs:
        if generator is None:
            matrices.append(draft_inverse_of_p01(parameters, mass_squared))
            norms.append(draft_normalization(parameters, mass_squared))
            has_alp.append(True)
        else:
            matrices.append(fa / (2 * B0 * fpi) * generator)
            norms.append(1.0)
            has_alp.append(False)
    trace = np.trace
    braces = 0
    for ordering in itertools.permutations(range(4)):
        braces += B0 ** 4 * trace(S[(0, 0)] @ ordered_product([matrices[n] for n in ordering]))
    for x1 in range(4):
        if has_alp[x1]:
            for ordering in itertools.permutations([y for y in range(4) if y != x1]):
                braces -= 2 * B0 ** 3 * trace(P[(0, 1)] @ ordered_product([matrices[n] for n in ordering]))
    for x1, x2 in itertools.combinations(range(4), 2):
        if has_alp[x1] and has_alp[x2]:
            x3, x4 = [z for z in range(4) if z not in (x1, x2)]
            braces -= 6 * B0 ** 2 * trace(S.get((0, 2), zero) @ anticommutator(matrices[x3], matrices[x4]))
    for x1 in range(4):
        if all(has_alp[y] for y in range(4) if y != x1):
            braces += 18 * B0 * trace(P.get((0, 3), zero) @ matrices[x1])
    if all(has_alp):
        braces += 36 * trace(S.get((0, 4), zero))
    for r, (u, v) in ((1, (2, 3)), (2, (1, 3)), (3, (1, 2))):
        braces -= B0 ** 3 * invariants[(u, v)] * trace(
            commutator(matrices[u], matrices[0]) @ commutator(matrices[v], matrices[r])
            + commutator(matrices[v], matrices[0]) @ commutator(matrices[u], matrices[r]))
    return complex(2 * B0 * fpi ** 2 / (3 * fa ** 4) * braces / np.prod(np.sqrt(norms)))


def draft_neutral_bracket(parameters, charges, mass_squared, flavour):
    """Bracket of eq:murayama-quartic-neutral: 2 B0 (M^{-1}_{m_x^2}(P_01))_qq + n_q."""
    return 2 * parameters['B_0'] * draft_inverse_of_p01(parameters, mass_squared)[flavour, flavour].real + charges[flavour]


def draft_neutral_quartic(parameters, charges, masses_squared):
    """eq:murayama-quartic-neutral, masses_squared = [m_a^2, m_i^2, m_j^2, m_k^2]:
    V4 = (B0 f_pi^2/f_a^4) prod_x N(m_x^2)^{-1/2} sum_q m_q prod_x [2 B0 (M^{-1}_{m_x^2}(P_01))_qq + n_q]."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    quark_masses = np.diag(parameters['sources']['S'][(0, 0)]).real
    total = sum(quark_masses[q] * np.prod([draft_neutral_bracket(parameters, charges, m2, q) for m2 in masses_squared])
                for q in range(3))
    return B0 * fpi ** 2 / fa ** 4 * total / np.prod([np.sqrt(draft_normalization(parameters, m2))
                                                     for m2 in masses_squared])


def draft_cosine_squared_mass(parameters, charges, flavour):
    """2 B0 m_q (1 + n_q^2 f_pi^2/(2 f_a^2))."""
    quark_masses = np.diag(parameters['sources']['S'][(0, 0)]).real
    return 2 * parameters['B_0'] * quark_masses[flavour] * (
        1 + charges[flavour] ** 2 * parameters['f_pi'] ** 2 / (2 * parameters['f_a'] ** 2))


def draft_single_flavour_quartic(parameters, charges, masses_squared):
    """eq:murayama-quartic-single-flavour:
    V4 = (B0 f_pi^2/f_a^4)(2 B0 m_1^2/3)^4 prod_x <M^{-1}_{m_x^2}(P_01)>/sqrt(N(m_x^2))
         sum_q m_q/prod_x [m_x^2 - 2 B0 m_q (1 + n_q^2 f_pi^2/(2 f_a^2))]."""
    B0, fpi, fa, singlet = parameters['B_0'], parameters['f_pi'], parameters['f_a'], parameters['m_1_squared']
    quark_masses = np.diag(parameters['sources']['S'][(0, 0)]).real
    total = sum(quark_masses[q] / np.prod([m2 - draft_cosine_squared_mass(parameters, charges, q)
                                          for m2 in masses_squared]) for q in range(3))
    traces = np.prod([np.trace(draft_inverse_of_p01(parameters, m2)).real / np.sqrt(draft_normalization(parameters, m2))
                      for m2 in masses_squared])
    return B0 * fpi ** 2 / fa ** 4 * (2 * B0 * singlet / 3) ** 4 * traces * total


def draft_charged_quartic(parameters, charges, pair, parent_mass_squared, neutral_mass_squared, invariants):
    """eq:murayama-quartic-charged, daughters (i, j, k) = (neutral, sqrt2|f><f'|, sqrt2|f'><f|),
    pair = (f, f'), invariants = (s_ij, s_ik, s_jk):
    V4 = B0/(f_a^2 sqrt(N(m_a^2) N(m_i^2))) { (2 B0^2/3)(m_f + m_f')[(M_a)_ff + (M_a)_f'f'][(M_i)_ff + (M_i)_f'f']
         + (4 B0^2/3) sum_{q=f,f'} m_q (M_a)_qq (M_i)_qq
         + (2 B0/3) sum_{(q,q')} m_q n_q [2 (M_a + M_i)_qq + (M_a + M_i)_q'q'] + sum_{q=f,f'} m_q n_q^2 }
       + B0^2 (2 s_jk - s_ij - s_ik)/(3 f_a^2 sqrt(N(m_a^2) N(m_i^2))) [(M_a)_ff - (M_a)_f'f'][(M_i)_ff - (M_i)_f'f'],
    M_x = M^{-1}_{m_x^2}(P_01) in this docstring only."""
    B0, fa = parameters['B_0'], parameters['f_a']
    m = np.diag(parameters['sources']['S'][(0, 0)]).real
    parent = np.diag(draft_inverse_of_p01(parameters, parent_mass_squared)).real
    daughter = np.diag(draft_inverse_of_p01(parameters, neutral_mass_squared)).real
    total = parent + daughter
    f, g = pair
    s_ij, s_ik, s_jk = invariants
    root = np.sqrt(draft_normalization(parameters, parent_mass_squared) * draft_normalization(parameters, neutral_mass_squared))
    braces = (2 * B0 ** 2 / 3 * (m[f] + m[g]) * (parent[f] + parent[g]) * (daughter[f] + daughter[g])
              + 4 * B0 ** 2 / 3 * sum(m[q] * parent[q] * daughter[q] for q in (f, g))
              + 2 * B0 / 3 * sum(m[q] * charges[q] * (2 * total[q] + total[p]) for q, p in ((f, g), (g, f)))
              + sum(m[q] * charges[q] ** 2 for q in (f, g)))
    return (B0 / (fa ** 2 * root) * braces
            + B0 ** 2 * (2 * s_jk - s_ij - s_ik) / (3 * fa ** 2 * root)
            * (parent[f] - parent[g]) * (daughter[f] - daughter[g]))


def draft_only_up_coupling(parameters, mass_x, mass_y, momentum_squared):
    """eq:murayama-only-up-offshell:
    g_xy(q^2) = [q^2 - m_x^2 - m_y^2 + m_x^2 m_y^2 f_a^2/(B0 m_u f_pi^2)]/(f_a sqrt(N(m_x^2) N(m_y^2)))."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    up = parameters['sources']['S'][(0, 0)][0, 0].real
    return ((momentum_squared - mass_x - mass_y + mass_x * mass_y * fa ** 2 / (B0 * up * fpi ** 2))
            / (fa * np.sqrt(draft_normalization(parameters, mass_x) * draft_normalization(parameters, mass_y))))


def draft_only_up_radial_mass_squared(parameters):
    """eq:murayama-only-up-vacuum: m_sigma^2 = B0 m_u f_pi^2/f_a^2 + 2 lambda f_a^2."""
    up = parameters['sources']['S'][(0, 0)][0, 0].real
    return parameters['B_0'] * up * parameters['f_pi'] ** 2 / parameters['f_a'] ** 2 + 2 * parameters['lambda'] * parameters['f_a'] ** 2


def draft_legs(daughters):
    """(m^2, generator or None) for draft_quartic_vertex from (direction, m^2) legs."""
    return [(mass_squared, None if abs(direction[1]) > 0 else column_from_direction(direction)[0])
            for direction, mass_squared in daughters]


def channel_invariants(invariants):
    return {(1, 2): invariants[frozenset({1, 2})], (1, 3): invariants[frozenset({1, 3})],
            (2, 3): invariants[frozenset({2, 3})]}


def off_diagonal_pair(generator):
    """(f, f') of sqrt2 |f><f'|."""
    return [(row, column) for row in range(3) for column in range(3) if abs(generator[row, column]) > 0][0]


PARAMETER_SETS = {'random': lambda: random_parameters(2026), 'only up': lambda: charge_model_parameters((1, 0, 0)),
                  'charges (2, 0, 1)': lambda: charge_model_parameters((2, 0, 1))}


class DraftPhysicalFieldsTest(unittest.TestCase):
    """eq:murayama-physical-fields against the eigenvectors of the polar Hessian (unit kinetic matrix,
    angular component positive); the entry-space M_s^{-1} against the generator-component solver."""

    def test_pole_directions(self):
        for name, make in PARAMETER_SETS.items():
            parameters = make()
            data = spectrum(parameters)
            with self.subTest(parameters=name):
                for n in data['mixed']:
                    mass_squared = data['masses_squared'][n]
                    meson_matrix, alp_component = draft_physical_fields(parameters, mass_squared)
                    self.assertLess(np.max(np.abs(meson_matrix - meson_matrix.conj().T)), 1e-10 * np.max(np.abs(meson_matrix)))
                    expected = flavour_direction(meson_matrix, alp_component)
                    self.assertLess(np.max(np.abs(data['vectors'][:, n] - expected)), 1e-9)
                    np.testing.assert_allclose(draft_inverse_of_p01(parameters, mass_squared),
                                               inverse_mass_operator_of_p01(parameters, mass_squared),
                                               rtol=0, atol=1e-9 * np.max(np.abs(draft_inverse_of_p01(parameters, mass_squared))))
                    self.assertAlmostEqual(draft_normalization(parameters, mass_squared),
                                           pole_normalization(parameters, mass_squared),
                                           delta=1e-9 * pole_normalization(parameters, mass_squared))


class DraftRadialCouplingTest(unittest.TestCase):
    """eq:murayama-polar-vertex, -cartesian-vertex, -polar-cartesian and -offshell-vertex against
    the direct sigma P_x P_y vertices of the polar and Cartesian Lagrangians."""

    MOMENTA_SQUARED = (0.0, 0.41, 2.5)

    def test_on_shell_polar_and_cartesian(self):
        for name, make in PARAMETER_SETS.items():
            parameters = make()
            data = spectrum(parameters)
            parent, others = neutral_legs(data)
            legs = [parent] + others
            radial_mass_squared = draft_radial_mass_squared(parameters)
            self.assertAlmostEqual(radial_mass_squared, data['mass_squared'][0, 0], places=12)
            for x in legs:
                for y in legs:
                    with self.subTest(parameters=name, masses=(x[1], y[1])):
                        polar = draft_polar_vertex(parameters, x[1], y[1])
                        cartesian = draft_cartesian_vertex(parameters, x[1], y[1])
                        direct_polar = direct_radial_vertex(parameters, 'polar', x, y, radial_mass_squared)
                        direct_cartesian = direct_radial_vertex(parameters, 'cartesian', x, y, radial_mass_squared)
                        scale = max(abs(polar), abs(direct_polar))
                        # precision of M^{-1}_{m^2} near meson-like poles (N(m^2) up to ~10^3)
                        self.assertLess(abs(polar - direct_polar), 1e-8 * scale)
                        self.assertLess(abs(cartesian - direct_cartesian), 1e-8 * scale)
                        self.assertLess(abs(polar - cartesian), 1e-8 * scale)
                        difference = draft_polar_cartesian_difference(parameters, x[1], y[1])
                        norms = draft_normalization(parameters, x[1]), draft_normalization(parameters, y[1])
                        left = parameters['f_a'] * np.sqrt(norms[0] * norms[1]) * (polar - cartesian)
                        # each pole-equation term is N(m_r^2) times the round-off of m_r^2 (N up to ~10^6)
                        tolerance = 1e-12 * radial_mass_squared + 1e-14 * radial_mass_squared * sum(norms)
                        self.assertLess(abs(difference - left), tolerance)
                        self.assertLess(abs(difference), tolerance)

    def test_off_shell_vertex(self):
        for name, make in PARAMETER_SETS.items():
            parameters = make()
            parent, others = neutral_legs(spectrum(parameters))
            legs = [parent] + others
            for x in legs:
                for y in legs:
                    for momentum_squared in self.MOMENTA_SQUARED:
                        with self.subTest(parameters=name, masses=(x[1], y[1]), momentum_squared=momentum_squared):
                            draft = draft_offshell_vertex(parameters, x[1], y[1], momentum_squared)
                            direct = direct_radial_vertex(parameters, 'polar', x, y, momentum_squared)
                            self.assertLess(abs(draft - direct), 1e-8 * max(abs(draft), abs(direct)))
            radial_mass_squared = draft_radial_mass_squared(parameters)
            for x in legs:  # m_sigma^2 -> q^2 in the polar vertex
                self.assertAlmostEqual(draft_offshell_vertex(parameters, x[1], parent[1], radial_mass_squared),
                                       draft_polar_vertex(parameters, x[1], parent[1]), delta=1e-12)

    def test_off_shell_vertex_charged_legs(self):
        """Charged and K^0-type legs in the PQ-charge models: all pairs of the eight final states."""
        for charges in ((1, 0, 0), (2, 0, 1)):
            parameters = charge_model_parameters(charges)
            data = spectrum(parameters)
            parent, _ = neutral_legs(data)
            legs = {(round(leg[1], 12), tuple(np.round(leg[0], 12))): leg
                    for daughters in final_states(data).values() for leg in [parent] + daughters}
            legs = list(legs.values())
            literal_failures = 0
            for x in legs:
                for y in legs:
                    if abs(x[0][1]) > 0 and abs(y[0][1]) > 0:
                        continue  # covered by test_off_shell_vertex
                    for momentum_squared in (0.2, 0.7):
                        direct = direct_radial_vertex(parameters, 'polar', x, y, momentum_squared)
                        legs_xy = draft_legs([x, y])
                        draft = draft_offshell_vertex_with_legs(parameters, *legs_xy, momentum_squared)
                        literal = draft_offshell_vertex_with_legs(parameters, *legs_xy, momentum_squared, True)
                        with self.subTest(charges=charges, masses=(x[1], y[1])):
                            self.assertLess(abs(draft - direct), 1e-12 * parameters['B_0'] * QUARK_MASSES[2]
                                            / parameters['f_a'])
                        literal_failures += abs(literal - direct) > 1e-3 * abs(direct)
            self.assertGreater(literal_failures, 0)  # literal reading of the sentence is not correct

    def test_field_change_without_mixing(self):
        """eq:murayama-T with w_k^{(a)} of eq:murayama-pole-direction (component form) against the
        same matrix written with M^{-1}_{m_a^2}(P_01): orthogonal, removes ALP-meson mixing."""
        for parameters in (only_up()[0], random_parameters(2026)):
            data = spectrum(parameters)
            parent_mass_squared = data['masses_squared'][data['parent']]
            invariant = no_mixing_field_change_invariant_form(parameters, parent_mass_squared)
            components = no_mixing_field_change_component_form(parameters, parent_mass_squared)
            np.testing.assert_allclose(invariant, components, rtol=0, atol=1e-13)
            np.testing.assert_allclose(invariant.T @ invariant, np.eye(10), rtol=0, atol=1e-12)
            order = list(range(2, NUMBER_OF_FIELDS)) + [1]  # (P_0..P_8, a)
            kinetic = data['kinetic'][np.ix_(order, order)]
            mass = data['mass_squared'][np.ix_(order, order)]
            np.testing.assert_allclose(invariant.T @ kinetic @ invariant, np.eye(10), rtol=0, atol=1e-12)
            transformed = invariant.T @ mass @ invariant
            self.assertLess(np.max(np.abs(transformed[:9, 9])), 1e-10 * abs(transformed[9, 9]) + 1e-14)
            self.assertLess(abs(transformed[9, 9] - parent_mass_squared), 1e-10 * parent_mass_squared)


class DraftQuarticVertexTest(unittest.TestCase):
    """eq:murayama-quartic with the charged-leg rule against the direct vertex of the original polar
    fields and of the basis without ALP-meson mixing (on shell)."""

    def compare(self, parameters, data, daughters):
        parent, _ = neutral_legs(data)
        legs = [parent] + daughters
        invariants = dalitz_invariants(parent[1], [leg[1] for leg in daughters], *POINT)
        gram = decay_gram(parent[1], [leg[1] for leg in daughters], invariants)
        direct = original_basis_direct_quartic(parameters, [leg[0] for leg in legs], gram)
        unmixed = unmixed_basis_direct_quartic(parameters, parent[1], [leg[0] for leg in legs], gram)
        draft = draft_quartic_vertex(parameters, draft_legs(legs), channel_invariants(invariants))
        self.assertGreater(abs(direct), 1e-8)
        self.assertLess(abs(draft - direct), 1e-10 * abs(direct))
        self.assertLess(abs(draft - unmixed), 1e-10 * abs(direct))

    def test_random_sources(self):
        parameters = random_parameters(2026)
        data = spectrum(parameters)
        _, others = neutral_legs(data)
        for selection in ((0, 1, 2), (3, 4, 5), (2, 2, 4), (1, 5, 5)):
            with self.subTest(selection=selection):
                self.compare(parameters, data, [others[n] for n in selection])

    def test_charge_models_all_final_states(self):
        for charges in ((1, 0, 0), (2, 0, 1)):
            parameters = charge_model_parameters(charges)
            data = spectrum(parameters)
            for name, daughters in final_states(data).items():
                if name == 'pi0 eta eta_prime' and charges == (1, 0, 0):
                    continue  # vanishes; covered in DraftOnlyUpAmplitudeTest
                with self.subTest(charges=charges, final_state=name):
                    self.compare(parameters, data, daughters)


class DraftDiagonalQuarticTest(unittest.TestCase):
    """eq:murayama-quartic-neutral and eq:murayama-quartic-charged against the direct amplitude minus
    sigma exchange (no cubic vertices for these sources) and against eq:murayama-quartic."""

    def test_neutral_and_charged_forms(self):
        for charges in ((1, 0, 0), (2, 0, 1), (1, 1, 0), (-1, 0, 0.5)):
            parameters = charge_model_parameters(charges)
            data = spectrum(parameters)
            parent, _ = neutral_legs(data)
            for name, daughters in final_states(data).items():
                pieces = direct_pieces(parameters, data, [leg[0] for leg in [parent] + daughters],
                                       [parent[1]] + [leg[1] for leg in daughters], *POINT)
                direct = pieces['total'] - pieces['radial_exchange']
                scale = abs(pieces['total'])
                masses = [leg[1] for leg in daughters]
                invariants = invariants_tuple(pieces['invariants'])
                general = draft_quartic_vertex(parameters, draft_legs([parent] + daughters),
                                               channel_invariants(pieces['invariants']))
                if name in ('pi0 eta eta_prime', 'pi0 pi0 eta', 'pi0 pi0 pi0'):
                    draft = draft_neutral_quartic(parameters, charges, [parent[1]] + masses)
                else:
                    pair = off_diagonal_pair(column_from_direction(daughters[1][0])[0])
                    draft = draft_charged_quartic(parameters, charges, pair, parent[1], masses[0], invariants)
                    swapped = draft_charged_quartic(parameters, charges, pair[::-1], parent[1], masses[0],
                                                    (invariants[1], invariants[0], invariants[2]))
                    self.assertLess(abs(draft - swapped), 1e-12 * scale)
                    # 2 s_jk - s_ij - s_ik = 3 s_jk - sum_x m_x^2 on shell
                    self.assertAlmostEqual(2 * invariants[2] - invariants[0] - invariants[1],
                                           3 * invariants[2] - parent[1] - sum(masses), places=12)
                with self.subTest(charges=charges, final_state=name):
                    self.assertLess(np.max(np.abs(pieces['non_radial_exchange'])), 1e-13 * scale)
                    self.assertLess(abs(draft - direct), 1e-10 * scale)
                    self.assertLess(abs(draft - general), 1e-10 * scale)

    def test_neutral_bracket_and_cosine_argument(self):
        """Bracket of eq:murayama-quartic-neutral against the cosine argument Phi_qq/f_pi + n_q a/f_a of
        the normalized pole vector: bracket = f_a sqrt(N(m_x^2)) (argument)."""
        for charges in ((1, 0, 0), (2, 0, 1)):
            parameters = charge_model_parameters(charges)
            data = spectrum(parameters)
            parent, others = neutral_legs(data)
            for direction, mass_squared in [parent] + others:
                meson_matrix, alp_component = column_from_direction(direction)
                root = np.sqrt(draft_normalization(parameters, mass_squared))
                for q in range(3):
                    argument = (meson_matrix[q, q].real / parameters['f_pi']
                                + charges[q] * alp_component.real / parameters['f_a'])
                    bracket = draft_neutral_bracket(parameters, charges, mass_squared, q)
                    self.assertLess(abs(bracket - parameters['f_a'] * root * argument), 1e-9 * (abs(bracket) + 1))


class DraftOnlyUpAmplitudeTest(unittest.TestCase):
    """eq:murayama-only-up-offshell, the q = u bracket, eq:murayama-only-up-three and the channel table."""

    @classmethod
    def setUpClass(cls):
        cls.parameters, cls.data = only_up()
        cls.parent, cls.neutral = neutral_legs(cls.data)

    def test_offshell_coupling_and_up_bracket(self):
        B0, fpi, fa = self.parameters['B_0'], self.parameters['f_pi'], self.parameters['f_a']
        up = QUARK_MASSES[0]
        legs = [self.parent] + self.neutral
        radial_mass_squared = draft_only_up_radial_mass_squared(self.parameters)
        self.assertAlmostEqual(radial_mass_squared, self.data['mass_squared'][0, 0], places=12)
        for x in legs:
            # 2 B0 (M^{-1}_{m_x^2}(P_01))_uu + 1 = m_x^2 f_a^2/(B0 m_u f_pi^2)
            bracket = draft_neutral_bracket(self.parameters, (1, 0, 0), x[1], 0)
            self.assertLess(abs(bracket - x[1] * fa ** 2 / (B0 * up * fpi ** 2)), 1e-9 * abs(bracket))
            for y in legs:
                for momentum_squared in (0.0, 0.41, radial_mass_squared):
                    direct = direct_radial_vertex(self.parameters, 'polar', x, y, momentum_squared)
                    draft = draft_only_up_coupling(self.parameters, x[1], y[1], momentum_squared)
                    general = draft_offshell_vertex(self.parameters, x[1], y[1], momentum_squared)
                    self.assertLess(abs(direct - draft), 1e-8 * abs(draft))
                    self.assertLess(abs(general - draft), 1e-8 * abs(draft))
                # eq:murayama-only-up-cartesian-vertex at q^2 = m_sigma^2
                inverse = [draft_inverse_of_p01(self.parameters, leg[1]) for leg in (x, y)]
                cartesian = ((2 * self.parameters['lambda'] * fa + 2 * B0 ** 3 * fpi ** 2 / fa ** 3 * np.trace(
                    self.parameters['sources']['S'][(1, 0)] @ anticommutator(*inverse)).real)
                    / np.sqrt(draft_normalization(self.parameters, x[1]) * draft_normalization(self.parameters, y[1])))
                draft = draft_only_up_coupling(self.parameters, x[1], y[1], radial_mass_squared)
                self.assertLess(abs(cartesian - draft), 1e-8 * abs(draft))

    def test_charged_couplings(self):
        """eq:murayama-charged, and g = 0 between a neutral leg and a charged or K^0-type leg."""
        states = final_states(self.data)
        up = QUARK_MASSES[0]
        expected = {'pi0 pi+ pi-': self.parameters['B_0'] * up / self.parameters['f_a'],
                    'pi0 K+ K-': self.parameters['B_0'] * up / self.parameters['f_a'], 'pi0 K0 K0bar': 0.0}
        for name, value in expected.items():
            neutral, plus, minus = states[name]
            for momentum_squared in (0.2, 0.7):
                coupling = direct_radial_vertex(self.parameters, 'polar', plus, minus, momentum_squared)
                self.assertLess(abs(coupling - value), 1e-12 * self.parameters['B_0'] * up / self.parameters['f_a'])
                for leg in (self.parent, neutral):
                    self.assertLess(abs(direct_radial_vertex(self.parameters, 'polar', leg, plus, momentum_squared)),
                                    1e-12 * self.parameters['B_0'] * up / self.parameters['f_a'])

    def draft_amplitude(self, daughters, invariants):
        """eq:murayama-only-up-three with the draft couplings; returns (V4, list of sigma channels)."""
        B0, fa = self.parameters['B_0'], self.parameters['f_a']
        up = QUARK_MASSES[0]
        legs = [self.parent] + daughters
        masses = [leg[1] for leg in daughters]
        neutral = [abs(leg[0][1]) > 0 for leg in legs]

        def coupling(x, y, momentum_squared):
            if neutral[x] and neutral[y]:
                return draft_only_up_coupling(self.parameters, legs[x][1], legs[y][1], momentum_squared)
            if neutral[x] or neutral[y]:
                return 0.0
            matrix = column_from_direction(legs[x][0])[0]
            contains_up = np.any(np.abs(matrix[0]) > 0) or np.any(np.abs(matrix[:, 0]) > 0)
            return B0 * up / fa if contains_up else 0.0  # pi+pi-, K+K-: B0 m_u/f_a; K0 K0bar: 0
        radial_mass_squared = draft_only_up_radial_mass_squared(self.parameters)
        channels = []
        for attached, (first, second) in PAIRINGS:
            momentum_squared = invariants[frozenset({first, second})]
            channels.append(coupling(0, attached, momentum_squared) * coupling(first, second, momentum_squared)
                            / (radial_mass_squared - momentum_squared))
        if all(neutral):
            contact = draft_neutral_quartic(self.parameters, (1, 0, 0), [self.parent[1]] + masses)
        else:
            pair = off_diagonal_pair(column_from_direction(daughters[1][0])[0])
            contact = draft_charged_quartic(self.parameters, (1, 0, 0), pair, self.parent[1], masses[0],
                                            invariants_tuple(invariants))
        return contact, channels

    def test_amplitude_and_channel_table(self):
        table = {'pi0 eta eta_prime': ('zero', (True, True, True)),
                 'pi0 pi0 eta': ('neutral', (True, True, True)), 'pi0 pi0 pi0': ('neutral', (True, True, True)),
                 'pi0 pi+ pi-': ('charged', (True, False, False)), 'eta pi+ pi-': ('charged', (True, False, False)),
                 'pi0 K+ K-': ('charged', (True, False, False)), 'eta K+ K-': ('charged', (True, False, False)),
                 'pi0 K0 K0bar': ('charged', (False, False, False))}
        states = final_states(self.data)
        self.assertEqual(set(states), set(table))
        for name, daughters in states.items():
            legs = [self.parent] + daughters
            pieces = direct_pieces(self.parameters, self.data, [leg[0] for leg in legs], [leg[1] for leg in legs],
                                   *POINT)
            contact, channels = self.draft_amplitude(daughters, pieces['invariants'])
            scale = abs(pieces['total'])
            contact_kind, open_channels = table[name]
            with self.subTest(final_state=name):
                self.assertLess(np.max(np.abs(pieces['non_radial_exchange'])), 1e-13 * scale)
                # Precision limited by N(m^2) near meson-like poles: 1e-8.
                self.assertLess(abs(contact + sum(channels) - pieces['total']), 1e-8 * scale)
                np.testing.assert_allclose(np.array(channels), pieces['radial_channels'].real, rtol=1e-8,
                                           atol=1e-14 * scale)
                for value, is_open in zip(channels, open_channels):
                    self.assertEqual(abs(value) > 1e-25, is_open)
                if contact_kind == 'zero':  # compared with the size of a non-vanishing neutral contact
                    reference = draft_neutral_quartic(self.parameters, (1, 0, 0),
                                                      [self.parent[1]] + [leg[1] for leg in states['pi0 pi0 eta']])
                    self.assertLess(abs(contact), 1e-12 * abs(reference))
                else:
                    self.assertGreater(abs(contact), 1e-3 * scale)


class DraftDividedDifferenceTest(unittest.TestCase):
    """eq:murayama-quartic-single-flavour, the pole relation and the product stated before
    eq:murayama-divided-difference, the factor -3/(2 B0 m_1^2) and eq:murayama-divided-difference."""

    def test_pole_relation_product_and_identity(self):
        for charges in ((1, 0, 0), (2, 0, 0), (0, 0, 1), (-1, 0, 0)):
            parameters = charge_model_parameters(charges)
            data = spectrum(parameters)
            parent, others = neutral_legs(data)
            poles = [parent[1]] + [leg[1] for leg in others]
            self.assertEqual(len(poles), 4)
            B0, singlet = parameters['B_0'], parameters['m_1_squared']
            quark_masses = np.array(QUARK_MASSES)
            nodes = [draft_cosine_squared_mass(parameters, charges, q) for q in range(3)]
            np.testing.assert_allclose(nodes, cosine_direction_masses_squared(parameters, charges), rtol=1e-15)
            with self.subTest(charges=charges):
                # [2 B0 (M^{-1}_{m_x^2}(P_01))_qq + n_q][m_x^2 - 2 B0 m_q (1 + n_q^2 f_pi^2/(2 f_a^2))]
                #   = (2 B0 m_1^2/3) <M^{-1}_{m_x^2}(P_01)>
                for mass_squared in poles:
                    right = 2 * B0 * singlet / 3 * np.trace(draft_inverse_of_p01(parameters, mass_squared)).real
                    for q in range(3):
                        left = draft_neutral_bracket(parameters, charges, mass_squared, q) * (mass_squared - nodes[q])
                        self.assertLess(abs(left - right), 1e-9 * abs(right))
                # prod_x [node_q - m_x^2] = -(2 B0 m_1^2/3) m_q prod_{p != q} [node_q - node_p]
                for q in range(3):
                    left = np.prod([nodes[q] - m2 for m2 in poles])
                    differences = np.prod([nodes[q] - nodes[p] for p in range(3) if p != q])
                    right = -2 * B0 * singlet / 3 * quark_masses[q] * differences
                    self.assertLess(abs(left - right), 1e-9 * abs(right))
                    # each term of the sum over q in eq:murayama-quartic-single-flavour is
                    # -3/(2 B0 m_1^2) times the term of eq:murayama-divided-difference
                    term = quark_masses[q] / np.prod([m2 - nodes[q] for m2 in poles])
                    self.assertLess(abs(term + 3 / (2 * B0 * singlet) / differences), 1e-9 * abs(term))
                terms = [1 / np.prod([nodes[q] - nodes[p] for p in range(3) if p != q]) for q in range(3)]
                self.assertLess(abs(sum(terms)), 1e-12 * max(abs(t) for t in terms))
                # closed form on four distinct poles: zero relative to a repeated-daughter value
                distinct = draft_single_flavour_quartic(parameters, charges, poles)
                repeated = draft_single_flavour_quartic(parameters, charges, [poles[0], poles[1], poles[1], poles[2]])
                self.assertLess(abs(distinct), 1e-10 * abs(repeated))
                self.assertLess(abs(repeated - draft_neutral_quartic(parameters, charges,
                                                                     [poles[0], poles[1], poles[1], poles[2]])),
                                1e-9 * abs(repeated))
                self.assertLess(abs(distinct - draft_neutral_quartic(parameters, charges, poles)),
                                1e-9 * abs(repeated))


class DraftParityTest(unittest.TestCase):
    """eq:murayama-parity and eq:murayama-cubic-parity-trace."""

    def test_cubic_source_terms(self):
        generator = np.random.default_rng(7)
        B0, fpi, fa = 2.6, 0.0924, 0.4

        def hermitian(scale):
            matrix = generator.normal(size=(3, 3)) + 1j * generator.normal(size=(3, 3))
            return scale * (matrix + matrix.conj().T) / 2
        scalar = [np.diag([0.004, 0.008, 0.1])] + [hermitian(0.02) for _ in range(4)]
        pseudoscalar = [hermitian(0.02) for _ in range(5)]
        meson = hermitian(0.05)
        angular = 0.07
        points = 64
        coefficient = 0
        for n in range(points):
            t = np.exp(2j * np.pi * n / points)
            x = t * angular / fa
            source = sum((scalar[s] + 1j * pseudoscalar[s]) * x ** s for s in range(5))
            source_bar = sum((scalar[s] - 1j * pseudoscalar[s]) * x ** s for s in range(5))
            field = expm(1j * t * meson / fpi)
            field_bar = expm(-1j * t * meson / fpi)
            value = fpi ** 2 / 4 * 2 * B0 * np.trace(source @ field_bar + source_bar @ field)
            coefficient += value * t ** -3 / points
        draft = B0 * fpi ** 2 * (angular ** 3 / fa ** 3 * np.trace(scalar[3])
                                 + angular ** 2 / (fa ** 2 * fpi) * np.trace(pseudoscalar[2] @ meson)
                                 - angular / (2 * fa * fpi ** 2) * np.trace(scalar[1] @ meson @ meson)
                                 - 1 / (6 * fpi ** 3) * np.trace(pseudoscalar[0] @ meson @ meson @ meson))
        self.assertLess(abs(coefficient - draft), 1e-12 * abs(draft))

    def test_parity_condition(self):
        def source(parameters, x):
            return sum((parameters['sources']['S'].get((0, s), 0) + 1j * parameters['sources']['P'].get((0, s), 0))
                       * x ** s for s in range(5))

        def coefficients_have_parity(parameters):
            return all(np.max(np.abs(matrix)) < 1e-14 for kind, table in parameters['sources'].items()
                       for (r, s), matrix in table.items() if (kind == 'S') == (s % 2 == 1))

        conserving = [charge_model_parameters((1, 0, 0)), charge_model_parameters((2, 0, 1)),
                      axial_rotation_parameters(0.5 * LAMBDA[6].real + 0.3 * LAMBDA[4].real)]
        base = charge_model_parameters((1, 0, 0))
        violating = [with_extra_sources(base, {('P', (0, 2)): 0.002 * np.diag([0.0, 1.0, -1.0])}),
                     with_extra_sources(base, {('S', (0, 1)): 0.002 * np.diag([1.0, -1.0, 0.0])})]
        for parameters, expected in [(p, True) for p in conserving] + [(p, False) for p in violating]:
            difference = max(np.max(np.abs(source(parameters, -x) - source(parameters, x).conj().T))
                             for x in (0.3, 0.7))
            self.assertEqual(difference < 1e-14, expected)
            self.assertEqual(coefficients_have_parity(parameters), expected)

    def test_chiral_shift_commutes_with_the_parity(self):
        """g(-a) = g(a)^dagger and U(-Phi', -a) = U(Phi', a)^dagger for U = g U' g (random sources)."""
        parameters = random_parameters(2026)
        data = spectrum(parameters)
        generator_matrix = alp_shift_matrix(parameters, data['masses_squared'][data['parent']])
        meson = np.einsum('k,kij->ij', np.random.default_rng(3).normal(size=9) * 0.05, LAMBDA)
        fpi = parameters['f_pi']
        for angular in (0.03, 0.11):
            rotation = lambda a: expm(1j * a * generator_matrix / (2 * fpi))
            field = lambda a, phi: rotation(a) @ expm(1j * phi / fpi) @ rotation(a)
            self.assertLess(np.max(np.abs(rotation(-angular) - rotation(angular).conj().T)), 1e-14)
            self.assertLess(np.max(np.abs(field(-angular, -meson) - field(angular, meson).conj().T)), 1e-13)


if __name__ == '__main__':
    unittest.main()
