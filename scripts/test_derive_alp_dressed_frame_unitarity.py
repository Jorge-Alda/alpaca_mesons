"""Unit tests for derive_alp_dressed_frame_unitarity.py (hypotheses H1-H4).

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_dressed_frame_unitarity.py
"""
import itertools
import random
import unittest

import mpmath as mp
import numpy as np
import sympy as sp

from derive_alp_dressed_frame_unitarity import (
    DIGITS, alp_like_and_pion_like, cayley_orthogonal, compare_frames, constant_vertex_closed_form,
    contact_unitarity_bounds, diagonal_linear_spurion_dressed_mass, dressed_mass_derivatives,
    dressed_sources, dressed_theta_derivatives, frame_theory, illustrative_inputs,
    kinetic_four_meson_vertex, linear_spurion_sources, on_shell_gram, pad_with_alp, pair_matrix,
    pole_amplitude, pullback_field_strength, rotate_tensor, s_wave_coefficient_tensor,
    singlet_s_wave, topological_vertex_closed_form, transformed_current_derivatives,
    unitarity_scale_squared)
from derive_alp_lo_joint_vertices import quadratic_blocks, vertex
from derive_alp_pole_projected_amplitudes import ccwz_vertex, projected_amplitude
from solve_alp_meson_spectrum import PoleSpectrum
from symbolic_utils import u3_generators

GENERATORS = [sp.Matrix(g) for g in u3_generators()]
TOLERANCE = sp.Float(10) ** (-(DIGITS - 10))


def symbolic_hermitian(name, size=2):
    entries = {}
    for i in range(size):
        entries[i, i] = sp.Symbol(f'{name}_{i}{i}', real=True)
        for j in range(i + 1, size):
            re, im = sp.symbols(f'{name}_re{i}{j} {name}_im{i}{j}', real=True)
            entries[i, j], entries[j, i] = re + sp.I * im, re - sp.I * im
    return sp.Matrix(size, size, lambda i, j: entries[i, j])


def rational_hermitian(generator, size=3):
    m = sp.Matrix(size, size, lambda i, j: sp.Rational(generator.randint(-4, 4), generator.randint(1, 3))
                  + sp.I * sp.Rational(generator.randint(-4, 4), generator.randint(1, 3)))
    return (m + m.H) / 2


def relative_difference(first, second):
    return abs(sp.N(first - second, DIGITS)) / max(abs(sp.N(first, DIGITS)), sp.Float(1))


class H1DressedFrame(unittest.TestCase):
    def test_connection_removed_symbolic_2x2(self):
        current = [symbolic_hermitian(f'C{n}') for n in range(4)]
        for value in transformed_current_derivatives(current, 3):
            self.assertEqual(value.applyfunc(sp.expand), sp.zeros(2))

    def test_connection_removed_rational_3x3_right_and_left(self):
        generator = random.Random(11)
        for _ in range(2):
            current = [rational_hermitian(generator) for _ in range(4)]
            for value in transformed_current_derivatives(current, 3):
                self.assertEqual(value, sp.zeros(3))

    def test_pullback_field_strength_vanishes(self):
        self.assertEqual(pullback_field_strength([symbolic_hermitian(f'D{n}') for n in range(3)]), sp.zeros(2))
        generator = random.Random(5)
        self.assertEqual(pullback_field_strength([rational_hermitian(generator) for _ in range(4)]), sp.zeros(3))

    def test_dressed_mass_closed_form_and_exponential(self):
        masses = sp.symbols('m_u m_d m_s', positive=True)
        axial = sp.symbols('a_u a_d a_s', real=True)
        pseudo = sp.symbols('p_u p_d p_s', real=True)
        chi, right, left, theta = linear_spurion_sources(sp.diag(*masses), sp.diag(*axial), sp.diag(*pseudo), 0, sp.Rational(1, 2))
        hat = dressed_mass_derivatives(chi, right, left, 4)
        for n in range(5):
            self.assertEqual((hat[n] - diagonal_linear_spurion_dressed_mass(masses, axial, pseudo, n)).applyfunc(sp.expand), sp.zeros(3))
        # non-diagonal: hat M = exp(-i x c_A) M exp(-i x c_A)
        x = sp.Symbol('x', real=True)
        data = illustrative_inputs(True)
        chi, right, left, theta = linear_spurion_sources(data['quark_mass'], data['axial'], data['pseudoscalar'], 0, sp.Rational(1, 2))
        rotation = sum(((-sp.I * x * data['axial'])**k / sp.factorial(k) for k in range(5)), sp.zeros(3))
        direct = (rotation * (chi[0] + x * chi[1]) * rotation).applyfunc(sp.expand)
        hat = dressed_mass_derivatives(chi, right, left, 4)
        for n in range(5):
            self.assertEqual(direct.applyfunc(lambda e: e.coeff(x, n)) * sp.factorial(n), hat[n])

    def test_dressed_theta_is_linear(self):
        data = illustrative_inputs(True)
        chi, right, left, theta = linear_spurion_sources(data['quark_mass'], data['axial'], data['pseudoscalar'], data['gluon'], 1)
        hat = dressed_theta_derivatives(theta, right, left, 4)
        self.assertEqual(hat, [0, data['gluon'] + data['axial'].trace(), 0, 0, 0])

    def test_dressed_kinetic_matrix_is_identity(self):
        data = illustrative_inputs(True)
        original = linear_spurion_sources(data['quark_mass'], data['axial'], data['pseudoscalar'], data['gluon'], data['B_0'])
        basis = [GENERATORS[i] for i in data['generator_indices']]
        for sources, expect_identity in ((original, False), (dressed_sources(*original), True)):
            kinetic, _ = quadratic_blocks(basis, *[array[:3] for array in sources], data['f_pi'], data['f_pi'],
                                          data['anomaly_mass_squared'], data['bare_alp_mass_squared'])
            self.assertEqual(kinetic == sp.eye(len(basis) + 1), expect_identity)


class H2OnShellEquivalence(unittest.TestCase):
    POINTS = [(sp.Integer(3), -sp.Rational(3, 5)), (sp.Integer(7), -sp.Rational(7, 4))]

    def check(self, ratio, off_diagonal, points):
        results = compare_frames(ratio, points, off_diagonal)
        original, dressed = results['original'], results['dressed']
        for first, second in zip(original['spectrum'].pole_squared_masses, dressed['spectrum'].pole_squared_masses):
            self.assertGreater(first, 0)  # stable vacuum, real poles
            self.assertLess(relative_difference(first, second), TOLERANCE)
        self.assertEqual((original['alp_like'], original['pion_like']), (dressed['alp_like'], dressed['pion_like']))
        # the ALP coordinate is not redefined: equal ALP components of the pole vectors
        for n in range(original['spectrum'].size):
            self.assertLess(relative_difference(original['spectrum'].pole_vectors[-1, n],
                                                dressed['spectrum'].pole_vectors[-1, n]), TOLERANCE)
        for key, value in original['amplitudes'].items():
            self.assertGreater(abs(value), 1e-6)
            self.assertLess(relative_difference(value, dressed['amplitudes'][key]), TOLERANCE, key)

    def test_diagonal_sources_three_ratios(self):
        for ratio in (sp.Rational(1, 10), sp.Integer(1), sp.Integer(3)):
            self.check(ratio, False, self.POINTS)

    def test_off_diagonal_axial_source(self):
        for ratio in (sp.Integer(1), sp.Integer(3)):
            self.check(ratio, True, self.POINTS)


class H3HighEnergy(unittest.TestCase):
    def test_kinetic_closed_form_matches_vertex(self):
        q = sp.symbols('q0:4', real=True)
        f_pi = sp.Rational(2, 3)
        zero = sp.zeros(3)
        for legs in [(1, 2, 1, 2), (6, 7, 3, 3), (1, 4, 6, 8), (3, 6, 0, 7)]:
            reference = vertex([GENERATORS[i] for i in legs], [0] * 4, q, [zero] * 5, [zero] * 5, [zero] * 5,
                               [0] * 5, f_pi, 1, 0, 0)
            closed = kinetic_four_meson_vertex(GENERATORS, f_pi, legs, lambda i, j: q[i] * q[j])
            self.assertEqual(sp.expand(reference - closed), 0)

    def test_dressed_alp_vertices_have_no_momentum_dependence(self):
        data = illustrative_inputs(True)
        dressed = dressed_sources(*linear_spurion_sources(data['quark_mass'], data['axial'], data['pseudoscalar'], data['gluon'], data['B_0']))
        q = sp.symbols('q0:4', real=True)
        zero = sp.zeros(3)
        cases = [([zero, GENERATORS[6], GENERATORS[7], GENERATORS[3]], [1, 0, 0, 0]),
                 ([zero, zero, GENERATORS[6], GENERATORS[3]], [1, 1, 0, 0]),
                 ([zero, zero, zero, GENERATORS[0]], [1, 1, 1, 0]),
                 ([zero, GENERATORS[7], GENERATORS[6]], [1, 0, 0])]
        for generators, alps in cases:
            result = vertex(generators, alps, q[:len(alps)], *dressed, data['f_pi'], data['f_pi'],
                            data['anomaly_mass_squared'], 0)
            self.assertFalse(result.has(*q))

    def test_largest_s_wave_eigenvalue(self):
        generators = [np.array(sp.N(g).tolist(), dtype=complex) for g in GENERATORS]
        for indices, expected in ((range(9), 3), (range(1, 9), 3), (range(1, 4), 2)):
            matrix, _ = pair_matrix(s_wave_coefficient_tensor([generators[i] for i in indices]))
            eigenvalues = np.linalg.eigvalsh(matrix)  # units of s / (32 pi f_pi^2)
            self.assertAlmostEqual(eigenvalues[-1], expected, places=10)

    def test_singlet_channel_with_common_mass(self):
        f_pi, mass_squared = sp.symbols('f m2', positive=True)
        for indices, flavours in (((1, 2, 3), 2), (tuple(range(1, 9)), 3)):
            wave, s = singlet_s_wave(indices, f_pi, mass_squared)
            expected = flavours * s / (32 * sp.pi * f_pi**2) - mass_squared / (16 * sp.pi * flavours * f_pi**2)
            self.assertEqual(sp.simplify(wave - expected), 0)

    def test_alp_states_do_not_change_eigenvalues(self):
        generators = [np.array(sp.N(g).tolist(), dtype=complex) for g in GENERATORS]
        padded = pad_with_alp(s_wave_coefficient_tensor(generators))
        reference = np.sort(np.linalg.eigvalsh(pair_matrix(padded)[0]))
        rotations = [np.array(sp.N(cayley_orthogonal(10, 3)).tolist(), dtype=float)]
        data = illustrative_inputs(False)
        original = linear_spurion_sources(data['quark_mass'], data['axial'], data['pseudoscalar'], data['gluon'], data['B_0'])
        for ratio in (sp.Integer(1), sp.Integer(3)):
            _, _, spectrum, _ = frame_theory(data['generator_indices'], dressed_sources(*original), data['f_pi'],
                                             data['f_pi'] / ratio, data['anomaly_mass_squared'], data['bare_alp_mass_squared'])
            neutral = np.array(sp.N(spectrum.pole_vectors, 20).tolist(), dtype=float)
            self.assertTrue(np.allclose(neutral.T @ neutral, np.eye(4)))
            full = np.eye(10)
            positions = list(data['generator_indices']) + [9]
            for row, i in enumerate(positions):
                for column, j in enumerate(positions):
                    full[i, j] = neutral[row, column]
            rotations.append(full)
        for rotation in rotations:
            rotated = np.sort(np.linalg.eigvalsh(pair_matrix(rotate_tensor(padded, rotation))[0]))
            self.assertTrue(np.allclose(rotated, reference, atol=1e-10))
            self.assertAlmostEqual(rotated[-1], 3, places=10)

    def test_full_amplitude_approaches_kinetic_contraction(self):
        data = illustrative_inputs(True)
        original = linear_spurion_sources(data['quark_mass'], data['axial'], data['pseudoscalar'], data['gluon'], data['B_0'])
        _, _, spectrum, function = frame_theory(data['generator_indices'], dressed_sources(*original), data['f_pi'],
                                                data['f_pi'], data['anomaly_mass_squared'], data['bare_alp_mass_squared'])
        size = spectrum.size
        vectors = spectrum.pole_vectors
        alp, pion = alp_like_and_pion_like(spectrum, 0)
        kaon = max(range(size), key=lambda n: abs(vectors[3, n]))
        basis = [GENERATORS[i] for i in data['generator_indices']]
        for indices in ([kaon, pion, kaon, alp], [kaon, pion, kaon, pion]):
            remainders = []
            for s, t in ((sp.Integer(2000), -sp.Integer(700)), (sp.Integer(8000), -sp.Integer(3000))):
                full = pole_amplitude(spectrum, function, indices, s, t)
                gram = on_shell_gram([spectrum.pole_squared_masses[i] for i in indices], s, t)
                generators = [sum((vectors[k, n] * basis[k] for k in range(size - 1)), sp.zeros(3)) for n in indices]
                asymptotic = sp.re(sp.N(kinetic_four_meson_vertex(generators, data['f_pi'], (0, 1, 2, 3),
                                                                  lambda i, j: gram[i, j])))
                self.assertGreater(abs(asymptotic), 1000)
                remainders.append(sp.re(sp.N(full - asymptotic)))
                self.assertLess(abs(remainders[-1]), 1e-2 * abs(asymptotic))
            # the remainder does not grow with s: it is the O(s^0) part
            self.assertLess(abs(remainders[1] - remainders[0]), 1e-2 * abs(remainders[0]))

    def test_multilinear_vertex_function_matches_ccwz_vertex(self):
        data = illustrative_inputs(True)
        original = linear_spurion_sources(data['quark_mass'], data['axial'], data['pseudoscalar'], data['gluon'], data['B_0'])
        basis = [GENERATORS[i] for i in data['generator_indices']]
        f_a = data['f_pi'] / 3
        _, _, _, function = frame_theory(data['generator_indices'], original, data['f_pi'], f_a,
                                         data['anomaly_mass_squared'], data['bare_alp_mass_squared'])
        generator = random.Random(2)
        directions = [sp.Matrix([sp.Rational(generator.randint(-3, 3), 2) for _ in range(5)]) for _ in range(4)]
        gram = sp.Matrix(4, 4, lambda i, j: sp.Rational(generator.randint(-9, 9), 4))
        gram = (gram + gram.T) / 2
        momenta = [sp.eye(4)[:, i] for i in range(4)]
        parameters = (*original, data['f_pi'], f_a, data['anomaly_mass_squared'], data['bare_alp_mass_squared'])
        for count in (3, 4):
            reference = ccwz_vertex(basis, parameters, directions[:count], momenta[:count], gram)
            self.assertLess(abs(sp.N(function(directions[:count], momenta[:count], gram) - reference, DIGITS)),
                            TOLERANCE * max(1, abs(sp.N(reference))))

    def test_lambda7_decouples_with_real_symmetric_sources(self):
        data = illustrative_inputs(True)
        original = linear_spurion_sources(data['quark_mass'], data['axial'], data['pseudoscalar'], data['gluon'], data['B_0'])
        zero = sp.zeros(3)
        q = sp.symbols('q0:4', real=True)
        even = [GENERATORS[i] for i in (0, 3, 8, 6)] + [zero]
        for sources in (original, dressed_sources(*original)):
            _, mass = quadratic_blocks([GENERATORS[i] for i in (3, 8, 0, 6, 7)], *[array[:3] for array in sources],
                                       data['f_pi'], data['f_pi'], data['anomaly_mass_squared'], data['bare_alp_mass_squared'])
            self.assertEqual([mass[4, j] for j in range(6) if j != 4], [0] * 5)
            for others in itertools.combinations_with_replacement(range(5), 2):
                generators = [GENERATORS[7]] + [even[k] for k in others]
                alps = [0] + [1 if k == 4 else 0 for k in others]
                self.assertEqual(sp.expand(vertex(generators, alps, q[:3], *sources, data['f_pi'], data['f_pi'],
                                                  data['anomaly_mass_squared'], 0)), 0)
            for others in ((1, 4, 3), (3, 3, 4), (1, 2, 3)):
                generators = [GENERATORS[7]] + [even[k] for k in others]
                alps = [0] + [1 if k == 4 else 0 for k in others]
                self.assertEqual(sp.expand(vertex(generators, alps, q, *sources, data['f_pi'], data['f_pi'],
                                                  data['anomaly_mass_squared'], 0)), 0)

    def test_unitarity_scale_numbers(self):
        f_pi = sp.Rational(9207, 100000)
        self.assertAlmostEqual(float(sp.sqrt(unitarity_scale_squared(3, f_pi, 0, sp.Rational(1, 2)))), 0.37690, places=4)
        self.assertAlmostEqual(float(sp.sqrt(unitarity_scale_squared(2, f_pi, 0, sp.Rational(1, 2)))), 0.46162, places=4)
        f, m2 = sp.symbols('f m2', positive=True)
        self.assertEqual(sp.simplify(unitarity_scale_squared(3, f, m2, sp.Rational(1, 2)) - 16 * sp.pi * f**2 / 3 - 2 * m2 / 9), 0)


class H4ContactTerms(unittest.TestCase):
    def setUp(self):
        data = illustrative_inputs(True)
        data['pseudoscalar'] = sp.Matrix([[sp.Rational(1, 5), 0, 0], [0, -sp.Rational(1, 4), sp.Rational(1, 3)], [0, sp.Rational(1, 3), sp.Rational(1, 2)]])
        self.data = data
        self.dressed = dressed_sources(*linear_spurion_sources(data['quark_mass'], data['axial'], data['pseudoscalar'], data['gluon'], data['B_0']))
        self.hat_mass = [value / (2 * data['B_0']) for value in self.dressed[0]]
        self.f_pi, self.f_a = data['f_pi'], sp.Rational(1, 7)

    def mass_vertex(self, generators, alps):
        q = sp.symbols(f'q0:{len(alps)}', real=True)
        result = vertex(generators, alps, q, self.dressed[0], self.dressed[1], self.dressed[2], [0] * 5, self.f_pi, self.f_a, 0, 0)
        return sp.expand(result.subs({value: 0 for value in q}))

    def test_closed_form_constant_vertices(self):
        zero = sp.zeros(3)
        for mesons in [[], [3, 3], [6, 8], [0, 6, 3], [3, 8, 6, 0], [7, 7], [6]]:
            for alps in range(0, 5 - len(mesons)):
                if len(mesons) + alps < 2:
                    continue
                generators = [zero] * alps + [GENERATORS[i] for i in mesons]
                expected = constant_vertex_closed_form(self.hat_mass, [GENERATORS[i] for i in mesons], alps,
                                                       self.data['B_0'], self.f_pi, self.f_a)
                self.assertEqual(sp.nsimplify(sp.expand(self.mass_vertex(generators, [1] * alps + [0] * len(mesons)) - expected)), 0)

    def test_topological_vertices(self):
        zero = sp.zeros(3)
        q = sp.symbols('q0:4', real=True)
        m0 = sp.Rational(4, 5)
        for theta in ([0, sp.Rational(2, 3), 0, 0, 0], [0, sp.Rational(2, 3), sp.Rational(1, 2), sp.Rational(-3, 2), 0]):
            for mesons, alps in (([], 4), ([0], 3), ([0, 0], 2), ([], 3), ([0], 2)):
                generators = [zero] * alps + [GENERATORS[i] for i in mesons]
                result = vertex(generators, [1] * alps + [0] * len(mesons), q[:len(generators)], [zero] * 5, [zero] * 5,
                                [zero] * 5, theta, self.f_pi, self.f_a, m0, 0)
                expected = topological_vertex_closed_form(theta, [GENERATORS[i] for i in mesons], alps, m0, self.f_pi, self.f_a)
                self.assertEqual(sp.expand(result - expected), 0)
                if theta[2] == theta[3] == 0 and alps + len(mesons) > 2:
                    self.assertEqual(result, 0)

    def test_cubic_vertices_vanish_with_cp_conserving_sources(self):
        zero = sp.zeros(3)
        cp_odd = [0, 3, 8, 6]  # real symmetric generators
        for count in range(4):
            for mesons in itertools.combinations_with_replacement(cp_odd, count):
                alps = 3 - count
                generators = [zero] * alps + [GENERATORS[i] for i in mesons]
                q = sp.symbols('q0:3', real=True)
                result = vertex(generators, [1] * alps + [0] * count, q, *self.dressed, self.f_pi, self.f_a,
                                self.data['anomaly_mass_squared'], 0)
                self.assertEqual(sp.nsimplify(result), 0, mesons)

    def test_contact_amplitude_without_mixing(self):
        """With mass mixing neglected (pure field poles), aa->aa and a P3 -> a P3 are the bare contact vertices."""
        zero = sp.zeros(3)
        basis = [GENERATORS[i] for i in (3, 8, 0)]
        dressed = dressed_sources(*linear_spurion_sources(sp.diag(sp.Rational(1, 500), sp.Rational(1, 200), sp.Rational(1, 10)),
                                                          sp.diag(1, sp.Rational(1, 2), -sp.Rational(1, 3)),
                                                          sp.diag(sp.Rational(1, 5), 0, 0), 0, sp.Rational(27, 10)))
        parameters = (*dressed, self.f_pi, self.f_a, 0, 0)
        masses = [sp.Rational(1, 50), sp.Rational(1, 4), sp.Rational(3, 4), sp.Integer(2)]
        spectrum = PoleSpectrum(field_labels=('P3', 'P8', 'P0', 'a'), kinetic=sp.eye(4), mass_squared=sp.diag(*masses),
                                pole_squared_masses=masses, pole_vectors=sp.eye(4), residues=[], exact=True)

        def function(directions, momenta, gram):
            return ccwz_vertex(basis, parameters, directions, momenta, gram)
        hat_mass = [value / sp.Rational(27, 5) for value in dressed[0]]
        for indices, mesons, alps in (([3, 3, 3, 3], [], 4), ([3, 0, 3, 0], [basis[0]] * 2, 2)):
            gram = on_shell_gram([masses[i] for i in indices], sp.Integer(5), -sp.Integer(1))
            amplitude = projected_amplitude(spectrum, indices, gram, function)
            expected = constant_vertex_closed_form(hat_mass, mesons, alps, sp.Rational(27, 10), self.f_pi, self.f_a)
            self.assertEqual(sp.nsimplify(sp.expand(amplitude - expected)), 0)

    def test_bounds_closed_form_and_illustration(self):
        B_0, f_pi = sp.symbols('B_0 f_pi', positive=True)
        masses = sp.symbols('m_u m_d m_s', positive=True)
        c = sp.symbols('c_u c_d c_s', real=True)
        hat = [diagonal_linear_spurion_dressed_mass(masses, c, (0, 0, 0), n) for n in range(5)]
        waves, bounds, r = contact_unitarity_bounds(hat, GENERATORS[3], B_0, f_pi)
        weighted4 = sum(m * x**4 for m, x in zip(masses, c))
        weighted2 = masses[0] * c[0]**2 + masses[1] * c[1]**2
        self.assertEqual(sp.simplify(waves['a a -> a a'] - B_0 * r**4 * weighted4 / (2 * sp.pi * f_pi**2)), 0)
        self.assertEqual(sp.simplify(waves['a P -> a P'] - B_0 * r**2 * weighted2 / (4 * sp.pi * f_pi**2)), 0)
        self.assertEqual(sp.simplify(waves['P P -> a a'] - B_0 * r**2 * weighted2 / (8 * sp.pi * f_pi**2)), 0)
        numbers = {masses[0]: sp.Rational(216, 100000), masses[1]: sp.Rational(470, 100000), masses[2]: sp.Rational(935, 10000),
                   c[0]: 1, c[1]: 0, c[2]: 0, f_pi: sp.Rational(9207, 100000)}
        numbers[B_0] = sp.Rational(134977, 1000000)**2 / (numbers[masses[0]] + numbers[masses[1]])
        values = [float(bounds[k].subs(numbers)) for k in ('a a -> a a', 'a P -> a P', 'P P -> a a')]
        for value, expected in zip(values, (1.468, 3.047, 4.309)):
            self.assertAlmostEqual(value, expected, places=3)

    def test_general_dressed_mass_depends_on_axial_and_pseudoscalar_separately(self):
        """Beyond first order hat M_n is not a function of c_q^{AP} = c_q^A - c_q^P alone."""
        masses = sp.symbols('m_u m_d m_s', positive=True)
        same_difference = [diagonal_linear_spurion_dressed_mass(masses, (1, 0, 0), (0, 0, 0), 2),
                           diagonal_linear_spurion_dressed_mass(masses, (2, 1, 1), (1, 1, 1), 2)]
        self.assertNotEqual(same_difference[0], same_difference[1])
        self.assertEqual(diagonal_linear_spurion_dressed_mass(masses, (1, 0, 0), (0, 0, 0), 1),
                         diagonal_linear_spurion_dressed_mass(masses, (2, 1, 1), (1, 1, 1), 1))


if __name__ == '__main__':
    unittest.main()
