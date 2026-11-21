"""Unit tests for derive_alp_loop_power_counting.py (units GeV).

Run: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_loop_power_counting.py
"""
import itertools
import unittest

import numpy as np
import sympy as sp

from derive_alp_loop_power_counting import (
    ILLUSTRATIVE, bubble_to_tree_ratios, bubble_unitarity_check, charged_loop_self_energy,
    coleman_weinberg_potential, dressed_frame_vertex, dressed_mass_derivative_closed_form,
    dressed_mass_function, flavour_basis_generators, illustrative_neutral_theory,
    illustrative_regime_numbers, loop_function_A0, loop_kinematics,
    mass_term_vertex_closed_form, neutral_derivative_tensor, neutral_dressed_lagrangian,
    pole_sum_decomposition, regime_map, singlet_derivative_terms, tadpole_mass_correction,
    toy_frame_comparison, toy_frame_series, unmixed_closed_forms)
from derive_alp_lo_joint_vertices import vertex


def hermitian_generators():
    return [sp.Matrix([[1, 2, 0], [2, -1, 1], [0, 1, 3]]),
            sp.Matrix([[0, sp.I, 1], [-sp.I, 2, 0], [1, 0, -2]]),
            sp.Matrix([[2, 0, sp.I], [0, 1, 1], [-sp.I, 1, 0]])]


class DressedVertexTests(unittest.TestCase):
    """L1: vertex scaling in the dressed frame."""

    def setUp(self):
        self.masses = (sp.Rational(1, 3), sp.Rational(1, 2), 2)
        self.axial = (sp.Rational(3, 2), sp.Rational(-1, 2), 1)
        self.pseudoscalar = (sp.Rational(1, 4), 1, sp.Rational(-2, 3))
        self.derivatives = [dressed_mass_derivative_closed_form(n, self.masses, self.axial,
                                                                self.pseudoscalar)
                            for n in range(5)]
        self.f, self.f_a, self.B_0 = sp.Rational(2, 3), sp.Rational(5, 4), sp.Rational(7, 3)

    def test_dressed_mass_taylor_derivatives(self):
        x = sp.Symbol('x', real=True)
        function = dressed_mass_function(x, self.masses, self.axial, self.pseudoscalar)
        for order in range(5):
            series = function.diff(x, order).subs(x, 0)
            self.assertEqual(sp.simplify(series-self.derivatives[order]), sp.zeros(3))

    def test_closed_form_matches_general_vertex(self):
        generators = hermitian_generators()
        zero = sp.zeros(3)
        for number_of_alps, number_of_mesons in ((1, 1), (2, 1), (1, 2), (2, 2), (1, 3), (0, 3)):
            mesons = generators[:number_of_mesons]
            legs = mesons+[zero]*number_of_alps
            alps = [0]*number_of_mesons+[1]*number_of_alps
            general = dressed_frame_vertex(legs, alps, [0]*len(legs), self.derivatives,
                                           self.B_0, 0, self.f, self.f_a, 0, 0)
            closed = mass_term_vertex_closed_form(mesons, number_of_alps,
                                                  self.derivatives[number_of_alps],
                                                  self.f, self.f_a, self.B_0)
            self.assertEqual(sp.simplify(general-closed), 0,
                             msg=f'n_a={number_of_alps}, n_P={number_of_mesons}')

    def test_alp_vertices_are_momentum_independent(self):
        generators = hermitian_generators()
        zero = sp.zeros(3)
        momenta = sp.symbols('q1:5', real=True)
        legs = generators[:2]+[zero, zero]
        with_momenta = dressed_frame_vertex(legs, [0, 0, 1, 1], momenta, self.derivatives,
                                            self.B_0, sp.Rational(1, 3), self.f, self.f_a, 1, 0)
        without = dressed_frame_vertex(legs, [0, 0, 1, 1], [0]*4, self.derivatives,
                                       self.B_0, sp.Rational(1, 3), self.f, self.f_a, 1, 0)
        self.assertEqual(sp.simplify(with_momenta-without), 0)

    def test_kinetic_term_has_no_alp_vertex(self):
        zero = sp.zeros(3)
        momenta = sp.symbols('q1:5', real=True)
        generators = hermitian_generators()
        no_mass = [zero]*5
        value = dressed_frame_vertex(generators[:2]+[zero, zero], [0, 0, 1, 1], momenta, no_mass,
                                     self.B_0, 0, self.f, self.f_a, 0, 0)
        self.assertEqual(value, 0)

    def test_linear_topological_source_is_only_quadratic(self):
        zero = sp.zeros(3)
        neutral, _ = flavour_basis_generators()
        no_mass = [zero]*5
        for legs, alps in (([neutral[0], zero, zero], [0, 1, 1]),
                           ([neutral[0], neutral[1], zero, zero], [0, 0, 1, 1]),
                           ([zero]*4, [1]*4)):
            value = dressed_frame_vertex(legs, alps, [0]*len(legs), no_mass, self.B_0,
                                         sp.Rational(5, 2), self.f, self.f_a, 1, 0)
            self.assertEqual(value, 0)
        quadratic = dressed_frame_vertex([zero, zero], [1, 1], [0, 0], no_mass, self.B_0,
                                         sp.Rational(5, 2), self.f, self.f_a, 1, 0)
        self.assertNotEqual(quadratic, 0)

    def test_neutral_closed_lagrangian_matches_general_vertices(self):
        fields = sp.symbols('P_u P_d P_s a', real=True)
        c_GA, m0sq, ma0sq = sp.Rational(4, 3), sp.Rational(3, 5), sp.Rational(1, 7)
        lagrangian = neutral_dressed_lagrangian(fields[:3], fields[3], self.masses, self.axial,
                                                self.pseudoscalar, self.B_0, c_GA, self.f,
                                                self.f_a, m0sq, ma0sq)
        neutral, _ = flavour_basis_generators()
        generators = neutral+[sp.zeros(3)]
        alps = [0, 0, 0, 1]
        for indices in ((0, 0, 0, 0), (0, 0, 3, 3), (2, 3, 3, 3), (3, 3, 3, 3), (1, 1, 1, 3),
                        (0, 1, 3, 3)):
            closed = sp.diff(lagrangian, *[fields[i] for i in indices]).subs(
                {field: 0 for field in fields})
            general = dressed_frame_vertex([generators[i] for i in indices],
                                           [alps[i] for i in indices], [0]*4, self.derivatives,
                                           self.B_0, c_GA, self.f, self.f_a, m0sq, ma0sq)
            self.assertEqual(sp.simplify(closed-general), 0, msg=str(indices))
        hessian = neutral_derivative_tensor(lagrangian, fields, 2)
        for i, j in itertools.combinations_with_replacement(range(4), 2):
            general = dressed_frame_vertex([generators[i], generators[j]], [alps[i], alps[j]],
                                           [0, 0], self.derivatives, self.B_0, c_GA, self.f,
                                           self.f_a, m0sq, ma0sq)
            self.assertEqual(sp.simplify(hessian[i, j]-general), 0)

    def test_neutral_kinetic_quartic_vanishes(self):
        neutral, _ = flavour_basis_generators()
        momenta = sp.symbols('q1:5', real=True)
        mass = [sp.zeros(3)]*5
        value = dressed_frame_vertex([neutral[0], neutral[0], neutral[1], neutral[1]], [0]*4,
                                     momenta, mass, self.B_0, 0, self.f, self.f_a, 0, 0)
        self.assertEqual(sp.expand(value), 0)


class TadpoleTests(unittest.TestCase):
    """L2: combinatorics, normalization and closed forms."""

    def test_factor_one_half_against_coleman_weinberg(self):
        rng = np.random.default_rng(3)
        size = 3
        mass_squared = np.diag([0.3, 0.7, 1.4]) + 0.05*np.ones((3, 3))
        quartic = rng.normal(size=(size,)*4)
        quartic = sum(np.transpose(quartic, p) for p in itertools.permutations(range(4)))/24
        scale = 0.8
        eigenvalues, vectors = np.linalg.eigh(mass_squared)
        lagrangian_vertex = -quartic  # V^L = -d^4 V
        predicted = tadpole_mass_correction(lagrangian_vertex, eigenvalues, vectors, scale)
        step = 1e-3
        numerical = np.zeros((size, size))
        for i in range(size):
            for j in range(size):
                def value(a, b):
                    field = np.zeros(size)
                    field[i] += a
                    field[j] += b
                    return coleman_weinberg_potential(field, mass_squared, quartic, scale)
                numerical[i, j] = (value(step, step)-value(step, -step)-value(-step, step)
                                   + value(-step, -step))/(4*step**2)
        np.testing.assert_allclose(numerical, predicted, rtol=1e-5, atol=1e-9)

    def test_su2_chiral_logarithm(self):
        f, B0, m, s, mu2 = sp.symbols('f B_0 m s mu2', positive=True)
        tau = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]),
               sp.Matrix([[1, 0], [0, -1]])]
        q = sp.symbols('q1:5', real=True)
        zero = sp.zeros(2)
        chi = [2*B0*m*sp.eye(2)]+[zero]*4
        mass_sq = 2*B0*m
        total = sum(loop_kinematics(vertex([tau[2], tau[2], tau[k], tau[k]], [0]*4, q, chi,
                                           [zero]*5, [zero]*5, [0]*5, f, 1, 0, 0),
                                    q, s, mass_sq) for k in range(3))
        sigma = -sp.Rational(1, 2)*total*loop_function_A0(mass_sq, mu2)
        pole = mass_sq-sigma.subs(s, mass_sq)
        log_coefficient = -sp.diff(pole, mu2)*mu2
        self.assertEqual(sp.simplify(log_coefficient-mass_sq**2/(32*sp.pi**2*f**2)), 0)

    def test_unmixed_closed_forms_against_tensor_contraction(self):
        forms = unmixed_closed_forms()
        sym = forms['symbols']
        f, f_a, B_0 = sym['f'], sym['f_a'], sym['B_0']
        masses, axial, pseudoscalar = sym['masses'], sym['axial'], sym['pseudoscalar']
        fields = sp.symbols('P_u P_d P_s a', real=True)
        lagrangian = neutral_dressed_lagrangian(fields[:3], fields[3], masses, axial,
                                                pseudoscalar, B_0, sp.Symbol('c_GA'), f, f_a,
                                                sp.Symbol('m0sq'), 0)
        zero = {field: 0 for field in fields}
        A0 = sym['A0']
        m_a_sq = sym['alp_mass_squared']
        # (a): ALP line, neutral meson entry
        for q in range(3):
            vertex_value = sp.diff(lagrangian, fields[q], fields[q], fields[3], fields[3]).subs(zero)
            correction = vertex_value/2*A0(m_a_sq)
            ratio = forms['alp_line_in_meson_entry_ratio'][(q, q)]
            self.assertEqual(sp.simplify(correction-ratio*2*B_0*masses[q]), 0)
        # (c): ALP line, ALP entry
        vertex_value = sp.diff(lagrangian, fields[3], 4).subs(zero)
        self.assertEqual(sp.simplify(vertex_value/2*A0(m_a_sq)
                                     - forms['alp_line_in_alp_entry']), 0)
        # tree quark part of M^2_aa (m0 -> 0)
        tree = -sp.diff(lagrangian, fields[3], 2).subs(zero).subs(sp.Symbol('m0sq'), 0)
        self.assertEqual(sp.simplify(tree-forms['tree_alp_entry_quark_part']), 0)
        # (b): neutral lines (flavour basis) plus charged lines
        neutral_part = sum(sp.diff(lagrangian, fields[3], fields[3], fields[q], fields[q]).subs(zero)
                           / 2*A0(2*B_0*masses[q]) for q in range(3))
        neutral_generators, _ = flavour_basis_generators()
        sigma = charged_loop_self_energy([sp.zeros(3)], [1], masses, axial, pseudoscalar, B_0,
                                         sp.Symbol('c_GA'), f, f_a, sp.Symbol('m0sq'),
                                         sp.Symbol('s'), 1, symbolic_loop=True)
        names = {'ud': (0, 1), 'us': (0, 2), 'ds': (1, 2)}
        charged_part = -sigma[0, 0].subs({sp.Symbol(f'A0_{k}'): A0(B_0*(masses[v[0]]+masses[v[1]]))
                                          for k, v in names.items()})
        self.assertEqual(sp.simplify(neutral_part+charged_part
                                     - forms['meson_lines_in_alp_entry']), 0)

    def test_pole_sum_approaches_unmixed_formula_at_large_ratio(self):
        data = pole_sum_decomposition(30.0, gluon_coefficient=-3.0)
        spectrum, n = data['spectrum'], data['alp_like_index']
        alp_mass_sq = spectrum.pole_squared_masses[n]
        f_a = data['info']['f_a']
        predicted = [2*loop_function_A0(alp_mass_sq, 1.0)/f_a**2*m for m in data['info']['two_B0_m']]
        np.testing.assert_allclose(np.diag(data['alp_like'])[:3], predicted, rtol=2e-3)
        self.assertAlmostEqual(alp_mass_sq/data['mass_squared'][3, 3], 1.0, delta=2e-3)

    def test_topological_term_gives_mass_but_no_self_coupling(self):
        base = illustrative_neutral_theory(2.0, gluon_coefficient=0.0)
        shifted = illustrative_neutral_theory(2.0, gluon_coefficient=4.0)
        self.assertNotAlmostEqual(base[1][3, 3], shifted[1][3, 3])
        np.testing.assert_allclose(base[2], shifted[2], rtol=1e-12, atol=1e-9)

    def test_shift_symmetric_point_has_massless_pole(self):
        _, mass_squared, _, info = illustrative_neutral_theory(2.0, gluon_coefficient=0.0)
        self.assertEqual(info['singlet'], 3.0)
        smallest = np.min(np.abs(np.linalg.eigvalsh(mass_squared)))
        self.assertLess(smallest, 1e-12)


class BubbleAndRegimeTests(unittest.TestCase):
    """L3 and L4."""

    def test_unitarity_of_alp_pair_bubble(self):
        residual, im_B0 = bubble_unitarity_check()
        self.assertEqual(residual, 0)
        w, m_a_sq = sp.symbols('w m_a^2', positive=True)
        self.assertEqual(sp.simplify(im_B0-sp.pi*sp.sqrt(1-4*m_a_sq/(4*m_a_sq+w))), 0)

    def test_bubble_ratios_match_regime_forms(self):
        ratios = bubble_to_tree_ratios()
        sym = ratios['symbols']
        regime, _, rsym = regime_map()
        substitutions = {rsym['r']: sym['f']/sym['f_a'], rsym['f']: sym['f'],
                         rsym['c_AP']: sym['c_AP'], rsym['kappa']: sym['kappa'],
                         rsym['m_qq']: sp.sqrt(2*sym['B_0']*sym['m_q']),
                         rsym['B0_loop']: sym['B0_loop']}
        for name in ('aP', 'aa'):
            self.assertEqual(sp.simplify(ratios[name]-regime[f'bubble_{name}'].subs(substitutions)), 0)
        self.assertEqual(sp.simplify(ratios['PP']-2*sym['B_0']*sym['m_q']*sym['B0_loop']
                                     / (16*sp.pi**2*sym['f']**2)), 0)

    def test_crossings_solve_ratio_equal_one(self):
        ratios, crossings, sym = regime_map()
        for name in ratios:
            self.assertEqual(sp.simplify(ratios[name].subs(sym['r'], crossings[name])-1), 0)

    def test_illustrative_numbers(self):
        numbers = illustrative_regime_numbers()
        self.assertAlmostEqual(numbers['alp_line_meson_entry_quark_mass'], 0.901, places=3)
        self.assertAlmostEqual(numbers['alp_line_meson_entry_topological'], 0.625, places=3)
        self.assertAlmostEqual(numbers['alp_line_meson_entry_fixed_mass'], 0.818, places=3)
        self.assertAlmostEqual(numbers['bubble_aa_uu'], 2.036, places=3)


class SingletDerivativeTermTests(unittest.TestCase):
    """L5."""

    def test_amplitudes_and_energy_growth(self):
        terms = singlet_derivative_terms()
        sym = terms['symbols']
        s = sym['s']
        self.assertEqual(sp.simplify(terms['pp_to_aa_kinetic']
                                     + 2*sym['uu_ratio']*sym['c_GA']**2
                                     * (s-2*sym['m_P_sq'])/sym['f_a']**2), 0)
        self.assertEqual(sp.simplify(terms['pp_to_a_eta0']
                                     + sp.sqrt(6)*sym['uu_ratio']*sym['c_GA']
                                     * (s-2*sym['m_P_sq'])/(sym['f']*sym['f_a'])), 0)
        self.assertEqual(sp.diff(terms['aa_to_aa'], s), 0)
        self.assertEqual(sp.diff(terms['aa_to_aa'], sym['t']), 0)
        self.assertNotEqual(sp.diff(terms['eta0_to_aa'], s), 0)
        self.assertEqual(sp.simplify(terms['alp_kinetic_shift']
                                     - 8*sym['V5_zero']*sym['c_GA']**2/sym['f_a']**2), 0)
        self.assertEqual(sp.simplify(terms['scale_pp_to_aa_nonet']*3
                                     - terms['scale_pp_to_aa_single']), 0)


class FrameComparisonTests(unittest.TestCase):
    """L6: toy model, derivative frame vs dressed frame."""

    def test_on_shell_shifts_agree_off_shell_differ(self):
        values = dict(f=0.092, m_P_sq=0.25, c=0.7, f_a=0.05, m_a0_sq=0.9, scale=1.0)
        output = toy_frame_comparison(values)
        np.testing.assert_allclose(output['derivative']['poles'], output['dressed']['poles'],
                                   rtol=1e-10)
        np.testing.assert_allclose(output['derivative']['shifts'], output['dressed']['shifts'],
                                   rtol=1e-9)
        self.assertGreater(abs(output['derivative']['off_shell'][1, 1]
                               - output['dressed']['off_shell'][1, 1]), 1e-3)

    def test_series_agree(self):
        series = toy_frame_series()
        self.assertEqual(sp.simplify(series['derivative']-series['dressed']), 0)


if __name__ == '__main__':
    unittest.main()
