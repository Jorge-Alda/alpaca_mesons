"""Unit tests: a_phys -> P_i P_j P_k with radial-mode exchange, polar coordinates (GeV units).

Run: .venv/bin/python -m unittest scripts.test_derive_alp_three_pseudoscalar_radial_exchange
"""
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

from scripts.derive_alp_three_pseudoscalar_radial_exchange import (
    B_0, LAMBDA, PAIRINGS, analytic_amplitude, compare_point, component_quadratic_blocks, dalitz_invariants,
    derive_cartesian_quartic_contacts, derive_off_shell_radial_vertex, derive_polar_radial_mixing,
    derive_record_mapping, derive_symmetric_map_second_order, direct_amplitude, f_a, f_pi, impose_stationarity,
    leg_column, only_up_parameters, polar_vertex_trace_form, pole_vectors, quadratic_forms, quartic,
    radial_exchange_trace_form, radial_mass_squared_trace_form, radial_vertex_trace_form, random_parameters,
    record_amplitude_without_radial, spectrum, trace, with_radial_mass)

SEEDS = (2026, 7, 11)


def daughter_choices(data):
    """Three daughter assignments: three distinct mixed poles, identical pair plus an exact meson
    pole, two exact meson poles plus a mixed pole."""
    mixed = [n for n in data['mixed'] if n != data['parent']]
    unmixed = data['unmixed']
    return ((mixed[0], mixed[1], mixed[2]), (mixed[0], mixed[0], unmixed[0]), (unmixed[0], unmixed[1], mixed[1]))


class SymbolicTest(unittest.TestCase):
    def test_symmetric_chiral_map_has_no_quadratic_term(self):
        self.assertEqual(derive_symmetric_map_second_order(), (0, 0))

    def test_polar_radial_mixing_needs_s11_and_p10(self):
        sigma_angular, sigma_meson = derive_polar_radial_mixing()
        self.assertEqual(sp.simplify(sigma_angular - B_0 * f_pi ** 2 * sp.Symbol('tr_S_11', real=True) / f_a ** 2), 0)
        self.assertEqual(sigma_meson, {(1, 0): B_0 * f_pi / f_a})

    def test_off_shell_vertex(self):
        _, conjecture_difference, on_shell_difference, parent_pair = derive_off_shell_radial_vertex()
        self.assertEqual((conjecture_difference, on_shell_difference, parent_pair), (0, 0, 0))

    def test_cartesian_quartic_contacts(self):
        result = derive_cartesian_quartic_contacts()
        t = {name: sp.Symbol(name, real=True) for name in ('tr_S_02', 'tr_S_10', 'tr_S_12', 'tr_S_20',
                                                           'tr_P_01_Phi', 'tr_P_11_Phi', 'tr_S_10_PhiPhi')}
        expected_p4 = (B_0 * f_pi ** 2 / (24 * f_a ** 4) * (-16 * t['tr_S_02'] - 3 * t['tr_S_10'] + 12 * t['tr_S_12']
                                                           + 6 * t['tr_S_20']) - quartic / 4)
        self.assertEqual(sp.simplify(result['p^4'] - expected_p4), 0)
        self.assertEqual(sp.simplify(result['p^3 Phi'] - B_0 * f_pi * (3 * t['tr_P_11_Phi'] - 2 * t['tr_P_01_Phi'])
                                     / (6 * f_a ** 3)), 0)
        self.assertEqual(sp.simplify(result['p^2 Phi^2'] + B_0 * t['tr_S_10_PhiPhi'] / (4 * f_a ** 2)), 0)
        self.assertEqual(result['p^2 (must vanish)'], 0)

    def test_record_mapping(self):
        coupling, alp_entry = derive_record_mapping()
        self.assertEqual(coupling, B_0)
        self.assertEqual(alp_entry, 0)


class MappingNumericalTest(unittest.TestCase):
    def test_quadratic_blocks_equal_existing_record_blocks(self):
        """Polar Hessian = component_quadratic_blocks of the record with C = B0 P_01, m_0^2 = m_1^2,
        m_{a,0} = 0, d_a^2 chi_+|_0 = 8 B0 S_02/f_a^2."""
        parameters = random_parameters(2026)
        data = spectrum(parameters)
        S, P = parameters['sources']['S'], parameters['sources']['P']
        B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
        to_sympy = lambda m: sp.Matrix(m.tolist())
        _, record_mass = component_quadratic_blocks(
            [to_sympy(LAMBDA[k]) for k in range(9)], to_sympy(S[(0, 0)]), parameters['m_1_squared'],
            to_sympy(B0 * P[(0, 1)]), -f_pi ** 2 * trace(8 * B0 * S[(0, 2)] / fa ** 2) / 4)
        record_mass = np.array(sp.N(record_mass.subs({f_pi: fpi, f_a: fa, B_0: B0})).tolist(), dtype=complex)
        order = list(range(2, 11)) + [1]
        polar_mass = data['mass_squared'][np.ix_(order, order)]
        self.assertLess(np.max(np.abs(record_mass.imag)), 1e-14)
        self.assertLess(np.max(np.abs(record_mass.real - polar_mass)), 1e-13)
        self.assertLess(np.max(np.abs(data['kinetic'] - np.eye(11))), 1e-14)
        self.assertLess(np.max(np.abs(data['gradient'][1:])), 1e-14)


class AmplitudeTest(unittest.TestCase):
    """(a) direct polar = analytic; (b) direct polar = direct Cartesian, for random sources."""

    def test_random_sources(self):
        for seed in SEEDS:
            parameters = random_parameters(seed)
            data = spectrum(parameters)
            self.assertTrue(np.all(data['masses_squared'] > 0))
            self.assertGreater(radial_mass_squared_trace_form(parameters), 0)
            self.assertLess(abs(data['gradient'][0]), 1e-14)
            self.assertLess(np.max(np.abs(data['mass_squared'][0, 1:])), 1e-14)
            for daughters in daughter_choices(data):
                for s12, s13 in ((0.9, 1.4), (0.37, 2.2)):
                    with self.subTest(seed=seed, daughters=daughters, s12=s12):
                        point = compare_point(parameters, data, daughters, s12, s13)
                        scale = abs(point['polar']['total'])
                        self.assertLess(point['cartesian_hessian_difference'], 1e-14)
                        self.assertLess(point['polar']['radial_offdiagonal'], 1e-14)
                        self.assertLess(abs(point['polar']['total'] - point['analytic']) / scale, 1e-10)
                        self.assertLess(abs(point['polar']['total'] - point['cartesian']['total']) / scale, 1e-12)
                        radial = point['polar']['radial_exchange']
                        self.assertLess(abs(radial - point['radial_trace_form']) / scale, 1e-10)
                        self.assertLess(abs(point['polar']['total'] - radial - point['record_without_radial']) / scale,
                                        1e-10)
                        # The sigma term is resolved: dropping it fails the comparison by far.
                        self.assertGreater(abs(radial), 1e3 * abs(point['polar']['total'] - point['analytic']))

    def test_stationarity_is_needed(self):
        """(c) shifting m_phi^2 off stationarity breaks polar = Cartesian."""
        parameters = random_parameters(2026)
        parameters['m_phi_squared'] += 0.05
        totals = {}
        for coordinates in ('polar', 'cartesian'):
            _, kinetic, mass_squared = quadratic_forms(parameters, coordinates)
            masses_squared, vectors = pole_vectors(mass_squared)
            parent = int(np.argmax(np.abs(vectors[1])))
            others = [n for n in range(10) if n != parent][:3]
            legs = (parent,) + tuple(others)
            invariants = dalitz_invariants(masses_squared[parent], [masses_squared[n] for n in others], 0.9, 1.4)
            totals[coordinates] = direct_amplitude(parameters, coordinates, mass_squared, kinetic, masses_squared,
                                                   vectors, legs, invariants)['total']
        self.assertGreater(abs(totals['polar'] - totals['cartesian']), 1e-3 * abs(totals['polar']))


class ResidueTest(unittest.TestCase):
    def test_residue_at_radial_pole(self):
        """(d) (m_sigma^2 - s_12) M -> g_{a3} g_{12} with the on-shell couplings of the sigma -> PP record."""
        parameters = random_parameters(2026)
        data = spectrum(parameters)
        mixed = [n for n in data['mixed'] if n != data['parent']]
        legs = (data['parent'], mixed[0], mixed[1], mixed[2])
        masses_squared = data['masses_squared']
        radial_mass_squared = radial_mass_squared_trace_form(parameters)
        relative_step = 1e-6
        estimates = []
        for sign in (1, -1):
            s12 = radial_mass_squared * (1 + sign * relative_step)
            invariants = dalitz_invariants(masses_squared[legs[0]], [masses_squared[n] for n in legs[1:]], s12, 1.1)
            total = direct_amplitude(parameters, 'polar', data['mass_squared'], data['kinetic'], masses_squared,
                                     data['vectors'], legs, invariants)['total']
            estimates.append((radial_mass_squared - s12) * total)
        residue = sum(estimates) / 2
        record = (polar_vertex_trace_form(parameters, radial_mass_squared, masses_squared[legs[0]], masses_squared[legs[3]])
                  * polar_vertex_trace_form(parameters, radial_mass_squared, masses_squared[legs[1]],
                                            masses_squared[legs[2]]))
        self.assertLess(abs(residue - record), 1e-6 * abs(record))
        columns = [(leg_column(parameters, masses_squared[n], data['vectors'][:, n]), masses_squared[n]) for n in legs]
        off_shell_form = (radial_vertex_trace_form(parameters, columns[0], columns[3], radial_mass_squared)
                          * radial_vertex_trace_form(parameters, columns[1], columns[2], radial_mass_squared))
        self.assertLess(abs(off_shell_form - record), 1e-12 * abs(record))


class LargeQuarticTest(unittest.TestCase):
    def test_polar_radial_term_decouples_and_cartesian_pieces_grow(self):
        """(e) polar sigma term ~ 1/m_sigma^2; Cartesian contact and s exchange ~ lambda; sums equal."""
        base = random_parameters(2026)
        data = spectrum(base)
        mixed = sorted((n for n in data['mixed'] if n != data['parent']), key=lambda n: -abs(data['vectors'][1, n]))
        daughters = tuple(mixed[:3])
        rows = []
        for quartic_value in (8.0, 80.0, 800.0, 8000.0):
            parameters = dict(base)
            parameters['lambda'] = quartic_value
            parameters = impose_stationarity(parameters)
            point = compare_point(parameters, spectrum(parameters), daughters, 0.9, 1.4)
            rows.append((quartic_value, radial_mass_squared_trace_form(parameters), point))
        for quartic_value, radial_mass_squared, point in rows:
            scale = abs(point['polar']['total'])
            self.assertLess(abs(point['polar']['total'] - point['cartesian']['total']) / scale, 1e-9)
            self.assertLess(abs(point['polar']['total'] - point['analytic']) / scale, 1e-9)
        (_, mass_3, point_3), (_, mass_4, point_4) = rows[2], rows[3]
        # Polar: m_sigma^2 x (sigma term) tends to a constant.
        product_3 = mass_3 * point_3['polar']['radial_exchange']
        product_4 = mass_4 * point_4['polar']['radial_exchange']
        self.assertLess(abs(product_3 - product_4), 0.02 * abs(product_4))
        # Polar contact does not depend on lambda.
        self.assertLess(abs(point_3['polar']['contact'] - point_4['polar']['contact']), 1e-12 * abs(point_4['polar']['contact']))
        # Cartesian contact and s exchange grow linearly in lambda and cancel.
        for piece in ('contact', 'radial_exchange'):
            slope_3 = (point_3['cartesian'][piece] - rows[1][2]['cartesian'][piece]) / (800.0 - 80.0)
            slope_4 = (point_4['cartesian'][piece] - point_3['cartesian'][piece]) / (8000.0 - 800.0)
            self.assertLess(abs(slope_3 - slope_4), 0.02 * abs(slope_4))
            self.assertGreater(abs(point_4['cartesian'][piece]), 50 * abs(point_4['polar']['total']))


class OnlyUpTest(unittest.TestCase):
    def test_only_up_model(self):
        """(f) only-up coefficients; direct polar = Cartesian = analytic; factorized sigma vertex."""
        parameters = only_up_parameters()
        S, P = parameters['sources']['S'], parameters['sources']['P']
        up = 0.0022
        e1 = np.diag([1.0, 0.0, 0.0])
        self.assertTrue(np.allclose(S[(1, 0)], up * e1) and np.allclose(P[(0, 1)], -up * e1)
                        and np.allclose(P[(1, 1)], -up * e1) and np.allclose(S[(0, 2)], -up / 2 * e1)
                        and np.allclose(S[(1, 2)], -up / 2 * e1) and np.allclose(P[(0, 3)], up / 6 * e1)
                        and np.allclose(S[(0, 4)], up / 24 * e1))
        data = spectrum(parameters)
        mixed = [n for n in data['mixed'] if n != data['parent']]
        self.assertEqual(len(data['mixed']), 4)  # pi0, eta, eta' and the ALP-like pole
        for daughters in ((mixed[0], mixed[1], mixed[2]), (mixed[0], data['unmixed'][0], data['unmixed'][0])):
            point = compare_point(parameters, data, daughters, 0.3, 0.5)
            scale = abs(point['polar']['total'])
            self.assertLess(abs(point['polar']['total'] - point['analytic']) / scale, 1e-9)
            self.assertLess(abs(point['polar']['total'] - point['cartesian']['total']) / scale, 1e-11)
            self.assertLess(abs(point['polar']['radial_exchange'] - point['radial_trace_form']),
                            1e-9 * abs(point['polar']['radial_exchange']))
        # Factorized only-up vertex for neutral (diagonal) columns:
        # f_a g = B0 m_u (Lambda_x,11 + f_pi xi_x/f_a)(Lambda_y,11 + f_pi xi_y/f_a) + xi_x xi_y (q^2 - m_x^2 - m_y^2).
        B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
        masses_squared = data['masses_squared']
        legs = [data['parent']] + mixed
        q_squared = 0.77
        for x in legs:
            for y in legs:
                column_x, column_y = (leg_column(parameters, masses_squared[n], data['vectors'][:, n]) for n in (x, y))
                general = radial_vertex_trace_form(parameters, (column_x, masses_squared[x]),
                                                   (column_y, masses_squared[y]), q_squared)
                factorized = (B0 * up * (column_x[0][0, 0].real + fpi * column_x[1] / fa)
                              * (column_y[0][0, 0].real + fpi * column_y[1] / fa)
                              + column_x[1] * column_y[1] * (q_squared - masses_squared[x] - masses_squared[y])) / fa
                self.assertLess(abs(general - factorized), 1e-12 * max(abs(general), 1e-12))

    def test_only_up_implementation_has_no_radial_exchange(self):
        """only_up_decays.three_pseudoscalar_width receives no radial inputs (lambda, m_sigma, S_10)."""
        source = (Path(__file__).resolve().parent / 'implementation' / 'only_up_decays.py').read_text()
        self.assertIn('radial scalar', source)
        self.assertNotIn('lambda', source.split('def three_pseudoscalar_width')[1])


class WeakMixingTest(unittest.TestCase):
    def test_radial_term_is_relative_order_f_pi_over_f_a_squared(self):
        """(g) fixed sources and m_sigma^2: f_a M_without_sigma -> constant, f_a^3 (sigma term) -> constant."""
        base = random_parameters(2026)
        rows = []
        for decay_constant in (4.0, 8.0, 16.0):
            parameters = dict(base)
            parameters['f_a'] = decay_constant
            parameters = with_radial_mass(parameters, 1.0)
            data = spectrum(parameters)
            mixed = sorted(n for n in data['mixed'] if n != data['parent'])
            legs = (data['parent'], mixed[0], mixed[1], mixed[2])
            masses_squared = data['masses_squared']
            invariants = dalitz_invariants(masses_squared[legs[0]], [masses_squared[n] for n in legs[1:]], 0.3, 0.45)
            without_radial, radial = analytic_amplitude(parameters, masses_squared, data['vectors'], legs, invariants)
            rows.append((decay_constant, without_radial, radial))
        for (fa_1, without_1, radial_1), (fa_2, without_2, radial_2) in zip(rows, rows[1:]):
            self.assertLess(abs(fa_1 * without_1 - fa_2 * without_2), 0.02 * abs(fa_2 * without_2))
            self.assertLess(abs(fa_1 ** 3 * radial_1 - fa_2 ** 3 * radial_2), 0.05 * abs(fa_2 ** 3 * radial_2))
        self.assertLess(abs(rows[-1][2] / rows[-1][1]), 1e-3)


if __name__ == '__main__':
    unittest.main()
