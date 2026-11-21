"""Unit tests of scripts/derive_murayama_alp_unitarity.py (natural units, GeV).

Run: .venv/bin/python -m unittest scripts.test_derive_murayama_alp_unitarity
"""
import unittest

import sympy as sp

from scripts import derive_murayama_alp_unitarity as model
from scripts.derive_murayama_alp_unitarity import (
    B_0, m_u, f_pi, f_a, lam, mandelstam_s, mandelstam_t, cos_theta)

SAMPLE = {B_0: sp.Rational(2656, 1000), m_u: sp.Rational(216, 100000),
          f_pi: sp.Rational(9207, 100000), f_a: sp.Rational(1, 2), lam: sp.Rational(13, 10)}
UP_SOURCE = B_0*m_u*f_pi**2


class VacuumAndVerticesTest(unittest.TestCase):
    def test_vacuum_masses_match_sigma_record(self):
        self.assertEqual(sp.simplify(model.VACUUM_PHI_MASS_SQUARED - (UP_SOURCE/f_a**2 - lam*f_a**2)), 0)
        self.assertEqual(sp.simplify(model.RADIAL_MASS_SQUARED - (UP_SOURCE/f_a**2 + 2*lam*f_a**2)), 0)

    def test_polar_vertices(self):
        alp_mass, contact, potential_vertex, kinetic = model.polar_vertices()
        self.assertEqual(sp.simplify(alp_mass - UP_SOURCE/f_a**2), 0)
        self.assertEqual(sp.simplify(contact - UP_SOURCE/f_a**4), 0)
        self.assertEqual(sp.simplify(potential_vertex + UP_SOURCE/f_a**3), 0)
        self.assertEqual(kinetic, 1/f_a)
        q2 = sp.Symbol('q2')
        vertex = model.polar_sigma_alp_alp_vertex(q2, alp_mass)
        self.assertEqual(sp.simplify(vertex + (q2 - UP_SOURCE/f_a**2)/f_a), 0)

    def test_quartic_term_gives_no_phase_contact(self):
        """lambda |Phi|^4 does not depend on a in polar fields (Part D)."""
        potential = model.polar_potential()
        self.assertEqual(sp.diff(potential.coeff(lam), model.phase_field), 0)

    def test_cartesian_vertices(self):
        vertices = model.cartesian_vertices()
        self.assertEqual(sp.simplify(vertices['mass_p'] - UP_SOURCE/f_a**2), 0)
        self.assertEqual(sp.simplify(vertices['spp'] + 2*lam*f_a), 0)
        self.assertEqual(sp.simplify(vertices['sss'] + 6*lam*f_a), 0)
        self.assertEqual(vertices['pppp'], -6*lam)
        self.assertEqual(vertices['ssss'], -6*lam)
        self.assertEqual(vertices['sspp'], -2*lam)
        self.assertEqual(vertices['ppp'], 0)
        self.assertEqual(vertices['spp_s'], 0)


class AlpAlpAmplitudeTest(unittest.TestCase):
    def test_polar_equals_closed_form_and_cartesian(self):
        alp_mass = model.polar_vertices()[0]
        closed = model.closed_form_alp_alp_amplitude()
        self.assertEqual(sp.simplify(model.polar_alp_alp_amplitude(alp_mass) - closed), 0)
        self.assertEqual(sp.simplify(model.hypothesis_polar_form() - closed), 0)
        self.assertEqual(sp.simplify(model.cartesian_alp_alp_amplitude() - closed), 0)

    def test_shifted_external_mass_spoils_equality(self):
        shift, _, asymptotic, residual = model.polar_minus_cartesian_with_external_mass()
        self.assertEqual(residual, 0)
        self.assertEqual(sp.simplify(asymptotic - 8*shift/f_a**2), 0)

    def test_external_phase_potential(self):
        value = model.asymptotic_constant_with_external_potential()
        expected = (-6*lam - model.external_potential_fourth
                    + 8*model.external_potential_second/f_a**2)
        self.assertEqual(sp.simplify(value - expected), 0)


class SWaveTest(unittest.TestCase):
    def test_integrated_wave_equals_closed_form(self):
        integrated, _ = model.alp_alp_s_wave()
        closed = model.alp_alp_s_wave_closed_form()
        for s_value in (sp.Rational(1, 5), sp.Rational(9, 10), 3, 50, 10**4):
            a = complex(sp.N(integrated.subs(SAMPLE).subs(mandelstam_s, s_value), 20))
            b = complex(sp.N(closed.subs(SAMPLE).subs(mandelstam_s, s_value), 20))
            self.assertAlmostEqual(a.real, b.real, places=12)
            self.assertAlmostEqual(a.imag, 0, places=12)

    def test_pole_term_helper(self):
        mass2, t_minus, t_plus = sp.Rational(7, 3), -5, -sp.Rational(1, 2)
        t_of_cos = (t_plus + t_minus)/2 + (t_plus - t_minus)/2*cos_theta
        direct = sp.integrate(1/(mass2 - t_of_cos), (cos_theta, -1, 1))
        self.assertAlmostEqual(float(direct), float(model.pole_term_s_wave(mass2, t_minus, t_plus)), 14)

    def test_high_energy_limit(self):
        limit = sp.limit(model.alp_alp_s_wave_closed_form(), mandelstam_s, sp.oo)
        self.assertEqual(sp.simplify(limit + 3*lam/(16*sp.pi)), 0)

    def test_intermediate_energy(self):
        a0 = model.intermediate_energy_s_wave()
        self.assertEqual(sp.simplify(a0 - 5*mandelstam_s**2/(96*sp.pi*f_a**2*model.m_sigma_squared)), 0)
        self.assertEqual(sp.simplify(model.intermediate_energy_scale()
                                     - sp.Rational(48, 5)*sp.pi*f_a**2*model.m_sigma_squared), 0)

    def test_intermediate_expansion_of_exact_amplitude(self):
        """Exact amplitude at m_a -> 0 (B_0 -> 0) expanded in s/m_sigma^2."""
        epsilon = sp.Symbol('epsilon', positive=True)
        closed = model.closed_form_alp_alp_amplitude().subs(B_0, 0)
        scaled = closed.subs({mandelstam_s: epsilon*mandelstam_s, mandelstam_t: epsilon*mandelstam_t})
        series = sp.series(scaled, epsilon, 0, 3).removeO()
        u = -mandelstam_s - mandelstam_t
        expected = epsilon**2*(mandelstam_s**2 + mandelstam_t**2 + u**2)/(2*lam*f_a**4)
        self.assertEqual(sp.simplify(series - expected), 0)


class CoupledChannelTest(unittest.TestCase):
    def test_contact_matrix_and_bounds(self):
        matrix, sigma_alp = model.asymptotic_contact_matrix()
        self.assertEqual(sp.simplify(matrix - lam/(16*sp.pi)*sp.Matrix([[-3, -1], [-1, -3]])), sp.zeros(2))
        self.assertEqual(sp.simplify(sigma_alp + lam/(8*sp.pi)), 0)
        bounds = model.coupling_bounds()
        self.assertEqual(bounds['alp_alp'], (8*sp.pi/3, 16*sp.pi/3))
        self.assertEqual(bounds['coupled'], (2*sp.pi, 4*sp.pi))
        self.assertEqual(bounds['sigma_alp'], (4*sp.pi, 8*sp.pi))
        self.assertEqual(model.bound_on_radial_mass(2*sp.pi), 4*sp.pi*f_a**2)

    def test_exchange_parts_vanish_asymptotically(self):
        contact, sigma_alp_contact = model.asymptotic_contact_matrix()
        previous = None
        for s_value in (10**4, 10**6, 10**8):
            matrix, sigma_alp, exchanges = model.cartesian_channel_s_waves(sp.Integer(s_value))
            sizes = [abs(complex(sp.N(v.subs(SAMPLE)))) for v in exchanges.values()]
            bound = 50*float(sp.log(s_value))/s_value
            for size in sizes:
                self.assertLess(size, bound)
            if previous is not None:
                for now, before in zip(sizes, previous):
                    self.assertLess(now, before)
            previous = sizes
        difference = (matrix - contact).subs(SAMPLE)
        self.assertLess(max(abs(complex(sp.N(x))) for x in difference), 1e-7)
        self.assertLess(abs(complex(sp.N((sigma_alp - sigma_alp_contact).subs(SAMPLE)))), 1e-7)


class MesonTest(unittest.TestCase):
    def test_alp_meson_scattering_polar_equals_cartesian(self):
        polar, cartesian = model.alp_meson_scattering_polar_cartesian()
        self.assertEqual(sp.simplify(polar - cartesian), 0)
        self.assertEqual(sp.limit(polar, mandelstam_t, -sp.oo), 0)

    def test_three_meson_vertices_and_bounds(self):
        v333, v311, v322, bounds = model.alp_meson_contact_bounds()
        self.assertEqual(sp.simplify(v333 - B_0*m_u/(f_a*f_pi)), 0)
        self.assertEqual(sp.simplify(v311 - v333/3), 0)
        self.assertEqual(sp.simplify(v322 - v311), 0)
        self.assertEqual(sp.simplify(bounds['P3P3P3'][1] - 8*sp.sqrt(2)*sp.pi*f_pi**2/(B_0*m_u)), 0)
        self.assertEqual(sp.simplify(bounds['P3pipi'][1] - 24*sp.pi*f_pi**2/(B_0*m_u)), 0)
        number = float(bounds['P3P3P3'][1].subs(model.NUMERICAL_INPUTS))
        self.assertAlmostEqual(number, 52.52, delta=0.01)


if __name__ == '__main__':
    unittest.main()
