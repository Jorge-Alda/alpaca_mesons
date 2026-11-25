"""Symbolic external-line loop vertex projections.

Run: PYTHONPATH=scripts .venv/bin/python scripts/derive_alp_loop_external_contractions.py
All mass parameters have GeV units, and all momenta have GeV units.
"""
import unittest

import sympy as sp
from derive_alp_dressed_vertex_power_counting import (
    anticommutator, closed_alp_line_potential_matrix_form,
    contract_potential_vertex, kinetic_closed_line_form,
)
from derive_alp_three_meson_loop_assembly import kinetic_vertex


def calculate():
    epsilon, decay_constant, low_energy_constant, anomaly_mass_squared = sp.symbols(
        r'\epsilon f_\pi B_0 m_0^2', real=True, nonzero=True)
    quark_u, quark_d, quark_s, scalar_ds, pseudoscalar_ds = sp.symbols(
        'm_u m_d m_s S_ds P_ds', real=True)
    scalar_second, pseudoscalar_second, theta_second = sp.symbols('S_2 P_2 Theta_2', real=True)
    overlap_u, overlap_d, overlap_s, overlap_ds = sp.symbols('z_u z_d z_s z_ds', real=True)
    quark_mass = sp.diag(quark_u, quark_d, quark_s)
    overlap = sp.Matrix([[overlap_u, 0, 0], [0, overlap_d, overlap_ds], [0, overlap_ds, overlap_s]])
    scalar_first = sp.Matrix([[0, 0, 0], [0, 0, scalar_ds], [0, scalar_ds, 0]])
    pseudoscalar_first = sp.Matrix([[0, 0, 0], [0, 0, pseudoscalar_ds], [0, pseudoscalar_ds, 0]])
    scalar_sources = [quark_mass, scalar_first, scalar_second * sp.eye(3)]
    pseudoscalar_sources = [sp.zeros(3), pseudoscalar_first, pseudoscalar_second * sp.eye(3)]
    invariants = {'f_pi': decay_constant, 'B_0': low_energy_constant,
                  'anomaly_mass_squared': anomaly_mass_squared,
                  'S': scalar_sources, 'P': pseudoscalar_sources,
                  'theta': [0, 0, theta_second]}
    meson_first = sp.diag(1, -1, 0)
    meson_second = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
    pole = (2 * epsilon * overlap, sp.S.One)
    first_leg, second_leg = (meson_first, 0), (meson_second, 0)
    first_cubic = contract_potential_vertex((first_leg, second_leg, pole), invariants, epsilon)
    expected_first = -epsilon * low_energy_constant / (2 * decay_constant) * (
        scalar_first * anticommutator(meson_first, meson_second)).trace()
    closed_cubic = contract_potential_vertex((first_leg, pole, pole), invariants, epsilon)
    expected_closed = epsilon**2 / decay_constant * (
        meson_first * (low_energy_constant * (pseudoscalar_sources[2]
                       - 2 * anticommutator(overlap, scalar_first))
                       - anomaly_mass_squared * theta_second * sp.eye(3) / 3)).trace()
    closed_quartic = contract_potential_vertex((first_leg, second_leg, pole, pole), invariants, epsilon)
    expected_quartic = epsilon**2 * closed_alp_line_potential_matrix_form(
        meson_first, meson_second, invariants, overlap)
    momentum_p, momentum_k = sp.symbols('p k', real=True)
    kinetic = kinetic_vertex(((meson_first, momentum_p), (meson_second, -momentum_p),
                              (2 * epsilon * overlap, momentum_k),
                              (2 * epsilon * overlap, -momentum_k)),
                             lambda left, right: left * right, decay_constant)
    expected_kinetic = epsilon**2 * kinetic_closed_line_form(
        meson_first, meson_second, invariants, overlap).subs({'p': momentum_p, 'k': momentum_k})
    return tuple(sp.simplify(value) for value in (
        first_cubic - expected_first, closed_cubic - expected_closed,
        closed_quartic - expected_quartic, kinetic - expected_kinetic))


class ExternalContractionTest(unittest.TestCase):
    def test_projected_vertices(self):
        self.assertEqual(calculate(), (0, 0, 0, 0))


if __name__ == '__main__':
    for value in calculate():
        print(sp.latex(value))
