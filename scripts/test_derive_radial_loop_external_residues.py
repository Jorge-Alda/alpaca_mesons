"""Exact checks of the external-leg self-energy slopes against the Lagrangian expansion.

Run:  PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_radial_loop_external_residues
"""
import unittest

import sympy as sp

import derive_radial_loop_external_residues as calculation
import derive_radial_loop_trace_contractions as common
import expand_radial_loop_lagrangian as lagrangian
from group_radial_loop_trace_permutations import Commutator, polynomial_difference, trace_polynomial

LEG_E = (calculation.lambda_e.name, calculation.xi_e, 0, 'e')
LEG_EBAR = (calculation.lambda_ebar.name, calculation.xi_ebar, 0, 'ebar')


def lagrangian_slope(loop_legs, scale):
    vertex = lagrangian.vertex_polynomial((LEG_E, LEG_EBAR) + loop_legs, scale, calculation.SELF_ENERGY_DOTS)
    return {key: sp.diff(value, calculation.p_squared) for key, value in vertex.items()}


class ResidueSlopeChecks(unittest.TestCase):
    def test_alp_line_slope_equals_lagrangian_expansion(self):
        expected = lagrangian_slope((lagrangian.INTERNAL_ALP,) * 2, common.ma**2)
        self.assertEqual(polynomial_difference(expected, trace_polynomial(calculation.alp_line_slope)), {})

    def test_radial_line_slope_equals_lagrangian_expansion(self):
        expected = lagrangian_slope(((None, 0, 1, 'loop'),) * 2, common.ms**2)
        self.assertEqual(polynomial_difference(expected, trace_polynomial(calculation.radial_line_slope)), {})

    def test_record_closed_forms(self):
        alp_display = common.tr(Commutator(common.la, calculation.lambda_e)
                                * Commutator(common.la, calculation.lambda_ebar)) / (12 * common.fpi**2)
        self.assertEqual(polynomial_difference(trace_polynomial(calculation.alp_line_slope),
                                               trace_polynomial(alp_display)), {})
        self.assertEqual(sp.expand(calculation.radial_line_slope
                                   - 2 * calculation.xi_e * calculation.xi_ebar / common.fa**2), 0)

    def test_alp_line_slope_vanishes_on_the_alp_leg(self):
        """For e = ebar = a the commutators [lambda_a, lambda_a] vanish."""
        on_alp_leg = calculation.alp_line_slope.xreplace({calculation.lambda_e: common.la,
                                                          calculation.lambda_ebar: common.la})
        self.assertEqual(on_alp_leg, 0)


if __name__ == '__main__':
    unittest.main()
