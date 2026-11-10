#!/usr/bin/env python3
"""Derive mixed UV-pole flavour coefficients and an ALP-dependent field shift.

Run: .venv/bin/python scripts/derive_mixed_field_redefinition.py
Units: GeV for masses and decay constants; s is in GeV^2.
The meson shift is local and linear in the ALP; terms quadratic in it are omitted.
"""
import unittest
import sympy as sp
from derive_nlo_loop_self_energies import (
    NONET_LABELS, GENERATORS, MANDELSTAM_S, leading_mass_matrix,
    mixed_tadpole_coefficients, mixed_local_pole_coefficients,
)


def mixed_poles():
    """Matrix coefficients inside tr(c_q^A ...), tr(c_q^P ...)."""
    axial_loop, pseudoscalar_loop = mixed_tadpole_coefficients(2 * leading_mass_matrix(NONET_LABELS))
    axial_local, pseudoscalar_local = mixed_local_pole_coefficients()
    return ([sp.simplify(a+b) for a, b in zip(axial_loop, axial_local)],
            [sp.simplify(a+b) for a, b in zip(pseudoscalar_loop, pseudoscalar_local)])


def field_shift():
    """Expand the quadratic form for old pi = new pi + epsilon a."""
    inverse_meson, mixed, inverse_alp, epsilon, order = sp.symbols('Q_ii Q_i_a Q_a_a epsilon t', real=True)
    inverse = sp.Matrix([[inverse_meson, mixed], [mixed, inverse_alp]])
    transform = sp.Matrix([[1, order * epsilon], [0, 1]])
    transformed = (transform.T * inverse * transform).applyfunc(lambda x: sp.expand(x).series(order, 0, 2).removeO().subs(order, 1))
    return inverse, transform, transformed


class Checks(unittest.TestCase):
    def test_shift(self):
        inverse, _, transformed = field_shift()
        epsilon = sp.Symbol('epsilon', real=True)
        self.assertEqual(sp.expand(transformed[0, 1] - inverse[0, 1]), epsilon * inverse[0, 0])

    def test_pseudoscalar_axial_relation(self):
        from derive_nlo_loop_self_energies import chi_zero
        axial, pseudoscalar = mixed_poles()
        for index in NONET_LABELS:
            self.assertEqual((pseudoscalar[index] + (chi_zero()*axial[index] + axial[index]*chi_zero())/2).applyfunc(sp.simplify), sp.zeros(3))

    def test_pion_axial(self):
        from derive_nlo_loop_self_energies import B_ZERO, F_PI, F_A, MASS_UP, MASS_DOWN, MASS_STRANGE
        axial, _ = mixed_poles()
        result = sp.simplify(sp.trace(GENERATORS[3] * axial[3]) * MANDELSTAM_S)
        expected = -2 * B_ZERO * (5*MASS_UP+5*MASS_DOWN+2*MASS_STRANGE) * MANDELSTAM_S / (3*F_A*F_PI)
        self.assertEqual(sp.simplify(result - expected), 0)


if __name__ == '__main__':
    print('Coefficients inside the traces with c_q^A and c_q^P; the axial coefficient is multiplied by s.')
    axial, pseudoscalar = mixed_poles()
    for label in NONET_LABELS:
        print('external generator', label, 'axial:', sp.latex(axial[label]), 'pseudoscalar:', sp.latex(pseudoscalar[label]))
    inverse, transformation, transformed = field_shift()
    for name, expression in [('original inverse propagator', inverse), ('field transformation', transformation), ('transformed inverse propagator', transformed)]:
        print(name + ':', sp.latex(expression))
    unittest.main(argv=[''])
