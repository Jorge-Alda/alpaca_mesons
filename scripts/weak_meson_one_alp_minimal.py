#!/usr/bin/env .venv/bin/python
"""Expand weak meson operators to one ALP and at most one meson.

Run: .venv/bin/python scripts/weak_meson_one_alp_minimal.py [--test]
Natural units. General ALP-quark coefficients, formal weak spurions w0,w1;
external sources are the ALP-induced sources of Alda:2026nrd. Flavor matrices
are abstract noncommuting symbols; each displayed product is inside a flavor
trace. The Lorentz index on the meson derivative is suppressed.
"""
import argparse
import unittest
import sympy as sp
from sympy.physics.quantum import TensorProduct

alp_order, meson_order = sp.symbols('a_order phi_order', commutative=True)
a_over_fa, partial_a_over_fa = sp.symbols('a_over_f_a partial_a_over_f_a', commutative=True)
f_pi = sp.symbols('f_pi', positive=True)
phi, partial_phi, w0, w1, chi0, cq_a, cq_v, cq_s, cq_p = sp.symbols(
    'Phi partial_Phi w_0 w_1 chi_0 c_q^A c_q^V c_q^S c_q^P', commutative=False)


def commutator(left, right):
    return left * right - right * left


def anticommutator(left, right):
    return left * right + right * left


def coefficient(expression, alp_power, meson_power):
    return sp.expand(expression).coeff(alp_order, alp_power).coeff(meson_order, meson_power).expand()


def calculate():
    identity = sp.S.One
    # Terms are sufficient for precisely one ALP and at most one meson.
    u = identity + sp.I * meson_order * phi / (2 * f_pi)
    u_dagger = identity - sp.I * meson_order * phi / (2 * f_pi)
    weak = w0 + alp_order * a_over_fa * w1
    u_mu = (-meson_order * partial_phi / f_pi
            + 2 * alp_order * partial_a_over_fa * cq_a
            - sp.I * alp_order * meson_order * partial_a_over_fa
            * commutator(phi, cq_v) / f_pi)
    chi_plus = (2 * chi0 + 2 * sp.I * alp_order * a_over_fa * commutator(cq_s, chi0)
                + alp_order * meson_order * a_over_fa
                * anticommutator(phi, anticommutator(cq_p, chi0)) / f_pi)
    dressed_weak = sp.expand(u * weak * u_dagger)
    dressed_current = sp.expand(u_dagger * u_mu * u)
    octet = coefficient(dressed_weak * u_mu * u_mu, 1, 1)
    mass_no_meson = coefficient(dressed_weak * chi_plus, 1, 0)
    mass_one_meson = coefficient(dressed_weak * chi_plus, 1, 1)
    current_no_alp_one_meson = coefficient(dressed_current, 0, 1)
    current_one_alp_no_meson = coefficient(dressed_current, 1, 0)
    # For the double-trace term, expand each trace argument independently.
    singlet_weak_alp = coefficient(dressed_weak * u_mu, 1, 0)
    singlet_weak_meson = coefficient(dressed_weak * u_mu, 0, 1)
    singlet_alp = coefficient(u_mu, 1, 0)
    singlet_meson = coefficient(u_mu, 0, 1)
    tensor_27 = (TensorProduct(current_one_alp_no_meson, current_no_alp_one_meson)
                 + TensorProduct(current_no_alp_one_meson, current_one_alp_no_meson))
    return locals()


def show(label, expression):
    print(f'{label}:\n{sp.latex(expression)}\n')


class WeakMesonExpansionTest(unittest.TestCase):
    def test_octet(self):
        result = calculate()
        expected = -2 * partial_a_over_fa / f_pi * w0 * anticommutator(cq_a, partial_phi)
        self.assertEqual(sp.expand(result['octet'] - expected), 0)

    def test_mass(self):
        result = calculate()
        tadpole = 2 * a_over_fa * (w1 * chi0 + sp.I * w0 * commutator(cq_s, chi0))
        mixed = a_over_fa / f_pi * (sp.I * commutator(phi, w1) * chi0
                 + w0 * anticommutator(phi, anticommutator(cq_p, chi0))
                 - commutator(phi, w0) * commutator(cq_s, chi0))
        self.assertEqual(sp.expand(result['mass_no_meson'] - tadpole), 0)
        self.assertEqual(sp.expand(result['mass_one_meson'] - mixed), 0)

    def test_singlet_and_27_inputs(self):
        result = calculate()
        self.assertEqual(result['singlet_weak_alp'], 2 * partial_a_over_fa * w0 * cq_a)
        self.assertEqual(result['singlet_weak_meson'], -w0 * partial_phi / f_pi)
        self.assertEqual(result['singlet_alp'], 2 * partial_a_over_fa * cq_a)
        self.assertEqual(result['singlet_meson'], -partial_phi / f_pi)
        self.assertEqual(result['current_no_alp_one_meson'], -partial_phi / f_pi)
        self.assertEqual(result['current_one_alp_no_meson'], 2 * partial_a_over_fa * cq_a)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test', action='store_true')
    args = parser.parse_args()
    if args.test:
        unittest.main(argv=[__file__])
        return
    result = calculate()
    print('Natural units; formal generic basis; one ALP and at most one meson. '
          'Multiply displayed flavor structures by f_pi^4 and their weak couplings; add h.c.')
    for key, label in [
        ('u_mu', 'Input u_mu'), ('chi_plus', 'Input chi_plus'),
        ('dressed_weak', 'Dressed weak matrix'), ('octet', 'Octet derivative structure inside trace'),
        ('singlet_weak_alp', 'Weak singlet first trace, one ALP'),
        ('singlet_meson', 'Weak singlet second trace, one meson'),
        ('singlet_weak_meson', 'Weak singlet first trace, one meson'),
        ('singlet_alp', 'Weak singlet second trace, one ALP'),
        ('mass_no_meson', 'Weak mass structure with no meson, inside trace'),
        ('mass_one_meson', 'Weak mass structure with one meson, inside trace'),
        ('current_no_alp_one_meson', 'Left current, one meson'),
        ('current_one_alp_no_meson', 'Left current, one ALP'),
        ('tensor_27', '27-plet tensor-product structure inside tensor trace'),
    ]:
        show(label, result[key])
    print('The ALP-dependent 27-plet coefficient multiplies two meson currents and starts at one ALP and two mesons.')


if __name__ == '__main__':
    main()
