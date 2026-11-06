#!/usr/bin/env python3
"""Derive rho--isoscalar pole and two-pion cut comparison identities.

Invocation: .venv/bin/python scripts/derive_alp_tff_rho_isoscalar_overlap.py
Unit tests: .venv/bin/python scripts/derive_alp_tff_rho_isoscalar_overlap.py --test

Natural units; squared invariants are in GeV^2 and the scalar TFF is in
GeV^-1. Use disc f = f(s+i0)-f(s-i0). The isoscalar vector is a stable
narrow pole in this calculation; the rho response has an arbitrary analytic
self-energy with upper and lower boundary values. Physical two-pion cut
inputs and their phases are not equated with the model parameters.
The rho cut projection uses an N-subtracted dispersion integral. With a
finite upper limit it reproduces the rho discontinuity only inside the
integrated interval. A residual above the cutoff, additional channels or
isolated poles cannot be absorbed into a subtraction polynomial without
additional assumptions. The full-support polynomial remainder requires a
dispersion representation with no other singularities and adequate
subtractions.

For the vector flavor matrices, lambda_i = sqrt(2) T_i, where the T_i
in the source TFF record obey tr(T_i T_j) = delta_ij. Thus the VVP trace
with two lambda matrices carries a factor 1/2 relative to the same trace
written with T_i. The photon residues R_i retain their original definition
tr(T_i S_gamma) = tr(lambda_i S_gamma)/sqrt(2).
"""

import sys
import unittest

import sympy as sp


def expressions():
    """Build ordered mixed terms and their analytic decompositions."""
    s1, s2, mass_I_squared = sp.symbols('s_1 s_2 M_I^2')
    R_I = sp.Symbol('R_I')
    d1, d2, d3, mass_a_squared = sp.symbols('d_1 d_2 d_3 m_a^2')
    lambda_rho, lambda_I, U_a, Y_a = (sp.MatrixSymbol(name, 3, 3) for name in
                                      ('lambda_rho', 'lambda_I', 'U_a', 'Y_a'))
    # Each lambda matrix is sqrt(2) times the unit-normalized T matrix.
    numerator = 2 * (sp.Trace(lambda_rho*lambda_I*U_a)
                     *((d1-d3)*mass_a_squared + d3*(s1+s2))
                     - 4*d2*sp.Trace(lambda_rho*lambda_I*Y_a))
    numerator_at_pole = numerator.subs(s2, mass_I_squared)
    pole_quotient = numerator/(mass_I_squared-s2)
    pole_part = numerator_at_pole/(mass_I_squared-s2)
    regular_part = sp.simplify(pole_quotient-pole_part)

    W_rho = sp.Function(r'W_\rho^{\mathrm{dress}}')
    # The indexed neutral generators commute. The numerator is symmetric in
    # the invariants and in the neutral-vector labels under that assumption.
    ordered_first = -2*R_I*W_rho(s1)*pole_quotient
    ordered_second = ordered_first.xreplace({s1: s2, s2: s1})
    mixed_tff_contribution = ordered_first+ordered_second

    # Compute the double subtraction directly from the two ordered terms.
    double_subtracted = sp.expand(
        mixed_tff_contribution - mixed_tff_contribution.subs(s1, 0)
        - mixed_tff_contribution.subs(s2, 0)
        + mixed_tff_contribution.subs({s1: 0, s2: 0}))
    slope = sp.diff(numerator, s1)
    constant = numerator.subs({s1: 0, s2: 0})
    # W_I is the selected narrow pole and its difference is kept explicit.
    pole = R_I/(mass_I_squared-s2)
    pole_at_zero = pole.subs(s2, 0)
    # An independently expanded expression for the first ordering.
    first_double = sp.expand(ordered_first - ordered_first.subs(s1, 0)
                             - ordered_first.subs(s2, 0)
                             + ordered_first.subs({s1: 0, s2: 0}))
    first_factored = -2*(W_rho(s1)-W_rho(0))*(pole-pole_at_zero)*(
        constant+slope*(s1+mass_I_squared))
    first_factored += -2*W_rho(0)*slope*s1*(pole-pole_at_zero)
    # The term regular in s2 cancels under the s2 subtraction.

    return locals()


def boundary_expressions():
    """Calculate jumps from boundary values with independent channels."""
    mass_rho_squared, s, self_real, self_imaginary, R_rho = sp.symbols(
        'M_0rho^2 s Re_Delta_rho Im_Delta_rho R_rho', real=True)
    rho_upper = R_rho/(mass_rho_squared-s-self_real-sp.I*self_imaginary)
    rho_lower = R_rho/(mass_rho_squared-s-self_real+sp.I*self_imaginary)
    rho_jump = sp.factor(rho_upper-rho_lower)
    rho_jump_1, rho_jump_2, isoscalar_jump_1, isoscalar_jump_2 = sp.symbols(
        'disc_W_rho_1 disc_W_rho_2 disc_W_I_1 disc_W_I_2')
    rho_lower_1, rho_lower_2, isoscalar_lower_1, isoscalar_lower_2 = sp.symbols(
        'W_rho_1_minus W_rho_2_minus W_I_1_minus W_I_2_minus')
    P = sp.Symbol('P_rhoI')
    mixed_boundary = -2*P*(rho_lower_1*isoscalar_lower_2
                            + isoscalar_lower_1*rho_lower_2)
    upper_1 = {rho_lower_1: rho_lower_1+rho_jump_1,
               isoscalar_lower_1: isoscalar_lower_1+isoscalar_jump_1}
    upper_2 = {rho_lower_2: rho_lower_2+rho_jump_2,
               isoscalar_lower_2: isoscalar_lower_2+isoscalar_jump_2}
    double_jump = sp.expand(mixed_boundary.subs(upper_1).subs(upper_2)
                            - mixed_boundary.subs(upper_1)
                            - mixed_boundary.subs(upper_2)+mixed_boundary)
    return locals()


def rho_cut_replacement_expressions():
    """Project both rho placements onto their subtracted two-pion cut."""
    values = expressions()
    s1, s2 = values['s1'], values['s2']
    spectral_variable = sp.Symbol('x', positive=True)
    threshold = sp.Symbol('s_threshold', positive=True)
    cutoff_squared = sp.Symbol('Lambda^2', positive=True)
    subtraction_count = sp.Symbol('N', positive=True, integer=True)
    subtraction_index = sp.Symbol('k', integer=True, nonnegative=True)
    subtraction_coefficient = sp.Function('c_rho')
    rho_discontinuity = sp.Function('disc_W_rho')
    rho_cut_response = sp.Function(r'W_\rho^{\pi\pi,[N]}')
    cut_variable = sp.Symbol('s', positive=True)
    dispersion_integral = s1**subtraction_count/(2*sp.pi*sp.I)*sp.Integral(
        rho_discontinuity(spectral_variable)/(
            spectral_variable**subtraction_count*(spectral_variable-s1)),
        (spectral_variable, threshold, cutoff_squared))
    # Applying Disc_s[1/(x-s)] = 2*pi*i*delta(x-s) to the dispersion
    # integral gives this jump. The subtraction powers cancel on the cut.
    dispersive_jump = sp.simplify(cut_variable**subtraction_count*sp.integrate(
        rho_discontinuity(spectral_variable)*sp.DiracDelta(spectral_variable-cut_variable)
        /spectral_variable**subtraction_count,
        (spectral_variable, threshold, cutoff_squared)))
    subtraction_polynomial_1 = sp.Sum(
        subtraction_coefficient(subtraction_index)*s1**subtraction_index,
        (subtraction_index, 0, subtraction_count-1))
    subtraction_polynomial_2 = subtraction_polynomial_1.subs(s1, s2)

    rho_response = values['W_rho']
    cut_tff_contribution = values['mixed_tff_contribution'].xreplace({
        rho_response(s1): rho_cut_response(s1),
        rho_response(s2): rho_cut_response(s2),
    })
    residual = values['mixed_tff_contribution']-cut_tff_contribution
    residual_with_subtraction_polynomial = residual.xreplace({
        rho_response(s1): rho_cut_response(s1)+subtraction_polynomial_1,
        rho_response(s2): rho_cut_response(s2)+subtraction_polynomial_2,
    }).expand()

    first_cut_regular = -2*values['R_I']*rho_cut_response(s1)*values['regular_part']
    first_cut_pole = -2*values['R_I']*rho_cut_response(s1)*values['pole_part']
    first_cut_ordering = values['ordered_first'].xreplace({
        rho_response(s1): rho_cut_response(s1)})

    # The jump is projected on the rho response in each ordering; an
    # independent discontinuity of the isoscalar pole is a separate term.
    rho_jump_1, rho_jump_2 = sp.symbols('disc_W_rho_1 disc_W_rho_2')
    support_1, support_2 = sp.symbols('cut_support_1 cut_support_2')
    first_rho_jump = -2*values['R_I']*values['pole_quotient']*rho_jump_1
    second_rho_jump = (-2*values['R_I']*values['pole_quotient']
                       .xreplace({s1: s2, s2: s1})*rho_jump_2)
    first_cut_jump = first_rho_jump*support_1
    second_cut_jump = second_rho_jump*support_2
    first_residual_jump = first_rho_jump-first_cut_jump
    second_residual_jump = second_rho_jump-second_cut_jump
    return locals()


class DerivationTests(unittest.TestCase):
    def test_affine_pole_and_regular_term(self):
        values = expressions()
        self.assertEqual(sp.simplify(values['regular_part']
                                   + sp.diff(values['numerator'], values['s2'])), 0)

    def test_rho_jump_sign(self):
        v = boundary_expressions()
        expected = 2*sp.I*v['R_rho']*v['self_imaginary']/(
            (v['mass_rho_squared']-v['s']-v['self_real'])**2+v['self_imaginary']**2)
        self.assertEqual(sp.simplify(v['rho_jump']-expected), 0)

    def test_both_ordered_double_jumps(self):
        v = boundary_expressions()
        expected = -2*v['P']*(v['rho_jump_1']*v['isoscalar_jump_2']
                               +v['isoscalar_jump_1']*v['rho_jump_2'])
        self.assertEqual(sp.expand(v['double_jump']-expected), 0)

    def test_double_subtraction(self):
        v = expressions()
        self.assertEqual(sp.simplify(v['first_double']-v['first_factored']), 0)
        self.assertEqual(v['double_subtracted'].subs(v['s1'], 0), 0)
        self.assertEqual(v['double_subtracted'].subs(v['s2'], 0), 0)

    def test_equal_pole_residue_does_not_fix_single_cut(self):
        s2, mass = sp.symbols('s_2 M_I^2')
        residue, regular_first, regular_second = sp.symbols('residue regular_1 regular_2')
        first = residue/(mass-s2)+regular_first
        second = residue/(mass-s2)+regular_second
        self.assertEqual(sp.limit((mass-s2)*(first-second), s2, mass), 0)
        self.assertEqual(sp.simplify(first-second), regular_first-regular_second)

    def test_vector_generator_and_photon_residue_normalization(self):
        light, strange, axial_up, axial_down, axial_strange = sp.symbols(
            'light strange axial_up axial_down axial_strange', real=True)
        normalized_rho = sp.diag(1, -1, 0)/sp.sqrt(2)
        normalized_isoscalar = sp.diag(light, light, strange)
        rho_flavor_matrix = sp.sqrt(2)*normalized_rho
        isoscalar_flavor_matrix = sp.sqrt(2)*normalized_isoscalar
        axial_source = sp.diag(axial_up, axial_down, axial_strange)
        photon_source = sp.diag(axial_up, axial_down, axial_strange)
        self.assertEqual(sp.simplify(
            4*sp.trace(normalized_rho*normalized_isoscalar*axial_source)
            - 2*sp.trace(rho_flavor_matrix*isoscalar_flavor_matrix*axial_source)), 0)
        self.assertEqual(sp.simplify(
            sp.trace(normalized_isoscalar*photon_source)
            - sp.trace(isoscalar_flavor_matrix*photon_source)/sp.sqrt(2)), 0)

    def test_rho_cut_replacement_and_full_support_polynomial(self):
        v = rho_cut_replacement_expressions()
        original = v['values']
        expected_jump = v['rho_discontinuity'](v['cut_variable'])*(
            sp.Heaviside(v['cutoff_squared']-v['cut_variable'])
            - sp.Heaviside(v['threshold']-v['cut_variable']))
        self.assertEqual(sp.simplify(v['dispersive_jump']-expected_jump), 0)
        expected_residual = -2*original['R_I']*(
            original['pole_quotient']*v['subtraction_polynomial_1']
            + original['pole_quotient'].xreplace({v['s1']: v['s2'], v['s2']: v['s1']})
            *v['subtraction_polynomial_2'])
        self.assertEqual(sp.simplify(
            v['residual_with_subtraction_polynomial']-expected_residual), 0)
        self.assertEqual(sp.simplify(
            v['first_residual_jump'].subs(v['support_1'], 1)), 0)
        self.assertEqual(sp.simplify(
            v['second_residual_jump'].subs(v['support_2'], 1)), 0)

    def test_rho_cut_projection_keeps_regular_single_cut(self):
        v = rho_cut_replacement_expressions()
        self.assertEqual(sp.simplify(
            v['first_cut_ordering']-v['first_cut_pole']-v['first_cut_regular']), 0)
        self.assertNotEqual(sp.simplify(v['first_cut_regular']), 0)
        self.assertEqual(sp.simplify(
            v['first_residual_jump']-(1-v['support_1'])*v['first_rho_jump']), 0)
        self.assertEqual(sp.simplify(
            v['second_residual_jump']-(1-v['support_2'])*v['second_rho_jump']), 0)


def show(label, expression):
    print(f'{label}:\n{sp.latex(expression)}\n')


def main():
    v = expressions()
    b = boundary_expressions()
    cut = rho_cut_replacement_expressions()
    print(r'Vector flavor matrices: \lambda_i=\sqrt{2}\,T_i, '
          r'\quad R_i=\operatorname{tr}(\lambda_i S_\gamma)/\sqrt{2}.')
    print(r'Alda:2026nrd source matrices: '
          r'U_a=-2\mathbb{M}_a^{-1}(\mathcal C)/f_a, '
          r'\quad Y_a=2\{c_q^{AP}+\mathbb{M}_a^{-1}(\mathcal C),\chi\}/f_a.')
    print(r'Dressed rho response: '
          r'W_\rho^{\mathrm{dress}}(s)=R_\rho/[M_{0\rho}^2-s-\Delta_\rho^{\mathrm{self}}(s)].')
    print(r'Narrow isoscalar pole: '
          r'W_I^{\mathrm{pole}}(s)=R_I/(M_I^2-s-i0).')
    print()
    show('Indexed mixed-vector numerator', v['numerator'])
    show('Numerator evaluated at the isoscalar pole', v['numerator_at_pole'])
    show('Regular part after extracting the isoscalar pole', v['regular_part'])
    show(r'Mixed TFF contribution F^{(\rho I),\mathrm{dress}}, both photon orderings',
         v['mixed_tff_contribution'])
    show('Rho boundary discontinuity', b['rho_jump'])
    show('Double discontinuity with both orderings', b['double_jump'])
    show('First ordering after one subtraction in each invariant', v['first_factored'])
    show('Residual of the direct double-subtraction identity',
         sp.simplify(v['first_double']-v['first_factored']))
    show('N-subtracted rho two-pion cut response with finite upper cutoff',
         cut['dispersion_integral'])
    show('Discontinuity of that response from the delta-function identity', cut['dispersive_jump'])
    show('Mixed TFF rho-cut projection, both photon orderings', cut['cut_tff_contribution'])
    show('First ordering: rho-cut projection regular at the isoscalar pole',
         cut['first_cut_regular'])
    show('Residual when the full-support rho response differs by its subtraction polynomial',
         cut['residual_with_subtraction_polynomial'])
    show('First ordered rho-cut discontinuity left in the residual', cut['first_residual_jump'])
    show('Second ordered rho-cut discontinuity left in the residual', cut['second_residual_jump'])
    print('Each cut_support_i is 1 on the integrated rho cut and 0 outside it.')
    print('At a finite cutoff the omitted rho cut above the cutoff remains in the residual.')
    print('A subtraction-polynomial-only residual requires full support and no other singularities.')
    print()

    # The physical two-pion discontinuity is independent input. C_I is the
    # pole coefficient in f_1, before conversion to Fcal = f_1/e.
    x, y, mass, charge = sp.symbols('x y M_I^2 e')
    sigma = sp.Function('sigma_pi')(x)
    pion_form_factor = sp.Symbol('conjugate_F_pi_V')
    C_I = sp.Function('C_I')
    physical_pole_coefficient = (2*sp.I*x*sigma**3/(96*sp.pi)
                                 *pion_form_factor*C_I(x, mass)/charge)
    physical_regular = sp.Function('disc_F_pipi_regular')(x, y)
    physical_single = physical_pole_coefficient/(mass-y)+physical_regular
    mixed_rho_jump = sp.Symbol('disc_W_rho')
    mixed_pole_coefficient = -2*v['R_I']*mixed_rho_jump*v['numerator_at_pole'].subs(
        {v['s1']: x, v['mass_I_squared']: mass})
    mixed_regular = -2*v['R_I']*mixed_rho_jump*v['regular_part'].subs(
        {v['s1']: x, v['s2']: y, v['mass_I_squared']: mass})
    show('Physical pion-cut pole coefficient from the production wave', physical_pole_coefficient)
    show('Mixed-model pion-cut pole coefficient', mixed_pole_coefficient)
    show('Physical minus model pion discontinuity, one ordering',
         sp.expand(physical_single - mixed_pole_coefficient/(mass-y)-mixed_regular))
    show('Narrow-pole jump factor', 2*sp.pi*sp.I*v['R_I']*sp.DiracDelta(mass-y))
    rho_jump_1 = sp.Symbol('disc_W_rho_1')
    rho_jump_2 = sp.Symbol('disc_W_rho_2')
    model_double_jump = -4*sp.pi*sp.I*v['R_I']*(
        rho_jump_1*sp.DiracDelta(v['mass_I_squared']-v['s2'])*v['numerator_at_pole']
        + rho_jump_2*sp.DiracDelta(v['mass_I_squared']-v['s1'])
        * v['numerator_at_pole'].subs(v['s1'], v['s2']))
    show('Mixed double discontinuity, both ordered narrow poles', model_double_jump)
    physical_double_jump = 2*sp.pi*sp.I*(
        physical_pole_coefficient.subs(x, v['s1'])
        *sp.DiracDelta(mass-v['s2'])
        +physical_pole_coefficient.subs(x, v['s2'])
        *sp.DiracDelta(mass-v['s1']))
    show('Physical production-pole double discontinuity, both orderings', physical_double_jump)
    print('The discontinuity comparison leaves physical residues, regular terms,')
    print('subtraction functions and relative phases independent until matched.')
    print('All expressions are analytic; no numerical input has been selected.')


if __name__ == '__main__':
    if '--test' in sys.argv:
        sys.argv.remove('--test')
        unittest.main()
    else:
        main()
