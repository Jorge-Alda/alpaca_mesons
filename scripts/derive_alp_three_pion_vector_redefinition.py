#!/usr/bin/env python3
"""Check an on-shell finite-mixing a -> pi0 pi+ pi- amplitude in two coordinates.

Run: .venv/bin/python scripts/derive_alp_three_pion_vector_redefinition.py
LO chiral theory, natural units, tree level; exact rational masses and mixing.
"""
import sympy as sp
from derive_alp_finite_three_pion import (
    f_pi, quartic_kinetic_vertex, quartic_mass_vertex,
)


def derive():
    """Evaluate the LO vertex and the generated contact and charged exchanges."""
    # Masses and pole columns are in units where f_pi = f_a = B_0 = 1.
    # The charged pion squared mass 2 B_0 m_q follows from the neutral
    # pi0 diagonal mass entry, not from an independently chosen spectrum.
    alp_mass_squared = sp.Rational(36)
    neutral_pion_mass_squared = sp.S.One
    alp_component = sp.Rational(12, 13)
    neutral_meson_component = sp.Rational(5, 13)
    pion_mass_squared = (neutral_meson_component**2*alp_mass_squared
                         +alp_component**2*neutral_pion_mass_squared)
    quark_mass = pion_mass_squared/2
    mass_mixing = (alp_mass_squared-neutral_pion_mass_squared)*alp_component*neutral_meson_component
    pseudoscalar_coupling = -mass_mixing/(4*quark_mass)
    assert -4*quark_mass*pseudoscalar_coupling == mass_mixing
    neutral_mass_matrix = sp.Matrix([
        [alp_component**2*alp_mass_squared+neutral_meson_component**2*neutral_pion_mass_squared,
         mass_mixing], [mass_mixing,pion_mass_squared]])
    alp_pole = sp.Matrix([alp_component,neutral_meson_component])
    pion_pole = sp.Matrix([-neutral_meson_component,alp_component])
    assert neutral_mass_matrix*alp_pole == alp_mass_squared*alp_pole
    assert neutral_mass_matrix*pion_pole == neutral_pion_mass_squared*pion_pole
    assert 6 > 1+2*sp.sqrt(pion_mass_squared)

    def matrix_unit(row,column):
        return sp.eye(3)[:,row]*sp.eye(3)[column,:]
    pion_neutral = sp.diag(1,-1,0)
    pion_plus = sp.sqrt(2)*matrix_unit(0,1)
    pion_minus = pion_plus.H
    zero = sp.zeros(3)
    poles = [(neutral_meson_component*pion_neutral,alp_component),
             (alp_component*pion_neutral,-neutral_meson_component),
             (pion_plus,0),(pion_minus,0)]
    chi_plus = sp.diag(4*quark_mass,4*quark_mass,4*quark_mass)
    chi_minus_first = 8*sp.I*quark_mass*pseudoscalar_coupling*pion_neutral

    # Physical Dalitz point: s>4 m_pi^2, s<(m_a-m_pi0)^2 and t=u.
    mandelstam_s = sp.Rational(124,5)
    mandelstam_t = (alp_mass_squared+neutral_pion_mass_squared+2*pion_mass_squared-mandelstam_s)/2
    mandelstam_u = mandelstam_t
    assert 4*pion_mass_squared < mandelstam_s < (6-1)**2
    scalar_products = {
        (0,1):(mandelstam_s-alp_mass_squared-neutral_pion_mass_squared)/2,
        (0,2):(mandelstam_t-alp_mass_squared-pion_mass_squared)/2,
        (0,3):(mandelstam_u-alp_mass_squared-pion_mass_squared)/2,
        (1,2):(mandelstam_u-neutral_pion_mass_squared-pion_mass_squared)/2,
        (1,3):(mandelstam_t-neutral_pion_mass_squared-pion_mass_squared)/2,
        (2,3):(mandelstam_s-2*pion_mass_squared)/2,
    }
    mass_vertex = quartic_mass_vertex(poles,chi_plus,chi_minus_first,zero,zero,zero).doit().subs(f_pi,1)
    kinetic_vertex = quartic_kinetic_vertex(poles,scalar_products).doit().subs(f_pi,1)
    original = sp.factor(mass_vertex+kinetic_vertex)

    # g_R=g_L=exp(i v a lambda_3/2) rotates pi+ by exp(i v a).
    # The old ALP field is the given combination of both physical neutral poles.
    rotation = sp.Symbol('v', real=True, nonzero=True)
    first_alp_component = alp_component
    second_alp_component = -neutral_meson_component
    # Extract cubic and quartic terms from |d pi + i v pi d a|^2.
    derivative_alp, derivative_pion0, pion, pion_conjugate, derivative_pion, derivative_pion_conjugate = sp.symbols(
        'dA dn pi_plus pi_minus dpi_plus dpi_minus')
    derivative_old_alp = first_alp_component*derivative_alp+second_alp_component*derivative_pion0
    charged_kinetic = sp.expand((derivative_pion+sp.I*rotation*pion*derivative_old_alp)
                                *(derivative_pion_conjugate-sp.I*rotation*pion_conjugate*derivative_old_alp))
    quartic_coefficient = sp.diff(charged_kinetic,derivative_alp,derivative_pion0,pion,pion_conjugate)
    contact = -quartic_coefficient*scalar_products[(0,1)]
    def exchange(invariant):
        """Differentiate both cubic vertices, then contract the internal pion."""
        first_external_dot = (invariant-alp_mass_squared-pion_mass_squared)/2
        second_external_dot = (invariant-neutral_pion_mass_squared-pion_mass_squared)/2
        first_internal_dot = -alp_mass_squared-first_external_dot
        second_internal_dot = -neutral_pion_mass_squared-second_external_dot
        first_vertex = (-sp.diff(charged_kinetic,derivative_alp,pion,derivative_pion_conjugate)*first_internal_dot
                        -sp.diff(charged_kinetic,derivative_alp,derivative_pion,pion_conjugate)*first_external_dot)
        second_vertex = (-sp.diff(charged_kinetic,derivative_pion0,pion,derivative_pion_conjugate)*second_external_dot
                         -sp.diff(charged_kinetic,derivative_pion0,derivative_pion,pion_conjugate)*second_internal_dot)
        return sp.factor(first_vertex*second_vertex/(pion_mass_squared-invariant))
    exchange_t = exchange(mandelstam_t)
    exchange_u = exchange(mandelstam_u)
    transformed = sp.factor(original+contact+exchange_t+exchange_u)
    assert original != 0 and contact != 0 and exchange_t != 0 and exchange_u != 0
    assert sp.simplify(contact+exchange_t+exchange_u) == 0
    assert sp.simplify(transformed-original) == 0
    assert sp.simplify(original+contact-original) != 0
    return original, contact, exchange_t, exchange_u, transformed


def main():
    original, contact, exchange_t, exchange_u, transformed = derive()
    print('Tree-level charge-conserving U(3) example; exact finite neutral mixing.')
    print('Original physical a -> pi0 pi+ pi- amplitude:',sp.latex(original))
    print('Generated quartic contribution:',sp.latex(contact))
    print('Charged pion exchange t:',sp.latex(exchange_t))
    print('Charged pion exchange u:',sp.latex(exchange_u))
    print('Transformed physical amplitude:',sp.latex(transformed))


if __name__ == '__main__':
    main()
