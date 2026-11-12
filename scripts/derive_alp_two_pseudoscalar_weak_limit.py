#!/usr/bin/env python3
"""Derive the leading dressed cubic vertex, natural units.

Run: .venv/bin/python scripts/derive_alp_two_pseudoscalar_weak_limit.py
The external meson matrices denote normalized unmixed eigenmatrices.
"""
import sympy as sp
from derive_alp_flavour_violating_two_pseudoscalars import (
    eom_identity_left, eom_identity_right,
)


def derive():
    epsilon, f_pi, B_0, m_0, N_f = sp.symbols('epsilon f_pi B_0 m_0 N_f', nonzero=True, real=True)
    theta_1, theta_2 = sp.symbols('hat_Theta_1 hat_Theta_2', real=True)
    parent, daughter_i, daughter_j = sp.symbols('t_a t_i t_j', real=True)
    alp_i, alp_j = sp.symbols('v_a_i v_a_j', real=True)
    parent_inverse = sp.MatrixSymbol('M_a_inverse_C', 3, 3)
    meson_i = sp.MatrixSymbol('Phi_i', 3, 3)
    meson_j = sp.MatrixSymbol('Phi_j', 3, 3)
    scalar_1 = sp.MatrixSymbol('hat_S_1', 3, 3)
    scalar_3 = sp.MatrixSymbol('hat_S_3', 3, 3)
    pseudoscalar_2 = sp.MatrixSymbol('hat_P_2', 3, 3)
    alp = parent + epsilon * (alp_i * daughter_i + alp_j * daughter_j)
    meson = 2 * epsilon * parent * parent_inverse + daughter_i * meson_i + daughter_j * meson_j
    # f_a=f_pi/epsilon. Trace linearity uses SymPy's matrix-expression expansion.
    directions = (2*epsilon*parent*parent_inverse, daughter_i*meson_i, daughter_j*meson_j)
    scalar_term = sum((sp.Trace(scalar_1 * left * right) for left in directions for right in directions), sp.S.Zero)
    lagrangian = (-B_0*epsilon*alp*scalar_term/(2*f_pi)
        + B_0*epsilon**2*alp**2*sp.Trace(meson*pseudoscalar_2)/(2*f_pi)
        - m_0**2*theta_2*epsilon**2*alp**2*sp.Trace(meson)/(2*N_f*f_pi)
        + B_0*epsilon**3*alp**3*sp.Trace(scalar_3)/(6*f_pi)
        - m_0**2*theta_1*theta_2*epsilon**3*alp**3/(N_f*f_pi))
    expanded = lagrangian.replace(lambda term: isinstance(term, sp.Trace),
                                   lambda term: sp.Trace(term.arg.expand()).doit()).expand()
    traces = sorted(expanded.atoms(sp.Trace), key=str)
    placeholders = sp.symbols(f'trace_0:{len(traces)}')
    derivative = expanded.xreplace(dict(zip(traces, placeholders))).diff(parent, daughter_i, daughter_j)
    derivative = derivative.subs({parent: 0, daughter_i: 0, daughter_j: 0}).xreplace(dict(zip(placeholders, traces))).expand()
    return lagrangian, derivative, derivative.coeff(epsilon, 1) * epsilon


def derive_projector_and_pole_shift():
    """Derive a spectral projection and a meson-pole shift, in natural units.

    A two-dimensional symmetric example makes the spectral identity explicit;
    its polynomial form extends to every distinct eigenvalue of the meson
    mass operator. The source components have mass dimension two.
    """
    meson_mass_1, meson_mass_2, alp_mass = sp.symbols('m_1_squared m_2_squared m_a_squared', real=True)
    mixing_1, mixing_2, epsilon, s = sp.symbols('c_1 c_2 epsilon s', real=True)
    rotation = sp.Matrix([[sp.Rational(3,5), -sp.Rational(4,5)], [sp.Rational(4,5), sp.Rational(3,5)]])
    meson_mass = rotation * sp.diag(meson_mass_1, meson_mass_2) * rotation.T
    source = rotation * sp.Matrix([mixing_1, mixing_2])
    projector = sp.simplify((meson_mass - meson_mass_2*sp.eye(2))/(meson_mass_1-meson_mass_2))
    residue = sp.simplify((meson_mass_1-s)*(meson_mass-s*sp.eye(2)).inv()*source)
    residue = sp.simplify(residue.subs(s, meson_mass_1))
    # The canonically normalized meson generator has trace norm two, hence
    # the off-diagonal quadratic matrix element is sqrt(2)*epsilon*c_1.
    shift = sp.symbols('delta', real=True)
    secular = (meson_mass_1+shift-alp_mass)*shift-2*epsilon**2*mixing_1**2
    shift_coefficient = sp.symbols('shift_coefficient', real=True)
    shift_equation = sp.expand(secular.subs(shift, epsilon**2*shift_coefficient)).coeff(epsilon, 2)
    shift_coefficient = sp.solve(shift_equation, shift_coefficient)[0]
    return projector, source, residue, sp.factor(shift_coefficient)


def derive_eigenbasis_pole_shift():
    """Resolve the inverse on normalized meson eigenvectors, in natural units.

    The two real components use the trace-orthonormal basis obtained by
    dividing each flavor matrix of trace norm two by sqrt(2).
    """
    mass_1, mass_2, alp_mass, s = sp.symbols('m_1_squared m_2_squared m_a_squared s', real=True)
    coupling_1, coupling_2, epsilon = sp.symbols('c_1 c_2 epsilon', real=True)
    mass_operator = sp.diag(mass_1, mass_2)
    source = sp.Matrix([coupling_1, coupling_2])
    inverse_source = (mass_operator-s*sp.eye(2)).inv()*source
    residue = sp.simplify(((mass_1-s)*inverse_source).subs(s, mass_1))
    shift, shift_coefficient = sp.symbols('delta shift_coefficient', real=True)
    # Pole equation after retaining the singular eigenvector contribution.
    secular = (mass_1+shift-alp_mass)*shift-2*epsilon**2*coupling_1**2
    equation = sp.expand(secular.subs(shift,epsilon**2*shift_coefficient)).coeff(epsilon,2)
    return inverse_source, residue, sp.solve(equation,shift_coefficient)[0]


def derive_kaon_pion_comparison():
    """Reduce the weak-mixing invariant for crossed K0 and pi0 generators."""
    m_d, m_s, B_0, f_a, real_c, imag_c = sp.symbols('m_d m_s B_0 f_a c_real c_imag', real=True)
    quark_mass = sp.diag(m_d, m_d, m_s)
    coupling = sp.Matrix([[0, 0, 0], [0, 0, real_c+sp.I*imag_c], [0, real_c-sp.I*imag_c, 0]])
    outgoing_kaon = sp.sqrt(2)*sp.Matrix([[0, 0, 0], [0, 0, 0], [0, 1, 0]])
    outgoing_pion = sp.diag(1, -1, 0)
    anticommutator = outgoing_kaon*outgoing_pion + outgoing_pion*outgoing_kaon
    mass_commutator = coupling*quark_mass - quark_mass*coupling
    reduced = sp.factor(-B_0*sp.trace(anticommutator*mass_commutator)/(2*f_a))
    generator_commutator = outgoing_kaon*outgoing_pion - outgoing_pion*outgoing_kaon
    earlier = sp.factor(-B_0*(m_s-m_d)*sp.trace(coupling*generator_commutator)/(2*f_a))
    return anticommutator, reduced, earlier


if __name__ == '__main__':
    for label, expression in zip(('Dressed cubic terms with pole directions', 'Differentiated amplitude', 'Leading amplitude'), derive()):
        print(f'{label}:\n{sp.latex(expression)}\n')
    for label, expression in zip(('Meson spectral projector', 'Mixing source', 'Inverse-operator residue', 'Meson pole shift coefficient'), derive_projector_and_pole_shift()):
        print(f'{label}:\n{sp.latex(expression)}\n')
    for label, expression in zip(('Inverse in normalized meson eigenbasis', 'First eigenmatrix residue', 'Eigenbasis meson pole shift coefficient'), derive_eigenbasis_pole_shift()):
        print(f'{label}:\n{sp.latex(expression)}\n')
    print('General meson mass-operator identity, difference:\n'
          f'{sp.latex(sp.expand(eom_identity_left-eom_identity_right))}\n')
    for label, expression in zip(('Crossed generator anticommutator', 'Leading K0 pi0 vertex iM', 'Earlier weak-mixing result iM'), derive_kaon_pion_comparison()):
        print(f'{label}:\n{sp.latex(expression)}\n')
