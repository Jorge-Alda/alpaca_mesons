#!/usr/bin/env python3
"""Derive pole-aligned cubic amplitude with invariant inverse inputs.

Run: .venv/bin/python scripts/derive_alp_invariant_pole_aligned_amplitude.py
Natural units; all momenta incoming; simple nondegenerate poles.
"""
import sympy as sp


def derive_dressed_contact(return_terms=False):
    """Differentiate the dressed CCWZ cubic mass and singlet terms.

    The dressed connections vanish; Hermitian source matrices need not commute.
    Symbolic inverse matrices denote the existing operator applied to C.
    """
    f_pi, f_a, B_0, m_0, N_f = sp.symbols('f_pi f_a B_0 m_0 N_f', nonzero=True, real=True)
    theta_1, theta_2 = sp.symbols('hat_Theta_1 hat_Theta_2', real=True)
    fields = sp.symbols('t_a t_1 t_2', real=True)
    inverses = [sp.MatrixSymbol(name,3,3) for name in ('M_a_inverse_C','M_1_inverse_C','M_2_inverse_C')]
    scalar_1=sp.MatrixSymbol('hat_S_1',3,3)
    scalar_3=sp.MatrixSymbol('hat_S_3',3,3)
    pseudoscalar_2=sp.MatrixSymbol('hat_P_2',3,3)
    alp=sum(fields)
    meson=sum((2*f_pi*fields[i]*inverses[i]/f_a for i in range(3)),sp.ZeroMatrix(3,3))
    terms = (
        -B_0*alp*(2*f_pi/f_a)**2*sum(fields[i]*fields[j]*sp.Trace(scalar_1*inverses[i]*inverses[j]) for i in range(3) for j in range(3))/(2*f_a),
        f_pi*B_0*alp**2*sp.Trace(meson*pseudoscalar_2)/(2*f_a**2),
        -f_pi*m_0**2*theta_2*alp**2*sp.Trace(meson)/(2*N_f*f_a**2),
        f_pi**2*B_0*alp**3*sp.Trace(scalar_3)/(6*f_a**3),
        -f_pi**2*m_0**2*theta_1*theta_2*alp**3/(N_f*f_a**3),
    )
    def differentiate(term):
        expanded=term.replace(lambda item: isinstance(item, sp.Trace),
                              lambda item: sp.Trace(item.arg.expand()).doit()).expand()
        traces=sorted(expanded.atoms(sp.Trace),key=str)
        placeholders=sp.symbols(f"trace_0:{len(traces)}")
        return expanded.xreplace(dict(zip(traces,placeholders))).diff(*fields).subs(dict.fromkeys(fields,0)).xreplace(dict(zip(placeholders,traces))).expand()
    contributions = tuple(differentiate(term) for term in terms)
    lagrangian = sum(terms)
    amplitude = sum(contributions)
    if return_terms:
        return lagrangian, amplitude, contributions
    return lagrangian, amplitude

if __name__ == '__main__':
    lagrangian, amplitude, contributions = derive_dressed_contact(return_terms=True)
    f_pi, f_a = sp.symbols('f_pi f_a', nonzero=True, real=True)
    inverses = [sp.MatrixSymbol(name, 3, 3) for name in ('M_a_inverse_C', 'M_1_inverse_C', 'M_2_inverse_C')]
    normalization = sp.prod((1+2*(f_pi/f_a)**2*sp.Trace(matrix*matrix))**sp.Rational(-1,2) for matrix in inverses)
    steps = [('Dressed cubic CCWZ terms', lagrangian)]
    steps.extend((f'Differentiated cubic term {index}', term) for index, term in enumerate(contributions, 1))
    steps.extend((('Differentiated amplitude', amplitude), ('Pole normalizations', normalization), ('Complete physical amplitude', normalization*amplitude)))
    for label, expression in steps:
        print(f'{label}:\n{sp.latex(expression)}\n')
