#!/usr/bin/env python3
"""Derive cubic and quartic vertices from the nonlinear topological source.

Run: .venv/bin/python scripts/derive_alp_three_pseudoscalar_topological_vertices.py
Natural units, all momenta incoming, with the overall Feynman-rule i omitted.
The singlet source is expanded about Theta(0)=0; theta_n are x derivatives.
"""
import sympy as sp


def topological_vertex(number_of_legs):
    """Differentiate the LO singlet Lagrangian for labeled field directions."""
    f_pi, f_a, m_0_squared, n_f = sp.symbols('f_pi f_a m_0_squared N_f', positive=True)
    leg = sp.symbols(f't_0:{number_of_legs}', real=True)
    alp_component = sp.symbols(f'xi_0:{number_of_legs}', real=True)
    meson_trace = sp.symbols(f'trace_lambda_0:{number_of_legs}', real=True)
    theta = sp.symbols('Theta_1:4', real=True)
    alp = sum(leg[i] * alp_component[i] for i in range(number_of_legs))
    trace_phi = sum(leg[i] * meson_trace[i] for i in range(number_of_legs))
    topological_source = trace_phi / f_pi + 2 * sum(
        theta[n-1] * (alp / f_a)**n / sp.factorial(n) for n in range(1, 4))
    lagrangian = -f_pi**2 * m_0_squared * topological_source**2 / (4 * n_f)
    return sp.factor(sp.diff(lagrangian, *leg).subs({value: 0 for value in leg}))


def main():
    print('LO singlet source, natural units, all momenta incoming; omitted overall i.')
    print('Theta(0)=0; coefficients are derivatives with respect to a/f_a.')
    for number_of_legs in (3, 4):
        print(f'{number_of_legs}-leg singlet vertex:\n{sp.latex(topological_vertex(number_of_legs))}\n')


if __name__ == '__main__':
    main()
