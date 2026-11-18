#!/usr/bin/env python3
"""Derive LO radial vertices from the dressed only-up source.

Run: .venv/bin/python scripts/derive_murayama_sigma_matrix_invariants.py
Natural units: fields and masses in GeV. Tree level, CP-conserving vacuum.
"""
import sympy as sp


def derive():
    B0, mu, fpi, fa = sp.symbols('B_0 m_u f_pi f_a', positive=True)
    sigma, a, pi, eta8, eta0 = sp.symbols('sigma a pi_0 eta_8 eta_0', real=True)
    meson_fields = sp.Matrix([pi, eta8, eta0])
    generators_up = sp.Matrix([1, 1/sp.sqrt(3), sp.sqrt(sp.Rational(2, 3))])
    phase = a/fa - (generators_up.dot(meson_fields))/fpi
    # This is f_pi^2 B_0 Re Tr(Mhat^dagger U^dagger), up term.
    mass_lagrangian = B0*fpi**2*mu*(1 + sigma/fa)*sp.cos(phase)
    pseudoscalars = sp.Matrix([pi, eta8, eta0, a])
    zero_fields = {field: 0 for field in pseudoscalars}
    zero_fields[sigma] = 0
    mass_matrix = -sp.hessian(mass_lagrangian, pseudoscalars).subs(zero_fields)
    radial_mass_derivative = -sp.hessian(sp.diff(mass_lagrangian, sigma), pseudoscalars).subs(zero_fields)
    radial_kinetic_derivative = sp.diag(0, 0, 0, 2/fa)
    msigma2, mi2, mj2 = sp.symbols('m_sigma_squared m_i_squared m_j_squared', real=True)
    vertex = radial_mass_derivative + (msigma2-mi2-mj2)*radial_kinetic_derivative/2
    return mass_matrix, radial_mass_derivative, radial_kinetic_derivative, vertex


def main():
    print('Tree level, CP-conserving vacuum; fields and masses in GeV.')
    for label, expression in zip(('Up-quark mass Hessian', 'Radial derivative of mass Hessian', 'Radial derivative of kinetic matrix', 'On-shell radial vertex'), derive()):
        print(f'{label}:\n{sp.latex(expression)}\n')


if __name__ == '__main__':
    main()
