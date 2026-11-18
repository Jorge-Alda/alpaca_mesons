#!/usr/bin/env python3
"""Derive only-up radial vertices from dressed sources; GeV units, LO.

Run: .venv/bin/python scripts/derive_murayama_sigma_source_coefficients.py
"""
import sympy as sp


def derive():
    up, down, strange, B0, fpi, fa, anomaly, scalar_mass = sp.symbols(
        'm_u m_d m_s B_0 f_pi f_a m_1 m_sigma', positive=True)
    x, sigma, s, mass_i, mass_j = sp.symbols('x sigma s m_i m_j', real=True)
    source = sp.diag(up*(1+sigma/fa)*sp.exp(-sp.I*x), down, strange)
    scalar = (source + source.conjugate().T)/2
    pseudoscalar = (source - source.conjugate().T)/(2*sp.I)
    scalar0 = scalar.subs({x: 0, sigma: 0})
    pseudoscalar1 = sp.diff(pseudoscalar, x).subs({x: 0, sigma: 0})
    scalar2 = sp.diff(scalar, x, 2).subs({x: 0, sigma: 0})
    scalar0_radial = sp.diff(scalar, sigma).subs({x: 0, sigma: 0})
    generators = [sp.diag(1, -1, 0), sp.diag(1, 1, -2)/sp.sqrt(3), sp.sqrt(sp.Rational(2,3))*sp.eye(3)]
    mesons = sp.Matrix(3, 3, lambda i, j: B0*sp.trace(scalar0*(generators[i]*generators[j]+generators[j]*generators[i]))/2 + (anomaly**2 if i == j == 2 else 0))
    mixing = sp.Matrix([B0*fpi/fa*sp.trace(pseudoscalar1*g) for g in generators])
    axion = -B0*fpi**2/fa**2*sp.trace(scalar2)
    full = mesons.row_join(mixing).col_join(mixing.T.row_join(sp.Matrix([[axion]])))
    meson_polynomial = sp.factor((s*sp.eye(3)-mesons).det())
    flavour = s*sp.eye(3)-2*B0*scalar0
    source_polynomial = sp.factor(flavour.det()-anomaly**2*sp.trace(flavour.adjugate())/3)
    pole_polynomial = sp.factor((s*sp.eye(4)-full).det())
    mass_bracket = sp.factor((scalar_mass**2-axion)/fa + (axion-mass_i**2)*(axion-mass_j**2)/(B0*fpi**2*sp.trace(scalar0_radial)))
    original_bracket = (scalar_mass**2-axion)/fa+fa/(B0*up*fpi**2)*(axion-mass_i**2)*(axion-mass_j**2)
    return locals()


def main():
    result = derive()
    for label, key in [('Dressed mass source', 'source'), ('Scalar source at zero phase', 'scalar0'),
                       ('Pseudoscalar first derivative', 'pseudoscalar1'), ('Scalar second derivative', 'scalar2'),
                       ('Radial scalar derivative', 'scalar0_radial'), ('Meson polynomial', 'meson_polynomial'),
                       ('Source form of meson polynomial', 'source_polynomial'), ('Full pole polynomial', 'pole_polynomial'),
                       ('On-shell coupling factor', 'mass_bracket')]:
        print(f'{label}:\n{sp.latex(result[key])}\n')


if __name__ == '__main__':
    main()
