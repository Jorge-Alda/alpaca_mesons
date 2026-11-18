#!/usr/bin/env python3
"""Derive the radial decay vertex with the meson mass operator.

Run: .venv/bin/python scripts/derive_murayama_sigma_operator_amplitude.py
Natural units, GeV; tree level, CP-conserving only-up polar coordinates.
"""
import sympy as sp
try:
    from .derive_murayama_sigma_matrix_invariants import derive as derive_radial_vertex
except ImportError:
    from derive_murayama_sigma_matrix_invariants import derive as derive_radial_vertex


def derive():
    B0, mu, fpi, fa = sp.symbols('B_0 m_u f_pi f_a', positive=True)
    mi2, mj2 = sp.symbols('m_i_squared m_j_squared', real=True)
    source_i = sp.symbols('A_i_u A_i_d A_i_s', real=True)
    source_j = sp.symbols('A_j_u A_j_d A_j_s', real=True)
    # These are diagonal entries of mathbb M_{m_i²}^{-1}(P_1).
    meson_i = sp.Matrix([-2*B0*fpi/fa*(source_i[0]-source_i[1])/2,
        -2*B0*fpi/fa*(source_i[0]+source_i[1]-2*source_i[2])/(2*sp.sqrt(3)),
        -2*B0*fpi/fa*sum(source_i)/sp.sqrt(6), 1])
    meson_j = sp.Matrix([-2*B0*fpi/fa*(source_j[0]-source_j[1])/2,
        -2*B0*fpi/fa*(source_j[0]+source_j[1]-2*source_j[2])/(2*sp.sqrt(3)),
        -2*B0*fpi/fa*sum(source_j)/sp.sqrt(6), 1])
    _, _, _, radial_vertex = derive_radial_vertex()
    projected = sp.factor((meson_i.T * radial_vertex * meson_j)[0])
    normal_i = sp.factor((meson_i.T * meson_i)[0])
    normal_j = sp.factor((meson_j.T * meson_j)[0])
    return meson_i, meson_j, projected, normal_i, normal_j


def derive_dimensionful_coefficients():
    """Derive ordinary coefficients in sigma and dimensionful a."""
    mu, md, ms, fa, sigma, angular = sp.symbols(
        'm_u m_d m_s f_a sigma a', real=True, positive=True
    )
    source = sp.diag(mu * (1 + sigma/fa) * sp.exp(-sp.I*angular/fa), md, ms)
    scalar = sp.simplify((source + source.conjugate().T) / 2)
    pseudoscalar = sp.simplify((source - source.conjugate().T) / (2*sp.I))
    coefficients = {}
    for name, expression, radial_order, phase_order in (
        ('S00', scalar, 0, 0), ('S10', scalar, 1, 0),
        ('P01', pseudoscalar, 0, 1), ('P11', pseudoscalar, 1, 1),
        ('S02', scalar, 0, 2), ('S12', scalar, 1, 2),
    ):
        coefficients[name] = sp.simplify(
            expression.diff(sigma, radial_order).diff(angular, phase_order)
            .subs({sigma: 0, angular: 0})
            / (sp.factorial(radial_order) * sp.factorial(phase_order))
        )
    return source, coefficients


def main():
    print('Tree level; no expansion in f_pi/f_a; masses and fields in GeV.')
    for label, expression in zip(('First pole direction', 'Second pole direction',
                                  'Unnormalized on-shell vertex', 'First norm', 'Second norm'), derive()):
        print(f'{label}:\n{sp.latex(expression)}\n')
    source, coefficients = derive_dimensionful_coefficients()
    print(f'Dressed mass source:\n{sp.latex(source)}\n')
    for name, coefficient in coefficients.items():
        print(f'Ordinary coefficient {name}:\n{sp.latex(coefficient)}\n')


if __name__ == '__main__':
    main()
