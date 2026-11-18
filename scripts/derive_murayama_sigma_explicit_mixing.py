#!/usr/bin/env python3
"""Derive the orthogonal neutral-pole transformation. GeV units; tree-level U(3).

Run: .venv/bin/python scripts/derive_murayama_sigma_explicit_mixing.py
"""
import sympy as sp


def derivation():
    B0, mu, md, ms, m1, fpi, fphi, s = sp.symbols('B_0 m_u m_d m_s m_1 f_pi f_varphi s', real=True, nonzero=True)
    # Canonical flavour-diagonal meson coordinates, obtained by an orthogonal
    # change of basis from pi0, eta8, eta0, with anomaly mass m1^2.
    meson_mass = sp.diag(2*B0*mu, 2*B0*md, 2*B0*ms) + m1**2*sp.ones(3)/3
    mixing = sp.Matrix([-sp.sqrt(2)*B0*mu*fpi/fphi, 0, 0])
    axion_mass = B0*mu*fpi**2/fphi**2
    mass_matrix = meson_mass.row_join(mixing).col_join(mixing.T.row_join(sp.Matrix([[axion_mass]])))
    meson_block = meson_mass - s*sp.eye(3)
    determinant = sp.factor(meson_block.det())
    cofactors = sp.simplify(meson_block.adjugate()[:, 0])
    unnormalized = (-mixing[0]*cofactors).col_join(sp.Matrix([determinant]))
    pole_polynomial = sp.factor((mass_matrix - s*sp.eye(4)).det())
    norm_squared = (unnormalized.T*unnormalized)[0]
    overlap_i, overlap_j, axion_i, axion_j, s_j = sp.symbols('overlap_i overlap_j axion_i axion_j s_j', real=True)
    quartic = sp.Symbol('lambda', real=True)
    # The a row of the mass eigenvalue equation fixes the up-flavour overlap.
    overlap_i_solution = sp.solve(sp.Eq(axion_mass*axion_i - sp.sqrt(2)*B0*mu*fpi/fphi*overlap_i, s*axion_i), overlap_i)[0]
    overlap_j_solution = overlap_i_solution.subs({s: s_j, axion_i: axion_j})
    coupling = sp.factor((2*quartic*fphi*axion_i*axion_j + 2*B0*mu/fphi*overlap_i*overlap_j).subs({overlap_i: overlap_i_solution, overlap_j: overlap_j_solution})/(axion_i*axion_j))
    return locals()


def main():
    values = derivation()
    print('Canonical neutral mass matrix, in GeV squared:')
    print(sp.latex(values['mass_matrix']))
    print('Meson block determinant:')
    print(sp.latex(values['determinant']))
    print('Meson block adjugate first column:')
    print(sp.latex(values['cofactors']))
    print('Unnormalized pole column:')
    print(sp.latex(values['unnormalized']))
    print('Pole polynomial:')
    print(sp.latex(values['pole_polynomial']))
    print('Squared pole-column norm:')
    print(sp.latex(values['norm_squared']))
    print('Neutral coupling divided by the two normalized axion components:')
    print(sp.latex(values['coupling']))


if __name__ == '__main__':
    main()
