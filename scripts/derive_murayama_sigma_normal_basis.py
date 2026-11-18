#!/usr/bin/env python3
"""Neutral radial widths in the pi0, eta8, eta0 basis; GeV units, LO.

Run: .venv/bin/python scripts/derive_murayama_sigma_normal_basis.py
"""
import sympy as sp


def derive():
    B0, mu, md, ms, m1, fpi, fphi, s = sp.symbols(
        'B_0 m_u m_d m_s m_1 f_pi f_varphi s', real=True, nonzero=True)
    meson = B0 * sp.Matrix([
        [mu+md, (mu-md)/sp.sqrt(3), sp.sqrt(sp.Rational(2, 3))*(mu-md)],
        [(mu-md)/sp.sqrt(3), (mu+md+4*ms)/3, sp.sqrt(2)*(mu+md-2*ms)/3],
        [sp.sqrt(sp.Rational(2, 3))*(mu-md), sp.sqrt(2)*(mu+md-2*ms)/3,
         2*(mu+md+ms)/3]]) + sp.diag(0, 0, m1**2)
    mixing = -B0*mu*fpi/fphi*sp.Matrix([1, 1/sp.sqrt(3), sp.sqrt(sp.Rational(2, 3))])
    mass = meson.row_join(mixing).col_join(mixing.T.row_join(sp.Matrix([[B0*mu*fpi**2/fphi**2]])))
    meson_polynomial = sp.factor((s*sp.eye(3)-meson).det())
    full_polynomial = sp.factor((s*sp.eye(4)-mass).det())
    shift = sp.simplify(-(meson-s*sp.eye(3)).inv()*mixing)
    # The residue of the a-a propagator is det(s-M_PP)/d det(s-M)/ds.
    residue = meson_polynomial/sp.diff(full_polynomial, s)
    mass_i, mass_j, quartic = sp.symbols('m_i m_j lambda', real=True)
    overlap_i, overlap_j, axion_i, axion_j = sp.symbols('P_i P_j a_i a_j', real=True)
    axion_mass = B0*mu*fpi**2/fphi**2
    # The a row of the pole equation fixes the up-meson overlap.
    overlap_i_solution = sp.solve(axion_mass*axion_i-B0*mu*fpi/fphi*overlap_i-mass_i**2*axion_i, overlap_i)[0]
    overlap_j_solution = overlap_i_solution.subs({mass_i: mass_j, axion_i: axion_j})
    coupling = (2*quartic*fphi*axion_i*axion_j+B0*mu/fphi*overlap_i*overlap_j).subs(
        {overlap_i: overlap_i_solution, overlap_j: overlap_j_solution})
    coupling = sp.factor(coupling.subs({axion_i: sp.sqrt(residue.subs(s, mass_i**2)),
                                         axion_j: sp.sqrt(residue.subs(s, mass_j**2))}))
    return locals()


def main():
    result = derive()
    for label, key in [('Neutral mass matrix', 'mass'), ('Meson polynomial', 'meson_polynomial'),
                       ('Full polynomial', 'full_polynomial'), ('ALP-dependent shift', 'shift'),
                       ('ALP propagator residue', 'residue'), ('Neutral coupling with residues substituted', 'coupling')]:
        print(label + ':\n' + sp.latex(result[key]) + '\n')


if __name__ == '__main__':
    main()
