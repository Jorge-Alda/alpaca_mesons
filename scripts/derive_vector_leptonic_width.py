#!/usr/bin/env python
"""Tree-level tensor-vector leptonic width and one-spurion source residues.
Natural units, metric (+---). Run .venv/bin/python scripts/derive_vector_leptonic_width.py.
"""
import sympy as sp
from sympy.physics.matrices import mgamma


def show(label, expression):
    print(label + ':\n' + sp.latex(expression) + '\n')


def main():
    M, ml, p, e, alpha, residue = sp.symbols('M m_l p e alpha_em R_V', positive=True)
    identity = sp.eye(4)
    # Explicit gamma matrices use the documented SymPy Dirac representation.
    slash_minus = M/2*mgamma(0) - p*mgamma(3)
    slash_plus = M/2*mgamma(0) + p*mgamma(3)
    spatial_trace = sum(sp.trace((slash_minus+ml*identity)*mgamma(i)*(slash_plus-ml*identity)*mgamma(i)) for i in (1,2,3))
    show('Rest-frame spatial lepton spin trace', sp.simplify(spatial_trace))
    coupling = sp.sqrt(2)*e**2*residue/M
    show('Effective conserved-current coupling from tensor photon transition', coupling)
    amplitude_squared = sp.simplify(coupling**2*spatial_trace/3)
    amplitude_squared = sp.simplify(amplitude_squared.subs(p**2, M**2/4-ml**2))
    show('Spin-summed and three-polarization-averaged amplitude squared', amplitude_squared)
    beta = sp.sqrt(1-4*ml**2/M**2)
    costheta = sp.symbols('z', real=True)
    # dPhi_2 after azimuthal integration = beta dz/(16 pi).
    phase_space = sp.integrate(beta/(16*sp.pi), (costheta, -1, 1))
    show('Integrated two-body phase space', phase_space)
    width = sp.simplify(amplitude_squared*phase_space/(2*M))
    width = sp.simplify(width.subs(e**4, (4*sp.pi*alpha)**2))
    show('Finite-lepton-mass width', width)
    show('Massless limit', sp.simplify(width.subs(ml, 0)))

    FV, z_V, xl, xs, theta = sp.symbols('F_V z_V x_l x_s theta_V', real=True)
    Q = sp.diag(sp.Rational(2,3), -sp.Rational(1,3), -sp.Rational(1,3))
    chi = sp.diag(xl, xl, xs)
    generators = [sp.diag(1,-1,0)/sp.sqrt(2), sp.diag(1,1,0)/sp.sqrt(2), sp.diag(0,0,1)]
    # L = FV/(2sqrt2) Tr(V f+) + z_V/(2sqrt2) Tr(V{chi+,f+}).
    residues = sp.Matrix([sp.simplify(FV*sp.trace(T*Q)+z_V*sp.trace(T*(chi*Q+Q*chi))) for T in generators])
    show('rho and ideal-light/strange current residues', residues)
    rotation = sp.Matrix([[sp.cos(theta),sp.sin(theta)],[-sp.sin(theta),sp.cos(theta)]])
    physical = sp.simplify(rotation.T*residues[1:3,0])
    show('Physical omega and phi current residues', physical)
    show('Singlet-only electromagnetic double trace', sp.trace(Q))
    show('Spurion double-trace source overlap Tr(chi+ Q)', sp.trace(chi*Q))
    show('Spurion double-trace ideal-light and strange generator weights', sp.Matrix([sp.trace(T)*sp.trace(chi*Q) for T in generators[1:]]))
    y_V = sp.symbols('y_V', real=True)
    full_residues = residues + y_V*sp.Matrix([sp.trace(T)*sp.trace(chi*Q) for T in generators])
    show('All ideal residues including spurion double trace', full_residues)
    source_map = full_residues.jacobian([FV,z_V,y_V])
    show('Three-source residue map determinant', sp.factor(source_map.det()))
    assert sp.simplify(width.subs(ml,0)-8*sp.pi*alpha**2*residue**2/(3*M)) == 0
    assert sp.simplify((physical.T*physical)[0]-(residues[1:3,0].T*residues[1:3,0])[0]) == 0


if __name__ == '__main__':
    main()
