#!/usr/bin/env python3
"""Derivatives with respect to p^2 of the external-leg self-energy vertices; GeV units.

For an external pole e = (lambda_e, xi_e) with incoming momentum p and its conjugate
ebar = (lambda_ebar, xi_ebar) with incoming momentum -p, so that k_e . k_ebar = -p^2:
    closed ALP line:    d V(e, ebar, a_loop, a_loop)/d p^2   (four-meson kinetic part only),
    closed radial line: d V(e, ebar, sigma, sigma)/d p^2      (ALP kinetic part only).
Products with the loop momentum are m_a^2 and do not depend on p^2.  Multiply by the tree
amplitude and the external residue factor 1/2; the tadpole symmetry factor 1/2 and the loop
factor are excluded.  No loop integral is evaluated.

Run:  PYTHONPATH=scripts .venv/bin/python scripts/derive_radial_loop_external_residues.py
"""
import sympy as sp

import derive_radial_loop_trace_contractions as common

lambda_e, lambda_ebar = sp.symbols('lambda_e lambda_ebar', commutative=False)
xi_e, xi_ebar = sp.symbols('xi_e xi_ebar', real=True)
p_squared = sp.Symbol('p^2', real=True)
leg_e = (lambda_e, xi_e, 'e')
leg_ebar = (lambda_ebar, xi_ebar, 'ebar')
SELF_ENERGY_DOTS = {frozenset(('e', 'ebar')): -p_squared}
NAMES = {lambda_e: r'\lambda_e', lambda_ebar: r'\lambda_{\bar e}', xi_e: r'\xi_e', xi_ebar: r'\xi_{\bar e}'}

alp_line_slope = sp.diff(common.vertex([leg_e, leg_ebar, common.internal, common.internal], common.ma**2,
                                       SELF_ENERGY_DOTS), p_squared)
radial_line_slope = sp.diff(common.vertex([leg_e, leg_ebar, None, None], common.ms**2, SELF_ENERGY_DOTS), p_squared)

if __name__ == '__main__':
    print(common.NOTATION + '\n')
    print('d V(e, ebar, a, a)/d p^2:\n' + common.to_latex(alp_line_slope, NAMES) + '\n')
    print('d V(e, ebar, sigma, sigma)/d p^2:\n' + common.to_latex(radial_line_slope, NAMES) + '\n')
