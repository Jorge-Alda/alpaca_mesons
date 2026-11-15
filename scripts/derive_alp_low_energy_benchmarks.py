#!/usr/bin/env python3
"""Conditional tree and double-spectral ALP TFF checks; natural units.

Run .venv/bin/python scripts/derive_alp_low_energy_benchmarks.py --test
"""
import sys
import unittest
import sympy as sp

s1, s2, ma2, mv, fa = sp.symbols('s_1 s_2 m_a2 M_V f_a', nonzero=True)
cg = sp.Symbol('c_gamma_hat')
c1, c2, c3, c5, c6, d1, d2, d3 = sp.symbols('c_1 c_2 c_3 c_5 c_6 d_1 d_2 d_3')
mass = sp.symbols('M_1:3', nonzero=True)
residue = sp.symbols('R_1:3')
source_u = sp.symbols('tr_TQU_1:3')
source_y = sp.symbols('tr_TQY_1:3')
pair_u = [[sp.Symbol(f'tr_TTU_{min(i,j)}{max(i,j)}') for j in range(2)] for i in range(2)]
pair_y = [[sp.Symbol(f'tr_TTY_{min(i,j)}{max(i,j)}') for j in range(2)] for i in range(2)]

def vector_vertex(i, first, second):
    return 2/mv*(source_u[i]*((c1+c2-c5)*ma2+(-c1+c2+c5-2*c6)*first+(c1-c2+c5)*second)-4*c3*source_y[i])

def pair_vertex(i, j, first, second):
    return 4*(pair_u[i][j]*((d1-d3)*ma2+d3*(first+second))-4*d2*pair_y[i][j])

def response(i, invariant):
    return residue[i]/(mass[i]**2-invariant)

tree = cg/(4*sp.pi**2*fa)+2*sp.sqrt(2)*sum(response(i,s1)*vector_vertex(i,s1,s2)+response(i,s2)*vector_vertex(i,s2,s1) for i in range(2))-2*sum(response(i,s1)*response(j,s2)*pair_vertex(i,j,s1,s2) for i in range(2) for j in range(2))
x, y, threshold, cutoff = sp.symbols('x y s_th Lambda2', positive=True)
rho = sp.Function('rho_a')
spectral_integrand = rho(x,y)*((x-s1)**-1*(y-s2)**-1+(x-s2)**-1*(y-s1)**-1)/sp.pi**2
dispersive = sp.Integral(spectral_integrand, (x,threshold,cutoff),(y,threshold,cutoff))

class Checks(unittest.TestCase):
    def test_tree_exchange(self):
        self.assertEqual(sp.simplify(tree-tree.xreplace({s1:s2,s2:s1})),0)
    def test_spectral_exchange(self):
        self.assertEqual(sp.simplify(spectral_integrand-spectral_integrand.xreplace({s1:s2,s2:s1})),0)
    def test_spectral_zero(self):
        self.assertEqual(sp.simplify(spectral_integrand.subs({s1:0,s2:0})-2*rho(x,y)/(sp.pi**2*x*y)),0)
    def test_contact(self):
        self.assertEqual(sp.simplify(tree.subs({r:0 for r in residue})-cg/(4*sp.pi**2*fa)),0)

if __name__ == '__main__':
    if '--test' in sys.argv:
        sys.argv.remove('--test')
        unittest.main()
    else:
        for label, value in [('ordered vector vertex',vector_vertex(0,s1,s2)),('two-vector vertex',pair_vertex(0,1,s1,s2)),('tree benchmark',tree),('symbolic dispersive benchmark',dispersive),('spectral zero-virtuality integrand',spectral_integrand.subs({s1:0,s2:0}))]:
            print(f'{label}:\n{sp.latex(value)}\n')
