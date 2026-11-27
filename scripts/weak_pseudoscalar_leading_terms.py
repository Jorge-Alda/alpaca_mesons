#!/usr/bin/env python3
"""Leading two-pseudoscalar terms of weak octet and 27-plet operators.

Example: .venv/bin/python scripts/weak_pseudoscalar_leading_terms.py --test
Natural units; external sources vanish; terms through two pseudoscalar fields.
Flavor matrices are abstract and need not commute.  The 27-plet matrix acts
on the tensor product of two flavor spaces.
"""
import sys
import unittest
import sympy as sp
from sympy.matrices.expressions.kronecker import KroneckerProduct

field_order = sp.Symbol('epsilon', real=True)
f_pi = sp.Symbol('f_pi', positive=True)
phi = sp.MatrixSymbol('Phi', 3, 3)
partial_mu_phi = sp.MatrixSymbol('partial_mu_Phi', 3, 3)
weak_matrix = sp.MatrixSymbol('w', 3, 3)
weak_tensor = sp.MatrixSymbol('T', 9, 9)
unit_matrix = sp.Identity(3)

# Series of u=exp(i Phi/2f_pi), with Phi -> epsilon Phi.
u = unit_matrix + sp.I*field_order*phi/(2*f_pi) - field_order**2*phi**2/(8*f_pi**2)
u_dagger = unit_matrix - sp.I*field_order*phi/(2*f_pi) - field_order**2*phi**2/(8*f_pi**2)
partial_mu_u = sp.I*field_order*partial_mu_phi/(2*f_pi) - field_order**2*(partial_mu_phi*phi + phi*partial_mu_phi)/(8*f_pi**2)
partial_mu_u_dagger = -sp.I*field_order*partial_mu_phi/(2*f_pi) - field_order**2*(partial_mu_phi*phi + phi*partial_mu_phi)/(8*f_pi**2)
# Alda:2026nrd vielbein, without external sources.
u_mu = sp.expand(sp.I*(u_dagger*partial_mu_u-u*partial_mu_u_dagger))
left_vielbein = sp.expand(u_dagger*u_mu*u)
left_vielbein_first = sp.diff(left_vielbein, field_order).subs(field_order, 0).doit()

# At quadratic field order, each of the two vielbeins contributes its
# linear term. SymPy has no derivative implementation for Trace of a
# KroneckerProduct, so use the product rule on the matrix factors first.
octet_quadratic_matrix = sp.expand(weak_matrix*left_vielbein_first*left_vielbein_first).doit()
twenty_seven_quadratic_matrix = weak_tensor*KroneckerProduct(left_vielbein_first, left_vielbein_first)

def show(label, expression):
    print(f'{label}:\n{sp.latex(expression)}\n')

def main():
    print('Natural units; no external sources; two pseudoscalar fields; one Lorentz contraction implicit.')
    show('u to quadratic order', u)
    show('u_mu through quadratic field order', u_mu.subs(field_order**3, 0))
    show('Coefficient of epsilon in u^dagger u_mu u', left_vielbein_first)
    show('Octet matrix in the trace at epsilon^2', octet_quadratic_matrix)
    show('27-plet matrix in the tensor-product trace at epsilon^2', twenty_seven_quadratic_matrix)
    show('Octet trace', sp.Trace(octet_quadratic_matrix))
    show('27-plet trace', sp.Trace(twenty_seven_quadratic_matrix))

class TestWeakMesonTerms(unittest.TestCase):
    def test_linear_vielbein(self):
        self.assertEqual(left_vielbein_first, -partial_mu_phi/f_pi)

    def test_octet_quadratic(self):
        self.assertEqual(octet_quadratic_matrix, weak_matrix*partial_mu_phi**2/f_pi**2)

    def test_27_quadratic(self):
        expected = weak_tensor*KroneckerProduct(partial_mu_phi, partial_mu_phi)/f_pi**2
        self.assertEqual(twenty_seven_quadratic_matrix, expected)

if __name__ == '__main__':
    if '--test' in sys.argv:
        sys.argv.remove('--test')
        unittest.main()
    else:
        main()
