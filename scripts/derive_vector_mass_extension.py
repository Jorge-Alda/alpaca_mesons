#!/usr/bin/env python
"""Assumed double-trace tensor mass extension; natural units.
Run .venv/bin/python scripts/derive_vector_mass_extension.py.
"""
import sympy as sp


def show(label, expression):
    print(f'{label}:\n{sp.latex(expression)}\n')


def derive():
    M, e_m, Delta, kap, xl, xs = sp.symbols('M_V e_m^V Delta_V kappa_V x_l x_s', real=True)
    w, p = sp.symbols('omega_I phi_I', real=True)
    chi = sp.diag(xl, xl, xs)
    Tw, Tp = sp.diag(1, 1, 0)/sp.sqrt(2), sp.diag(0, 0, 1)
    V = w*Tw+p*Tp
    extra = kap*sp.trace(V)*sp.trace(V*chi)/4
    action = M**2*sp.trace(V*V)/4-e_m*sp.trace(V*V*chi)+Delta*sp.trace(V)**2/12+extra
    matrix = sp.simplify(2*sp.hessian(action, (w,p)))
    r, k, W, P = sp.symbols('r k W P', real=True)
    charged = []
    for T in (sp.Matrix([[0,1,0],[1,0,0],[0,0,0]])/sp.sqrt(2),sp.Matrix([[0,0,1],[0,0,0],[1,0,0]])/sp.sqrt(2)):
        field = w*T
        L = action.subs({w:0,p:0})+M**2*sp.trace(field*field)/4-e_m*sp.trace(field*field*chi)+Delta*sp.trace(field)**2/12+kap*sp.trace(field)*sp.trace(field*chi)/4
        charged.append(sp.simplify(2*sp.diff(L,w,2)))
    base = sp.solve([charged[0]-r,charged[1]-k], [M**2,e_m])
    calibrated = sp.simplify(matrix.subs(base))
    delta = sp.solve(sp.trace(calibrated)-W-P,Delta)[0]
    reduced = sp.simplify(calibrated.subs(Delta,delta))
    polynomial = sp.Poly(sp.expand(reduced.det()-W*P),kap)
    discr = sp.factor(sp.discriminant(polynomial.as_expr(),kap))
    return locals()


def main():
    z=derive()
    for label,key in [('Extra mass invariant','extra'),('Full neutral mass matrix','matrix'),('Charged masses squared','charged'),('Charged calibration','base'),('Trace solution for Delta_V','delta'),('Neutral matrix after trace','reduced'),('Determinant equation polynomial','polynomial'),('Quadratic discriminant','discr')]:
        value=z[key]
        if key=='charged': value=sp.Tuple(*value)
        if key=='base': value=sp.Dict(value)
        if key=='polynomial': value=value.as_expr()
        show(label,value)
    roots=sp.solve(z['polynomial'].as_expr(),z['kap'])
    show('Mass-compatible kappa branches',sp.Tuple(*roots))
    assert sp.simplify(sp.trace(z['reduced'])-z['W']-z['P'])==0
    assert sp.factor(z['polynomial'].coeff_monomial(z['kap']**2)) == -(z['xl']-z['xs'])**2/2
    for root in roots:
        assert sp.simplify(z['polynomial'].as_expr().subs(z['kap'],root))==0
    show('SU(3)-symmetric reduced determinant equation',sp.factor(z['polynomial'].as_expr().subs(z['xs'],z['xl'])))
    print('Assertions passed: trace, determinant roots and quadratic coefficient.')


if __name__=='__main__':
    main()
