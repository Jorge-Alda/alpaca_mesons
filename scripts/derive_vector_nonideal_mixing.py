#!/usr/bin/env python
"""Tree-level tensor mass extension. Natural units; run with .venv/bin/python."""
import sympy as sp


def show(label, expression):
    print(f'{label}:\n{sp.latex(expression)}\n')


def main():
    M_V, e_m_V, Delta_V, x_l, x_s, s = sp.symbols('M_V e_m^V Delta_V x_l x_s s', real=True)
    w, p = sp.symbols('omega_I phi_I', real=True)
    theta = sp.symbols('theta_V', real=True)
    T_w = sp.diag(1, 1, 0) / sp.sqrt(2)
    T_p = sp.diag(0, 0, 1)
    T_0 = sp.eye(3) / sp.sqrt(3)
    Q = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
    chi_plus = sp.diag(x_l, x_l, x_s)
    V = T_w*w + T_p*p
    L_mass = M_V**2*sp.trace(V*V)/4 - e_m_V*sp.trace(V*V*chi_plus) + Delta_V*sp.trace(V)**2/12
    mass = sp.simplify(2*sp.hessian(L_mass, (w, p)))
    show('Input vacuum spurion', chi_plus)
    show('Mass action with Lorentz contraction suppressed', L_mass)
    show('Derived mass-squared matrix', mass)
    show('Singlet projector', sp.Matrix([sp.trace(T_w*T_0), sp.trace(T_p*T_0)]) * sp.Matrix([[sp.trace(T_w*T_0), sp.trace(T_p*T_0)]]))
    A, B, D, lam = sp.symbols('A B D lambda', real=True)
    generic = sp.Matrix([[A, B], [B, D]])
    polynomial = sp.expand((generic-lam*sp.eye(2)).det())
    roots = sp.solve(polynomial, lam)
    show('Characteristic polynomial', polynomial)
    show('Squared-mass eigenvalues', sp.Tuple(*roots))
    R = sp.Matrix([[sp.cos(theta), sp.sin(theta)], [-sp.sin(theta), sp.cos(theta)]])
    rotated = sp.simplify(sp.trigsimp(R.T*generic*R))
    show('Rotated off-diagonal entry; set to zero', sp.expand_trig(rotated[0, 1]))
    t = sp.symbols('t', real=True)
    angle_condition = sp.factor(sp.expand_trig(rotated[0,1]).subs({sp.sin(theta):t/sp.sqrt(1+t**2), sp.cos(theta):1/sp.sqrt(1+t**2)}))
    show('Equation for t=tan(theta_V)', angle_condition)
    show('Solutions for tan(theta_V); branch fixed by light eigenstate', sp.Tuple(*sp.solve(angle_condition, t)))
    gap = sp.sqrt((D-A)**2+4*B**2)
    # atan2 selects cos(2 theta)=(D-A)/gap, sin(2 theta)=2B/gap.
    C2, S2 = sp.symbols('C_2 S_2', real=True)
    doubled = rotated.applyfunc(lambda item: sp.expand_trig(sp.expand(item)).subs({sp.sin(theta)**2:(1-C2)/2, sp.cos(theta)**2:(1+C2)/2, sp.sin(theta)*sp.cos(theta):S2/2}))
    diagonal = sp.simplify(doubled.subs({C2:(D-A)/gap, S2:2*B/gap}))
    show('Rotated mass matrix on selected nondegenerate branch', diagonal)
    assert sp.simplify(diagonal-sp.diag(*roots)) == sp.zeros(2)
    show('Branch omega connected to light ideal state when D>A', sp.atan2(2*B, D-A)/2)
    show('Physical omega generator', sp.simplify(T_w*R[0, 0]+T_p*R[1, 0]))
    show('Physical phi generator; positive strange ideal phase', sp.simplify(T_w*R[0, 1]+T_p*R[1, 1]))
    c = sp.Matrix([sp.trace(Q*T_w), sp.trace(Q*T_p)])
    show('Ideal photon overlaps', c)
    show('Physical photon overlaps', sp.simplify(R.T*c))
    show('Isoscalar overlap norm', sp.simplify((c.T*c)[0]))
    e, F_V = sp.symbols('e F_V', real=True)
    g = sp.sqrt(2)*e*F_V*c
    r_w, r_p = sp.symbols('r_omega r_phi', real=True)
    r = sp.Matrix([r_w, r_p])
    # Restrict V_i^{mu nu}=r_i F^{mu nu}; suppress the common F^{mu nu}F_mu nu.
    # The tensor kinetic term on this transverse electromagnetic subspace gives -s I.
    L_em = (r.T*(generic-s*sp.eye(2))*r)[0]/4-(g.T*r)[0]/2
    eom = sp.Matrix([2*sp.diff(L_em, coordinate) for coordinate in r])
    solutions = sp.solve(list(eom), list(r))
    r_solution = sp.Matrix([sp.factor(solutions[coordinate]) for coordinate in r])
    show('Quadratic action restricted to electromagnetic tensor subspace', L_em)
    show('Derived electromagnetic-subspace equations of motion', eom)
    show('Derived response vector r=(mass-sI)^(-1)g', r_solution)
    assert sp.simplify((generic-s*sp.eye(2))*r_solution-g) == sp.zeros(2,1)
    ideal_response = r_solution.subs(B,0)
    assert sp.simplify(ideal_response-sp.Matrix([g[0]/(A-s),g[1]/(D-s)])) == sp.zeros(2,1)
    response = sp.factor((c.T*(generic-s*sp.eye(2)).inv()*c)[0])
    show('Isoscalar flavor response before multiplying F_V factors', response)
    j_w, j_p, M_common = sp.symbols('j_omega j_phi M_common', real=True)
    assert sp.simplify(r_solution.subs({A:M_common**2,D:M_common**2,B:0})-g/(M_common**2-s)) == sp.zeros(2,1)
    assert sp.simplify(rotated.subs({theta:0,B:0})-sp.diag(A,D)) == sp.zeros(2)
    j = sp.Matrix([j_w, j_p])
    common = sp.trigsimp(((R.T*c).T*(R.T*j))[0]/(M_common**2-s))
    show('Common-pole exchange with consistently rotated vertex', common)
    assert sp.simplify(common - (c.T*j)[0]/(M_common**2-s)) == 0
    assert sp.simplify(mass[0,1]-sp.sqrt(2)*Delta_V/3) == 0
    assert sp.simplify((R.T*R)-sp.eye(2)) == sp.zeros(2)
    # Derive the charged K* and rho masses from normalized off-diagonal generators.
    rho_generator = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])/sp.sqrt(2)
    k_generator = sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]])/sp.sqrt(2)
    for name, T in [('rho', rho_generator), ('K_star', k_generator)]:
        show(f'Derived {name} squared mass', sp.simplify(M_V**2-4*e_m_V*sp.trace(T*T*chi_plus)))
    print('Assertions passed: canonical rotation, positive-strange phase, common-pole invariance.')


if __name__ == '__main__':
    main()
