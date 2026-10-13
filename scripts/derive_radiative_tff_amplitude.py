"""Derive the neutral-P two-photon amplitude from the actual RChT tensors.

Run: .venv/bin/python scripts/derive_radiative_tff_amplitude.py
Natural units, GeV; epsilon_0123=+1, (+---), all momenta incoming.
Tree level, common vector mass, leading electromagnetic order; diagonal
neutral-P and Q flavor matrices. No short-distance substitution until the end.
The finite component sums implement the explicit Lorentz contractions of the
published operators; SymPy performs every algebraic reduction.
"""
from itertools import permutations
import sympy as sp

m, s1, s2, b = sp.symbols("m_P s_1 s_2 b", real=True)
F, FV, MV, NC = sp.symbols("F F_V M_V N_C", positive=True)
FQ, FH = sp.symbols("F_Q F_H", real=True)
c = sp.symbols("c_1:8", real=True)
d = sp.symbols("d_1:5", real=True)
g = (1, -1, -1, -1)
a = (m**2 + s1 - s2)/(2*m)
q1, q2, p = sp.Matrix([a,0,0,b]), sp.Matrix([m-a,0,0,-b]), sp.Matrix([-m,0,0,0])
e1, e2 = sp.Matrix([0,1,0,0]), sp.Matrix([0,0,1,0])
EPS = [(t, sp.LeviCivita(*t)) for t in permutations(range(4))]

def reduce(x):
    return sp.factor(sp.cancel(sp.expand(x).subs(b**2, a**2-s1)))

def wedge(q,e):
    return q*e.T-e*q.T

def reference(q,r,ev,er):
    return sum(sign*ev[mu]*er[nu]*q[rho]*r[sigma]
               for (mu,nu,rho,sigma),sign in EPS)

def contractions(q,r,ev,er, normalize=True):
    """Strip building-block coefficients, retain Fourier i factors."""
    v, f = sp.I*wedge(q,ev), sp.I*wedge(r,er)
    u = sp.I*p
    result = [0]*7
    for (mu,nu,rho,sigma),sign in EPS:
        for alpha in range(4):
            result[0] += sign*v[mu,nu]*f[rho,alpha]*sp.I*g[alpha]*p[alpha]*u[sigma]
            result[1] += sign*v[mu,alpha]*f[rho,sigma]*sp.I*g[alpha]*p[alpha]*u[nu]
            result[4] += sign*sp.I*g[alpha]*q[alpha]*v[mu,nu]*f[rho,alpha]*u[sigma]
            result[5] += sign*sp.I*g[alpha]*q[alpha]*v[mu,alpha]*f[rho,sigma]*u[nu]
            result[6] += sign*sp.I*q[sigma]*v[mu,nu]*f[rho,alpha]*g[alpha]*u[alpha]
        result[2] += sign*v[mu,nu]*f[rho,sigma]
    S = reference(q,r,ev,er) if normalize else 1
    return [reduce(x/S) for x in result]

def show(label,expr):
    print(label+":\n"+sp.latex(expr)+"\n")

def derive():
    first = contractions(q1,q2,e1,e2)
    second = contractions(q2,q1,e2,e1)
    # u=-sqrt(2)/F dPhi; f+=-2eQ dA; i chi-=sqrt(2)/F {Phi,chi}.
    # Derivative coefficient is +2sqrt(2)e/F, mass coefficient its negative.
    # Vector solution: V=sqrt(2)eFV/(MV^2-s) Q dA, from quadratic EOM.
    # FQ=Tr({Q,Q}TP), FH=Tr({Q,Q}{TP,chi}), chi=2B0 mq.
    direct = sum(c[i]*(-FH if i==2 else FQ)*first[i] for i in range(7))/(MV**2-s1)
    direct += sum(c[i]*(-FH if i==2 else FQ)*second[i] for i in range(7))/(MV**2-s2)
    direct = reduce(4*FV/(MV*F)*direct)
    # VVP requires both Wick assignments of its two vector fields.
    exchange = d[0]*FQ*(first[0]+second[0])-d[1]*FH*(first[2]+second[2])
    exchange += d[2]*FQ*(first[4]+second[4])+d[3]*FQ*(first[6]+second[6])
    exchange = reduce(-2*sp.sqrt(2)*FV**2/F*exchange/((MV**2-s1)*(MV**2-s2)))
    # WZW: two assignments of dA dA in epsilon(mu,nu,rho,sigma).
    wzw_l = sum(sign*((sp.I*q1[mu])*e1[nu]*(sp.I*q2[rho])*e2[sigma]
                     +(sp.I*q2[mu])*e2[nu]*(sp.I*q1[rho])*e1[sigma])
                for (mu,nu,rho,sigma),sign in EPS)
    wzw = reduce(-sp.sqrt(2)*NC/(8*sp.pi**2*F)*(FQ/2)*wzw_l/reference(q1,q2,e1,e2))
    # Global meson phase chosen to make F_pi(0,0)>0 for the anomaly alone.
    return first,second,-wzw,-direct,-exchange

def pole_response():
    """Solve quadratic action in the exact tensor subspace V=r Fgamma."""
    r,gamma,s = sp.symbols('r g_gamma s', real=True)
    # Factor Fgamma_mu_nu(q) Fgamma^mu_nu(-q)/4 from the quadratic action.
    action = (MV**2-s)*r**2-2*gamma*r
    return action,sp.solve(sp.diff(action,r),r)[0]

def main():
    action,response=pole_response()
    show("Quadratic action after factoring Fgamma squared/4",action)
    show("Photon-induced tensor response V/Fgamma",response)
    first,second,wzw,direct,exchange = derive()
    show("All-incoming p+q1+q2", p+q1+q2)
    show("q1^2", reduce(sum(g[i]*q1[i]**2 for i in range(4))))
    show("q2^2", reduce(sum(g[i]*q2[i]**2 for i in range(4))))
    for i in (0,1,2,4,5,6):
        show(f"VJP O{i+1}: photon 1 on vector line", first[i])
        show(f"VJP O{i+1}: photon 2 on vector line", second[i])
    for label,expr in (("WZW F",wzw),("VJP F",direct),("VVP F",exchange)):
        show(label,expr)
    pi = {FQ:sp.sqrt(2)/3,FH:2*m**2*sp.sqrt(2)/3, NC:3}
    onshell = {s1:0,s2:0}
    show("pi VJP both-real",reduce(direct.subs(pi).subs(onshell)))
    show("pi VVP both-real",reduce(exchange.subs(pi).subs(onshell)))
    constraints = {c[2]:-c[0]/4,c[1]:c[0]+c[4],c[5]:c[4]-NC*MV/(64*sp.pi**2*sp.sqrt(2)*FV),
                   d[0]:-8*d[1]-NC*MV**2/(64*sp.pi**2*FV**2)+F**2/(4*FV**2),
                   d[2]:-NC*MV**2/(64*sp.pi**2*FV**2)+F**2/(8*FV**2)}
    constrained = reduce((wzw+direct+exchange).subs({FH:2*m**2*FQ}).subs(constraints, simultaneous=True))
    show("OPE-constrained full F, equal-mass spurion",constrained)
    show("Spacelike pi full F",reduce((wzw+direct+exchange).subs(pi).subs({s1:-sp.Symbol('Q^2',nonnegative=True),s2:0})))
    parameters = c+d
    ope = sp.Matrix([[sp.diff(expr,x) for x in parameters] for expr in (
        c[0]+4*c[2],c[0]-c[1]+c[4],c[4]-c[5],d[0]+8*d[1],d[2])])
    sensitivity = sp.Matrix([[sp.diff(reduce((direct+exchange).subs(FH,2*m**2*FQ)),x)
                              for x in parameters]])
    # All six homogeneous OPE null directions annihilate the complete TFF.
    residual=sp.Matrix([[reduce((sensitivity*v)[0]) for v in ope.nullspace()]])
    show("Complete equal-spurion TFF along OPE null directions",residual)
    assert residual==sp.zeros(1,6)
    show("OPE rank with complete TFF rows",sp.Integer(ope.rank()))
    generic_sensitivity = sp.Matrix([[sp.diff(direct+exchange,x) for x in parameters]])
    invisible = [sp.eye(11)[:,j] for j in (3,6,10)]
    invisible += [sp.Matrix([0,1,0,0,1,1,0,0,0,0,0])]
    generic_residual = sp.Matrix([[reduce((generic_sensitivity*v)[0]) for v in invisible]])
    show("Generic-spurion TFF along c4,c7,d4 and common c2,c5,c6 shift",generic_residual)
    assert generic_residual==sp.zeros(1,4)
    assert reduce(direct-direct.xreplace({s1:s2,s2:s1}))==0
    assert reduce(exchange-exchange.xreplace({s1:s2,s2:s1}))==0
    assert reduce(direct.subs(pi).subs(onshell)-8*sp.sqrt(2)*FV*m**2/(3*F*MV**3)*(c[0]+c[1]+8*c[2]-c[4]))==0
    assert reduce(exchange.subs(pi).subs(onshell)+8*FV**2*m**2/(3*F*MV**4)*(d[0]+8*d[1]-d[2]))==0
    assert reduce(constrained.subs(pi).subs(onshell)-(1/(4*sp.pi**2*F)-F*m**2/(3*MV**4)))==0
    # Replace either photon polarization by its momentum before contraction.
    for ev,er in ((q1,e2),(e1,q2)):
        assert contractions(q1,q2,ev,er,False)==[0]*7
        assert contractions(q2,q1,er,ev,False)==[0]*7

if __name__=="__main__":
    main()
