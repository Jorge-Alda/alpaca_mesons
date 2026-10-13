"""Stage-4 action derivation; .venv/bin/python scripts/derive_vector_stage4_amplitudes.py.
Natural units GeV, tree level, canonical tensors, real signed source residues.
Explicit finite Lorentz contractions reuse the published operator implementation;
SymPy has no single built-in contraction for this indexed operator list.
"""
import sympy as sp
import derive_radiative_tff_amplitude as t
import derive_radiative_sm_amplitude as r


def ac(a,b): return a*b+b*a

def main():
    show=t.show
    s=sp.symbols('s',real=True)
    masses=sp.symbols('M_rho M_omega M_phi',positive=True)
    residues=sp.symbols('R_rho R_omega R_phi',real=True)
    theta=sp.symbols('theta_V',real=True)
    xl,xs,zv,yv=sp.symbols('x_l x_s z_V y_V',real=True)
    U=sp.Matrix([[sp.cos(theta),-sp.sin(theta)],[sp.sin(theta),sp.cos(theta)]])
    Ts=(r.RHO,sp.cos(theta)*r.OMEGA-sp.sin(theta)*r.PHI,
        sp.sin(theta)*r.OMEGA+sp.cos(theta)*r.PHI)
    X=sp.diag(xl,xl,xs)
    Rs=tuple(sp.simplify(t.FV*sp.trace(T*r.Q)+zv*sp.trace(T*ac(X,r.Q))
                         +yv*sp.trace(T)*sp.trace(X*r.Q)) for T in Ts)
    show('Source-action signed residues (rho, omega, phi)',sp.Matrix(Rs))
    # Differentiate the tensor-subspace action, rather than postulate pole factors.
    rr=sp.symbols('r_rho r_omega r_phi',real=True)
    action=sum((M**2-s)*a**2-2*sp.sqrt(2)*r.e*R*a
               for M,R,a in zip(masses,residues,rr))
    response=sp.solve([sp.diff(action,a) for a in rr],rr)
    W=sum((response[a]*T/(sp.sqrt(2)*r.e) for a,T in zip(rr,Ts)),sp.zeros(3))
    show('Tensor-subspace quadratic action (Fgamma squared / 4 removed)',action)
    show('W(s) from stationary action; V/Fgamma=sqrt(2)e W',W)
    # Covariance under physical versus ideal-coordinate tensor descriptions.
    K=U.T*sp.diag(masses[1]**2-s,masses[2]**2-s)*U
    R=sp.Matrix(residues[1:])
    physical=U.T*sp.diag(1/(masses[1]**2-s),1/(masses[2]**2-s))*R
    assert sp.simplify(K*physical-U.T*R)==sp.zeros(2,1)
    show('Ideal-basis response K(s)^(-1) R_ideal verified by multiplication',physical)
    # Derive a generic pair from actual Fourier tensors before assembling W slots.
    aq,ah=sp.symbols('A_Q A_H',real=True)
    first=t.contractions(t.q1,t.q2,t.e1,t.e2)
    second=t.contractions(t.q2,t.q1,t.e2,t.e1)
    C=t.reduce(-sum(t.c[i]*(-ah if i==2 else aq)*first[i] for i in range(7)))
    H=t.reduce(-(t.d[0]*aq*(first[0]+second[0])-t.d[1]*ah*(first[2]+second[2])
                 +t.d[2]*aq*(first[4]+second[4])+t.d[3]*aq*(first[6]+second[6]))/2)
    show('Generic VJP pair C(A,B;s1,s2) from contractions',C)
    show('Generic VVP pair H(A,B;s1,s2), both assignments',H)
    # Physical neutral TP is independent of vector generators and masses.
    pu,pd,ps=sp.symbols('t_Pu t_Pd t_Ps',real=True)
    TP=sp.diag(pu,pd,ps)
    def traces(A,B):
        return {aq:sp.trace(ac(A,B)*TP),ah:sp.trace(ac(A,B)*ac(TP,r.CHI))}
    # Independent symbolic W entries keep the amplitude manageable and allow arbitrary residues.
    w1=sp.diag(*sp.symbols('w_1u w_1d w_1s',real=True))
    w2=sp.diag(*sp.symbols('w_2u w_2d w_2s',real=True))
    direct=4/(t.F*t.MV)*(C.subs(traces(w1,r.Q))+C.xreplace({t.s1:t.s2,t.s2:t.s1}).subs(traces(w2,r.Q)))
    exchange=-4*sp.sqrt(2)/t.F*H.subs(traces(w1,w2))
    _,_,wzw,old_direct,old_exchange=t.derive()
    anomaly=wzw.subs(t.FQ,traces(r.Q,r.Q)[aq])
    show('General neutral-P TFF WZW term',anomaly)
    show('General neutral-P TFF VJP, W1=W(s1), W2=W(s2)',sp.factor(direct))
    show('General neutral-P TFF VVP',sp.factor(exchange))
    swap={t.s1:t.s2,t.s2:t.s1}
    swap.update(dict(zip(list(w1.diagonal())+list(w2.diagonal()),list(w2.diagonal())+list(w1.diagonal()))))
    assert sp.simplify(direct-direct.xreplace(swap))==0
    assert sp.simplify(exchange-exchange.xreplace(swap))==0
    baseline=dict(zip(w1.diagonal(),(t.FV*r.Q/(t.MV**2-t.s1)).diagonal()))
    baseline.update(dict(zip(w2.diagonal(),(t.FV*r.Q/(t.MV**2-t.s2)).diagonal())))
    old={t.FQ:traces(r.Q,r.Q)[aq],t.FH:traces(r.Q,r.Q)[ah]}
    assert sp.simplify(direct.subs(baseline)-old_direct.subs(old))==0
    assert sp.simplify(exchange.subs(baseline)-old_exchange.subs(old))==0
    # Full neutral basis with universal source, common poles gives FV Q / D.
    equal=sum((t.FV*sp.trace(T*r.Q)*T for T in Ts),sp.zeros(3))
    assert sp.simplify(equal-t.FV*r.Q)==sp.zeros(3)
    # Radiative amplitude retains external tensor wavefunction 1/Mext, not pole 1/Mi.
    C_rad,H_rad=r.derive()
    TV=Ts[1]
    W0=w1
    fq,fh=r.flavor(TV,TP)
    fk=sp.trace(TV*r.comm(r.comm(TP,r.Q),r.CHI))
    gd=2*sp.sqrt(2)*r.e/(r.F*r.Mext*r.MV)*C_rad.subs({r.fq:fq,r.fh:fh,r.fk:fk,t.s1:r.Mext**2,t.s2:0})
    eq,eh=r.flavor(TV,TP,W0)
    gx=-4*r.e/(r.F*r.Mext)*H_rad.subs({r.fq:eq,r.fh:eh,t.s1:r.Mext**2,t.s2:0})
    show('omega-physical direct G; TV can be any physical external vector',sp.factor(gd))
    show('omega-physical exchange G; W0=W(0)',sp.factor(gx))
    for ev,er in ((t.q1,t.e2),(t.e1,t.q2)):
        assert t.contractions(t.q1,t.q2,ev,er,False)==[0]*7
        assert t.contractions(t.q2,t.q1,er,ev,False)==[0]*7
    # Even source/mass insertions have no term linear in one SM pseudoscalar.
    book,F=sp.symbols('lambda_P F',real=True)
    uu=sp.eye(3)+sp.I*book*TP/(sp.sqrt(2)*F)
    ud=sp.eye(3)-sp.I*book*TP/(sp.sqrt(2)*F)
    chip=ud*r.CHI*ud+uu*r.CHI*uu
    fp=-r.e*(uu*r.Q*ud+ud*r.Q*uu)
    assert sp.diff(chip,book).subs(book,0)==sp.zeros(3)
    assert sp.diff(fp,book).subs(book,0)==sp.zeros(3)
    show('One-SM-P coefficient of chi+',sp.diff(chip,book).subs(book,0))
    show('One-SM-P coefficient of f+',sp.diff(fp,book).subs(book,0))
    print('Verified: response covariance, common-action baseline, Bose symmetry, Ward identities, no linear SM-P even insertion.')

if __name__=='__main__': main()
