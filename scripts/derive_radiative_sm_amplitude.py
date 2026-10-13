"""Action-based SM V-P-photon amplitude and complete neutral-vector flavor sum.

Run .venv/bin/python scripts/derive_radiative_sm_amplitude.py.
Uses RChT fields Phi_R and their explicit conversion to Alda's Phi_A;
all masses are GeV, chi entries GeV^2. No fitted parameter is introduced.
"""
import sympy as sp
from derive_radiative_tff_amplitude import contractions, reference, derive as derive_tff, q1, q2, e1, e2, m, s1, s2, reduce, show

F,FV,MV,Mext,e = sp.symbols('F F_V M_V M_ext e',positive=True)
fq,fh,fk=sp.symbols('F_Q F_H F_K',real=True)
c=sp.symbols('c_1:8',real=True)
d=sp.symbols('d_1:5',real=True)
chi_u,chi_d,chi_s=sp.symbols('chi_u chi_d chi_s',real=True)
CHI=sp.diag(chi_u,chi_d,chi_s)
Q=sp.diag(sp.Rational(2,3),-sp.Rational(1,3),-sp.Rational(1,3))
RHO=sp.diag(1,-1,0)/sp.sqrt(2)
OMEGA=sp.diag(1,1,0)/sp.sqrt(2)
PHI=sp.diag(0,0,1)
NEUTRAL=(('rho',RHO),('omega',OMEGA),('phi',PHI))

def ac(a,b): return a*b+b*a
def comm(a,b): return a*b-b*a

def flavor(v,p,charge=Q):
    return sp.simplify(sp.trace(ac(v,charge)*p)),sp.simplify(sp.trace(ac(v,charge)*ac(p,CHI)))

def derive():
    first=contractions(q1,q2,e1,e2)
    second=contractions(q2,q1,e2,e1)
    # Remove common -2sqrt2 e/(FMext) building blocks and external tensor factor.
    C=reduce(-sum(c[i]*(-fh if i==2 else fq)*first[i] for i in range(7)))
    # O4 has i[f_-,chi+] and therefore the derivative-block prefactor sign.
    C=sp.expand(C-c[3]*fk*first[2])
    H=reduce(-(d[0]*fq*(first[0]+second[0])-d[1]*fh*(first[2]+second[2])
             +d[2]*fq*(first[4]+second[4])+d[3]*fq*(first[6]+second[6]))/2)
    return C,H

def derive_widths():
    G=sp.symbols('G',real=True)
    photons=[sp.Matrix([0,1,0,0]),sp.Matrix([0,0,1,0])]
    axes=photons+[sp.Matrix([0,0,0,1])]
    kv=(Mext**2-m**2)/(2*Mext)
    qv=sp.Matrix([Mext,0,0,0])
    pp=sp.Matrix([(Mext**2+m**2)/(2*Mext),0,0,-kv])
    sumv=sp.simplify(sum(reference(qv,pp,v,g)**2 for v in axes for g in photons))
    widthv=sp.factor(kv/(8*sp.pi*Mext**2)*G**2*sumv/3)
    kp=(m**2-Mext**2)/(2*m)
    Ev=(m**2+Mext**2)/(2*m)
    qv=sp.Matrix([Ev,0,0,kp]);pp=sp.Matrix([m,0,0,0])
    vectors=photons+[sp.Matrix([kp/Mext,0,0,Ev/Mext])]
    sump=sp.simplify(sum(reference(qv,pp,v,g)**2 for v in vectors for g in photons))
    widthp=sp.factor(kp/(8*sp.pi*m**2)*G**2*sump)
    show('V decay polarization sum',sumv)
    show('V decay width, vector spin average 1/3',widthv)
    show('P decay polarization sum',sump)
    show('P decay width, spin-zero initial state',widthp)
    assert sp.simplify(widthv-G**2*(Mext**2-m**2)**3/(96*sp.pi*Mext**3))==0
    assert sp.simplify(widthp-G**2*(m**2-Mext**2)**3/(32*sp.pi*m**3))==0

def rank_probe():
    """Illustrative action rows: explicit old pseudoscalar inputs, canonical vectors.

    This supplies no physical eta fit: the existing /sqrt(12) pseudoscalar
    matrices are retained as explicit exploratory inputs, not attributed to
    the canonical RChT octet. chi is a stated leading-order isospin assumption.
    """
    import numpy as np
    C,H=derive()
    pars=c+d
    A=sp.Matrix([[sp.diff(x,y) for y in pars] for x in (
        c[0]+4*c[2],c[0]-c[1]+c[4],c[4]-c[5],d[0]+8*d[1],d[2])])
    eta8=sp.diag(1,1,-2)/sp.sqrt(12)
    eta1=sp.eye(3)/sp.sqrt(3)
    eta=2*sp.sqrt(2)/3*eta8-eta1/3
    etap=eta8/3+2*sp.sqrt(2)/3*eta1
    def E(i,j):
        x=sp.zeros(3);x[i,j]=1;return x
    masses={name:sp.Rational(value) for name,value in (
        ('rho','.7753'),('omega','.7827'),('phi','1.0195'),('eta','.5479'),
        ('etap','.9578'),('pi','.1396'),('K0','.4976'),('K+','.4937'),
        ('Kstar0','.8956'),('Kstar+','.8922'))}
    channels=[('omega pi',OMEGA,RHO,'omega','pi'),('rho pi',RHO,RHO,'rho','pi'),
              ('rho eta',RHO,eta,'rho','eta'),('omega eta',OMEGA,eta,'omega','eta'),
              ('etap rho',RHO,etap,'rho','etap'),('etap omega',OMEGA,etap,'omega','etap'),
              ('phi eta',PHI,eta,'phi','eta'),('Kstar0 K0',E(2,1),E(1,2),'Kstar0','K0'),
              ('Kstar+ K+',E(2,0),E(0,2),'Kstar+','K+'),('rho+ pi+',E(1,0),E(0,1),'rho','pi')]
    values={F:sp.Rational('.092'),FV:sp.Rational('.13'),MV:sp.Rational('.770'),e:sp.Rational('.3028'),
            chi_u:masses['pi']**2,chi_d:masses['pi']**2,
            chi_s:2*masses['K+']**2-masses['pi']**2}
    vectors=[sp.eye(11)[:,j] for j in (6,10)]
    vectors += [sp.Matrix([0,1,0,0,1,1,0,0,0,0,0])]
    for scheme in ('common','split'):
        rows=[]
        print('Action-based rank probe: '+scheme+' internal masses')
        for label,v,p,vname,pname in channels:
            fqv,fhv=flavor(v,p)
            fkv=sp.simplify(sp.trace(v*comm(comm(p,Q),CHI)))
            Cvp=C.subs({fq:fqv,fh:fhv,fk:fkv,s1:masses[vname]**2,s2:0,m:masses[pname]})
            exchange=0
            for name,T in NEUTRAL:
                aq,ah=flavor(v,p,T)
                D=MV**2 if scheme=='common' else masses[name]**2
                exchange+=sp.trace(T*Q)*H.subs({fq:aq,fh:ah,s1:masses[vname]**2,s2:0,m:masses[pname]})/D
            G=2*sp.sqrt(2)*e/(F*MV*masses[vname])*Cvp-4*e*FV/(F*masses[vname])*exchange
            row=sp.Matrix([[sp.simplify(sp.diff(G,x).subs(values)) for x in pars]])
            for null in vectors:
                assert sp.simplify((row*null)[0])==0
            numeric=np.array(row,dtype=float).reshape(11)
            rows.append(numeric/np.linalg.norm(numeric))
            print('  '+label+' c4='+str(sp.simplify(row[3])))
        matrix=np.vstack([np.array(A,dtype=float)]+rows)
        sv=np.linalg.svd(matrix,compute_uv=False)
        rank=int((sv>1e-11).sum())
        print('  OPE rank 5 -> '+str(rank)+' of 11; singular values '+str(sv))
        assert rank==(7 if scheme=='common' else 8)
        # TFF comparison uses the same explicit exploratory pseudoscalars and
        # chi input, not FH=2*mP^2*FQ for eta/eta'. Its poles remain common.
        import derive_radiative_tff_amplitude as tff
        _,_,_,direct_tff,exchange_tff=derive_tff()
        for p,pname in ((RHO,'pi'),(eta,'eta'),(etap,'etap')):
            fq_t,fh_t=flavor(Q,p)
            sub={tff.F:values[F],tff.FV:values[FV],tff.MV:values[MV],
                 tff.FQ:fq_t,tff.FH:fh_t.subs(values),tff.m:masses[pname],
                 tff.s1:-5,tff.s2:0}
            r=np.array([float(sp.diff(direct_tff+exchange_tff,x).subs(sub)) for x in pars])
            rows.append(r/np.linalg.norm(r))
        together=np.vstack([np.array(A,dtype=float)]+rows)
        combined_rank=int((np.linalg.svd(together,compute_uv=False)>1e-11).sum())
        print('  Plus full TFF at Q2=5, same independent chi inputs: rank '+str(combined_rank))
        assert combined_rank==8

def main():
    C,H=derive()
    show('VJP numerator from tensor contractions',C)
    show('VVP numerator, both vector Wick assignments',H)
    # f_- one-field expansion obtained by differentiation with a bookkeeper.
    z,Fchi=sp.symbols('z F',real=True)
    TP=sp.MatrixSymbol('T_P',3,3)
    QQ=sp.MatrixSymbol('Q',3,3)
    uu=sp.Identity(3)+sp.I*z*TP/(sp.sqrt(2)*Fchi)
    ud=sp.Identity(3)-sp.I*z*TP/(sp.sqrt(2)*Fchi)
    fm=sp.expand(-e*(uu*QQ*ud-ud*QQ*uu))
    show('Linear f_- / Fgamma',sp.diff(fm,z).subs(z,0))
    response=sum((sp.trace(T*Q)*T for _,T in NEUTRAL),sp.zeros(3))
    show('Complete equal-pole photon flavor response',sp.simplify(response))
    assert sp.simplify(response-Q)==sp.zeros(3)
    for label,v,p in (('omega pi',OMEGA,RHO),('rho pi',RHO,RHO),
                      ('Kstar+ K+',sp.eye(3)[:,0]*sp.eye(3)[2,:],sp.eye(3)[:,2]*sp.eye(3)[0,:])):
        # Charged entries are field matrices appropriate to conjugate external legs.
        FQ,FH=flavor(v,p)
        FK=sp.simplify(sp.trace(v*comm(comm(p,Q),CHI)))
        show(label+' direct (FQ,FH,FK)',sp.Matrix([FQ,FH,FK]))
        eq=0;eh=0
        for name,T in NEUTRAL:
            w=sp.trace(T*Q)
            aq,ah=flavor(v,p,T)
            show(label+' '+name+' exchange weights',sp.Matrix([sp.simplify(w*aq),sp.simplify(w*ah)]))
            eq+=w*aq;eh+=w*ah
        assert sp.simplify(eq-FQ)==0 and sp.simplify(eh-FH)==0
    omega_sub={fq:1,fh:2*m**2,fk:0,s1:Mext**2,s2:0}
    Gdir=reduce(2*sp.sqrt(2)*e/(F*Mext*MV)*C.subs(omega_sub))
    Gex=reduce(-4*e*FV/(F*Mext*MV**2)*H.subs(omega_sub))
    show('Canonical omega pi direct Eq27',Gdir)
    show('Canonical omega pi exchange Eq27',Gex)
    assert reduce(Gdir-2*sp.sqrt(2)*e/(F*Mext*MV)*((c[1]-c[0]+c[4]-2*c[5])*Mext**2+(c[0]+c[1]+8*c[2]-c[4])*m**2))==0
    assert reduce(Gex+4*e*FV/(F*Mext*MV**2)*(d[2]*Mext**2+(d[0]+8*d[1]-d[2])*m**2))==0
    derive_widths()
    rank_probe()

if __name__=='__main__': main()
