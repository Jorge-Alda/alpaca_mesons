"""Exact covariant closure identities and explicit-vector EOM convention audit.

Run .venv/bin/python scripts/derive_odd_covariant_closure.py
Abstract 3x3 flavor matrices; no full-action operator count inferred.
"""
from pathlib import Path
import contextlib
import hashlib
import json
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]

def comm(A,B):return A*B-B*A
def ac(A,B):return A*B+B*A
def simplify(A):return sp.expand(A).doit()
def zero(A):return simplify(A)==sp.ZeroMatrix(3,3)
def show(label,A):print(label+':\n$$\n'+sp.latex(A)+'\n$$\n')

def chiral():
    u1,u2,R,L,K=tuple(sp.MatrixSymbol(n,3,3) for n in ('u_mu','u_nu','R_munu','L_munu','K'))
    fp=L+R;fm=L-R
    # Curvatures of dressed right/left connections Omega_R/L=Gamma +/- i u/2
    # are -i R/L. Sum and difference solve for Gamma curvature and curl u.
    gamma=simplify((-sp.I*(R+L)+comm(u1,u2)/2)/2)
    curl=simplify((-sp.I*R+sp.I*L)/(-sp.I))
    show('Connection curvature from sum of dressed curvatures',gamma)
    show('Covariant curl u from their difference',curl)
    assert zero(gamma-comm(u1,u2)/4+sp.I*fp/2)
    assert zero(curl+fm)
    # Dressed source Bianchi identities D_OmegaR R=D_OmegaL L=0.
    bR=sp.I*comm(u1,R)/2; bL=-sp.I*comm(u1,L)/2
    plus=simplify(bL+bR);minus=simplify(bL-bR)
    show('Cyclic nabla f_plus (one representative cyclic slot)',plus)
    show('Cyclic nabla f_minus (one representative cyclic slot)',minus)
    assert zero(plus+sp.I*comm(u1,fm)/2)
    assert zero(minus+sp.I*comm(u1,fp)/2)
    # Explicit two-covariant-derivative expansion; commuting partials cancel.
    g1,g2,dg12,dg21,dK1,dK2=tuple(sp.MatrixSymbol(n,3,3) for n in ('Gamma_mu','Gamma_nu','d_mu_Gamma_nu','d_nu_Gamma_mu','d_mu_K','d_nu_K'))
    difference=comm(dg12-dg21,K)+comm(g2,dK1)+comm(g1,dK2)
    difference+=comm(g1,comm(g2,K))-comm(g1,dK2)-comm(g2,dK1)-comm(g2,comm(g1,K))
    assert zero(difference-comm(dg12-dg21+comm(g1,g2),K))
    curvature_remainder=simplify(comm(gamma,K))
    show('Commutator of covariant derivatives acting on K',curvature_remainder)
    # Dressed bifundamental spurions S1=u^dag chi u^dag, S2=u chi^dag u.
    s1,s2,c1,c2=tuple(sp.MatrixSymbol(n,3,3) for n in ('S1','S2','C_mu_1','C_mu_2'))
    ds1=c1+sp.I*ac(u1,s1)/2; ds2=c2-sp.I*ac(u1,s2)/2
    X=s1+s2;Y=sp.I*(s1-s2)
    dx=simplify(ds1+ds2);dy=simplify(sp.I*(ds1-ds2))
    show('nabla chi_plus from dressed derivatives',dx)
    show('nabla (i chi_minus) from dressed derivatives',dy)
    assert zero(dx-(c1+c2)-ac(u1,Y)/2)
    assert zero(dy-sp.I*(c1-c2)+ac(u1,X)/2)
    # Differentiate dy once more, then insert dx; no spurion source is zeroed.
    du,cplus,dcm=tuple(sp.MatrixSymbol(n,3,3) for n in ('nabla_mu_u_nu','chi_mu_plus','nabla_mu_chi_nu_minus'))
    uNu,uMu,Xb,Yb=tuple(sp.MatrixSymbol(n,3,3) for n in ('u_nu','u_mu','chi_plus','i_chi_minus'))
    ddy=simplify(sp.I*dcm-ac(du,Xb)/2-ac(uNu,cplus+ac(uMu,Yb)/2)/2)
    show('Second covariant derivative of i chi_minus: source and three-odd-block remainder',ddy)
    return {'identities_verified':True,
        'curvature_latex':sp.latex(gamma),'curl_u_latex':sp.latex(curl),
        'bianchi_plus_latex':sp.latex(plus),'bianchi_minus_latex':sp.latex(minus),
        'second_chi_minus_derivative_latex':sp.latex(ddy),
        'remainders':['u f_minus','u f_plus','[u_mu,u_nu] K','f_plus K',
                      'covariant scalar/pseudoscalar-source derivatives','u u (i chi_minus)']}

def vector_redefinition():
    # Single normalized flavor component suffices to exhibit the nonclosure.
    # Keep all six independent components of a 4D antisymmetric tensor.
    pairs=[(i,j) for i in range(4) for j in range(i+1,4)];g=(1,-1,-1,-1)
    v=sp.Matrix(sp.symbols('V0:6'));J=sp.Matrix(sp.symbols('J0:6'))
    odd=sp.Matrix(sp.symbols('O0:6'));q=sp.symbols('q0:4')
    M,e=sp.symbols('M epsilon',nonzero=True)
    jet=sp.Matrix(4,6,lambda mu,a:sp.Symbol(f'd{mu}V{a}'))
    V=sp.zeros(4)
    for (i,j),value in zip(pairs,v):V[i,j]=value;V[j,i]=-value
    def component(mu,i,j):
        if i==j:return 0
        sign=1 if i<j else -1
        index=pairs.index((min(i,j),max(i,j)))
        return sign*jet[mu,index]
    divergence=[g[nu]*sum(component(mu,mu,nu) for mu in range(4)) for nu in range(4)]
    kinetic=-sum(g[nu]*divergence[nu]**2 for nu in range(4))/2
    H=sp.diag(*(g[i]*g[j] for i,j in pairs))
    mass=(M**2*(v.T*H*v)[0])/2
    # Euler derivative, then Fourier substitute commuting second derivatives.
    euler=[]
    for a in range(6):
        expression=sp.diff(mass,v[a])
        for mu in range(4):
            first=sp.diff(kinetic,jet[mu,a])
            # q carries upper indices; derivative jets carry lower indices.
            expression-=sum(sp.diff(first,jet[nu,b])*(-g[mu]*q[mu]*g[nu]*q[nu]*v[b]) for nu in range(4) for b in range(6))
        euler.append(sp.expand(expression/H[a,a]))
    E=sp.Matrix(euler).jacobian(v);T=E-M**2*sp.eye(6)
    assert sp.simplify(H*E-E.T*H)==sp.zeros(6)
    pol=sp.symbols('pol0:4');wedge=sp.Matrix([q[i]*pol[j]-q[j]*pol[i] for i,j in pairs])
    qSquared=sum(g[i]*q[i]**2 for i in range(4))
    assert sp.simplify(E*wedge-(M**2-qSquared)*wedge)==sp.zeros(6,1)
    action=(v.T*H*E*v)[0]/2+(v.T*H*J)[0]+e*(v.T*H*odd)[0]
    shifted=action.subs(dict(zip(v,v-e*odd/M**2)),simultaneous=True)
    order=sp.expand(sp.diff(shifted,e).subs(e,0))
    expected=-(v.T*H*T*odd)[0]/M**2-(J.T*H*odd)[0]/M**2
    assert sp.simplify(order-expected)==0
    show('Tensor quadratic Euler operator E(q)',E)
    show('Odd term after field redefinition at first order',sp.factor(order))
    show('Canonical wedge pole eigenvalue',M**2-qSquared)
    return {'explicit_vector_EOM_checked':True,
        'induced_terms':['two-extra-derivative vector interaction','source-only contact'],
        'selected_convention':'retain explicit vector fields; no resonance-EOM elimination',
        'full_action_operator_count':None}

def main():
    transcript=ROOT/'data/vector_sector/odd_covariant_closure.md'
    with transcript.open('w') as handle,contextlib.redirect_stdout(handle):
        print('# Executed chiral identity and explicit-vector closure audit\n')
        identities=chiral();vector=vector_redefinition()
    paths=['scripts/derive_odd_covariant_closure.py','memory/rcht/foundations/rcht-basics.md',
           'memory/rcht/odd-parity/odd-trace-source-eom-audit.md']
    output=ROOT/'data/vector_sector/odd_covariant_closure.json'
    output.write_text(json.dumps({'chiral':identities,'vector':vector,
        'sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}},indent=2)+'\n')
    print(transcript.relative_to(ROOT));print(output.relative_to(ROOT))

if __name__=='__main__':main()
