"""Tree-action short-distance matching; no fit or implicit OPE prior.

Run .venv/bin/python scripts/match_vector_short_distance.py
GeV, positive-strange phi, fixed central 64 scenarios; chiral octet OPE
and physical TFF necessary conditions are separate calculations.
"""
from pathlib import Path
import contextlib
import hashlib
import json
import sympy as sp
import numpy as np
import diagnose_vector_stage4 as stage4

ROOT=Path(__file__).resolve().parents[1]
F,FV,MV,M0,NC=sp.symbols('F F_V M_V M_0 N_C',positive=True)
c=sp.symbols('c1:8',real=True); d=sp.symbols('d1:5',real=True)
theta=c+d
x=sp.symbols('Q2',positive=True)

def show(label,expr):
    print(label+'\n$$\n'+sp.latex(expr)+'\n$$\n')

def octet_ope():
    # Literature Eq.(18), including independent pseudoscalar-density contacts.
    # M0 is the octet pole in the chiral limit; MV is operator normalization.
    p,q,r,z=sp.symbols('p2 q2 r2 z',nonzero=True)
    A=-c[0]+c[1]+c[4]-2*c[5]; B=c[0]-c[1]+c[4]
    C=c[0]+c[1]-c[4]
    Dp=M0**2-p; Dq=M0**2-q
    bracket=4*FV**2*((d[0]-d[2])*r+d[2]*(p+q))/(Dp*Dq*r)
    bracket-=2*sp.sqrt(2)*FV/MV*((r*C+p*A+q*B)/(Dp*r)+(r*C+q*A+p*B)/(Dq*r))
    bracket+=32*FV**2*d[1]/(Dp*Dq)-16*sp.sqrt(2)*FV*c[2]/MV*(1/Dp+1/Dq)-NC/(8*sp.pi**2*r)
    show('Chiral octet correlator bracket: Pi=-condensate/F^2 times bracket',bracket)
    expanded=sp.series(bracket.subs({p:p/z,q:q/z,r:r/z}),z,0,3).removeO().expand()
    leading=sp.factor(expanded.coeff(z,1))
    next_term=sp.factor(expanded.coeff(z,2)-F**2*(p+q+r)/(2*p*q*r))
    show('Order lambda^-2, z=lambda^-2',leading)
    show('Order lambda^-4 minus QCD target',next_term)
    equations=[]
    for term in (leading,next_term):
        equations.extend(sp.Poly(sp.together(term).as_numer_denom()[0],p,q,r).coeffs())
    solved=sp.solve(equations,[c[2],c[1],c[5],d[0],d[2]],dict=True)
    assert len(solved)==1
    show('Solved chiral-octet relations',sp.Dict(solved[0]))
    assert sp.simplify(leading.subs(solved[0],simultaneous=True))==0
    assert sp.simplify(next_term.subs(solved[0],simultaneous=True))==0
    return solved[0]

def chiral_action_projection():
    delta,s=sp.symbols('Delta_V s',real=True)
    singlet=sp.Matrix([sp.sqrt(2),1])/sp.sqrt(3)
    K=(M0**2-s)*sp.eye(2)+delta*singlet*singlet.T
    source=FV*sp.Matrix([1/(3*sp.sqrt(2)),-sp.Rational(1,3)])
    response=sp.simplify(K.inv()*source)
    show('Chiral neutral isoscalar inverse operator',K)
    show('Chiral electromagnetic response: singlet mass drops out',response)
    assert sp.simplify(response-source/(M0**2-s))==sp.zeros(2,1)
    assert (singlet.T*source)[0]==0

def tff_coefficients():
    # Independent scalar flavor-pair traces, as in stage-5 preflight.
    # Derive Laurent coefficients by inserting spectral moment expansions.
    a,b=sp.symbols('alpha beta',positive=True)
    mp=sp.symbols('m_P',real=True)
    qQQ,qSQ,qJQ,q0Q,qSS,qS0,qJ0=sp.symbols('q_QQ q_SQ q_JQ q_0Q q_SS q_S0 q_J0')
    hSQ,h0Q,hS0=sp.symbols('h_SQ h_0Q h_S0')
    A=-c[0]+c[1]+c[4]-2*c[5]; B=c[0]-c[1]+c[4]
    C=lambda qt,ht,s,t: qt*((c[0]+c[1]-c[4])*mp**2+A*s+B*t)+4*c[2]*ht
    H=lambda qt,ht,s,t: qt*((d[0]-d[2])*mp**2+d[2]*(s+t))+4*d[1]*ht
    anomaly=sp.sqrt(2)*NC*qQQ/(8*sp.pi**2*F)
    # q(W(-x),B)=q(S,B)/x-q(J,B)/x^2+... . Higher
    # h moments cannot affect coefficients through 1/x and are omitted.
    single=anomaly+4/(F*MV)*(C(qSQ/x-qJQ/x**2,hSQ/x,-x,0)
                             +C(q0Q,h0Q,0,-x))
    single-=4*sp.sqrt(2)/F*H(qS0/x-qJ0/x**2,hS0/x,-x,0)
    single_coeff=[sp.expand(single).coeff(x,k) for k in (1,0,-1)]
    # For s1=-alpha*x,s2=-beta*x, C and W expansions are bilinear.
    double=anomaly+4/(F*MV)*(C(qSQ/(a*x)-qJQ/(a*x)**2,hSQ/(a*x),-a*x,-b*x)
                             +C(qSQ/(b*x)-qJQ/(b*x)**2,hSQ/(b*x),-b*x,-a*x))
    hSS=sp.symbols('h_SS')
    double-=4*sp.sqrt(2)/F*H(qSS/(a*b*x**2),hSS/(a*b*x**2),-a*x,-b*x)
    double_coeff=[sp.expand(double).coeff(x,k) for k in (0,-1)]
    assert sp.simplify((double_coeff[0]-double_coeff[0].subs({a:1,b:1})).subs(c[1],c[0]+c[4]))==0
    for label,value in zip(('single growing','single constant','single inverse','double constant','double inverse'),single_coeff+double_coeff):
        show(label,sp.factor(value))
    return dict(single=single_coeff,double=double_coeff,
                inputs=[mp,qQQ,qSQ,qJQ,q0Q,qSS,qS0,qJ0,hSQ,h0Q,hS0,a,b])

def common_limit(ope):
    import derive_radiative_tff_amplitude as old
    _,_,wzw,vjp,vvp=old.derive()
    mapping={old.F:F,old.FV:FV,old.MV:MV,old.NC:NC,old.m:0,old.FH:0}
    mapping.update(dict(zip(old.c,c))); mapping.update(dict(zip(old.d,d)))
    form=(wzw+vjp+vvp).subs(mapping,simultaneous=True)
    single=form.subs({old.s1:-x,old.s2:0})
    double=form.subs({old.s1:-x,old.s2:-x})
    s_a=sp.simplify(sp.limit(single/x,x,sp.oo))
    s_b=sp.simplify(sp.limit(single-s_a*x,x,sp.oo))
    db=sp.simplify(sp.limit(double,x,sp.oo))
    solution=sp.solve([s_a,s_b,db],[c[1],c[5],d[2]],dict=True)[0]
    show('Common-pole massless TFF: simultaneous single/double vanishing solution',sp.Dict(solution))
    ds=sp.factor(sp.limit(x*double.subs(solution,simultaneous=True),x,sp.oo))
    show('Double inverse coefficient after simultaneous vanishing',ds)
    show('Difference from the pion symmetric OPE target, positive anomaly phase',sp.simplify(ds-2*F/3))
    assert ds==0
    old_ope={key:value.subs(M0,MV) for key,value in ope.items()}
    mismatch=sp.factor(sp.limit(single.subs(old_ope,simultaneous=True),x,sp.oo))
    show('Single constant under chiral-octet OPE (M0=MV)',mismatch)
    assert sp.simplify(mismatch-sp.sqrt(2)*old.FQ*F/(2*MV**2))==0
    return {'common_double_inverse_after_vanishing':str(ds),'ope_single_constant_latex':sp.latex(mismatch)}

def source_matching():
    # Actual diagonal electromagnetic source from the declared action.
    xl,xs,zv,yv,Lambda=sp.symbols('x_l x_s z_V y_V Lambda',real=True)
    Q=sp.diag(sp.Rational(2,3),-sp.Rational(1,3),-sp.Rational(1,3))
    X=sp.diag(xl,xl,xs)
    S=FV*Q+zv*(X*Q+Q*X)+yv*sp.trace(X*Q)*sp.eye(3)
    show('UV source moment S from action',S)
    solved=sp.solve(list((S-Lambda*Q).diagonal()),[zv,yv,Lambda],dict=True)
    assert len(solved)==1
    show('S proportional to Q, nondegenerate spurion x_l != x_s',sp.Dict(solved[0]))
    assert solved[0][zv]==0 and solved[0][yv]==0
    show('SU(3)-symmetric spurion limit: double trace vanishes',S.subs(xs,xl))

def numeric(coefficients):
    sectors=stage4.module('vector_sector'); alda=stage4.module('alda_neutral_inputs')
    api=stage4.module('vector_sector_amplitudes'); Q=stage4.module('radiative_sm').Q
    inputs=alda.alda_neutral_inputs(mass_scheme='physical',F=.093,fpi=.093)
    compiled=sp.lambdify([*theta,*coefficients['inputs'],F,MV,NC],coefficients['single']+coefficients['double'],'numpy')
    previous=json.loads((ROOT/'data/vector_sector/stage4_rank_results.json').read_text())
    for path,digest in previous['sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    output=[]
    for case in previous['scenarios']:
        selection={k:case[k] for k in ('rho_input','kstar_input','phi_input','branch','source_hypothesis')}
        sector=sectors.VectorSector.from_stage3(**selection)
        S=sum((r*t for r,t in zip(sector.residues,sector.generators)),np.zeros((3,3),complex))
        J=sum((r*m*m*t for r,m,t in zip(sector.residues,sector.masses,sector.generators)),np.zeros((3,3),complex))
        W0=sector.response(0)
        rows=[]; offsets=[]; labels=[]; moments={}
        for name,entry in inputs.items():
            model=api.VectorAmplitudeInputs(.093,.770,entry['mp'],entry['TP'],entry['chi'],sector)
            q=lambda A,B: float(api.flavor(model,A,B)[0].real)
            h=lambda A,B: float(api.flavor(model,A,B)[1].real)
            args=[entry['mp'],q(Q,Q),q(S,Q),q(J,Q),q(W0,Q),q(S,S),q(S,W0),q(J,W0),
                  h(S,Q),h(W0,Q),h(S,W0),1.,1.,.093,.770,3]
            base=np.asarray(compiled(*([0]*11),*args),dtype=float)
            jac=np.column_stack([np.asarray(compiled(*v,*args),dtype=float)-base for v in np.eye(11)])
            # Necessary zero coefficients: single growing/constant and double constant.
            for index in (0,1,3):
                rows.append(jac[index]); offsets.append(base[index]); labels.append(f'{name}:{index}')
            moments[name]={'q_QQ':q(Q,Q),'q_SQ':q(S,Q),'ratio_qSQ_qQQ':q(S,Q)/q(Q,Q)}
            # Independently compare Laurent coefficients to the actual kernel at large x.
            arbitrary=np.arange(1,12)*.0007
            no_growth=arbitrary.copy(); no_growth[1]=no_growth[0]+no_growth[4]
            for vector in (arbitrary,no_growth):
                coupling=stage4.odd.OddParityCouplings(stage4.odd.VJPCouplings(*vector[:7]),stage4.odd.VVPCouplings(*vector[7:]))
                values=np.asarray(compiled(*vector,*args),dtype=float)
                for scale in (1e5,2e5):
                    actual=sum(api.tff_components(coupling,model,-scale,0)).real
                    predicted=values[0]*scale+values[1]+values[2]/scale
                    assert abs(actual-predicted)<1e-8*max(1,abs(actual))
                    for ratio in (1.,2.):
                        ray_args=list(args); ray_args[12]=ratio
                        ray_values=np.asarray(compiled(*vector,*ray_args),dtype=float)
                        actual=sum(api.tff_components(coupling,model,-scale,-ratio*scale)).real
                        predicted=ray_values[3]+ray_values[4]/scale
                        assert abs(actual-predicted)<1e-8*max(1,abs(actual))
        # Row normalization gives dimensionless equation residuals, NOT chi-square.
        matrix=np.asarray(rows); target=-np.asarray(offsets)
        scales=np.maximum(np.linalg.norm(matrix,axis=1),np.abs(target))
        A=matrix/scales[:,None]; b=target/scales
        solution=np.linalg.lstsq(A,b,rcond=1e-11)[0]
        residual=A@solution-b
        # Double-only equations have one common A-combination after B=0.
        ds=np.array([moments[n]['q_SQ'] for n in inputs]); dq=np.array([moments[n]['q_QQ'] for n in inputs])
        proportional=ds@dq/(dq@dq)
        source_residual=np.linalg.norm(ds-proportional*dq)/np.linalg.norm(ds)
        output.append({**selection,'equation_labels':labels,'matrix':matrix.tolist(),'target':target.tolist(),
            'rank':int(np.linalg.matrix_rank(A,tol=1e-11)),
            'normalized_residual_norm':float(np.linalg.norm(residual)),
            'zero_constraints_consistent_at_1e_minus_9':bool(np.linalg.norm(residual)<1e-9),
            'source_proportionality_relative_residual':float(source_residual),'moments':moments,
            'coefficient_solution_minimum_norm':solution.tolist()})
    return output

def main():
    transcript=ROOT/'data/vector_sector/short_distance_derivation.md'
    with transcript.open('w') as handle,contextlib.redirect_stdout(handle):
        print('# Executed short-distance algebra\n\nTree level; GeV. Scalar traces have the stage-4 definitions.\n')
        chiral_action_projection(); ope=octet_ope(); coefficients=tff_coefficients(); control=common_limit(ope); source_matching()
    scenarios=numeric(coefficients)
    paths=['scripts/match_vector_short_distance.py','scripts/derive_radiative_tff_amplitude.py',
           'data/vector_sector/stage4_rank_results.json','memory/sources/Roig:2013baa/AcceptedVersion.tex',
           'memory/sources/Ruiz-Femenia:2003jdx/vvp.tex']
    paths.append('memory/sources/Knecht:2001qf/pionpole.tex')
    result={'scope':'chiral-octet correlator matching and physical TFF necessary zero constraints; no physical fit',
        'common_control':control,'scenarios':scenarios,'sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        'coupling_order':stage4.ORDER,'ope_solution_latex':sp.latex(sp.Dict(ope))}
    target=ROOT/'data/vector_sector/short_distance_results.json'
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(transcript.relative_to(ROOT)); print(target.relative_to(ROOT))

if __name__=='__main__':main()
