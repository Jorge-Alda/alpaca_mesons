"""Neutral linear U(3) EOM/source-contact audit and baseline-trace contractions.

Run .venv/bin/python scripts/audit_odd_trace_sources.py
F=f_pi convention; independent chi_P, axial and topological sources.
No EOM-minimal full-action count or physical fit is inferred.
"""
from pathlib import Path
import contextlib
import hashlib
import json
import sympy as sp
import derive_radiative_tff_amplitude as t

ROOT=Path(__file__).resolve().parents[1]

def show(label,value):print(label+':\n$$\n'+sp.latex(value)+'\n$$\n')

def source_eom():
    F,m0=sp.symbols('F m_0',positive=True)
    theta,z=sp.symbols('theta_TC z',real=True)
    phi=sp.Matrix(sp.symbols('phi_3 phi_8 phi_0',real=True))
    T=(sp.diag(1,-1,0)/sp.sqrt(2),sp.diag(1,1,-2)/sp.sqrt(6),sp.eye(3)/sp.sqrt(3))
    Phi=sum((v*M for v,M in zip(phi,T)),sp.zeros(3))
    chi=sp.diag(*sp.symbols('chi_u chi_d chi_s',real=True))
    cp=sp.diag(*sp.symbols('chi_Pu chi_Pd chi_Ps',real=True))
    da=sp.diag(*sp.symbols('div_a_u div_a_d div_a_s',real=True))
    uu=sp.eye(3)+sp.I*z*Phi/(sp.sqrt(2)*F)-z**2*Phi**2/(4*F**2)
    ud=sp.eye(3)-sp.I*z*Phi/(sp.sqrt(2)*F)-z**2*Phi**2/(4*F**2)
    source=chi+sp.I*z*cp;dagger=chi-sp.I*z*cp
    chip=sp.expand(ud*source*ud+uu*dagger*uu)
    chim=sp.expand(ud*source*ud-uu*dagger*uu)
    Y=sp.simplify(sp.diff(sp.I*chim,z).subs(z,0))
    # Linear X_TC follows Alda's invariant; theta_TC is an independent phase
    # source, specializing to c_G a/f_a only after variation.
    XTC=-sp.sqrt(2)*sp.trace(Phi)/F-2*theta
    mass=(F**2*sp.trace(chip)/4).expand().coeff(z,2)
    anomaly=-F**2*m0**2*XTC**2/12  # N_f=3, Alda LO coefficient.
    # The cross term from F^2 Tr(u_mu u^mu)/4 gives, after IBP,
    # +sqrt(2) F Tr(Phi div_a), since u_mu=2a_mu-sqrt(2)/F d_mu Phi.
    action=sp.expand(mass+anomaly+sp.sqrt(2)*F*sp.trace(Phi*da))
    box_components=sp.Matrix([sp.diff(action,v) for v in phi])
    boxPhi=sum((v*M for v,M in zip(box_components,T)),sp.zeros(3))
    divu=sp.simplify(2*da-sp.sqrt(2)*boxPhi/F)
    residual=sp.simplify(divu-Y/2+m0**2*XTC*sp.eye(3)/3)
    assert residual==sp.zeros(3)
    show('Neutral canonical Phi_R in the 3,8,0 basis',Phi)
    show('Linear independent pseudoscalar block Y=i chi_minus',Y)
    show('Linear Alda X_TC; theta_TC independent before ALP substitution',XTC)
    show('Executed quadratic meson potential and integrated axial-source cross term',action)
    show('Euler-Lagrange box Phi_R',boxPhi)
    show('Derived divergence u_mu',divu)
    show('EOM identity residual: divu-Y/2+m0^2 X_TC I/3',residual)
    # Four JP divergence tensors and their four existing Y tensors.
    # Induced topological traces span Tr(AB) and Tr(A)Tr(B).
    a,b=sp.symbols('Tr_AB Tr_A_Tr_B')
    qI=sp.Matrix([2*a,2*b,2*b,6*a])
    contact=sp.Matrix([[sp.diff(expr,p) for expr in qI] for p in (a,b)])
    null=contact.nullspace()
    show('JP divergence-to-topological-contact matrix',contact)
    show('Combinations canceling the induced JP contact',sp.Matrix.hstack(*null))
    vv=contact[:,[0,1,3]]
    show('Symmetric VV divergence-to-contact matrix',vv)
    show('Combinations canceling the induced VV contact',sp.Matrix.hstack(*vv.nullspace()))
    return {'scope':'neutral linear U(3) EOM, with independent sources and LO meson action',
        'JP_divergence_contact_rank':contact.rank(),'JP_relations_in_divergence_mass_subsector':len(null),
        'VV_divergence_contact_rank':vv.rank(),'VV_relations_in_divergence_mass_subsector':len(vv.nullspace()),
        'full_action_new_operator_count':None,
        'count_limit':'not a complete full-action quotient; vector EOM and nonlinear contacts are not reduced'}

def contractions():
    first=t.contractions(t.q1,t.q2,t.e1,t.e2)
    second=t.contractions(t.q2,t.q1,t.e2,t.e1)
    q,h,qa=sp.symbols('q_l h_l q_minus',real=True)
    C=t.reduce(-sum(t.c[i]*(-h if i==2 else q)*first[i] for i in (0,1,2,4,5,6)))
    H=t.reduce(-(t.d[0]*q*(first[0]+second[0])-t.d[1]*h*(first[2]+second[2])
                +t.d[2]*q*(first[4]+second[4])+t.d[3]*q*(first[6]+second[6]))/2)
    dm=sp.symbols('d_minus',real=True)
    # Antisymmetric trace placement on the c6 Lorentz pattern, with both
    # vector Wick assignments: F1-F2 changes sign on exchange.
    HA=t.reduce(-dm*qa*(first[5]-second[5])/2)
    show('C_l: normalized VJP trace-family contraction',C)
    show('H_l: normalized symmetric VVP trace-family contraction',H)
    show('H_minus: antisymmetric trace and Lorentz contraction',HA)
    assert sp.simplify(HA+dm*qa*(t.s1-t.s2))==0
    assert sp.simplify(HA-HA.xreplace({qa:-qa,t.s1:t.s2,t.s2:t.s1}))==0
    for ev,er in ((t.q1,t.e2),(t.e1,t.q2)):
        assert t.contractions(t.q1,t.q2,ev,er,False)==[0]*7
    # Generic scalar EM kernel-space count, not a full-action operator count.
    ctable=sp.symbols('C0:24');dtable=sp.symbols('D0:12')
    c4=sp.symbols('c4'); fk=sp.symbols('F_K')
    qj=sp.symbols('qJ0:4');hj=sp.symbols('hJ0:4')
    qv=sp.symbols('qV0:3');hv=sp.symbols('hV0:3')
    values=[];symbols=[]
    for ell in range(4):
        cs=ctable[6*ell:6*(ell+1)]
        repl=dict(zip((t.c[0],t.c[1],t.c[2],t.c[4],t.c[5],t.c[6]),cs))
        values.append(C.subs(repl).subs({q:0 if ell==2 else qj[ell],h:0 if ell==2 else hj[ell]}))
        symbols.extend(cs)
    for ell in range(3):
        ds=dtable[4*ell:4*(ell+1)]
        values.append(H.subs(dict(zip(t.d,ds))).subs({q:qv[ell],h:hv[ell]}))
        symbols.extend(ds)
    symbols.extend((c4,dm))
    generic=sp.expand(sum(values)-4*c4*fk+HA)
    monomials=(t.m,t.s1,t.s2,*qj,*hj,*qv,*hv,fk,qa)
    polynomial=sp.Poly(generic,*monomials)
    design=sp.Matrix([[sp.diff(coef,p) for p in symbols] for coef in polynomial.coeffs()])
    baseline_indices=list(range(6))+list(range(24,28))+[36]
    rank=design.rank();baseline_rank=design[:,baseline_indices].rank()
    assert rank==23 and baseline_rank==8
    show('Generic scalar EM design rank (includes all virtualities, independent flavor traces)',sp.Integer(rank))
    show('Single-trace scalar EM design rank',sp.Integer(baseline_rank))
    return {'C_l_latex':sp.latex(C),'H_l_latex':sp.latex(H),'H_minus_latex':sp.latex(HA),
        'scalar_EM_kernel_rank':rank,'baseline_scalar_EM_kernel_rank':baseline_rank,
        'additional_scalar_EM_kernel_directions':rank-baseline_rank,
        'interpretation':'formal flavor/virtuality polynomial directions, NOT full-action operators or measurement rank'}

def main():
    path=ROOT/'data/vector_sector/odd_trace_source_audit.md'
    with path.open('w') as handle,contextlib.redirect_stdout(handle):
        print('# Executed neutral source EOM and trace-family contraction audit\n')
        eom=source_eom();kernels=contractions()
    inputs=['scripts/audit_odd_trace_sources.py','scripts/derive_radiative_tff_amplitude.py',
            'memory/sources/Alda:2026nrd/main_arXiv.tex']
    output=ROOT/'data/vector_sector/odd_trace_source_audit.json'
    output.write_text(json.dumps({'source_eom':eom,'normalized_kernels':kernels,
        'sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs}},indent=2)+'\n')
    print(path.relative_to(ROOT));print(output.relative_to(ROOT))

if __name__=='__main__':main()
