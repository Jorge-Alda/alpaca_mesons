"""Stage-5 preflight, without inherited OPE rows or a physical fit.

Run: .venv/bin/python scripts/diagnose_vector_stage5.py
Natural units, GeV; tree-level single-nonet action, fixed central inputs.
"""
from pathlib import Path
import hashlib
import json
import sympy as sp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def asymptotics():
    # q_i=q(T_i,Q), h_i=h(T_i,Q), q_ij/h_ij likewise.
    # Independent scalar trace coefficients preserve the input matrix algebra.
    x = sp.symbols('Q2', positive=True)
    F, MV, mp, NC, qQQ = sp.symbols('F M_V m_P N_C q_QQ', nonzero=True)
    c1,c2,c3,c5,c6,d1,d2,d3 = sp.symbols('c1 c2 c3 c5 c6 d1 d2 d3')
    M = sp.symbols('M_rho M_omega M_phi', positive=True)
    R = sp.symbols('R_rho R_omega R_phi', real=True)
    qi = sp.symbols('q_rho q_omega q_phi')
    hi = sp.symbols('h_rho h_omega h_phi')
    qij = sp.Matrix(3,3,lambda i,j: sp.Symbol(f'q_{i}{j}'))
    hij = sp.Matrix(3,3,lambda i,j: sp.Symbol(f'h_{i}{j}'))
    def C(q,h,s,t):
        return q*((c1+c2-c5)*mp**2+(-c1+c2+c5-2*c6)*s+(c1-c2+c5)*t)+4*c3*h
    def H(q,h,s,t):
        return q*((d1-d3)*mp**2+d3*(s+t))+4*d2*h
    wzw=sp.sqrt(2)*NC*qQQ/(8*sp.pi**2*F)
    vjp=4/(F*MV)*sum(R[i]*(C(qi[i],hi[i],-x,0)/(M[i]**2+x)
                         +C(qi[i],hi[i],0,-x)/M[i]**2) for i in range(3))
    vvp=-4*sp.sqrt(2)/F*sum(R[i]*R[j]*H(qij[i,j],hij[i,j],-x,0)
                   /((M[i]**2+x)*M[j]**2) for i in range(3) for j in range(3))
    expression=wzw+vjp+vvp
    slope=sp.simplify(sp.limit(expression/x,x,sp.oo))
    constant=sp.simplify(sp.limit(expression-slope*x,x,sp.oo))
    inverse=sp.simplify(sp.limit(x*(expression-slope*x-constant),x,sp.oo))
    # Compare coefficients extracted independently by a Taylor expansion at 1/Q2=0.
    z=sp.symbols('z',positive=True)
    series=sp.series((expression*x**-1).subs(x,1/z),z,0,3).removeO().expand()
    assert sp.simplify(series.coeff(z,0)-slope)==0
    assert sp.simplify(series.coeff(z,1)-constant)==0
    assert sp.simplify(series.coeff(z,2)-inverse)==0
    return {name:sp.latex(value) for name,value in
            [('input_F',expression),('Q2_coefficient',slope),
             ('constant_coefficient',constant),('inverse_Q2_coefficient',inverse)]}


def rank(rows):
    matrix=np.asarray(rows,dtype=float)
    norms=np.linalg.norm(matrix,axis=1)
    normalized=matrix[norms>0]/norms[norms>0,None]
    singular=np.linalg.svd(normalized,compute_uv=False)
    return {'rank':int(np.count_nonzero(singular>1e-11)),
            'singular_values':singular.tolist()}


def main():
    path=ROOT/'data/vector_sector/stage4_rank_results.json'
    previous=json.loads(path.read_text())
    # Refuse stale amplitude results rather than silently reusing them.
    for name,digest in previous['sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    scenarios=[]
    for item in previous['scenarios']:
        rad=item['radiative_rows']; real=item['real_tff_rows']; virtual=item['virtual_tff_rows']
        scenarios.append({**{k:item[k] for k in ('rho_input','kstar_input','phi_input','branch','source_hypothesis')},
            'radiative_only':rank(rad),'radiative_real':rank(rad+real),
            'radiative_real_virtual':rank(rad+real+virtual)})
    result={'scope':'preflight only: amplitude sensitivities are not measurement Jacobians',
            'ope_rows_used':0,'self_energies':'deferred by user: tree-level ALP accuracy',
            'asymptotics':asymptotics(),'scenarios':scenarios,
            'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in (path,Path(__file__).resolve())}}
    target=ROOT/'data/vector_sector/stage5_preflight.json'
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(target.relative_to(ROOT))


if __name__=='__main__':
    main()
