"""Lift the leading divergence relation through the nonlinear f-plus Bianchi.

Run .venv/bin/python scripts/lift_odd_leading_covariant.py.
Flavor words label abstract 3x3 matrices with Lorentz-component indices;
cyclic trace canonicalization reuses the documented finite-orbit routine.
No generic flavor-component replacement is made.
"""
from pathlib import Path
import json
import contextlib
import hashlib
import sympy as sp
import enumerate_odd_breaking_vertices as b
from reduce_odd_goldstone_eom import image, rank

ROOT=Path(__file__).resolve().parents[1]

def show(label,expr):
    print(label+':\n$$\n'+sp.latex(expr)+'\n$$\n')

def antisymmetric_label(name,i,j):
    if i==j:return 0,None
    return (1 if i<j else -1),f'{name}{min(i,j)}{max(i,j)}'

def bianchi_remainder(family):
    # N=epsilon_(i j k l) A^(i a) u_a nabla^l f_plus^(j k).
    # The full Bianchi gives N=-i/2 epsilon A u [u^l,f_minus^(j k)].
    out={}
    for (i,j,k,l),eps in b.EPS:
        sf,flabel=antisymmetric_label('fminus',j,k)
        for a in range(4):
            sa,alabel=antisymmetric_label('V',i,a)
            if not sa or not sf:continue
            for traces,coefficient in family.items():
                for sign,replacement in ((1,(f'u{l}',flabel)),(-1,(flabel,f'u{l}'))):
                    words=[]
                    for word in traces:
                        expanded=[]
                        for block in word:
                            expanded.extend(replacement if block=='B' else
                                            (alabel,) if block=='A' else
                                            ('chi_plus',) if block=='X' else (f'u{a}',))
                        words.append(tuple(expanded))
                    key=b.canonical(words)
                    out[key]=out.get(key,0)-sp.I*sp.Rational(1,2)*eps*b.G[a]*sa*sf*coefficient*sign
    return {key:sp.simplify(value) for key,value in out.items() if sp.simplify(value)!=0}

def main():
    data=json.loads((ROOT/'data/vector_sector/odd_vertex_enumeration.json').read_text())
    leading=data['lorentz']['VJP']
    generic=[]
    for label in leading['basis_labels']:
        D=b.qa if label['derivative_on']=='A' else b.qb
        generic.append(b.contraction(tuple(label['metric_pair']),b.AT,b.BT,b.u,D))
    N=b.contraction((1,4),b.AT,b.BT,b.u,b.qb)
    qU=tuple(-b.qa[i]-b.qb[i] for i in range(4))
    D=b.contraction((4,5),b.AT,b.BT,b.u,qU)
    pivots,coords=b.reduction(generic+[N,D])
    assert pivots==tuple(range(6))
    Dcoordinates=coords[:,6]
    # Independent confirmation in the closed source projection.
    closed=[expr.subs(dict(zip(b.Bvar,[b.BF[i,j] for i,j in b.AB])),simultaneous=True) for expr in generic+[N,D]]
    assert sp.expand(closed[5])==0
    assert sp.expand(closed[6]-sum(Dcoordinates[i]*closed[i] for i in range(5)))==0
    families=b.trace_catalog(('A','B','U'))
    contacts=[image(f,None) for f in families]
    remainders=[bianchi_remainder(f) for f in families]
    rr,_,_,R=rank(remainders)
    cr,_,_,C=rank(contacts)
    # The generic divergence relation is D=sum c_i O_i + c_N N.
    augmented=[{**{('contact',key):value for key,value in c.items()},
                **{('bianchi',key):Dcoordinates[5]*value for key,value in r.items()}}
               for c,r in zip(contacts,remainders)]
    ar,null,keys,M=rank(augmented)
    assert all(M*v==sp.zeros(M.rows,1) for v in null)
    vvfamilies=b.flavor_exchange_basis(families,1)
    vr,vnull,_,VC=rank([image(f,None) for f in vvfamilies])
    output={'leading_VJP':{'generic_divergence_to_selected_plus_bianchi':[str(x) for x in Dcoordinates],
                          'contact_rank':cr,'bianchi_remainder_rank':rr,
                          'combined_remainder_rank':ar,'exact_contact_free_relations':len(null),
                          'null_vectors':[[str(x) for x in v] for v in null]},
            'leading_VVP':{'contact_rank':vr,'exact_contact_free_relations':len(vnull),
                          'null_vectors':[[str(x) for x in v] for v in vnull]},
            'scope':'leading divergence/mass subsector only; fixed explicit-vector convention',
            'full_new_operator_count':None}
    transcript=ROOT/'data/vector_sector/odd_leading_covariant_lift.md'
    with transcript.open('w') as handle,contextlib.redirect_stdout(handle):
        print('# Leading covariant divergence lift\n')
        show('Generic divergence coordinates in five saved tensors plus N',Dcoordinates)
        show('Nonlinear Bianchi image rank',sp.Integer(rr))
        show('Topological contact image',C)
        show('Combined contact and Bianchi image rank',sp.Integer(ar))
        show('Exact leading VJP null vectors',sp.Matrix.hstack(*null) if null else sp.zeros(4,0))
        show('Exact leading VVP null vectors',sp.Matrix.hstack(*vnull))
    paths=['scripts/lift_odd_leading_covariant.py','data/vector_sector/odd_vertex_enumeration.json',
           'scripts/derive_odd_covariant_closure.py','scripts/reduce_odd_goldstone_eom.py']
    output['sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    result=ROOT/'data/vector_sector/odd_leading_covariant_lift.json'
    result.write_text(json.dumps(output,indent=2)+'\n')
    print(transcript.relative_to(ROOT));print(result.relative_to(ROOT))

if __name__=='__main__':main()
