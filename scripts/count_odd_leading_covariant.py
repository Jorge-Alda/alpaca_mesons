"""Exact complete leading covariant seed quotient, with original-action embedding.

Run .venv/bin/python scripts/count_odd_leading_covariant.py.
Uses formal derivative symbols (no Fourier factors of i), the source-preserving
LO Goldstone EOM, and no resonance EOM. Nonlinear Bianchi images are evaluated
as abstract flavor trace words by lift_odd_leading_covariant.py.
"""
from pathlib import Path
import contextlib
import hashlib
import json
import sympy as sp
import enumerate_odd_breaking_vertices as b
from reduce_odd_goldstone_eom import image
from lift_odd_leading_covariant import bianchi_remainder

ROOT=Path(__file__).resolve().parents[1]

def show(label,expr):
    print(label+':\n$$\n'+sp.latex(expr)+'\n$$\n')

def tensor_column(polynomial,family,sector):
    out={}
    polynomial=sp.Poly(sp.expand(polynomial),*VARIABLES)
    for monomial,coefficient in polynomial.terms():
        for word,value in family.items():
            out[(sector,word,monomial)]=coefficient*value
    return out

def add(*columns):
    out={}
    for column in columns:
        for key,value in column.items():out[key]=out.get(key,0)+value
    return {key:sp.expand(value) for key,value in out.items() if sp.expand(value)!=0}

def scaled(column,coefficient):return {key:value*coefficient for key,value in column.items()}

Z=sp.Symbol('Z_TC')
VARIABLES=(*b.Avar,*b.Bvar,*b.u,*b.qa,*b.qb,*b.ev,Z)

def main():
    source=ROOT/'data/vector_sector/odd_vertex_enumeration.json'
    data=json.loads(source.read_text())
    qU=tuple(-b.qa[i]-b.qb[i] for i in range(4))
    generic=[]
    for label in data['lorentz']['VJP']['basis_labels']:
        generic.append(b.contraction(tuple(label['metric_pair']),b.AT,b.BT,b.u,
                      b.qa if label['derivative_on']=='A' else b.qb))
    N=b.contraction((1,4),b.AT,b.BT,b.u,b.qb)
    generic.append(N)
    D=b.contraction((4,5),b.AT,b.BT,b.u,qU)
    pivots,coordinates=b.reduction(generic+[D])
    assert pivots==tuple(range(6))
    dc=coordinates[:,6]
    substitution=dict(zip(b.Bvar,[b.BF[i,j] for i,j in b.AB]))
    closed=[sp.expand(p.subs(substitution,simultaneous=True)) for p in generic]
    closedD=sp.expand(D.subs(substitution,simultaneous=True))
    wedge=sum(sign*b.AT[i,j]*b.BT[k,l] for (i,j,k,l),sign in b.EPS)
    closedWedge=sp.expand(wedge.subs(substitution,simultaneous=True))
    families=b.trace_catalog(('A','B','U'))
    rem=[bianchi_remainder(f) for f in families]
    columns=[];identifiers=[]
    def append(identifier,column):identifiers.append(identifier);columns.append(column)
    for f,family in enumerate(families):
        for l,polynomial in enumerate(closed):
            nonlinear={('Bianchi',word,()):value for word,value in rem[f].items()} if l==5 else {}
            append(f'VJP:F{f}:L{l}',add(tensor_column(polynomial,family,'JP'),nonlinear))
    for f,family in enumerate(families):
        nonlinear={('Bianchi',word,()):2*dc[5]*value for word,value in rem[f].items()}
        append(f'VJP:Y:F{f}',add(tensor_column(2*closedD,family,'JP'),
              tensor_column(2*Z*closedWedge,image(family,None),'JP'),nonlinear))
    vvfamilies=b.flavor_exchange_basis(families,1)
    exchange={**dict(zip(b.Avar,b.Bvar)),**dict(zip(b.Bvar,b.Avar)),
              **dict(zip(b.qa,b.qb)),**dict(zip(b.qb,b.qa))}
    vv=[]
    for entry in data['lorentz']['VVP']['symmetric']['basis_labels']:
        label=entry['parent']
        polynomial=b.contraction(tuple(label['metric_pair']),b.AT,b.BT,b.u,
                                 b.qa if label['derivative_on']=='A' else b.qb)
        vv.append(sp.expand(polynomial+polynomial.xreplace(exchange)))
    # Antisymmetric flavor and Lorentz sectors are separately Bose-even.
    afamilies=b.flavor_exchange_basis(families,-1)
    avv=[]
    for entry in data['lorentz']['VVP']['antisymmetric']['basis_labels']:
        label=entry['parent']
        polynomial=b.contraction(tuple(label['metric_pair']),b.AT,b.BT,b.u,
                                 b.qa if label['derivative_on']=='A' else b.qb)
        avv.append(sp.expand(polynomial-polynomial.xreplace(exchange)))
    for group,flavors,polynomials in (('S',vvfamilies,vv),('A',afamilies,avv)):
        for f,family in enumerate(flavors):
            for l,polynomial in enumerate(polynomials):
                append(f'VVP:{group}:F{f}:L{l}',tensor_column(polynomial,family,'VV'))
    for f,family in enumerate(vvfamilies):
        append(f'VVP:Y:F{f}',add(tensor_column(2*D,family,'VV'),
                               tensor_column(2*Z*wedge,image(family,None),'VV')))
    # f-minus chi-plus is in a distinct independent source sector here.
    append('VJP:fminus', {('fminus_X',(),()):sp.Integer(1)})
    assert len(columns)==44
    keys=sorted(set().union(*(set(c) for c in columns)),key=repr)
    matrix=sp.SparseMatrix(len(keys),len(columns),
                          {(i,j):c[key] for j,c in enumerate(columns) for i,key in enumerate(keys) if key in c})
    reduced,pivots=matrix.rref()
    null=matrix.nullspace()
    assert matrix[:,list(pivots)]*reduced[:len(pivots),:]==matrix
    # Embed the actual five original derivative seeds before their closed-
    # current projection, rather than assume vertex recovery implies action recovery.
    original=[b.contraction(pair,b.AT,b.BT,b.u,momentum) for pair,momentum in
              (((3,5),qU),((1,5),qU),((3,5),b.qa),((1,5),b.qa),((3,4),b.qa))]
    _,seedcoords=b.reduction(generic+original)
    seedcoords=seedcoords[:,6:]
    single=next(i for i,f in enumerate(families) if len(next(iter(f)))==1)
    vsingle=next(i for i,f in enumerate(vvfamilies) if len(next(iter(f)))==1)
    embedded=[]
    for j in range(5):
        embedded.append(add(*(scaled(columns[single*6+l],seedcoords[l,j]) for l in range(6))))
    embedded.append(columns[identifiers.index(f'VJP:Y:F{single}')])
    embedded.append(columns[-1])
    # Original d1,d3,d4 span precisely the symmetric single-trace Lorentz space.
    original_vv=[b.contraction(pair,b.AT,b.BT,b.u,momentum) for pair,momentum in
                 (((3,5),qU),((3,5),b.qa),((3,4),b.qa))]
    original_vv=[sp.expand(p+p.xreplace(exchange)) for p in original_vv]
    vp,vc=b.reduction(vv+original_vv)
    assert vp==tuple(range(3))
    vc=vc[:,3:]
    for j in range(3):
        embedded.append(add(*(scaled(columns[identifiers.index(f'VVP:S:F{vsingle}:L{l}')],vc[l,j])
                              for l in range(3))))
    embedded.append(columns[identifiers.index(f'VVP:Y:F{vsingle}')])
    base=sp.Matrix([[c.get(key,0) for c in embedded] for key in keys])
    ordered=base.row_join(matrix)
    _,selected=ordered.rref()
    assert selected[:11]==tuple(range(11))
    assert len(selected)==len(pivots)
    result={'scope':'complete leading covariant seeds only; higher-derivative/one-X derivative extension not counted',
            'raw_Lorentz_flavor_basis':len(columns),'independent_leading_action_operators':len(pivots),
            'original_action_rank':base.rank(),'new_leading_independent_operators':len(pivots)-base.rank(),
            'relations':[{identifiers[i]:str(v) for i,v in enumerate(vector) if v!=0} for vector in null],
            'original_VJP_derivatives_to_generic_basis':[[str(v) for v in seedcoords[:,j]] for j in range(5)],
            'original_VVP_derivatives_to_symmetric_basis':[[str(v) for v in vc[:,j]] for j in range(3)],
            'basis_preserving_original_action':['original:'+name for name in
                ('c1','c2','c5','c6','c7','c3','c4','d1','d3','d4','d2')]+
               [identifiers[j-11] for j in selected[11:]],
            'row_count':len(keys),'global_extension_count':None}
    transcript=ROOT/'data/vector_sector/odd_leading_action_count.md'
    with transcript.open('w') as handle,contextlib.redirect_stdout(handle):
        print('# Exact complete leading covariant quotient\n')
        show('Five original VJP seeds in generic six-dimensional space',seedcoords)
        show('Generic divergence coordinates',dc)
        show('Leading raw, quotient and original ranks',sp.Matrix([len(columns),len(pivots),base.rank()]))
        show('Kernel relations',sp.Matrix.hstack(*null))
        show('Number of new leading operators',sp.Integer(len(pivots)-base.rank()))
    paths=['scripts/count_odd_leading_covariant.py','scripts/lift_odd_leading_covariant.py',
           'scripts/reduce_odd_goldstone_eom.py','scripts/enumerate_odd_breaking_vertices.py',
           'data/vector_sector/odd_vertex_enumeration.json',
           'memory/extracts/Alda:2026nrd/lo-topological-action.md',
           'memory/rcht/odd-parity/odd-covariant-closure.md']
    result['sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    output=ROOT/'data/vector_sector/odd_leading_action_count.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(transcript.relative_to(ROOT));print(output.relative_to(ROOT))

if __name__=='__main__':main()
