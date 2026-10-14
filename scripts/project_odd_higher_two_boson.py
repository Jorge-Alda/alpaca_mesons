"""Exact stripped cubic two-boson projection; not normalized physical predictors.
Run: .venv/bin/python scripts/project_odd_higher_two_boson.py
Abstract noncommuting U(3) flavour matrices; constant scalar vacuum X.
Metric (+---), epsilon0123=+1, incoming momenta, canonical wedge tensors.
"""
from pathlib import Path
from functools import lru_cache
import contextlib
import json
import sympy as sp
import enumerate_odd_breaking_vertices as b
from derive_radiative_sm_amplitude import Q,CHI,RHO,OMEGA,PHI

ROOT=Path(__file__).resolve().parents[1]
s,t,h=sp.symbols('s t h')
AP=sp.Matrix(4,4,lambda i,j:b.qa[i]*b.ea[j]-b.qa[j]*b.ea[i])
BP=sp.Matrix(4,4,lambda i,j:b.qb[i]*b.ev[j]-b.qb[j]*b.ev[i])
p=tuple(-b.qa[i]-b.qb[i] for i in range(4))
E=sp.expand(sum(sign*b.qa[ix[0]]*b.ea[ix[1]]*b.qb[ix[2]]*b.ev[ix[3]] for ix,sign in b.EPS))
dot=lambda a,c:sum(b.G[i]*a[i]*c[i] for i in range(4))
refs=[sp.expand(E*dot(b.qa,b.qa)),sp.expand(E*dot(b.qb,b.qb)),sp.expand(E*dot(b.qa,b.qb))]
exchange={**dict(zip(b.qa,b.qb)),**dict(zip(b.qb,b.qa)),**dict(zip(b.ea,b.ev)),**dict(zip(b.ev,b.ea))}

def show(label,e):print(label+':\n$$\n'+sp.latex(e)+'\n$$\n')

@lru_cache(None)
def project(mode,where,pair,group):
    if mode=='X':
        D={'A':b.qa,'B':b.qb,'X':(0,0,0,0)}[where]
        raw=b.contraction(pair,AP,BP,p,D)
    else:
        U,D={'AA':(b.qa,b.qa),'AB':(b.qa,b.qb),'BB':(b.qb,b.qb)}[where]
        raw=b.contraction(pair,AP,BP,U,D)
    if group!='all':raw=sp.expand(raw+(1 if group=='symmetric' else -1)*raw.xreplace(exchange))
    pivots,coords=b.reduction(refs+[raw])
    assert pivots==(0,1,2), (mode,where,pair,group,pivots)
    c=coords[:,3]
    assert sp.expand(raw-sum(c[i]*refs[i] for i in range(3)))==0
    return tuple(c)

def family(entry):
    return {tuple(tuple(word) for word in k.split(';')):sp.Rational(v)
            for k,v in entry['flavor_trace_words'][0].items()}

def flavour_latex(words):
    mats={name:sp.MatrixSymbol(name,3,3) for name in 'ABUXY'}
    return sp.Mul(*(sp.Trace(sp.MatMul(*(mats[c] for c in word))) for word in words))

def neutral_projection(selected,images):
    """Restrict to the saved physical neutral generator span, not generic matrices."""
    generators=(RHO,OMEGA,PHI)
    def weight(prefix):
        return sum((c*T for c,T in zip(sp.symbols(prefix+'0:3'),generators)),sp.zeros(3))
    A,B,U,Y=(weight(prefix) for prefix in ('vA','vB','uA','yA'))
    commutators=[sp.simplify(W*Q-Q*W) for W in (U,Y)]
    assert all(C==sp.zeros(3) for C in commutators)
    polynomials=[]
    for entry in selected:
        matrices={'A':A,'B':Q if entry['block'].endswith('VJP') else B,'U':U,'Y':Y,'X':2*CHI}
        flavor=sum(value*sp.Mul(*(sp.trace(sp.MatMul(*(matrices[n] for n in word)).doit())
                       for word in words)) for words,value in family(entry).items())
        # Distinct VJP/VVP kernels and U/Y source weights must not be merged.
        sector=sp.Symbol('J' if entry['block'].endswith('VJP') else 'V')
        polynomials.append(sp.expand(sector*flavor*sp.sympify(images[entry['id']]['kinematic_polynomial'])))
    # All sqrt(2) factors come solely from the pinned neutral generators.
    # Scale rows over Q(sqrt(2)); no numerical generator or mass values are used.
    matrix=b.coefficient_matrix(polynomials)
    reduced,pivots=matrix.rref();rank=len(pivots)
    assert matrix[:,list(pivots)]*reduced[:rank,:]==matrix
    combinations=[{'representative':selected[p]['id'],'coefficients':{
        entry['id']:str(reduced[i,j]) for j,entry in enumerate(selected) if reduced[i,j]}}
        for i,p in enumerate(pivots)]
    blocks={}
    for block in sorted({e['block'] for e in selected}):
        indices=[i for i,e in enumerate(selected) if e['block']==block]
        blocks[block]={'rank':matrix[:,indices].rank(),
                       'zero_images':sum(polynomials[i]==0 for i in indices)}
    return {'rank':rank,'blocks':blocks,'independent_combinations':combinations,
            'neutral_source_charge_commutators':[sp.latex(C) for C in commutators],
            'scope':'neutral generator span, B=Q only for direct VJP; independent diagonal U,Y; symbolic quark masses; not a fit rank'}

def main():
    source=ROOT/'data/vector_sector/odd_covariant_action_basis.json'
    inventory=json.loads(source.read_text())['higher_basis']
    rows=[];columns=[];selected=[];images={};blocks={}
    for entry in inventory:
        if entry['block'].startswith('fmdd'):continue
        mode,kind=entry['block'].split(':');label=entry['lorentz_label']
        coords=project(mode,label['derivative_on'],tuple(label['metric_pair']),entry['bose_sector'])
        column={}
        for words,value in family(entry).items():
            for i,c in enumerate(coords):
                if c:column[(kind,words,i)]=value*c
        columns.append(column);selected.append(entry)
        images[entry['id']]={'kinematic_polynomial':str(sum(c*z for c,z in zip(coords,(s,t,h)))),
                            'block':entry['block'],'zero':not column}
    keys=sorted(set().union(*(set(c) for c in columns)),key=repr)
    matrix=sp.SparseMatrix(len(keys),len(columns),{(i,j):column[key] for i,key in enumerate(keys)
                         for j,column in enumerate(columns) if key in column})
    reduced,pivots=matrix.rref();rank=len(pivots)
    assert matrix[:,list(pivots)]*reduced[:rank,:]==matrix
    combinations=[]
    for i,pivot in enumerate(pivots):
        combinations.append({'representative':selected[pivot]['id'],
            'coefficients':{entry['id']:str(reduced[i,j]) for j,entry in enumerate(selected) if reduced[i,j]}})
    for block in sorted({e['block'] for e in selected}):
        indices=[i for i,e in enumerate(selected) if e['block']==block]
        blocks[block]={'operators':len(indices),'rank':matrix[:,indices].rank(),
                       'zero_images':sum(not columns[i] for i in indices)}
    output={'scope':'stripped cubic canonical-wedge projection with constant X and independent U,Y flavour weights; NOT complete physical ALP rank',
            'invariants':{'s':'q_A^2','t':'q_B^2','h':'q_A dot q_B; m_a^2=s+t+2h'},
            'projected_operator_count':len(selected),'formal_amplitude_rank':rank,
            'zero_images':sum(not c for c in columns),'blocks':blocks,'images':images,
            'independent_combinations':combinations,
            'deferred':['seven fmdd operators','leading/higher joint physical flavour projection','source/contact and LSZ normalization','measurement identifiability']}
    output['neutral_projection']=neutral_projection(selected,images)
    output['neutral_fminus_zero_ids']=[entry['id'] for entry in inventory if entry['block'].startswith('fmdd')]
    output['neutral_projection']['fminus_scope']='electromagnetic-compatible diagonal ALP/meson sources only; f_- at scalar-photon order is a vanishing charge commutator, including its ordered derivatives'
    path=ROOT/'data/vector_sector/odd_higher_two_boson_projection.json'
    path.write_text(json.dumps(output,indent=2)+'\n')
    transcript=ROOT/'data/vector_sector/odd_higher_two_boson_projection.md'
    with transcript.open('w') as handle,contextlib.redirect_stdout(handle):
        print('# Stripped higher-operator two-boson combinations\n')
        show('Canonical epsilon invariant',E)
        show('Reference invariants',sp.Matrix([dot(b.qa,b.qa),dot(b.qb,b.qb),dot(b.qa,b.qb)]))
        print('The following images omit Fourier, source, electric-charge and canonical LSZ factors.\n')
        for entry in selected:
            image=images[entry['id']]
            flavor=sum(v*flavour_latex(w) for w,v in family(entry).items())
            show(entry['id'],flavor*sp.sympify(image['kinematic_polynomial']))
        show('Formal amplitude rank',sp.Integer(rank))
        for block,values in blocks.items():show(block+' operators/rank/zero',sp.Matrix(list(values.values())))
        print('## Independent coefficient combinations\n')
        for combination in combinations:
            show(combination['representative'],sum(sp.Rational(c)*sp.Symbol('g_{'+identifier+'}')
                 for identifier,c in combination['coefficients'].items()))
        print('## Neutral electromagnetic/neutral-vector restriction\n')
        show('Neutral formal rank',sp.Integer(output['neutral_projection']['rank']))
        for C in output['neutral_projection']['neutral_source_charge_commutators']:
            print('Neutral source charge commutator:\n$$\n'+C+'\n$$\n')
        for combination in output['neutral_projection']['independent_combinations']:
            show(combination['representative'],sum(sp.sympify(c)*sp.Symbol('g_{'+identifier+'}')
                 for identifier,c in combination['coefficients'].items()))
    print(path.relative_to(ROOT));print(transcript.relative_to(ROOT));print(json.dumps(blocks));print('rank',rank)
    print('neutral',json.dumps(output['neutral_projection']['blocks']))

if __name__=='__main__':main()
