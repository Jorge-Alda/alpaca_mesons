"""Charged and neutral cubic SM kernels of all 241 bounded raw parents.

Run .venv/bin/python scripts/derive_odd_charged_kernels.py.
GeV units, (+---), epsilon0123=+1; all momenta incoming; tree level.
Derivative-on-fminus acts on the entire P-gamma composite, with momentum
p_P+q_gamma=-q_V. No constant/neutral projection is applied to its commutator.
"""
from pathlib import Path
from functools import lru_cache
import contextlib
import hashlib
import json
import sympy as sp
import derive_odd_power_counted_channels as d
import enumerate_odd_breaking_vertices as b
import project_odd_higher_two_boson as hp

ROOT = Path(__file__).resolve().parents[1]


@lru_cache(None)
def fminus_projection(where,pair):
    negative = tuple(-q for q in b.qa)
    left,right = {'AA':(b.qa,b.qa),'AB':(b.qa,negative),'BB':(negative,negative)}[where]
    raw = b.contraction(pair,hp.AP,hp.BP,left,right)
    pivots,coordinates = b.reduction(hp.refs+[raw])
    assert pivots==(0,1,2)
    result = coordinates[:,3]
    assert sp.expand(raw-sum(result[i]*hp.refs[i] for i in range(3)))==0
    return tuple(result)


def parents():
    records = []
    for spec in d.specifications():
        if spec['mode']=='fminus':
            polynomial = sp.Integer(4)
        elif spec['mode']=='fmdd':
            polynomial = sum(c*x for c,x in zip(fminus_projection(spec['label']['derivative_on'],
                                  tuple(spec['label']['metric_pair'])),(d.s,d.t,d.h)))
        else:
            label = spec['label']
            polynomial = d.lorentz(spec['mode'],label['derivative_on'] if label else '',
                                   tuple(label['metric_pair']) if label else (),spec['group'])
        coefficients = [sp.expand(polynomial).coeff(x) for x in (d.s,d.t,d.h)]
        constant = sp.expand(polynomial).subs({d.s:0,d.t:0,d.h:0})
        assert sp.expand(polynomial-constant-sum(c*x for c,x in zip(coefficients,(d.s,d.t,d.h))))==0
        records.append(dict(id=spec['id'],kind=spec['kind'],mode=spec['mode'],
            level=spec['level'],traces=spec['traces'],group=spec['group'],label=spec['label'],
            flavor_words=[dict(words=[list(w) for w in words],coefficient=str(value))
                          for words,value in spec['family'].items()],
            kinematic_coefficients=[str(constant),*[str(c) for c in coefficients]],
            operator_MV_power=1+2*spec['level'] if spec['kind']=='VJP' else 2*spec['level'],
            phase=-1 if spec['mode'] in ('Y','fminus') else 1,
            bose_factor=2 if spec['kind']=='VVP' and spec['level']==1 else 1))
    return records


def sm_kernel(spec,A,B,TP):
    U = -sp.sqrt(2)*TP/d.F
    Y = sp.sqrt(2)*(TP*d.CHI+d.CHI*TP)/d.F
    if spec['mode'] not in ('fminus','fmdd'):
        return -2*d.e*d.degree_normalization(spec)*d.kernel(spec,A,B,U,Y)/d.M1
    flavor = d.flavor(spec['family'],dict(A=A,B=TP*d.Q-d.Q*TP,X=2*d.CHI))
    if spec['mode']=='fminus':
        polynomial,phase = sp.Integer(4),-1
    else:
        polynomial = sum(c*x for c,x in zip(fminus_projection(spec['label']['derivative_on'],
                          tuple(spec['label']['metric_pair'])),(d.s,d.t,d.h)))
        phase = 1
    return sp.expand(phase*sp.sqrt(2)*d.e*flavor*polynomial*d.degree_normalization(spec)/(d.F*d.M1))


def main():
    specs = d.specifications()
    records = parents()
    by_id = {spec['id']:spec for spec in specs}
    leading = json.loads((ROOT/'data/vector_sector/odd_leading_action_count.json').read_text())
    higher = json.loads((ROOT/'data/vector_sector/odd_three_flavor_action_quotient.json').read_text())
    relations = leading['relations']+higher['all_higher_relations']
    def unit(i,j):
        matrix = sp.zeros(3);matrix[i,j]=1;return matrix
    W0 = sp.diag(*sp.symbols('W_0u W_0d W_0s'))
    channels = {}
    for name,i,j in [('rho_charged_pi',0,1),('Kstar_charged_K',0,2),('Kstar_neutral_K',1,2)]:
        TP,A = unit(i,j),unit(j,i)
        U = -sp.sqrt(2)*TP/d.F
        Y = sp.sqrt(2)*(TP*d.CHI+d.CHI*TP)/d.F
        values = []
        for spec in specs:
            value = (sm_kernel(spec,A,d.Q,TP) if spec['kind']=='VJP' else
                     sp.sqrt(2)*d.e*d.degree_normalization(spec)*d.kernel(spec,A,W0,U,Y)/d.M1)
            values.append(sp.expand(value.subs({d.t:0,d.h:(d.w-d.s)/2})))
        raw_residuals = [sp.expand(sum(sp.sympify(c)*values[[s['id'] for s in specs].index(identifier)]
                                  for identifier,c in relation.items())) for relation in relations]
        lo_mass2 = (d.CHI[i,i]+d.CHI[j,j])/2
        residuals = [sp.expand(value.subs(d.w,lo_mass2)) for value in raw_residuals]
        assert all(value==0 for value in residuals),(name,{k:str(v) for k,v in enumerate(residuals) if v})
        fk = sp.trace(A*((TP*d.Q-d.Q*TP)*d.CHI-d.CHI*(TP*d.Q-d.Q*TP)))
        original_c4 = 2*sm_kernel(by_id['VJP:fminus'],A,d.Q,TP)
        assert sp.expand(original_c4+8*sp.sqrt(2)*d.e*fk/(d.F*d.M1*d.MV))==0
        channels[name] = dict(tree_mass2=str(lo_mass2),fminus_commutator_trace=str(fk),
            original_c4_gradient=str(original_c4),all_13_tree_relation_images_zero=True,
            raw_parent_gradients={spec['id']:str(v) for spec,v in zip(specs,values)},
            unmatched_mass_relation_residuals={str(k):str(v) for k,v in enumerate(raw_residuals) if v})
    for matrix in (hp.AP.subs(dict(zip(b.ea,b.qa)),simultaneous=True),
                   hp.BP.subs(dict(zip(b.ev,b.qb)),simultaneous=True)):
        assert matrix==sp.zeros(4)
    fm = [record for record in records if record['mode']=='fmdd']
    polynomials = [sp.sympify(record['kinematic_coefficients'][0])+sum(sp.sympify(c)*x
                   for c,x in zip(record['kinematic_coefficients'][1:],(d.s,d.t,d.h))) for record in fm]
    matrix = b.coefficient_matrix(polynomials)
    _,pivots = matrix.rref()
    paths = ['scripts/derive_odd_charged_kernels.py','scripts/derive_odd_power_counted_channels.py',
             'data/vector_sector/odd_vertex_enumeration.json','data/vector_sector/odd_covariant_block_preparation.json',
             'data/vector_sector/odd_leading_action_count.json','data/vector_sector/odd_three_flavor_action_quotient.json']
    output = ROOT/'data/vector_sector/odd_cubic_parent_kernels.json'
    output.write_text(json.dumps(dict(scope='all 241 raw SM cubic parents; neutral and charged sources; no fit',
        parents=records,charged_channels=channels,higher_fminus_cubic_rank=len(pivots),
        higher_fminus_zero_images=sum(p==0 for p in polynomials),
        fminus_derivative_momentum='p_P+q_gamma=-q_V',
        checks=['exact wedge reconstruction','photon Ward tensors','charged c4 normalization',
                'all thirteen tree-mass action relations in three charged channels'],
        sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}),indent=2)+'\n')
    transcript = ROOT/'data/vector_sector/odd_charged_cubic_kernels.md'
    with transcript.open('w') as handle,contextlib.redirect_stdout(handle):
        print('# Charged SM cubic kernels\n')
        d.show('Composite fminus derivative momentum',sp.Matrix([hp.p[k]+b.qb[k] for k in range(4)]))
        d.show('Higher fminus cubic kinematic rank',sp.Integer(len(pivots)))
        for record,polynomial in zip(fm,polynomials):d.show(record['id'],polynomial)
        for name,record in channels.items():
            d.show(name+' LO source mass squared',sp.sympify(record['tree_mass2']))
            d.show(name+' charged c4 gradient',sp.sympify(record['original_c4_gradient']))
        print('\nAll 13 action-relation images vanish at the LO source mass in each charged channel.\n')
    print(output.relative_to(ROOT));print(transcript.relative_to(ROOT))


if __name__=='__main__':main()
