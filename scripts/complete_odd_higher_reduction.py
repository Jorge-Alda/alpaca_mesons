"""Complete the joint projection by matching f-minus DD remainders to X seeds.

Run after reduce_odd_covariant_higher.py:
.venv/bin/python scripts/complete_odd_higher_reduction.py

A seed-dependent linearization is not a common projection: a nonlinear
f-minus DD remainder can already lie in the leading source grade of a
V f-plus u X derivative seed. This script removes that overlap before
computing the covariant quotient rank.
"""
from pathlib import Path
import json
import pickle
import hashlib
import sympy as sp
import reduce_odd_covariant_higher as r

ROOT=Path(__file__).resolve().parents[1]
CACHE=Path('/private/tmp/alpaca-odd-covariant-cache')

def legs(key):return sum(len(word) for kind,word in key[0])

def linear_X_operator(spec):
    terms=[];pair=tuple(spec['label']['metric_pair']);remaining=[i for i in range(6) if i not in pair]
    for epsilon,sign in r.b.EPS:
        for metric in range(4):
            slots=[0]*6
            for index,value in zip(remaining,epsilon):slots[index]=value
            slots[pair[0]]=slots[pair[1]]=metric
            A=r.nc_V(slots[0],slots[1])
            B=2*(r.nc_partial(r.nc_jet(f'v{slots[3]}'),slots[2])
                 -r.nc_partial(r.nc_jet(f'v{slots[2]}'),slots[3]))
            if A==0 or B==0:continue
            matrices={'A':A,'B':B,'U':r.nc_u(slots[4]),'X':r.nc_jet('X')}
            name=spec['label']['derivative_on']
            matrices[name]=r.nc_partial(matrices[name],slots[5])
            terms.append(sign*r.b.G[metric]*r.flavor_expression(spec['family'],matrices))
    return sp.expand(sp.Add(*terms))

def add_scaled(output,column,coefficient):
    for key,value in column.items():output[key]=output.get(key,0)+coefficient*value

def main():
    source=ROOT/'data/vector_sector/odd_covariant_higher_kernel.json'
    prepared=json.loads(source.read_text())
    digest=prepared['script_sha256'][:16]
    assert hashlib.sha256((ROOT/'scripts/reduce_odd_covariant_higher.py').read_bytes()).hexdigest()==prepared['script_sha256']
    assert prepared['source_grade_ceiling'] is None
    specs=r.specifications();ids={spec['id']:i for i,spec in enumerate(specs)}
    vectors=[sp.Matrix([sp.sympify(vector.get(spec['id'],'0')) for spec in specs])
             for vector in prepared['kernel_vectors']]
    kernels=[pickle.loads((CACHE/f'{digest}-kernel-{j}-legs-None.pickle').read_bytes())
             for j in range(len(vectors))]
    fm=[j for j,vector in enumerate(prepared['kernel_vectors']) if all(name.startswith('fmdd:') for name in vector)]
    assert len(fm)==4
    xindices=[i for i,spec in enumerate(specs) if spec['mode']=='X' and spec['kind']=='VJP']
    phi=[]
    for i in xindices:
        path=CACHE/f'{digest}-linear-X-{i}.pickle'
        if path.exists():column=pickle.loads(path.read_bytes())
        else:
            column=r.normal_form(linear_X_operator(specs[i]));path.write_bytes(pickle.dumps(column))
        assert all(legs(key)==4 for key in column)
        phi.append(column)
    targets=[{key:value for key,value in kernels[j].items() if legs(key)==4} for j in fm]
    assert all(not any(legs(key)<4 for key in kernels[j]) for j in fm)
    keys,matrix=r.matrix_of(phi+targets)
    reduced,pivots=matrix.to_DM().to_field().rref();reduced=reduced.to_Matrix()
    assert len(pivots)==80 and all(j<len(phi) for j in pivots), \
           'The mixed source grade contains directions outside the X-derivative span.'
    coefficients=reduced[:80,len(phi):]
    correction=[]
    for k,j in enumerate(fm):
        column=dict(kernels[j]);vector=vectors[j].copy();terms={}
        for row,pivot in enumerate(pivots):
            coefficient=coefficients[row,k]
            if coefficient==0:continue
            i=xindices[pivot];terms[specs[i]['id']]=str(coefficient)
            path=CACHE/f'{digest}-operator-{i}.pickle'
            if path.exists():expression=pickle.loads(path.read_bytes());r.registry_restore(expression)
            else:expression=r.local_operator(specs[i]);path.write_bytes(pickle.dumps(expression))
            operator=r.normal_form(expression)
            add_scaled(column,operator,-coefficient);vector[i]-=coefficient
        column={key:sp.expand(value) for key,value in column.items() if sp.expand(value)!=0}
        assert all(legs(key)>4 for key in column)
        kernels[j]=column;vectors[j]=vector
        correction.append({'candidate':j,'X_derivative_subtraction':terms,'remaining_rows':len(column)})
    for j,column in enumerate(kernels):
        (CACHE/f'{digest}-joint-kernel-{j}.pickle').write_bytes(pickle.dumps(column))
        if j not in fm:
            mode=specs[next(i for i,value in enumerate(vectors[j]) if value)]['mode']
            assert all(legs(key)>=(5 if mode=='X' else 4) for key in column)
    remainder_keys,remainder=r.matrix_of(kernels)
    reduced,rp=remainder.to_DM().to_field().rref();reduced=reduced.to_Matrix()
    relations=[]
    for free in (j for j in range(len(kernels)) if j not in rp):
        vector=sp.zeros(len(kernels),1);vector[free]=1
        for row,pivot in enumerate(rp):vector[pivot]=-reduced[row,free]
        raw=sum((vectors[j]*value for j,value in enumerate(vector) if value),sp.zeros(len(specs),1))
        relations.append({specs[i]['id']:str(value) for i,value in enumerate(raw) if value})
    result={'scope':'joint higher covariant seeds; independent sources and LO Goldstone EOM; no vector EOM',
            'ordered_covariant_seed_count':len(specs),'common_source_projection_rank':prepared['flat_source_EOM_rank'],
            'mixed_grade_X_projection_rank':80,'cross_block_corrections':correction,
            'remainder_rank':len(rp),'formal_trace_action_rank_upper_bound':prepared['flat_source_EOM_rank']+len(rp),
            'source_grade_ceiling':None,'kernel_dimension':len(kernels),'kernel_digest':digest,
            'exact_formal_relations':relations,'remainder_rows':len(remainder_keys),
            'joint_action_count':None,
            'raw_ids':[spec['id'] for spec in specs],
            'corrected_kernel_vectors':[{specs[i]['id']:str(value) for i,value in enumerate(vector) if value}
                                        for vector in vectors],
            'sha256':{'scripts/complete_odd_higher_reduction.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      'scripts/reduce_odd_covariant_higher.py':prepared['script_sha256'],
                      'data/vector_sector/odd_covariant_higher_kernel.json':hashlib.sha256(source.read_bytes()).hexdigest()}}
    output=ROOT/'data/vector_sector/odd_covariant_higher_reduction.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(output.relative_to(ROOT))

if __name__=='__main__':main()
