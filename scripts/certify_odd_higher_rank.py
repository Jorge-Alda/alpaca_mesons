"""Exact three-flavor lower-bound certificate for the covariant remainder rank.

Run after reduce_odd_covariant_higher.py:
.venv/bin/python scripts/certify_odd_higher_rank.py

The primary basis is defined by symbolic reduction. Deterministic polarized
U(3) Fourier-source evaluations supply an exact lower bound, to compare with
that symbolic upper bound. No floating-point or random-rank criterion is used.
"""
from pathlib import Path
from itertools import permutations, product
from functools import lru_cache
from collections import Counter, defaultdict
import hashlib
import json
import pickle
import sympy as sp
import reduce_odd_covariant_higher as r

ROOT=Path(__file__).resolve().parents[1]
E=lambda i,j:sp.eye(3)[:,i]*sp.eye(3)[j,:]
DIRECTIONS=[sp.eye(3),E(0,1)+E(1,0),-sp.I*E(0,1)+sp.I*E(1,0),
            sp.diag(1,-1,0),E(0,2)+E(2,0),-sp.I*E(0,2)+sp.I*E(2,0),
            E(1,2)+E(2,1),-sp.I*E(1,2)+sp.I*E(2,1),sp.diag(1,1,-2)]

@lru_cache(None)
def color_trace(colors):
    matrix=sp.eye(3)
    for color in colors:matrix=(matrix*DIRECTIONS[color]).applyfunc(sp.expand)
    return sp.trace(matrix)

def coefficients(seed,count,modulus):
    return [hashlib.sha256((seed+':'+str(j)).encode()).digest()[0]%modulus for j in range(count)]

@lru_cache(None)
def leg_assignments(root):
    bases=tuple(base for kind,word in root for base in word)
    physical=tuple(sorted(bases))
    choices=[]
    for base in sorted(set(physical)):
        destinations=[i for i,name in enumerate(bases) if name==base]
        origins=[i for i,name in enumerate(physical) if name==base]
        choices.append([tuple(zip(destinations,ordering)) for ordering in permutations(origins)])
    assignments=[]
    for selections in product(*choices):
        mapping=[None]*len(bases)
        for pairs in selections:
            for destination,origin in pairs:mapping[destination]=origin
        assignments.append(tuple(mapping))
    return physical,tuple(assignments)

def evaluate_row(root,monomial,fixture):
    physical,assignments=leg_assignments(root)
    seed=repr(physical)+':'+str(fixture)
    colors=coefficients(seed+':color',len(physical),len(DIRECTIONS))
    values=coefficients(seed+':momentum',4*(len(physical)-1),5)
    momenta=[None]+[tuple(value-2 for value in values[4*j:4*j+4]) for j in range(len(physical)-1)]
    momenta[0]=tuple(-sum(p[mu] for p in momenta[1:]) for mu in range(4))
    result=0
    for mapping in assignments:
        p=sp.Integer(1)
        for index,power in enumerate(monomial):
            if power:p*=momenta[mapping[1+index//4]][index%4]**power
        if p==0:continue
        offset=0
        for kind,word in root:
            if kind=='trace':p*=color_trace(tuple(colors[mapping[offset+j]] for j in range(len(word))))
            offset+=len(word)
        result+=p
    return sp.expand(result)

def main():
    path=ROOT/'data/vector_sector/odd_covariant_higher_reduction.json'
    data=json.loads(path.read_text())
    if data.get('source_grade_ceiling') is not None:
        raise ValueError('A truncated source-grade run cannot certify the full count.')
    target=data['remainder_rank']
    digest=data['kernel_digest']
    assert hashlib.sha256((ROOT/'scripts/reduce_odd_covariant_higher.py').read_bytes()).hexdigest()==data['sha256']['scripts/reduce_odd_covariant_higher.py']
    cache=Path('/private/tmp/alpaca-odd-covariant-cache')
    columns=[pickle.loads((cache/f'{digest}-joint-kernel-{j}.pickle').read_bytes())
             for j in range(data['kernel_dimension'])]
    groups=defaultdict(dict)
    for j,column in enumerate(columns):
        for (root,monomial),coefficient in column.items():
            fields=tuple(sorted(base for kind,word in root for base in word))
            row=groups[fields].setdefault((root,monomial),{})
            row[j]=coefficient
    ordered=sorted(groups,key=lambda key:hashlib.sha256(repr(key).encode()).hexdigest())
    selected=[];certificates=[];rank=0;attempted=0
    # Revalidate saved independent probes against the current exact source
    # polynomials. This replays the certificate, rather than repeating every
    # unsuccessful probe whenever the reporting layer is revised.
    output=ROOT/'data/vector_sector/odd_higher_rank_certificate.json'
    if output.exists():
        previous=json.loads(output.read_text())
        if previous.get('sha256',{}).get('data/vector_sector/odd_covariant_higher_reduction.json')==hashlib.sha256(path.read_bytes()).hexdigest():
            for certificate in previous['polarized_source_probes']:
                vector=[sp.Integer(0)]*len(columns)
                for (root,monomial),row in groups[tuple(certificate['fields'])].items():
                    value=evaluate_row(root,monomial,certificate['fixture'])
                    if value:
                        for j,coefficient in row.items():vector[j]+=value*coefficient
                vector=[sp.expand(value) for value in vector]
                assert vector==[sp.sympify(value) for value in certificate['values']]
                selected.append(vector);certificates.append(certificate)
            _,pivots=sp.Matrix(selected).T.to_DM().to_field().rref()
            assert len(pivots)==len(selected)
            rank=len(pivots);attempted=previous['attempted_nonzero_probes']
    for fixture in range(12) if not selected else ():
        for fields in ordered:
            vector=[sp.Integer(0)]*len(columns)
            for (root,monomial),row in groups[fields].items():
                value=evaluate_row(root,monomial,fixture)
                if value:
                    for j,coefficient in row.items():vector[j]+=value*coefficient
            vector=[sp.expand(value) for value in vector]
            if not any(vector):continue
            attempted+=1
            trial=sp.Matrix(selected+[vector]).T
            _,pivots=trial.to_DM().to_field().rref()
            if len(pivots)>rank:
                selected.append(vector);rank=len(pivots)
                certificates.append({'fields':list(fields),'fixture':fixture,
                                     'values':[str(value) for value in vector]})
                progress=ROOT/'data/vector_sector/odd_higher_rank_certificate_progress.json'
                progress.write_text(json.dumps({'exact_probe_rank':rank,'symbolic_remainder_rank':target,
                                               'attempted_nonzero_probes':attempted})+'\n')
            if rank==target:break
        if rank==target:break
    result={'three_flavor_remainder_rank_lower_bound':rank,'symbolic_remainder_rank_upper_bound':target,
            'certified':rank==target,'polarized_source_probes':certificates,
            'attempted_nonzero_probes':attempted,'probe_matrix_rank':sp.Matrix(selected).rank(),
            'joint_higher_action_count':data['common_source_projection_rank']+rank if rank==target else None,
            'full_action_count':40+data['common_source_projection_rank']+rank if rank==target else None,
            'new_independent_operators':29+data['common_source_projection_rank']+rank if rank==target else None,
            'sha256':{'scripts/certify_odd_higher_rank.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      'data/vector_sector/odd_covariant_higher_reduction.json':hashlib.sha256(path.read_bytes()).hexdigest(),
                      'scripts/reduce_odd_covariant_higher.py':data['sha256']['scripts/reduce_odd_covariant_higher.py']}}
    result['role']='finite-source lower bound; use close_odd_three_flavor_quotient.py for complete closure'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(output.relative_to(ROOT))

if __name__=='__main__':main()
