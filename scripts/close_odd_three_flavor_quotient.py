"""Close any gap between exact source probes and the abstract trace-word bound.

Run .venv/bin/python scripts/close_odd_three_flavor_quotient.py.
Uses the exact U(3) generator completeness relation and the permutation
Gram matrix of polarized traces. This includes all finite-flavor trace
identities, without generic flavor-component matrices or floating point.
"""
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
import hashlib
import json
import pickle
import sympy as sp
from sympy.combinatorics import Permutation
import reduce_odd_covariant_higher as r
from certify_odd_higher_rank import DIRECTIONS, leg_assignments

ROOT=Path(__file__).resolve().parents[1]
CACHE=Path('/private/tmp/alpaca-odd-covariant-cache')

@lru_cache(None)
def polarized_monomials(monomial,mapping):
    n=len(mapping);momenta=sp.symbols(f'k0:{4*(n-1)}')
    polynomial=sp.Integer(1)
    for index,power in enumerate(monomial):
        if not power:continue
        physical=mapping[1+index//4];mu=index%4
        value=-sum(momenta[4*(j-1)+mu] for j in range(1,n)) if physical==0 \
              else momenta[4*(physical-1)+mu]
        polynomial*=value**power
    return sp.Poly(sp.expand(polynomial),*momenta).as_dict()

@lru_cache(None)
def trace_permutation(root,mapping):
    physical,assignments=leg_assignments(root)
    color_legs=[i for i,name in enumerate(physical) if name!='Z']
    rank={physical:i for i,physical in enumerate(color_legs)}
    sigma=list(range(len(color_legs)));offset=0
    for kind,word in root:
        if kind=='trace':
            cycle=[rank[mapping[offset+j]] for j in range(len(word))]
            for a,b in zip(cycle,cycle[1:]+cycle[:1]):sigma[a]=b
        offset+=len(word)
    return tuple(sigma)

def polarize(column):
    result=defaultdict(dict)
    for (root,monomial),coefficient in column.items():
        physical,assignments=leg_assignments(root)
        for mapping in assignments:
            sigma=trace_permutation(root,mapping)
            for power,value in polarized_monomials(monomial,mapping).items():
                row=result[physical,power]
                row[sigma]=row.get(sigma,0)+coefficient*value
    return {key:{sigma:sp.expand(value) for sigma,value in row.items() if sp.expand(value)!=0}
            for key,row in result.items() if any(sp.expand(value)!=0 for value in row.values())}

@lru_cache(None)
def loops(sigma,tau):
    return (Permutation(list(sigma),size=len(sigma))**-1*Permutation(list(tau),size=len(tau))).cycles

def gram_polynomial(columns):
    n=sp.Symbol('N_f',integer=True,positive=True)
    groups=defaultdict(dict)
    for i,column in enumerate(columns):
        for key,row in polarize(column).items():groups[key][i]=row
    entries=defaultdict(lambda:defaultdict(lambda:sp.Integer(0)))
    for key,rows in groups.items():
        for i,row in rows.items():
            for j,other in rows.items():
                if j<i:continue
                for sigma,c in row.items():
                    for tau,d in other.items():
                        entries[i,j][loops(sigma,tau)]+=sp.conjugate(c)*d
    G=sp.zeros(len(columns))
    for (i,j),coefficients in entries.items():
        value=sp.expand(sum(sp.expand(c)*n**power for power,c in coefficients.items()))
        G[i,j]=value;G[j,i]=sp.conjugate(value)
    return n,G,len(groups)

def verify_completeness():
    T=[direction/sp.sqrt(sp.trace(direction*direction)) for direction in DIRECTIONS]
    # This finite SymPy contraction derives the normalized U(3) Fierz tensor.
    contractions=[]
    for a in range(3):
        for b in range(3):
            for c in range(3):
                for d in range(3):
                    value=sp.simplify(sum(matrix[a,b]*matrix[c,d] for matrix in T))
                    assert value==sp.KroneckerDelta(a,d)*sp.KroneckerDelta(b,c)
                    contractions.append(value)
    return True

def main():
    assert verify_completeness()
    reduction_path=ROOT/'data/vector_sector/odd_covariant_higher_reduction.json'
    probe_path=ROOT/'data/vector_sector/odd_higher_rank_certificate.json'
    reduction=json.loads(reduction_path.read_text());probes=json.loads(probe_path.read_text())
    digest=reduction['kernel_digest']
    columns=[pickle.loads((CACHE/f'{digest}-joint-kernel-{j}.pickle').read_bytes())
             for j in range(reduction['kernel_dimension'])]
    keys,R=r.matrix_of(columns)
    _,pivots=R.to_DM().to_field().rref()
    assert len(pivots)==reduction['remainder_rank']
    P=sp.Matrix([[sp.sympify(value) for value in probe['values']] for probe in probes['polarized_source_probes']])
    P=P[:,list(pivots)]
    dm,pp=P.to_DM().to_field().rref();echelon=dm.to_Matrix()
    candidates=[];coordinates=[]
    for free in (j for j in range(P.cols) if j not in pp):
        vector=sp.zeros(P.cols,1);vector[free]=1
        for row,pivot in enumerate(pp):vector[pivot]=-echelon[row,free]
        assert P*vector==sp.zeros(P.rows,1)
        output={}
        for j,value in enumerate(vector):
            if not value:continue
            for key,coefficient in columns[pivots[j]].items():output[key]=output.get(key,0)+value*coefficient
        output={key:sp.expand(value) for key,value in output.items() if sp.expand(value)!=0}
        coordinates.append(vector);candidates.append(output)
    n,G,groups=gram_polynomial(candidates)
    G3=G.subs(n,3)
    _,gp=G3.to_DM().to_field().rref() if G3.rows else (None,())
    rank=len(pp)+len(gp)
    assert rank<=len(pivots)
    null=G3.nullspace() if G3.rows else []
    raw_vectors=[sp.Matrix([sp.sympify(vector.get(identifier,'0')) for identifier in reduction['raw_ids']])
                 for vector in reduction['corrected_kernel_vectors']]
    extra_relations=[]
    for vector in null:
        raw=sp.zeros(len(reduction['raw_ids']),1)
        for i,value in enumerate(vector):
            if not value:continue
            for j,coefficient in enumerate(coordinates[i]):
                if coefficient:raw+=value*coefficient*raw_vectors[pivots[j]]
        extra_relations.append({identifier:str(value) for identifier,value in zip(reduction['raw_ids'],raw) if value})
    # Preserve early derivative blocks when selecting representatives: eliminate
    # later columns first, using the complete exact relation space.
    relations=reduction['exact_formal_relations']+extra_relations
    ids=reduction['raw_ids']
    N=sp.Matrix([[sp.sympify(relation.get(identifier,'0')) for identifier in reversed(ids)] for relation in relations])
    _,deleted=N.to_DM().to_field().rref()
    removed={len(ids)-1-j for j in deleted}
    higher_basis=[identifier for j,identifier in enumerate(ids) if j not in removed]
    higher=reduction['common_source_projection_rank']+rank
    assert len(higher_basis)==higher
    counts={}
    for identifier in higher_basis:
        mode,kind,*rest=identifier.split(':')
        name=mode+':'+kind;counts[name]=counts.get(name,0)+1
    result={'certified':True,'exact_source_probe_rank':len(pp),
            'abstract_trace_remainder_upper_bound':len(pivots),'unresolved_directions':len(candidates),
            'candidate_Gram_polynomial':[[str(value) for value in G.row(i)] for i in range(G.rows)],
            'candidate_Gram_at_Nf3':[[str(value) for value in G3.row(i)] for i in range(G3.rows)],
            'candidate_Gram_at_Nf7':[[str(value) for value in G.subs(n,7).row(i)] for i in range(G.rows)],
            'candidate_polarized_polynomial_groups':groups,'additional_true_rank':len(gp),
            'three_flavor_remainder_rank':rank,'additional_exact_relations':extra_relations,
            'joint_higher_action_count':higher,'full_action_count':40+higher,
            'original_action_rank':11,'new_independent_operators':29+higher,
            'higher_basis_ids':higher_basis,'basis_counts':counts,
            'all_higher_relations':relations,
            'sha256':{'scripts/close_odd_three_flavor_quotient.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      'data/vector_sector/odd_covariant_higher_reduction.json':hashlib.sha256(reduction_path.read_bytes()).hexdigest(),
                      'data/vector_sector/odd_higher_rank_certificate.json':hashlib.sha256(probe_path.read_bytes()).hexdigest()}}
    output=ROOT/'data/vector_sector/odd_three_flavor_action_quotient.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(output.relative_to(ROOT))

if __name__=='__main__':main()
