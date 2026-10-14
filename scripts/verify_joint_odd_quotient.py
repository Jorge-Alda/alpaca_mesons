"""Verify complete counts, source-functional leading rank and input certificates.

Run .venv/bin/python scripts/verify_joint_odd_quotient.py.
"""
from pathlib import Path
from itertools import permutations
import hashlib
import json
import sympy as sp
from sympy.combinatorics import Permutation
import reduce_odd_covariant_higher as r
from close_odd_three_flavor_quotient import gram_polynomial, loops, verify_completeness

ROOT=Path(__file__).resolve().parents[1]

def main():
    assert verify_completeness()
    color_ranks={}
    for n in (2,3,4):
        words=[p for p in permutations(range(n)) if Permutation(list(p),size=n).cycles<=2]
        matrix=sp.Matrix([[3**loops(a,b) for b in words] for a in words])
        _,pivots=matrix.to_DM().to_field().rref()
        assert len(pivots)==len(words)
        color_ranks[str(n)]=len(pivots)
    # Independent functional check of the leading nonlinear Bianchi sector.
    leading=[]
    for family in r.b.trace_catalog(('A','B','U')):
        expression=r.local_operator({'family':family,'mode':'X','kind':'VJP','group':'all',
                                      'label':{'derivative_on':'B','metric_pair':[0,4]}})
        leading.append(r.normal_form(expression))
    Nf,G,groups=gram_polynomial(leading)
    leading_gram=G.subs(Nf,3)
    _,pivots=leading_gram.to_DM().to_field().rref()
    assert len(pivots)==3
    paths=['data/vector_sector/odd_leading_action_count.json',
           'data/vector_sector/odd_source_gauge_verification.json',
           'data/vector_sector/odd_covariant_higher_reduction.json',
           'data/vector_sector/odd_three_flavor_action_quotient.json']
    documents={path:json.loads((ROOT/path).read_text()) for path in paths}
    for document in documents.values():
        for path,digest in document.get('sha256',{}).items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    result=documents[paths[-1]]
    assert result['certified']
    assert (result['joint_higher_action_count'],result['full_action_count'],
            result['original_action_rank'],result['new_independent_operators'])==(188,228,11,217)
    assert sum(result['basis_counts'].values())==188
    assert len(result['higher_basis_ids'])==188
    assert len(result['all_higher_relations'])==197-188
    output=ROOT/'data/vector_sector/odd_joint_quotient_verification.json'
    output.write_text(json.dumps({'passed':True,'three_flavor_low_grade_color_ranks':color_ranks,
        'leading_source_functional_Bianchi_rank':len(pivots),
        'leading_Bianchi_Gram_at_Nf3':[[str(x) for x in leading_gram.row(i)] for i in range(4)],
        'checks':['complete U(3) Fierz tensor','injectivity of at-most-double-trace low-grade tensors',
                  'leading source-functional Bianchi rank','closed higher rank certificate',
                  'basis cardinalities and original rank','all provenance hashes'],
        'sha256':{**{path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
                  'scripts/verify_joint_odd_quotient.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}},indent=2)+'\n')
    print(output.relative_to(ROOT))

if __name__=='__main__':main()
