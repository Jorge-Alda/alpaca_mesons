"""Concrete U(3) probes of the nonlinear Bianchi rank and count certificates.

Run .venv/bin/python scripts/verify_odd_leading_action_count.py.
These deterministic physical flavor directions verify the abstract calculation;
they do not replace its trace-word derivation or define a fitted model.
"""
from pathlib import Path
import hashlib
import json
import sympy as sp
import enumerate_odd_breaking_vertices as b
from lift_odd_leading_covariant import bianchi_remainder

ROOT=Path(__file__).resolve().parents[1]

def main():
    # Rationally scaled Gell-Mann/identity directions, to keep arithmetic exact.
    E=lambda i,j:sp.eye(3)[:,i]*sp.eye(3)[j,:]
    directions=[sp.eye(3),E(0,1)+E(1,0),-sp.I*E(0,1)+sp.I*E(1,0),
                sp.diag(1,-1,0),E(0,2)+E(2,0),-sp.I*E(0,2)+sp.I*E(2,0),
                E(1,2)+E(2,1),-sp.I*E(1,2)+sp.I*E(2,1),sp.diag(1,1,-2)]
    families=b.trace_catalog(('A','B','U'))
    remainder=[bianchi_remainder(f) for f in families]
    names=sorted({name for column in remainder for traces in column for word in traces for name in word})
    rows=[]
    for fixture in range(8):
        # Explicit finite source/tensor fields in named U(3) directions.
        matrices={name:sum(((((index+2)*(fixture+3)+(j+1)**2)%7-3)*T
                          for j,T in enumerate(directions)),start=sp.zeros(3))
                  for index,name in enumerate(names)}
        cache={}
        for column in remainder:
            for traces in column:
                for word in traces:
                    if word not in cache:
                        product=sp.eye(3)
                        for name in word:product=product*matrices[name]
                        cache[word]=sp.trace(product)
        rows.append([sp.expand(sum(coef*sp.prod(cache[word] for word in traces)
                                   for traces,coef in column.items())) for column in remainder])
    probes=sp.Matrix(rows)
    assert probes.rank()==3
    assert all(probes[i,3]==0 for i in range(probes.rows))
    path=ROOT/'data/vector_sector/odd_leading_action_count.json'
    data=json.loads(path.read_text())
    assert (data['raw_Lorentz_flavor_basis'],data['independent_leading_action_operators'],
            data['original_action_rank'],data['new_leading_independent_operators'])==(44,40,11,29)
    assert len(data['relations'])==4
    assert len(data['basis_preserving_original_action'])==40
    for source,digest in data['sha256'].items():
        assert hashlib.sha256((ROOT/source).read_bytes()).hexdigest()==digest,source
    output=ROOT/'data/vector_sector/odd_leading_action_verification.json'
    output.write_text(json.dumps({'passed':True,'concrete_U3_Bianchi_probe_rank':probes.rank(),
              'probe_matrix':[[str(x) for x in probes.row(i)] for i in range(probes.rows)],
              'checked':['three independent nonlinear Bianchi directions',
                         'exact traced-Bianchi zero','count and basis dimensions','input hashes'],
              'count_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
              'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
    print(output.relative_to(ROOT))

if __name__=='__main__':main()
