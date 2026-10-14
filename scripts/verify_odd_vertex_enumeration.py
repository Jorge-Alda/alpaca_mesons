"""Independent physical invariance checks for the stage-A tensor catalogue.

Run .venv/bin/python scripts/verify_odd_vertex_enumeration.py
Uses exact rational boosts, concrete U(3) flavor generators and photon Ward
replacement; does not run a calibration or replace symbolic rank certificates.
"""
from pathlib import Path
import contextlib
import hashlib
import io
import json
import sympy as sp
import enumerate_odd_breaking_vertices as basis

ROOT=Path(__file__).resolve().parents[1]

def main():
    eta=sp.diag(1,-1,-1,-1)
    boost=sp.Matrix([[sp.Rational(5,3),sp.Rational(4,3),0,0],
                     [sp.Rational(4,3),sp.Rational(5,3),0,0],[0,0,1,0],[0,0,0,1]])
    assert boost.T*eta*boost==eta and boost.det()==1
    A=basis.antisymmetric([1,2,-1,3,4,2]); B=basis.antisymmetric([-2,1,4,-1,3,5])
    U=sp.Matrix([1,-2,3,2]); D=sp.Matrix([2,1,-3,4])
    for pair in basis.PAIR:
        assert sp.expand(basis.contraction(pair,A,B,U,D)
                 -basis.contraction(pair,boost*A*boost.T,boost*B*boost.T,boost*U,boost*D))==0
        assert basis.contraction(pair,A,sp.zeros(4),U,D)==0
    # Photon field strength, before any on-shell specialization.
    k=sp.Matrix([3,-2,1,4]); e=sp.Matrix([1,2,-1,3])
    F=k*e.T-e*k.T
    for pair in basis.PAIR:
        actual=basis.contraction(pair,A,F,U,D)
        transformed=basis.contraction(pair,boost*A*boost.T,boost*F*boost.T,boost*U,boost*D)
        assert sp.expand(actual-transformed)==0
        assert basis.contraction(pair,A,k*k.T-k*k.T,U,D)==0
    # Concrete normalized-flavor directions rather than a generic matrix
    # replacement for the abstract flavor calculation.
    def E(i,j):
        result=sp.zeros(3);result[i,j]=1;return result
    matrices={'A':E(0,1)+E(1,0),'B':E(1,2)+E(2,1),
              'U':E(0,2)+E(2,0),'X':sp.diag(2,-1,3)}
    with contextlib.redirect_stdout(io.StringIO()):identity=basis.chiral_trace_identity()
    value=sum(coef*sp.Mul(*(sp.trace(sp.MatMul(*(matrices[n] for n in word)).doit())
                       for word in traces)) for traces,coef in identity.items())
    assert sp.simplify(value)==0
    # Reflection changes the epsilon sign. This checks pseudoscalar parity,
    # separately from proper-Lorentz covariance.
    reflection=sp.diag(1,-1,1,1)
    for pair in basis.PAIR:
        assert sp.expand(basis.contraction(pair,reflection*A*reflection.T,
             reflection*B*reflection.T,reflection*U,reflection*D)
             +basis.contraction(pair,A,B,U,D))==0
    path=ROOT/'data/vector_sector/odd_vertex_enumeration.json'
    data=json.loads(path.read_text())
    for name,digest in data['sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    assert len(data['operator_inventory'])==sum(block['independent_source_resolved_vertex_count'] for block in data['summary'].values())
    output=ROOT/'data/vector_sector/odd_vertex_verification.json'
    output.write_text(json.dumps({'checks':['exact rational proper-Lorentz boost','orientation/parity reflection',
        'off-shell photon Ward replacement','concrete U(3) polarized trace identity','input hashes and inventory consistency'],
        'passed':True,'enumeration_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
    print(output.relative_to(ROOT))

if __name__=='__main__':main()
