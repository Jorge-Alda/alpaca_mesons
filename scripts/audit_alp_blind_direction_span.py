"""Compare ALP tensor projections with the current SM sensitivity row space.

Run .venv/bin/python scripts/audit_alp_blind_direction_span.py.
No physical fit has run. Independent formal ALP spurion traces test a span,
not amplitude normalization, UV coefficients or marginal statistical errors.
"""
from contextlib import redirect_stdout
from pathlib import Path
import hashlib
import io
import json
import numpy as np
import sympy as sp
from audit_alp_vector_blindness import audit as vector_audit
from audit_alp_gamma_gamma_blindness import audit as photon_audit
from diagnose_vector_stage4 import OPE

ROOT=Path(__file__).resolve().parents[1]
TOL=1e-11


def row_space(rows):
    matrix=np.asarray(rows,dtype=float)
    norms=np.linalg.norm(matrix,axis=1)
    matrix=matrix[norms>0]/norms[norms>0,None]
    _,singular,vt=np.linalg.svd(matrix,full_matrices=False)
    rank=int(np.count_nonzero(singular>TOL))
    return rank,vt[:rank],matrix


def main():
    with redirect_stdout(io.StringIO()):
        vector=vector_audit()
    photon=photon_audit()
    assert sp.simplify(vector['rows']*vector['nulls'])==sp.zeros(2,3)
    assert photon['blind_projection']==sp.zeros(2,3)
    assert vector['nulls'].rank()==3
    # Deliberately independent derivative and mass weights. The exact
    # symbolic zero is primary; samples only compare stored numeric spans.
    examples=[{'D':1.,'M':0.,'m_P':1.8,'s_1':.16,'s_2':.25},
              {'D':0.,'M':1.,'m_P':1.8,'s_1':.16,'s_2':.25},
              {'D':1.,'M':-.3,'m_P':2.1,'s_1':1.,'s_2':.25},
              {'D':.6,'M':.8,'m_P':1.2,'s_1':.49,'s_2':0.}]
    targets=[]
    for example in examples:
        values={symbol:example[str(symbol)] for symbol in vector['rows'].free_symbols}
        targets.extend(np.array(vector['rows'].subs(values),dtype=float))
    targets=np.asarray(targets)
    target_norm=np.linalg.norm(targets,axis=1)
    targets=targets[target_norm>0]/target_norm[target_norm>0,None]
    current=json.loads((ROOT/'data/vector_sector/stage4_rank_results.json').read_text())
    results=[]
    for scenario in current['scenarios']:
        rows=np.vstack([scenario[k] for k in ('radiative_rows','real_tff_rows','virtual_tff_rows')])
        blocks={}
        for label,matrix in (('sm_observable_rows_only',rows),('with_inherited_ope_rows',np.vstack([OPE,rows]))):
            rank,basis,normalized=row_space(matrix)
            residual=targets-targets@basis.T@basis
            error=float(np.max(np.abs(residual)))
            blocks[label]=dict(rank=rank,max_target_span_residual=error,
                              target_rows_in_span=error<1e-9)
            np.testing.assert_allclose(normalized@np.array(vector['nulls'],dtype=float),0,atol=1e-10)
            if label=='with_inherited_ope_rows':
                assert rank==8 and error<1e-9
        results.append(dict(hypothesis={k:scenario[k] for k in
            ('rho_input','kstar_input','phi_input','branch','source_hypothesis')},**blocks))
    paths=['scripts/audit_alp_blind_direction_span.py','scripts/audit_alp_vector_blindness.py',
           'scripts/audit_alp_gamma_gamma_blindness.py','scripts/audit_alp_three_p_odd_parity.py',
           'data/vector_sector/stage4_rank_results.json']
    output=ROOT/'data/vector_sector/alp_blind_direction_audit.json'
    output.write_text(json.dumps(dict(
        scope='projected one-ALP/two-boson operator sensitivity; not a complete normalized ALP predictor',
        assumptions='tree retained single-trace action; real photon/external-vector wedges; independent ALP derivative/mass traces; EM-compatible f_minus=0',
        coupling_order=[str(x) for x in vector['parameters']],
        null_columns=np.array(vector['nulls'],dtype=int).tolist(),
        exact_projected_zero=True,formal_probe_inputs=examples,
        homogeneous_ope_status='inherited comparison only, not rematched',
        sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        scenarios=results),indent=2,allow_nan=False)+'\n')
    print(output.relative_to(ROOT))


if __name__=='__main__':main()
