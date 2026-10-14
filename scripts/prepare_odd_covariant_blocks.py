"""Generic Lorentz bases before chiral Bianchi or derivative-commutator reduction.

Run .venv/bin/python scripts/prepare_odd_covariant_blocks.py.
Explicit Lorentz components; flavor matrices remain abstract. The two
ordered derivative slots are polarized independently, so their antisymmetric
part is retained rather than silently declared zero.
"""
from pathlib import Path
import contextlib
import hashlib
import json
import sympy as sp
import enumerate_odd_breaking_vertices as b

ROOT=Path(__file__).resolve().parents[1]
ra=sp.symbols('ra0:4');rb=sp.symbols('rb0:4')
exchange={**dict(zip(b.Avar,b.Bvar)),**dict(zip(b.Bvar,b.Avar)),
          **dict(zip(b.qa,b.qb)),**dict(zip(b.qb,b.qa)),
          **dict(zip(ra,rb)),**dict(zip(rb,ra))}

def show(label,expression):print(label+':\n$$\n'+sp.latex(expression)+'\n$$\n')

def lorentz(mode,kind):
    labels=[];polynomials=[]
    slots=(('A',b.qa),('B',b.qb),('X',b.qx)) if mode=='X' else \
          (('AA',(b.qa,ra)),('AB',(b.qa,b.qb)),('BB',(b.qb,rb)))
    for name,momenta in slots:
        for pair in b.PAIR:
            labels.append({'derivative_on':name,'metric_pair':list(pair)})
            polynomials.append(b.contraction(pair,b.AT,b.BT,b.u,momenta) if mode=='X'
                               else b.contraction(pair,b.AT,b.BT,*momenta))
    selected,coords=b.reduction(polynomials)
    result={'rank':len(selected),'basis_labels':[labels[i] for i in selected],
            'raw_labels':labels,'raw_to_basis':[[str(x) for x in coords[:,j]] for j in range(len(labels))]}
    if kind=='VVP':
        for parity,name in ((1,'symmetric'),(-1,'antisymmetric')):
            projected=[sp.expand(polynomials[j]+parity*polynomials[j].xreplace(exchange)) for j in selected]
            pivots,_=b.reduction(projected)
            result[name]={'rank':len(pivots),'basis_labels':[labels[selected[j]] for j in pivots]}
    return result

def main():
    transcript=ROOT/'data/vector_sector/odd_covariant_block_preparation.md'
    with transcript.open('w') as handle,contextlib.redirect_stdout(handle):
        print('# Ordered-derivative Lorentz bases before chiral reduction\n')
        result={}
        for mode in ('X','dd'):
            for kind in ('VJP','VVP'):
                name=kind+'_'+mode
                result[name]=lorentz(mode,kind)
                show(name+' generic rank',sp.Integer(result[name]['rank']))
                if kind=='VVP':
                    show(name+' symmetric/antisymmetric ranks',sp.Matrix([
                        result[name]['symmetric']['rank'],result[name]['antisymmetric']['rank']]))
    result['scope']='generic Lorentz preparation only; NOT independent chiral-action counts'
    paths=['scripts/prepare_odd_covariant_blocks.py','scripts/enumerate_odd_breaking_vertices.py']
    result['sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    output=ROOT/'data/vector_sector/odd_covariant_block_preparation.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(transcript.relative_to(ROOT));print(output.relative_to(ROOT))

if __name__=='__main__':main()
