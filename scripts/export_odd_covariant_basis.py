"""Export the certified covariant basis and executed derivation transcript.

Run .venv/bin/python scripts/export_odd_covariant_basis.py.
Representative normalization is specified here; no amplitude fit or matching
to physical new couplings is inferred from an operator-count certificate.
"""
from pathlib import Path
import contextlib
import hashlib
import json
import sympy as sp
import reduce_odd_covariant_higher as r

ROOT=Path(__file__).resolve().parents[1]

def show(label,expression):print(label+':\n$$\n'+sp.latex(expression)+'\n$$\n')

def representative(spec):
    pair=tuple(spec['label']['metric_pair']);remaining=[i for i in range(6) if i not in pair]
    indices=['']*6
    for slot,name in zip(remaining,(r'\mu',r'\nu',r'\rho',r'\sigma')):indices[slot]=name
    indices[pair[0]]=r'\alpha';indices[pair[1]]=r'\beta'
    names={'A':rf'V^{{{indices[0]}{indices[1]}}}',
           'B':rf'{"V" if spec["kind"]=="VVP" else "f_-" if spec["mode"]=="fmdd" else "f_+"}^{{{indices[2]}{indices[3]}}}'}
    mode=spec['mode'];where=spec['label']['derivative_on']
    if mode=='X':
        names.update(U=rf'u^{{{indices[4]}}}',X=r'\chi_+')
        names[where]=rf'\left(\nabla^{{{indices[5]}}} '+names[where]+r'\right)'
    else:
        names['X' if mode=='fmdd' else 'Y']=r'\chi_+' if mode=='fmdd' else r'i\chi_-'
        if where=='AB':
            names['A']=rf'\left(\nabla^{{{indices[4]}}} '+names['A']+r'\right)'
            names['B']=rf'\left(\nabla^{{{indices[5]}}} '+names['B']+r'\right)'
        else:
            block=where[0]
            names[block]=rf'\left(\nabla^{{{indices[4]}}}\nabla^{{{indices[5]}}} '+names[block]+r'\right)'
    matrices={name:sp.MatrixSymbol(latex,3,3) for name,latex in names.items()}
    flavor=sp.Add(*(value*sp.Mul(*(sp.Trace(sp.MatMul(*(matrices[name] for name in word)))
                                 for word in traces)) for traces,value in spec['family'].items()))
    phase=2 if spec['kind']=='VVP' else sp.I if mode=='fmdd' else 1
    return phase*sp.Symbol(r'\eta_{\alpha\beta}')*sp.Symbol(r'\epsilon_{\mu\nu\rho\sigma}')*flavor

def main():
    path=ROOT/'data/vector_sector/odd_three_flavor_action_quotient.json'
    final=json.loads(path.read_text());assert final['certified']
    leading=json.loads((ROOT/'data/vector_sector/odd_leading_action_count.json').read_text())
    specs={spec['id']:spec for spec in r.specifications()}
    inventory=[]
    for identifier in final['higher_basis_ids']:
        spec=specs[identifier]
        inventory.append({'id':identifier,'block':spec['mode']+':'+spec['kind'],
            'flavor_trace_words':[{';'.join(''.join(word) for word in traces):str(value)
                                  for traces,value in spec['family'].items()}],
            'flavor_trace_count':len(next(iter(spec['family']))),'spurion_count':1,
            'lorentz_label':spec['label'],'derivative_order':'one covariant derivative' if spec['mode']=='X'
                else 'ordered outer slot 4, inner slot 5 for AA/BB',
            'operator_mass_dimension':7 if spec['kind']=='VJP' else 6,
            'coefficient_normalization':'1/M_V^3' if spec['kind']=='VJP' else '1/M_V^2',
            'bose_sector':spec['group'],'hermitian_prefactor':'i' if spec['mode']=='fmdd' else '1',
            'representative_latex':sp.latex(representative(spec))})
    basis={'scope':'declared VJP/VVP covariant seed space only; not arbitrary multi-Goldstone U(3) odd action',
           'full_action_count':final['full_action_count'],'original_action_rank':11,
           'new_independent_operators':final['new_independent_operators'],
           'leading_basis':leading['basis_preserving_original_action'],
           'higher_basis':inventory,'higher_relations':final['all_higher_relations'],
           'normalization_note':'explicit new representatives; physical SM/ALP kernel normalization remains a separate calculation',
           'intermediate_projection_note':'odd_covariant_higher_kernel.json supplies a seed-dependent pre-completion bound; the final count is certified by odd_three_flavor_action_quotient.json',
           'sha256':{'data/vector_sector/odd_three_flavor_action_quotient.json':hashlib.sha256(path.read_bytes()).hexdigest(),
                     'scripts/export_odd_covariant_basis.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}
    output=ROOT/'data/vector_sector/odd_covariant_action_basis.json'
    output.write_text(json.dumps(basis,indent=2)+'\n')
    transcript=ROOT/'data/vector_sector/odd_covariant_action_derivation.md'
    with transcript.open('w') as handle,contextlib.redirect_stdout(handle):
        print('# Executed joint covariant reduction\n')
        show('Independent-source u in the U=I frame',r.matrix_from_polynomial(r.nc_u(0)))
        show('Independent-source f-plus component',r.matrix_from_polynomial(r.nc_fplus(0,1)))
        show('Independent-source f-minus component',r.matrix_from_polynomial(r.nc_fminus(0,1)))
        show('Source-preserving Goldstone EOM image of Y',r.matrix_from_polynomial(r.nc_Y()))
        show('Common source projection and retained remainder ranks',sp.Matrix([163,final['three_flavor_remainder_rank']]))
        n=sp.Symbol('N_f',integer=True,positive=True)
        polynomial=sp.sympify(final['candidate_Gram_polynomial'][0][0],locals={'N_f':n})
        show('Exact unresolved-direction squared norm, factored',sp.factor(polynomial))
        show('Squared norm at three flavors',polynomial.subs(n,3))
        show('Higher quotient dimensions',sp.Matrix([197,final['joint_higher_action_count'],197-final['joint_higher_action_count']]))
        show('Full and new operator counts',sp.Matrix([final['full_action_count'],final['original_action_rank'],final['new_independent_operators']]))
        symbols=sp.symbols(f'O0:{len(specs)}');identifiers=list(specs)
        for j,relation in enumerate(final['all_higher_relations']):
            expression=sum(sp.sympify(value)*symbols[identifiers.index(identifier)] for identifier,value in relation.items())
            show('Exact higher relation '+str(j)+' (equals zero)',expression)
    readable=ROOT/'data/vector_sector/odd_covariant_higher_basis.md'
    with readable.open('w') as handle:
        handle.write('# Certified higher covariant representatives\n\n')
        handle.write('Ordered covariant derivatives; metric $(+---)$ and $\\epsilon_{0123}=+1$. '
                     'Coefficients are normalized by the displayed power of $M_V$.\n\n')
        for item in inventory:
            handle.write('## '+item['id']+'\n\n$'+item['coefficient_normalization']+'$\n\n$$\n'+item['representative_latex']+'\n$$\n\n')
    print(output.relative_to(ROOT));print(transcript.relative_to(ROOT));print(readable.relative_to(ROOT))

if __name__=='__main__':main()
