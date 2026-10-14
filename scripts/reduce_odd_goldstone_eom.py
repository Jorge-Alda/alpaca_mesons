"""Covariant LO Goldstone EOM and exact flavor-multiplier image ranks.

Run .venv/bin/python scripts/reduce_odd_goldstone_eom.py.
Natural units; three flavors; independent sources; no resonance EOM.
Finite cyclic trace orbits reuse the catalogue implementation because SymPy
does not canonicalize cyclic words or products of flavor traces.
"""
from pathlib import Path
import contextlib
import hashlib
import json
import sympy as sp
from enumerate_odd_breaking_vertices import canonical, trace_catalog, flavor_exchange_basis

ROOT=Path(__file__).resolve().parents[1]

def show(label, expr):
    print(label+':\n$$\n'+sp.latex(expr)+'\n$$\n')

def eom():
    eta, divergence, Y = (sp.MatrixSymbol(n,3,3) for n in ('eta','nabla_dot_u','i_chi_minus'))
    F,m0,XTC=sp.symbols('F m_0 X_TC',real=True)
    # Covariant axial tangent delta U=i u eta u. Integration by parts
    # of delta u_mu=-nabla_mu eta gives the first expression.
    kinetic=F**2*sp.Trace(eta*divergence)/2
    # delta chi_plus=-{eta,Y}/2; trace cyclicity gives this variation.
    scalar=-F**2*sp.Trace(eta*Y)/4
    anomaly=-F**2*m0**2*XTC**2/12
    topological=sp.diff(anomaly,XTC)*(-sp.Trace(eta))
    for label,expr in [('Integrated kinetic variation',kinetic),('Scalar-source variation',scalar),('Topological variation',topological)]:show(label,expr)
    tangent_d,tangent_y,tangent_i=sp.symbols('tangent_D tangent_Y tangent_I')
    variation=(kinetic+scalar+topological).xreplace({
        sp.Trace(eta*divergence):tangent_d,sp.Trace(eta*Y):tangent_y,
        sp.Trace(eta):tangent_i})
    normalization=sp.diff(variation,tangent_d)
    residual=sp.Add(*(sp.diff(variation,t)/normalization*block
                     for t,block in ((tangent_d,divergence),(tangent_y,Y),
                                     (tangent_i,sp.Identity(3)))))
    show('Euler matrix after division by F squared over two',residual)
    return {'EOM_residual_latex':sp.latex(residual),
            'scope':'LO Goldstone action only; independent sources; no even resonance interaction variation'}

def image(terms, replacement):
    result={}
    for traces,coef in terms.items():
        kept=[]
        for word in traces:
            mapped=tuple(replacement if n=='U' else n for n in word)
            mapped=tuple(n for n in mapped if n is not None)
            if mapped:kept.append(mapped)
            else:coef*=3
        key=canonical(kept)
        result[key]=result.get(key,0)+coef
    return {k:v for k,v in result.items() if v}

def rank(columns):
    keys=sorted(set().union(*(set(c) for c in columns)))
    matrix=sp.Matrix([[c.get(k,0) for c in columns] for k in keys])
    reduced,pivots=matrix.rref()
    null=matrix.nullspace()
    assert all(matrix*v==sp.zeros(len(keys),1) for v in null)
    return len(pivots),null,keys,matrix

def multipliers():
    out={}
    for spurion in (False,True):
        base=trace_catalog(('A','B','U','X') if spurion else ('A','B','U'))
        for kind in ('VJP','VVP'):
            family=base if kind=='VJP' else flavor_exchange_basis(base,1)
            mass=[image(c,'Y') for c in family]
            contact=[image(c,None) for c in family]
            mr,_,_,_=rank(mass);cr,null,keys,matrix=rank(contact)
            # Tagged rows keep distinct mass and topological-source monomials.
            augmented=[{**{('mass',k):v/2 for k,v in a.items()},
                        **{('topological',k):-v/3 for k,v in b.items()}}
                       for a,b in zip(mass,contact)]
            ar,_,_,_=rank(augmented)
            name=kind+('_X' if spurion else '_leading')
            show(name+' topological-contact image matrix',matrix)
            show(name+' contact-free multiplier basis',sp.Matrix.hstack(*null) if null else sp.zeros(len(family),0))
            assert ar==mr==len(family)
            out[name]={'multipliers':len(family),'mass_image_rank':mr,
                       'topological_contact_rank':cr,'augmented_image_rank':ar,
                       'contact_free_dimension':len(null),
                       'mass_spurion_degree':2 if spurion else 1,
                       'mass_image_within_one_spurion_ceiling':not spurion,
                       'contact_free_null_vectors':[[str(x) for x in v] for v in null],
                       'contact_matrix':[[str(x) for x in matrix.row(i)] for i in range(matrix.rows)],
                       'contact_trace_rows':[repr(k) for k in keys]}
    return out

def main():
    transcript=ROOT/'data/vector_sector/odd_goldstone_eom_reduction.md'
    with transcript.open('w') as handle,contextlib.redirect_stdout(handle):
        print('# Covariant LO Goldstone EOM multiplier calculation\n')
        result={'goldstone':eom(),'divergence_multiplier_sectors':multipliers(),
                'full_operator_count':None}
    paths=['scripts/reduce_odd_goldstone_eom.py','scripts/enumerate_odd_breaking_vertices.py',
           'memory/rcht/odd-parity/odd-covariant-closure.md']
    result['sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    output=ROOT/'data/vector_sector/odd_goldstone_eom_reduction.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(transcript.relative_to(ROOT));print(output.relative_to(ROOT))

if __name__=='__main__':main()
