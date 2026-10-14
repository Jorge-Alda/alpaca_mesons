"""Exact source-resolved vertex enumeration for the bounded stage-A space.

Run: .venv/bin/python scripts/enumerate_odd_breaking_vertices.py
This reduces the specified one-u/one-chi_minus vertex sector, not a complete
nonlinear U(3) generating functional. No EOM or OPE constraint is imposed.
"""
from pathlib import Path
from itertools import permutations, combinations, product
import contextlib
import hashlib
import json
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
G=(1,-1,-1,-1)
EPS=[(p,sp.LeviCivita(*p)) for p in permutations(range(4))]
PAIR=list(combinations(range(6),2))
AB=list(combinations(range(4),2))
Avar=sp.symbols('A01 A02 A03 A12 A13 A23')
Bvar=sp.symbols('B01 B02 B03 B12 B13 B23')
u=sp.symbols('u0:4'); qa=sp.symbols('ka0:4'); qb=sp.symbols('kb0:4'); qx=sp.symbols('kx0:4')
ev=sp.symbols('eb0:4')
ea=sp.symbols('ea0:4')

def antisymmetric(values):
    M=sp.zeros(4)
    for (i,j),v in zip(AB,values): M[i,j]=v; M[j,i]=-v
    return M

AT=antisymmetric(Avar); BT=antisymmetric(Bvar)
BF=sp.Matrix(4,4,lambda i,j:qb[i]*ev[j]-qb[j]*ev[i])

def contraction(pair,A,B,U,D):
    """All epsilon metric contractions of two 2-forms and two vectors.

    Explicit Lorentz components only. Flavor matrices remain abstract.
    SymPy supplies LeviCivita, expansion, Poly and exact row reduction.
    """
    rest=[i for i in range(6) if i not in pair]
    out=0
    for ep,sign in EPS:
        for k in range(4):
            indices=[0]*6
            for slot,value in zip(rest,ep):indices[slot]=value
            indices[pair[0]]=indices[pair[1]]=k
            out+=sign*G[k]*A[indices[0],indices[1]]*B[indices[2],indices[3]]*U[indices[4]]*D[indices[5]]
    return sp.expand(out)

def coefficient_matrix(expressions):
    symbols=sorted(set().union(*(e.free_symbols for e in expressions)),key=str)
    if not symbols:
        return sp.Matrix([expressions])
    polynomials=[sp.Poly(e,*symbols) for e in expressions]
    monomials=sorted(set().union(*(set(p.monoms()) for p in polynomials)))
    dictionaries=[p.as_dict() for p in polynomials]
    return sp.SparseMatrix(len(monomials),len(expressions),
        {(i,j):d[m] for j,d in enumerate(dictionaries) for i,m in enumerate(monomials) if m in d})

def reduction(expressions):
    matrix=coefficient_matrix(expressions)
    echelon,pivots=matrix.rref()
    basis=matrix[:,list(pivots)]
    coordinates=echelon[:len(pivots),:]
    assert basis*coordinates==matrix
    return pivots,coordinates

def cyclic(word):
    return min(word[i:]+word[:i] for i in range(len(word)))

def canonical(traces):
    # Trace cyclic orbit and commuting scalar-trace order are not automatically
    # canonicalized by SymPy Trace. Only these finite label orbits are custom.
    return tuple(sorted(cyclic(tuple(word)) for word in traces))

def reverse(traces):return canonical(tuple(tuple(reversed(w)) for w in traces))

def trace_catalog(labels,charge_parity=1):
    catalog=set()
    for word in permutations(labels):
        catalog.add(canonical((word,)))
        for split in range(1,len(labels)):
            catalog.add(canonical((word[:split],word[split:])))
    orbits=[]
    for traces in sorted(catalog):
        rev=reverse(traces)
        if traces>rev:continue
        if traces==rev and charge_parity==-1:continue
        # Mean is C-even/Hermitian; odd case uses i times the antisymmetric
        # combination. The final imaginary unit is carried in the action.
        terms={traces:sp.Rational(1,2)}
        terms[rev]=terms.get(rev,0)+sp.Rational(charge_parity,2)
        terms={k:v for k,v in terms.items() if v}
        orbits.append(terms)
    return orbits

def trace_expression(terms):
    matrices={n:sp.MatrixSymbol(n,3,3) for n in ('A','B','U','X','Y')}
    return sp.Add(*(coef*sp.Mul(*(sp.Trace(sp.MatMul(*(matrices[n] for n in word))) for word in traces))
                    for traces,coef in terms.items()))

def swap_flavor(terms):
    return {canonical(tuple(tuple('B' if s=='A' else 'A' if s=='B' else s for s in w) for w in traces)):v
            for traces,v in terms.items()}

def flavor_exchange_basis(catalog,parity):
    keys=sorted(set().union(*(set(c)|set(swap_flavor(c)) for c in catalog)))
    candidates=[]
    for terms in catalog:
        switched=swap_flavor(terms)
        candidates.append({k:terms.get(k,0)+parity*switched.get(k,0) for k in keys})
    matrix=sp.Matrix([[col.get(k,0) for col in candidates] for k in keys])
    echelon,pivots=matrix.rref()
    assert matrix[:,list(pivots)]*echelon[:len(pivots),:]==matrix
    return [{k:v for k,v in candidates[j].items() if v} for j in pivots]

def chiral_trace_identity():
    # Polarized CH for 3x3 matrices: the fundamental degree-four trace identity.
    # All four variables remain abstract. Each permutation gives cycles/traces.
    labels=('A','B','U','X'); identity={}
    for perm in permutations(range(4)):
        seen=set(); cycles=[]
        for i in range(4):
            if i in seen:continue
            cycle=[]; j=i
            while j not in seen:
                seen.add(j); cycle.append(labels[j]); j=perm[j]
            cycles.append(tuple(cycle))
        key=canonical(cycles)
        identity[key]=identity.get(key,0)+(-1)**(4-len(cycles))
    identity={k:sp.Integer(v) for k,v in identity.items() if v}
    assert any(len(k)==4 for k in identity)
    # The sole multilinear degree-four FTI has a nonzero four-trace term,
    # hence cannot give a relation wholly inside our <=2-trace subspace.
    print('Polarized 3x3 Cayley-Hamilton trace identity (equals zero):\n$$\n'+sp.latex(trace_expression(identity))+'\n$$\n')
    return identity

def singlet_spurion_projection(catalog):
    base=trace_catalog(('A','B','U'))
    projected=[]
    for terms in catalog:
        result={}
        for traces,coef in terms.items():
            kept=[];factor=coef
            for word in traces:
                reduced=tuple(n for n in word if n!='X')
                if reduced:kept.append(reduced)
                else:factor*=3
            key=canonical(kept);result[key]=result.get(key,0)+factor
        projected.append(result)
    keys=sorted(set().union(*(set(c) for c in base+projected)))
    matrix=sp.Matrix([[c.get(k,0) for c in base+projected] for k in keys])
    echelon,pivots=matrix.rref()
    assert pivots==tuple(range(len(base)))
    return {'rank':len(base),'X_identity_to_zero_spurion_basis':
          [[str(v) for v in echelon[:len(base),len(base)+j]] for j in range(len(projected))]}

def lorentz_basis(kind,spurion):
    B=BF if kind=='VJP' else BT
    momenta=[('A',qa),('B',qb)]+([('X',qx)] if spurion else [])
    raw=[]; labels=[]
    for derivative,D in momenta:
        for pair in PAIR:
            raw.append(contraction(pair,AT,B,u,D));labels.append({'derivative_on':derivative,'metric_pair':list(pair)})
    # derivative on U eliminated by total momentum conservation: IBP.
    qU=tuple(-sum(momentum[i] for _,momentum in momenta) for i in range(4))
    for pair_index,pair in enumerate(PAIR):
        expression=contraction(pair,AT,B,u,qU)
        assert sp.expand(expression+sum(raw[j*15+pair_index] for j in range(len(momenta))))==0
    pivots,coordinates=reduction(raw)
    selected=[raw[i] for i in pivots]
    exchange={**dict(zip(Avar,Bvar)),**dict(zip(Bvar,Avar)),**dict(zip(qa,qb)),**dict(zip(qb,qa))}
    result={'raw_count':len(raw),'rank':len(pivots),'basis_labels':[labels[i] for i in pivots],
            'basis_polynomial_latex':[sp.latex(e) for e in selected],
            'all_raw_labels':labels,'all_raw_to_basis':[[str(v) for v in coordinates[:,j]] for j in range(len(raw))]}
    if kind=='VVP':
        for parity,key in ((1,'symmetric'),(-1,'antisymmetric')):
            projected=[sp.expand(e+parity*e.xreplace(exchange)) for e in selected]
            p,_=reduction(projected)
            result[key]={'rank':len(p),'basis_polynomial_latex':[sp.latex(projected[i]) for i in p],
                         'basis_labels':[{'parent':labels[pivots[i]],'exchange_parity':parity} for i in p]}
        assert result['symmetric']['rank']+result['antisymmetric']['rank']==len(pivots)
    if not spurion:
        # Explicit Fourier contractions of the published derivative seeds.
        seed=[contraction((3,5),AT,B,u,qU),contraction((1,5),AT,B,u,qU),
              contraction((3,5),AT,B,u,qa),contraction((1,5),AT,B,u,qa),
              contraction((3,4),AT,B,u,qa)]
        if kind=='VVP':
            seed=[sp.expand(seed[j]+seed[j].xreplace(exchange)) for j in (0,2,4)]
        _,combined=reduction(selected+seed)
        assert combined.rows==len(pivots)
        result['published_seed_to_basis']=[[str(v) for v in combined[:,len(selected)+j]] for j in range(len(seed))]
        result['published_seed_rank']=len(reduction(seed)[0])
        assert result['published_seed_rank']==(5 if kind=='VJP' else 3)
    # Distinguish vacuum spurions and source-resolved action before scalar P.
    if spurion:
        frozen=[sp.expand(e.subs(dict(zip(qx,[0]*4)))) for e in selected]
        p,coords=reduction(frozen)
        result['constant_spurion_projection']={'rank':len(p),'basis_indices':list(p),
                        'relations_to_selected':[[str(v) for v in coords[:,j]] for j in range(len(selected))]}
    # Canonical tensor legs and a scalar P: u is proportional to p_P, not
    # an independent axial source. These relations are NOT action reductions.
    AP=sp.Matrix(4,4,lambda i,j:qa[i]*ea[j]-qa[j]*ea[i])
    BP=sp.Matrix(4,4,lambda i,j:qb[i]*ev[j]-qb[j]*ev[i])
    sub={**dict(zip(Avar,[AP[i,j] for i,j in AB])),
         **dict(zip(Bvar,[BP[i,j] for i,j in AB])),
         **dict(zip(u,[-qa[i]-qb[i] for i in range(4)])),
         **dict(zip(qx,[0]*4))}
    projected=[sp.expand(e.subs(sub,simultaneous=True)) for e in selected]
    p,coords=reduction(projected)
    result['canonical_tensor_scalar_projection']={'rank':len(p),'basis_indices':list(p),
        'relations_to_selected':[[str(v) for v in coords[:,j]] for j in range(len(selected))]}
    if kind=='VVP':
        for parity,key in ((1,'symmetric'),(-1,'antisymmetric')):
            rawproj=[sp.expand((e+parity*e.xreplace(exchange)).subs(sub,simultaneous=True)) for e in selected]
            result[key]['canonical_tensor_scalar_rank']=len(reduction(rawproj)[0])
    return result

def projection_flavor(catalog,kind='VJP'):
    # Actual repository neutral/EM matrices: no generic component replacement.
    from derive_radiative_sm_amplitude import Q,CHI
    matrices={'A':sp.diag(*sp.symbols('a_u a_d a_s')),
              'B':Q if kind=='VJP' else sp.diag(*sp.symbols('b_u b_d b_s')),
              'U':sp.diag(*sp.symbols('p_u p_d p_s')),
              'Y':sp.diag(*sp.symbols('y_u y_d y_s')),'X':2*CHI}
    expressions=[]
    for item in catalog:
        expr=0
        for traces,coef in item.items():
            expr+=coef*sp.Mul(*(sp.trace(sp.MatMul(*(matrices[n] for n in w)).doit()) for w in traces))
        expressions.append(sp.expand(expr))
    pivots,coordinates=reduction(expressions)
    return {'rank':len(pivots),'basis_indices':list(pivots),
            'relations_to_selected':[[str(v) for v in coordinates[:,j]] for j in range(len(expressions))]}

def scalar_derivative_basis(kind,fminus=False):
    """Two additional derivatives on a one-spurion scalar block.

    chi_minus or chi_plus is scalar; its momentum eliminated by IBP.
    f_minus is not assumed closed: its chiral Bianchi relation contains u f_plus
    terms in the retained one-P photon sector and must not be discarded.
    """
    B=BT if kind=='VVP' or fminus else BF
    raw=[]; labels=[]
    for derivative_pair,U,D in (('AA',qa,qa),('AB',qa,qb),('BB',qb,qb)):
        for pair in PAIR:
            raw.append(contraction(pair,AT,B,U,D))
            labels.append({'derivative_pair':derivative_pair,'metric_pair':list(pair)})
    pivots,coords=reduction(raw);selected=[raw[i] for i in pivots]
    result={'raw_count':len(raw),'rank':len(pivots),'basis_labels':[labels[i] for i in pivots],
        'basis_polynomial_latex':[sp.latex(e) for e in selected],
        'all_raw_labels':labels,'all_raw_to_basis':[[str(v) for v in coords[:,j]] for j in range(len(raw))]}
    exchange={**dict(zip(Avar,Bvar)),**dict(zip(Bvar,Avar)),**dict(zip(qa,qb)),**dict(zip(qb,qa))}
    if kind=='VVP':
        for parity,key in ((1,'symmetric'),(-1,'antisymmetric')):
            projected=[sp.expand(e+parity*e.xreplace(exchange)) for e in selected]
            p,_=reduction(projected)
            result[key]={'rank':len(p),'basis_polynomial_latex':[sp.latex(projected[i]) for i in p],
                 'basis_labels':[{'parent':labels[pivots[i]],'exchange_parity':parity} for i in p]}
        assert sum(result[key]['rank'] for key in ('symmetric','antisymmetric'))==len(pivots)
    return result

def dictionary_audit():
    nf,MV=sp.symbols('n_f M_V',positive=True)
    d1,d2,c2,c3,c6=sp.symbols('d_1 d_2 c_2 c_3 c_6')
    kvv1,kvv2,kv13,kv14=sp.symbols('kappa_1_VV kappa_2_VV kappa_13_V kappa_14_V')
    equations=[kvv1+d1/(8*nf),kvv2-d1/8-d2,
               nf*MV*kv13-(-c2+c6)/4,MV*kv14-(c2+4*c3-c6)/4]
    solved=sp.solve(equations,[d1,d2,c2,c3],dict=True)[0]
    print('Literature basis dictionary, not an imposed EOM reduction:\n$$\n'+sp.latex(sp.Dict(solved))+'\n$$\n')
    return {str(k):sp.latex(v) for k,v in solved.items()}

def materialize(summary,catalogs,lorentz):
    result=[]
    for block,meta in summary.items():
        flavor_name='fminus' if 'fminus' in block else 'mass' if 'mass' in block else 'derivative_spurion' if 'spurion' in block else 'derivative'
        flavor=catalogs[flavor_name]
        groups=('symmetric','antisymmetric') if block.startswith('VVP') else ('all',)
        for group in groups:
            traces=flavor['trace_basis_latex'] if group=='all' else flavor[group]['trace_basis_latex']
            if block in lorentz:
                L=lorentz[block] if group=='all' else lorentz[block][group]
                labels=L['basis_labels']; expressions=L['basis_polynomial_latex']
            else:
                if group=='antisymmetric':continue
                labels=[{'epsilon_slots':'A(mu,nu) B(rho,sigma)'}]; expressions=['\\epsilon_{\\mu\\nu\\rho\\sigma}A^{\\mu\\nu}B^{\\rho\\sigma}']
            for fi,trace in enumerate(traces):
                for li,(label,expression) in enumerate(zip(labels,expressions)):
                    result.append({'id':f'{block}:{group}:F{fi}:L{li}','block':block,'exchange_parity':group,
                       'flavor_index':fi,'lorentz_index':li,'flavor_latex':trace,
                       'trace_count':(flavor['trace_counts'] if group=='all' else flavor[group]['trace_counts'])[fi],
                       'lorentz_label':label,'lorentz_reference':f'{block}/{group}/L{li}',
                       'hermitian_prefactor':'i' if 'fminus' in block else '1',
                       **{k:v for k,v in meta.items() if k!='independent_source_resolved_vertex_count'}})
    for block,meta in summary.items():
        assert sum(item['block']==block for item in result)==meta['independent_source_resolved_vertex_count']
    return result

def main():
    output=ROOT/'data/vector_sector/odd_vertex_enumeration.json'
    transcript=ROOT/'data/vector_sector/odd_vertex_enumeration.md'
    catalogs={}
    with transcript.open('w') as f,contextlib.redirect_stdout(f):
        print('# Executed flavor-trace and Lorentz enumeration\n')
        print('Scope: one-u/source or one-i-chi_minus sector; <=2 flavor traces; one spurion. No EOM or OPE.\n')
        for name,labels,parity in (('derivative',('A','B','U'),1),
                ('derivative_spurion',('A','B','U','X'),1),
                ('mass',('A','B','Y'),1),('fminus',('A','B','X'),-1)):
            catalog=trace_catalog(labels,parity)
            data={'count':len(catalog),'trace_basis_latex':[sp.latex(trace_expression(c)) for c in catalog],
                  'trace_counts':[len(next(iter(c))) for c in catalog],
                  'trace_words':[{';'.join(''.join(w) for w in key):str(value) for key,value in c.items()} for c in catalog]}
            print(name+' flavor basis:')
            for i,item in enumerate(catalog):print(str(i)+':\n$$\n'+sp.latex((sp.I if parity==-1 else 1)*trace_expression(item))+'\n$$\n')
            if name!='fminus':
                data['neutral_EM_projection']=projection_flavor(catalog)
                data['two_neutral_vectors_projection']=projection_flavor(catalog,'VVP')
                for sign,key in ((1,'symmetric'),(-1,'antisymmetric')):
                    exchange=flavor_exchange_basis(catalog,sign)
                    data[key]={'count':len(exchange),'trace_basis_latex':[sp.latex(trace_expression(c)) for c in exchange],
                               'trace_counts':[len(next(iter(c))) for c in exchange]}
                assert data['symmetric']['count']+data['antisymmetric']['count']==len(catalog)
                if name=='derivative_spurion':data['singlet_spurion_projection']=singlet_spurion_projection(catalog)
            catalogs[name]=data
        chiral_trace_identity()
        dictionary=dictionary_audit()
        lorentz={}
        for kind in ('VJP','VVP'):
            for spurion in (False,True):
                name=kind+('_spurion' if spurion else '')
                data=lorentz_basis(kind,spurion);lorentz[name]=data
                print(name+' Lorentz rank: '+str(data['rank']))
                for lab,expr in zip(data['basis_labels'],data['basis_polynomial_latex']):
                    print(str(lab)+':\n$$\n'+expr+'\n$$\n')
        for name,kind,fminus in (('VJP_mass_dd','VJP',False),('VVP_mass_dd','VVP',False),('VJP_fminus_dd','VJP',True)):
            data=scalar_derivative_basis(kind,fminus);lorentz[name]=data
            print(name+' Lorentz rank: '+str(data['rank']))
            for lab,expr in zip(data['basis_labels'],data['basis_polynomial_latex']):
                print(str(lab)+':\n$$\n'+expr+'\n$$\n')
        summary={}
        for kind in ('VJP','VVP'):
            for spurion in (False,True):
                name=kind+('_spurion' if spurion else '')
                flavor=catalogs['derivative_spurion' if spurion else 'derivative']; L=lorentz[name]
                count=flavor['count']*L['rank'] if kind=='VJP' else sum(flavor[key]['count']*L[key]['rank'] for key in ('symmetric','antisymmetric'))
                summary[name]={'independent_source_resolved_vertex_count':count,
                    'operator_mass_dimension':(5 if kind=='VJP' else 4)+(2 if spurion else 0),
                    'coefficient_normalization':'1/M_V^'+str((1 if kind=='VJP' else 0)+(2 if spurion else 0)),
                    'sources':'u_mu plus f_plus' if kind=='VJP' else 'u_mu',
                    'trace_ceiling':2,'spurion_count':int(spurion)}
        summary['VJP_mass']={'independent_source_resolved_vertex_count':catalogs['mass']['count'],'operator_mass_dimension':5,'coefficient_normalization':'1/M_V','spurion_count':1,'sources':'i chi_minus, f_plus'}
        summary['VVP_mass']={'independent_source_resolved_vertex_count':catalogs['mass']['symmetric']['count'],'operator_mass_dimension':4,'coefficient_normalization':'1','spurion_count':1,'sources':'i chi_minus'}
        summary['VJP_fminus']={'independent_source_resolved_vertex_count':catalogs['fminus']['count'],'operator_mass_dimension':5,'coefficient_normalization':'1/M_V','spurion_count':1,'sources':'f_minus, chi_plus'}
        for name,flavor_name,kind in (('VJP_mass_dd','mass','VJP'),('VVP_mass_dd','mass','VVP'),('VJP_fminus_dd','fminus','VJP')):
            flavor=catalogs[flavor_name];L=lorentz[name]
            count=flavor['count']*L['rank'] if kind=='VJP' else sum(flavor[key]['count']*L[key]['rank'] for key in ('symmetric','antisymmetric'))
            summary[name]={'independent_source_resolved_vertex_count':count,
                'operator_mass_dimension':7 if kind=='VJP' else 6,
                'coefficient_normalization':'1/M_V^3' if kind=='VJP' else '1/M_V^2',
                'spurion_count':1,'sources':'f_minus, chi_plus' if flavor_name=='fminus' else 'i chi_minus',
                'additional_derivatives':2}
        print('Source-resolved block counts:',summary)
    operators=materialize(summary,catalogs,lorentz)
    inputs=['scripts/enumerate_odd_breaking_vertices.py','scripts/derive_radiative_sm_amplitude.py',
            'memory/sources/Roig:2013baa/AcceptedVersion.tex','memory/sources/Ruiz-Femenia:2003jdx/vvp.tex']
    output.write_text(json.dumps({'scope':'source-resolved specified local vertex sector; NOT a complete nonlinear or EOM-minimal U(3) action',
        'flavor':catalogs,'lorentz':lorentz,'summary':summary,
        'operator_inventory':operators,'literature_dictionary_inverse_latex':dictionary,
        'sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs}},indent=2)+'\n')
    print(transcript.relative_to(ROOT));print(output.relative_to(ROOT))

if __name__=='__main__':main()
