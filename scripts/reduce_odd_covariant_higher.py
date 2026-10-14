"""Joint reduction with independent sources in the local U=I source gauge.

Run .venv/bin/python scripts/reduce_odd_covariant_higher.py.
All flavor objects are abstract 3x3 MatrixSymbols. SymPy supplies matrix
algebra, Frechet directional differentiation, Trace, Poly and exact rref.
Custom finite trace-orbit metadata is needed because Trace does not implement
cyclic-word canonicalization or IBP on matrix-valued field jets.
"""
from pathlib import Path
from functools import lru_cache
from itertools import product, permutations
import argparse
import contextlib
import hashlib
import json
import pickle
import sympy as sp
import enumerate_odd_breaking_vertices as b
from reduce_odd_goldstone_eom import image

ROOT=Path(__file__).resolve().parents[1]
REGISTRY={}
JET_CACHE={}
NC_REGISTRY={}
ZERO=sp.ZeroMatrix(3,3)
Z=sp.Symbol('Z_TC')
REGISTRY[Z]=('Z',())

def jet(base,indices=()):
    key=(base,tuple(sorted(indices)))
    if key not in JET_CACHE:
        name=base+('_d'+''.join(map(str,key[1])) if key[1] else '')
        result=sp.Symbol(name) if base=='Z' else sp.MatrixSymbol(name,3,3)
        REGISTRY[result]=key;JET_CACHE[key]=result
    return JET_CACHE[key]

JET_CACHE['Z',()]=Z

@lru_cache(None)
def matrix_polynomial(expression):
    """Adapter to SymPy's noncommutative polynomial differentiation.

    MatrixSymbol coordinate dependence and general MatPow Frechet derivatives
    are unavailable in SymPy. This adapter uses its existing NC Add/Mul/Pow
    algebra, not a replacement algebra or generic flavor components.
    """
    if isinstance(expression,sp.MatrixSymbol):
        atom=sp.Symbol(str(expression),commutative=False)
        NC_REGISTRY[atom]=expression
        return atom
    if isinstance(expression,sp.Identity):return sp.Integer(1)
    if isinstance(expression,sp.ZeroMatrix):return sp.Integer(0)
    if isinstance(expression,sp.MatAdd):return sp.Add(*(matrix_polynomial(x) for x in expression.args))
    if isinstance(expression,sp.MatMul):return sp.Mul(*(matrix_polynomial(x) for x in expression.args))
    if isinstance(expression,sp.MatPow):return matrix_polynomial(expression.base)**expression.exp
    assert not expression.is_Matrix,expression
    return expression

def matrix_from_polynomial(expression):
    terms=[]
    for term in sp.Add.make_args(sp.expand(expression)):
        coefficient,factors=term.args_cnc()
        matrices=[]
        for factor in factors:
            if isinstance(factor,sp.Pow):
                assert factor.exp.is_Integer and factor.exp>0
                matrices.extend([NC_REGISTRY[factor.base]]*int(factor.exp))
            else:matrices.append(NC_REGISTRY[factor])
        terms.append(sp.Mul(*coefficient)*sp.MatMul(*matrices) if matrices else
                     sp.Mul(*coefficient)*sp.Identity(3))
    return sp.expand(sp.MatAdd(*terms)) if terms else ZERO

@lru_cache(None)
def partial(expression,mu):
    if expression==ZERO or expression==0:return expression
    tangent=sp.Dummy('t')
    if expression.is_Matrix:
        polynomial=matrix_polynomial(expression)
        atoms=polynomial.atoms(sp.Symbol)
        mapping={}
        for atom in atoms:
            original=NC_REGISTRY.get(atom,atom)
            if original not in REGISTRY:continue
            base,indices=REGISTRY[original]
            derivative=jet(base,indices+(mu,))
            mapping[atom]=atom+tangent*(matrix_polynomial(derivative) if derivative.is_Matrix else derivative)
        derivative=sp.diff(sp.expand(polynomial.xreplace(mapping)),tangent).subs(tangent,0)
        return matrix_from_polynomial(derivative)
    # General Trace.diff is unavailable; use its linear Frechet variation,
    # then let SymPy differentiate the scalar polynomial in the tangent.
    mapping={tr:tr+tangent*sp.Trace(partial(tr.arg,mu)).doit() for tr in expression.atoms(sp.Trace)}
    mapping.update({s:s+tangent*jet(REGISTRY[s][0],REGISTRY[s][1]+(mu,))
                    for s in expression.free_symbols if s in REGISTRY and not s.is_Matrix})
    return sp.expand(expression.xreplace(mapping).diff(tangent).subs(tangent,0))

@lru_cache(None)
def covariant(expression,mu):
    v=jet(f'v{mu}')
    return sp.expand(partial(expression,mu)-sp.I*(v*expression-expression*v))

def V(i,j):
    if i==j:return ZERO
    return (1 if i<j else -1)*jet(f'V{min(i,j)}{max(i,j)}')

def u(mu):return 2*jet(f'a{mu}')

@lru_cache(None)
def fplus(i,j):
    if i==j:return ZERO
    vi,vj,ai,aj=(jet(f'{name}{mu}') for name,mu in (('v',i),('v',j),('a',i),('a',j)))
    return sp.expand(2*(partial(vj,i)-partial(vi,j)-sp.I*(vi*vj-vj*vi)-sp.I*(ai*aj-aj*ai)))

@lru_cache(None)
def fminus(i,j):
    if i==j:return ZERO
    vi,vj,ai,aj=(jet(f'{name}{mu}') for name,mu in (('v',i),('v',j),('a',i),('a',j)))
    return sp.expand(-2*(partial(aj,i)-partial(ai,j)-sp.I*(vi*aj-aj*vi)-sp.I*(ai*vj-vj*ai)))

@lru_cache(None)
def Y():
    return sp.expand(2*sum((b.G[mu]*covariant(u(mu),mu) for mu in range(4)),ZERO)+2*Z*sp.Identity(3))

@lru_cache(None)
def trace_expanded(matrices):
    polynomial=sp.expand(sp.Mul(*(matrix_polynomial(matrix) for matrix in matrices)))
    terms=[]
    for term in sp.Add.make_args(polynomial):
        coefficient,factors=term.args_cnc();word=[]
        for factor in factors:
            if isinstance(factor,sp.Pow):
                assert factor.exp.is_Integer and factor.exp>0
                word.extend([NC_REGISTRY[factor.base]]*int(factor.exp))
            else:word.append(NC_REGISTRY[factor])
        if word:
            # Finite cyclic orbit only; multiplication/expansion remain SymPy.
            word=min((word[i:]+word[:i] for i in range(len(word))),key=lambda w:tuple(map(str,w)))
            terms.append(sp.Mul(*coefficient)*sp.Trace(sp.MatMul(*word)))
        else:terms.append(3*sp.Mul(*coefficient))
    return sp.Add(*terms)

def flavor_expression(family,matrices):
    return sp.Add(*(coefficient*sp.Mul(*(trace_expanded(tuple(matrices[name] for name in word))
                                       for word in traces)) for traces,coefficient in family.items()))

def nc_jet(base,indices=()):return matrix_polynomial(jet(base,indices)) if base!='Z' else jet(base,indices)

@lru_cache(None)
def nc_partial(expression,mu):
    tangent=sp.Dummy('t');mapping={}
    for atom in expression.atoms(sp.Symbol):
        original=NC_REGISTRY.get(atom,atom)
        if original not in REGISTRY:continue
        base,indices=REGISTRY[original]
        mapping[atom]=atom+tangent*nc_jet(base,indices+(mu,))
    return sp.expand(sp.diff(sp.expand(expression.xreplace(mapping)),tangent).subs(tangent,0))

@lru_cache(None)
def nc_covariant(expression,mu):
    v=nc_jet(f'v{mu}')
    return sp.expand(nc_partial(expression,mu)-sp.I*(v*expression-expression*v))

def nc_V(i,j):return 0 if i==j else (1 if i<j else -1)*nc_jet(f'V{min(i,j)}{max(i,j)}')
def nc_u(mu):return 2*nc_jet(f'a{mu}')

@lru_cache(None)
def nc_fplus(i,j):
    if i==j:return sp.Integer(0)
    vi,vj,ai,aj=(nc_jet(f'{name}{mu}') for name,mu in (('v',i),('v',j),('a',i),('a',j)))
    return sp.expand(2*(nc_partial(vj,i)-nc_partial(vi,j)-sp.I*(vi*vj-vj*vi)-sp.I*(ai*aj-aj*ai)))

@lru_cache(None)
def nc_fminus(i,j):
    if i==j:return sp.Integer(0)
    vi,vj,ai,aj=(nc_jet(f'{name}{mu}') for name,mu in (('v',i),('v',j),('a',i),('a',j)))
    return sp.expand(-2*(nc_partial(aj,i)-nc_partial(ai,j)-sp.I*(vi*aj-aj*vi)-sp.I*(ai*vj-vj*ai)))

@lru_cache(None)
def nc_Y():return sp.expand(2*sum(b.G[mu]*nc_covariant(nc_u(mu),mu) for mu in range(4))+2*Z)

def local_operator(spec):
    family=spec['family'];mode=spec['mode'];kind=spec['kind'];label=spec['label']
    pair=tuple(label['metric_pair'])
    remaining=[i for i in range(6) if i not in pair]
    terms=[]
    for eps,sign in b.EPS:
        for metric in range(4):
            slots=[0]*6
            for index,value in zip(remaining,eps):slots[index]=value
            slots[pair[0]]=slots[pair[1]]=metric
            A=nc_V(slots[0],slots[1])
            B=nc_V(slots[2],slots[3]) if kind=='VVP' else \
              nc_fminus(slots[2],slots[3]) if mode=='fmdd' else nc_fplus(slots[2],slots[3])
            if A==0 or B==0:continue
            matrices={'A':A,'B':B}
            if mode=='X':
                matrices.update(U=nc_u(slots[4]),X=nc_jet('X'))
                where=label['derivative_on']
                if where=='B' and pair in ((0,4),(1,4)):
                    # All three Bianchi indices occupy epsilon slots. Their
                    # cyclic sum is three copies of this exact contraction.
                    matrices['B']=-sp.I*(nc_u(slots[5])*nc_fminus(slots[2],slots[3])
                                        -nc_fminus(slots[2],slots[3])*nc_u(slots[5]))/2
                else:matrices[where]=nc_covariant(matrices[where],slots[5])
            else:
                matrices['X' if mode=='fmdd' else 'Y']=nc_jet('X') if mode=='fmdd' else nc_Y()
                where=label['derivative_on']
                if where=='AB':
                    matrices['A']=nc_covariant(matrices['A'],slots[4])
                    matrices['B']=nc_covariant(matrices['B'],slots[5])
                else:
                    name=where[0]
                    matrices[name]=nc_covariant(nc_covariant(matrices[name],slots[5]),slots[4])
            terms.append(sign*b.G[metric]*flavor_expression(family,matrices))
    return sp.expand(sp.Add(*terms)*(2 if kind=='VVP' else sp.I if mode=='fmdd' else 1))

def parse_matrix_word(expression):
    coefficient=sp.Integer(1);word=[]
    for factor in expression.args if isinstance(expression,sp.MatMul) else (expression,):
        if factor in REGISTRY:word.append(REGISTRY[factor])
        elif isinstance(factor,sp.MatPow):
            assert factor.exp.is_Integer and factor.exp>0
            word.extend([REGISTRY[factor.base]]*int(factor.exp))
        elif isinstance(factor,sp.Identity):pass
        elif factor.is_commutative:coefficient*=factor
        else:raise ValueError(('unparsed matrix factor',factor))
    return coefficient,tuple(word)

def parse_term(term):
    coefficient=sp.Integer(1);words=[]
    for factor in sp.Mul.make_args(term):
        if isinstance(factor,sp.Trace):
            c,word=parse_matrix_word(factor.arg);coefficient*=c
            if word:words.append(('trace',word))
            else:coefficient*=3
        elif isinstance(factor,sp.Pow) and isinstance(factor.base,sp.Trace):
            assert factor.exp.is_Integer and factor.exp>0
            c,word=parse_matrix_word(factor.base.arg)
            coefficient*=c**factor.exp
            words.extend([('trace',word)]*int(factor.exp))
        elif factor in REGISTRY:
            assert REGISTRY[factor][0]=='Z'
            words.append(('scalar',(REGISTRY[factor],)))
        elif isinstance(factor,sp.Pow) and factor.base in REGISTRY:
            assert factor.exp.is_Integer and factor.exp>0
            words.extend([('scalar',(REGISTRY[factor.base],))]*int(factor.exp))
        else:coefficient*=factor
    assert not coefficient.free_symbols,coefficient
    return sp.expand(coefficient),tuple(words)

@lru_cache(None)
def orbit_normal_form(words):
    """Finite Reynolds orbit followed by total-momentum elimination.

    Cyclic rotations and scalar-trace permutations act on derivative momenta.
    Averaging their canonicalizing orbit before eliminating one momentum
    handles stabilizers, e.g. repeated undifferentiated fields in a trace.
    """
    options=[]
    for kind,word in words:
        rotations=[word[i:]+word[:i] for i in range(len(word))] if kind=='trace' else [word]
        minimum=min(tuple(base for base,indices in rotated) for rotated in rotations)
        options.append([(kind,minimum,rotated) for rotated in rotations
                        if tuple(base for base,indices in rotated)==minimum])
    mappings=set();root=None
    for rotated in product(*options):
        # At most double flavor traces; scalar topological factors may add
        # singleton entries. Enumerate tied trace orderings only.
        for ordering in permutations(range(len(rotated))):
            arranged=[rotated[i] for i in ordering]
            signature=tuple((kind,base) for kind,base,word in arranged)
            if signature!=tuple(sorted(signature)):continue
            if root is None:root=signature
            assert root==signature
            mappings.add(tuple(indices for kind,base,word in arranged for name,indices in word))
    legs=sum(len(word) for kind,word in root)
    if legs==1:
        return root,{():sp.Integer(1)} if not next(iter(mappings))[0] else {}
    momenta=sp.symbols(f'p0:{4*(legs-1)}')
    polynomial=0
    for mapping in mappings:
        term=sp.Integer(1)
        for leg,derivatives in enumerate(mapping):
            for mu in derivatives:
                term*=(-sum(momenta[4*(j-1)+mu] for j in range(1,legs))) if leg==0 \
                       else momenta[4*(leg-1)+mu]
        polynomial+=term/sp.Integer(len(mappings))
    return root,sp.Poly(sp.expand(polynomial),*momenta).as_dict()

def normal_form(expression,max_legs=None):
    output={}
    for term in sp.Add.make_args(expression):
        coefficient,words=parse_term(term)
        if coefficient==0:continue
        if max_legs is not None and sum(len(word) for kind,word in words)>max_legs:continue
        root,polynomial=orbit_normal_form(words)
        for monomial,value in polynomial.items():
            key=(root,monomial)
            output[key]=output.get(key,0)+coefficient*value
    return {key:sp.expand(value) for key,value in output.items() if sp.expand(value)!=0}

FLAT_VARIABLES=(*b.Avar,*b.Bvar,*b.u,*b.qa,*b.qb,*b.qx,*b.ev,Z)

def flat_column(spec):
    family=spec['family'];mode=spec['mode'];kind=spec['kind'];label=spec['label']
    B=b.BT if kind=='VVP' else b.BF
    if mode=='fmdd':B=-sp.Matrix(4,4,lambda i,j:b.qb[i]*b.u[j]-b.qb[j]*b.u[i])
    pair=tuple(label['metric_pair'])
    if mode=='X':
        polynomial=b.contraction(pair,b.AT,B,b.u,{'A':b.qa,'B':b.qb,'X':b.qx}[label['derivative_on']])
        if kind=='VVP':
            exchange={**dict(zip(b.Avar,b.Bvar)),**dict(zip(b.Bvar,b.Avar)),
                      **dict(zip(b.qa,b.qb)),**dict(zip(b.qb,b.qa))}
            polynomial=sp.expand(polynomial+(1 if spec['group']=='symmetric' else -1)*polynomial.xreplace(exchange))
        sector='X:'+kind
        pieces=[(polynomial,family)]
    else:
        where=label['derivative_on']
        U,D={'AA':(b.qa,b.qa),'AB':(b.qa,b.qb),'BB':(b.qb,b.qb)}[where]
        polynomial=b.contraction(pair,b.AT,B,U,D)
        if kind=='VVP':
            exchange={**dict(zip(b.Avar,b.Bvar)),**dict(zip(b.Bvar,b.Avar)),
                      **dict(zip(b.qa,b.qb)),**dict(zip(b.qb,b.qa))}
            polynomial=sp.expand(polynomial+(1 if spec['group']=='symmetric' else -1)*polynomial.xreplace(exchange))
        if mode=='fmdd':
            sector='X:fminus'
            pieces=[(sp.I*polynomial,family)]
        else:
            sector='Y:'+kind
            qY=tuple(-b.qa[i]-b.qb[i] for i in range(4))
            divergence=sum(b.G[i]*qY[i]*b.u[i] for i in range(4))
            replaced={b.canonical(tuple(tuple('U' if n=='Y' else n for n in word) for word in traces)):c
                      for traces,c in family.items()}
            contact=image({b.canonical(tuple(tuple('U' if n=='Y' else n for n in word) for word in traces)):c
                           for traces,c in family.items()},None)
            pieces=[(2*polynomial*divergence,replaced),(2*polynomial*Z,contact)]
    output={}
    for polynomial,flavor in pieces:
        for monomial,coefficient in sp.Poly(sp.expand(polynomial),*FLAT_VARIABLES).terms():
            for word,value in flavor.items():
                key=(sector,word,monomial)
                output[key]=output.get(key,0)+coefficient*value
    return {key:value for key,value in output.items() if value!=0}

def specifications():
    prepared=json.loads((ROOT/'data/vector_sector/odd_covariant_block_preparation.json').read_text())
    output=[]
    for mode,labels,parity in (('X',('A','B','U','X'),1),('Ydd',('A','B','Y'),1),('fmdd',('A','B','X'),-1)):
        family=b.trace_catalog(labels,parity)
        for kind in ('VJP','VVP') if mode!='fmdd' else ('VJP',):
            L=prepared[kind+('_X' if mode=='X' else '_dd')]
            for group,flavors in [('all',family)] if kind=='VJP' else \
                [('symmetric',b.flavor_exchange_basis(family,1)),('antisymmetric',b.flavor_exchange_basis(family,-1))]:
                entries=L['basis_labels'] if group=='all' else L[group]['basis_labels']
                for f,flavor in enumerate(flavors):
                    for l,label in enumerate(entries):
                        output.append({'id':f'{mode}:{kind}:{group}:F{f}:L{l}',
                                       'mode':mode,'kind':kind,'group':group,'family':flavor,'label':label})
    return output

def matrix_of(columns):
    keys=sorted(set().union(*(set(c) for c in columns)),key=repr)
    indices={key:i for i,key in enumerate(keys)}
    return keys,sp.SparseMatrix(len(keys),len(columns),
                  {(indices[key],j):value for j,column in enumerate(columns) for key,value in column.items()})

def registry_restore(expression):
    for atom in expression.atoms(sp.MatrixSymbol)|expression.free_symbols:
        name=str(atom)
        if name=='Z_TC':REGISTRY[atom]=('Z',());JET_CACHE['Z',()]=atom
        elif name.startswith(('V','v','a','X','Z')):
            base,separator,indices=name.partition('_d')
            if isinstance(atom,sp.MatrixSymbol) or base=='Z':
                key=(base,tuple(map(int,indices)) if separator else ())
                REGISTRY[atom]=key;JET_CACHE[key]=atom

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--kernel',type=int)
    parser.add_argument('--max-legs',type=int)
    args=parser.parse_args()
    specs=specifications()
    columns=[flat_column(spec) for spec in specs]
    keys,matrix=matrix_of(columns)
    domain_reduced,pivots=matrix.to_DM().to_field().rref()
    reduced=domain_reduced.to_Matrix()
    null=[]
    for free in (j for j in range(matrix.cols) if j not in pivots):
        vector=sp.zeros(matrix.cols,1);vector[free]=1
        for row,pivot in enumerate(pivots):vector[pivot]=-reduced[row,free]
        assert matrix*vector==sp.zeros(matrix.rows,1)
        null.append(vector)
    assert matrix[:,list(pivots)]*reduced[:len(pivots),:]==matrix
    cache=Path('/private/tmp/alpaca-odd-covariant-cache')
    cache.mkdir(exist_ok=True)
    digest=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]
    metadata={'ordered_covariant_seed_count':len(specs),'flat_source_EOM_rank':len(pivots),
              'kernel_dimension':len(null),'kernel_vectors':[{specs[i]['id']:str(value) for i,value in enumerate(v) if value}
                                                          for v in null],
              'joint_action_count':None,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    preparation=ROOT/'data/vector_sector/odd_covariant_higher_kernel.json'
    preparation.write_text(json.dumps(metadata,indent=2)+'\n')
    if args.prepare_only:
        print(preparation.relative_to(ROOT));return
    selected=range(len(null)) if args.kernel is None else [args.kernel]
    results=[]
    for j in selected:
        path=cache/f'{digest}-kernel-{j}-legs-{args.max_legs}.pickle'
        if path.exists():
            results.append(pickle.loads(path.read_bytes()));continue
        column={}
        for i,coefficient in enumerate(null[j]):
            if coefficient==0:continue
            path=cache/f'{digest}-operator-{i}.pickle'
            if path.exists():
                expression=pickle.loads(path.read_bytes());registry_restore(expression)
            else:
                expression=local_operator(specs[i]);path.write_bytes(pickle.dumps(expression))
            nf=normal_form(expression,args.max_legs)
            for key,value in nf.items():column[key]=column.get(key,0)+coefficient*value
        column={key:sp.expand(value) for key,value in column.items() if sp.expand(value)!=0}
        path=cache/f'{digest}-kernel-{j}-legs-{args.max_legs}.pickle'
        path.write_bytes(pickle.dumps(column));results.append(column)
        progress=ROOT/'data/vector_sector/odd_covariant_higher_progress.json'
        progress.write_text(json.dumps({'completed_kernel':j,'normal_form_rows':len(column),
                                      'kernel_dimension':len(null),'script_digest':digest})+'\n')
    if args.kernel is None:
        remainder_keys,remainder=matrix_of(results)
        _,rp=remainder.to_DM().to_field().rref()
        metadata.update(remainder_rank=len(rp),formal_trace_quotient_rank=len(pivots)+len(rp),
                        source_grade_ceiling=args.max_legs,remainder_rows=len(remainder_keys))
        preparation.write_text(json.dumps(metadata,indent=2)+'\n')
    print(preparation.relative_to(ROOT))

if __name__=='__main__':main()
