"""Joint neutral cubic parent maps and the on-shell derivative-counting audit.

Run .venv/bin/python scripts/derive_odd_power_counted_channels.py.
Natural units; (+---), epsilon0123=+1, all momenta incoming, tree level.
All 241 raw parents are filtered before projection. U,Y are independent
until the separately labelled LO EOM-compatible source specialization.
No data fit, vector EOM, numerical mass substitution or ALP UV choice.
"""
from pathlib import Path
from functools import lru_cache
import contextlib
import hashlib
import json
import pickle
import sympy as sp
import enumerate_odd_breaking_vertices as b
import reduce_odd_covariant_higher as r
import project_odd_higher_two_boson as hp
from derive_radiative_sm_amplitude import RHO, OMEGA, PHI, Q, CHI
from derive_odd_bounded_hilbert_series import relation_matrix

ROOT = Path(__file__).resolve().parents[1]
s,t,h = hp.s,hp.t,hp.h
w,Z = sp.symbols('m_a2 Z_a', real=True)
e,MV,M1,M2 = sp.symbols('e M_V M_1 M_2', positive=True)
F,fa = sp.symbols('F f_a', positive=True)
J,V = sp.symbols('J V')


def show(label,expression):
    print(label+':\n$$\n'+sp.latex(expression)+'\n$$\n')


def diagonal_weight(prefix):
    coefficients = sp.symbols(prefix+'_0:3')
    assert len(coefficients)==3
    return sum((c*T for c,T in zip(coefficients,(RHO,OMEGA,PHI))),sp.zeros(3))


@lru_cache(None)
def lorentz(mode,where,pair,group):
    if mode in ('leading','X'):
        return sp.expand(sum(c*x for c,x in zip(hp.project('X',where,pair,group),(s,t,h))))
    if mode == 'Y':
        raw = sp.expand(sum(sign*hp.AP[i,j]*hp.BP[k,l] for (i,j,k,l),sign in b.EPS))
        pivots,coords = b.reduction([hp.E,raw])
        assert pivots == (0,)
        return coords[0,1]
    return sp.expand(sum(c*x for c,x in zip(hp.project('Ydd',where,pair,group),(s,t,h))))


def specifications():
    saved = json.loads((ROOT/'data/vector_sector/odd_vertex_enumeration.json').read_text())
    f3 = b.trace_catalog(('A','B','U'))
    specs = []
    def add(identifier,kind,mode,family,label=None,group='all'):
        specs.append(dict(id=identifier,kind=kind,mode=mode,family=family,
                          label=label,group=group,level=0,
                          traces=len(next(iter(family)))))
    labels = saved['lorentz']['VJP']['basis_labels'] + [dict(derivative_on='B',metric_pair=[1,4])]
    for f,family in enumerate(f3):
        for ell,label in enumerate(labels):
            add(f'VJP:F{f}:L{ell}','VJP','leading',family,label)
    for f,family in enumerate(f3):
        add(f'VJP:Y:F{f}','VJP','Y',family)
    for key,group,sign in [('S','symmetric',1),('A','antisymmetric',-1)]:
        for f,family in enumerate(b.flavor_exchange_basis(f3,sign)):
            for ell,entry in enumerate(saved['lorentz']['VVP'][group]['basis_labels']):
                add(f'VVP:{key}:F{f}:L{ell}','VVP','leading',family,entry['parent'],group)
    for f,family in enumerate(b.flavor_exchange_basis(f3,1)):
        add(f'VVP:Y:F{f}','VVP','Y',family,group='symmetric')
    add('VJP:fminus','VJP','fminus',b.trace_catalog(('A','B','X'),-1)[0])
    for spec in r.specifications():
        specs.append(dict(spec,level=1,traces=len(next(iter(spec['family'])))))
    return specs


def flavor(family,matrices):
    return sp.expand(sum(value*sp.Mul(*(sp.trace(sp.MatMul(*(matrices[name] for name in word)).doit())
                                       for word in words)) for words,value in family.items()))


def kernel(spec,A,B,U,Y):
    if spec['mode'] in ('fmdd','fminus'):
        # Independently checked diagonal electromagnetic compatibility below.
        return sp.Integer(0)
    mode = spec['mode']
    family = spec['family']
    if mode in ('Y','Ydd'):
        family = {tuple(tuple('Y' if n=='U' else n for n in word) for word in words):value
                  for words,value in family.items()}
    matrix_flavor = flavor(family,dict(A=A,B=B,U=U,Y=Y,X=2*CHI))
    label = spec['label']
    polynomial = lorentz(mode,label['derivative_on'] if label else '',
                         tuple(label['metric_pair']) if label else (),spec['group'])
    # All tensors carry i wedge. u=i p U; each explicit derivative carries i q.
    # Derivative seeds: i**4=+1; Y seeds: i**2=-1; Ydd: i**4=+1.
    phase = -1 if mode=='Y' else 1
    # The leading quotient's scalar column uses a normalized Bose average,
    # whereas derivative columns already store their symmetrized polynomial.
    # The published d2 maps to twice this raw scalar column. Higher exported
    # representatives carry an explicit action factor two before both placements.
    bose = 2 if spec['kind']=='VVP' and spec['level']==1 else 1
    return sp.expand(phase*bose*matrix_flavor*polynomial)


def canonical_source(U,Y,compatible):
    if compatible:
        return U,-2*w*U+2*Z*sp.eye(3)
    return U,Y


def degree_normalization(spec):
    return MV**(-(1+2*spec['level']) if spec['kind']=='VJP' else -2*spec['level'])


def exact_map(polynomials,identifiers,variables=None):
    cache = Path('/private/tmp/alpaca-power-counted-channel-cache')
    cache.mkdir(exist_ok=True)
    key = hashlib.sha256(Path(__file__).read_bytes()+pickle.dumps((polynomials,identifiers,variables))).hexdigest()
    cached = cache/(key+'.pickle')
    if cached.exists():
        return pickle.loads(cached.read_bytes())
    if variables is None:
        matrix = b.coefficient_matrix(polynomials)
    else:
        expanded = [sp.Poly(value,*variables) for value in polynomials]
        monomials = sorted(set().union(*(set(poly.monoms()) for poly in expanded)))
        dictionaries = [poly.as_dict() for poly in expanded]
        matrix = sp.SparseMatrix(len(monomials),len(polynomials),{
            (i,k):poly[monomial] for k,poly in enumerate(dictionaries)
            for i,monomial in enumerate(monomials) if monomial in poly})
    coefficient_domain = (sp.QQ.algebraic_field(sp.sqrt(2)) if
                          any(value.has(sp.sqrt(2)) for value in matrix) else sp.QQ)
    parameters = tuple(sorted(matrix.free_symbols,key=str))
    if parameters:
        coefficient_domain = coefficient_domain.poly_ring(*parameters)
    domain = matrix.to_DM(domain=coefficient_domain)
    if parameters:
        # Exact rational evaluation only chooses candidate independent rows.
        # Symbolic polynomial reconstruction below proves the upper bound,
        # while the evaluation proves the matching lower bound. No sampled
        # saturation is accepted without this symbolic certificate.
        probe = matrix.subs({parameter:sp.Rational(i+2) for i,parameter in enumerate(parameters)})
        _,rows = probe.T.to_DM().to_field().rref()
        small = domain.extract(list(rows),list(range(matrix.cols)))
        numerator,denominator,pivots = small.rref_den(method='FF')
        rank = len(pivots)
        assert rank == len(rows)
        reconstructed = domain.extract(list(range(matrix.rows)),list(pivots))*numerator.extract(
            list(range(rank)),list(range(matrix.cols)))
        difference = reconstructed-domain.scalarmul(denominator)
        if not difference.is_zero_matrix:
            # A special evaluation can miss directions. Fall back to the full
            # symbolic calculation rather than treating the probe as a proof.
            numerator,denominator,pivots = domain.rref_den(method='FF')
            rank = len(pivots)
            reconstructed = domain.extract(list(range(matrix.rows)),list(pivots))*numerator.extract(
                list(range(rank)),list(range(matrix.cols)))
            assert (reconstructed-domain.scalarmul(denominator)).is_zero_matrix
        reduced = numerator.to_Matrix()/coefficient_domain.to_sympy(denominator)
    else:
        reduced_domain,pivots = domain.rref()
        reduced = reduced_domain.to_Matrix()
        rank = len(pivots)
        reconstructed = domain.extract(list(range(matrix.rows)),list(pivots))*reduced_domain.extract(
            list(range(rank)),list(range(matrix.cols)))
        assert (reconstructed-domain).is_zero_matrix
    rank = len(pivots)
    result = matrix,dict(rank=rank,zero_images=sum(p==0 for p in polynomials),
                      combinations=[dict(representative=identifiers[p],coefficients={
                          identifiers[k]:str(reduced[i,k]) for k in range(len(identifiers)) if reduced[i,k]})
                          for i,p in enumerate(pivots)])
    cached.write_bytes(pickle.dumps(result))
    return result


def baseline_check(specs,A,TP):
    """Independent comparison with the established C,H physical formulae."""
    leading = json.loads((ROOT/'data/vector_sector/odd_leading_action_count.json').read_text())
    by_id = {spec['id']:spec for spec in specs}
    f3 = b.trace_catalog(('A','B','U'))
    single = next(k for k,family in enumerate(f3) if len(next(iter(family)))==1)
    sf = b.flavor_exchange_basis(f3,1)
    vsingle = next(k for k,family in enumerate(sf) if len(next(iter(family)))==1)
    U = -sp.sqrt(2)*TP/F
    Y = sp.sqrt(2)*(TP*CHI+CHI*TP)/F
    c = sp.symbols('c_1:8')
    d = sp.symbols('d_1:5')
    JP = sp.Integer(0)
    for coefficient,coordinates in zip((c[0],c[1],c[4],c[5],c[6]),leading['original_VJP_derivatives_to_generic_basis']):
        JP += 2*coefficient*sum(sp.Rational(value)*kernel(by_id[f'VJP:F{single}:L{k}'],A,Q,U,Y)
                               for k,value in enumerate(coordinates))
    JP += 2*c[2]*kernel(by_id[f'VJP:Y:F{single}'],A,Q,U,Y)
    VV = sp.Integer(0)
    for coefficient,coordinates in zip((d[0],d[2],d[3]),leading['original_VVP_derivatives_to_symmetric_basis']):
        VV += coefficient*sum(sp.Rational(value)*kernel(by_id[f'VVP:S:F{vsingle}:L{k}'],A,Q,U,Y)
                               for k,value in enumerate(coordinates))
    VV += 2*d[1]*kernel(by_id[f'VVP:Y:F{vsingle}'],A,Q,U,Y)
    JP = sp.expand((-2*e*JP/(M1*MV)).subs(h,(w-s-t)/2))
    fq = sp.trace((A*Q+Q*A)*TP)
    fh = sp.trace((A*Q+Q*A)*(TP*CHI+CHI*TP))
    C = fq*((c[0]+c[1]-c[4])*w+(-c[0]+c[1]+c[4]-2*c[5])*s
             +(c[0]-c[1]+c[4])*t)+4*c[2]*fh
    assert sp.expand(JP-2*sp.sqrt(2)*e*C/(F*M1*MV)) == 0
    B = diagonal_weight('baselineB')
    # Rebuild VV with B instead of Q; JP deliberately retained Q above.
    VV = sp.Integer(0)
    for coefficient,coordinates in zip((d[0],d[2],d[3]),leading['original_VVP_derivatives_to_symmetric_basis']):
        VV += coefficient*sum(sp.Rational(value)*kernel(by_id[f'VVP:S:F{vsingle}:L{k}'],A,B,U,Y)
                               for k,value in enumerate(coordinates))
    VV += 2*d[1]*kernel(by_id[f'VVP:Y:F{vsingle}'],A,B,U,Y)
    VV = sp.expand((sp.sqrt(2)*e*VV/M1).subs(h,(w-s-t)/2))
    fq = sp.trace((A*B+B*A)*TP)
    fh = sp.trace((A*B+B*A)*(TP*CHI+CHI*TP))
    H = fq*((d[0]-d[2])*w+d[2]*(s+t))+4*d[1]*fh
    assert sp.expand(VV+4*e*H/(F*M1)) == 0
    return dict(original_VJP_recovered=True,original_VVP_recovered=True,
                scalar_VJP_mass_column_factor=2,scalar_VVP_mass_column_factor=2)


def main():
    specs = specifications()
    ids = [spec['id'] for spec in specs]
    assert len(ids)==241 and len(set(ids))==241
    A,B,U,Y = (diagonal_weight(prefix) for prefix in ('vA','vB','uA','yA'))
    for source in (U,Y):
        assert source*Q-Q*source == sp.zeros(3)
    leading = json.loads((ROOT/'data/vector_sector/odd_leading_action_count.json').read_text())
    higher = json.loads((ROOT/'data/vector_sector/odd_three_flavor_action_quotient.json').read_text())
    relations = leading['relations']+higher['all_higher_relations']
    R = relation_matrix(ids,relations)
    images = {}
    maps = {}
    residuals = {}
    for compatible,label in [(False,'independent_sources'),(True,'LO_EOM_compatible_sources')]:
        uu,yy = canonical_source(U,Y,compatible)
        values = [sp.expand((J if spec['kind']=='VJP' else V)*kernel(spec,A,Q if spec['kind']=='VJP' else B,uu,yy)
                            .subs(h,(w-s-t)/2)) for spec in specs]
        images[label] = {identifier:str(value) for identifier,value in zip(ids,values)}
        matrix,result = exact_map(values,ids)
        relation_residuals = [sp.expand(sum(sp.sympify(value)*values[ids.index(identifier)]
                                         for identifier,value in relation.items())) for relation in relations]
        residuals[label] = {str(i):str(value) for i,value in enumerate(relation_residuals) if value!=0}
        if compatible:
            assert not residuals[label], residuals[label]
            assert matrix*R.T == sp.zeros(matrix.rows,R.rows)
        result['nested'] = {}
        for scheme in ('all_derivatives_soft','hard_vector_derivatives_promoted'):
            result['nested'][scheme] = {}
            for ceiling in (0,1,2):
                selected = [k for k,spec in enumerate(specs) if
                            ((0 if scheme=='hard_vector_derivatives_promoted' and spec['mode'] in ('Ydd','fmdd')
                              else spec['level'])+spec['traces']-1)<=ceiling]
                if len(selected)==len(ids):
                    filtered = {key:result[key] for key in ('rank','zero_images','combinations')}
                else:
                    _,filtered = exact_map([values[k] for k in selected],[ids[k] for k in selected])
                result['nested'][scheme][str(ceiling)] = dict(raw_parents=len(selected),**filtered)
        maps[label] = result
    # Scalar-source normalization is derived from the independent-source expansion.
    TP = sp.diag(*sp.symbols('t_Pu t_Pd t_Ps'))
    S = sp.diag(*sp.symbols('S_au S_ad S_as'))
    AP = sp.diag(*sp.symbols('c_APu c_APd c_APs'))
    sm_U = -sp.sqrt(2)*TP/F
    sm_Y = sp.sqrt(2)*(TP*CHI+CHI*TP)/F
    alp_U = -2*S/fa
    alp_Y = 2*((AP+S)*CHI+CHI*(AP+S))/fa
    baseline = baseline_check(specs,A,TP)
    alp_source_defect = sp.expand(alp_Y+2*w*alp_U-2*Z*sp.eye(3))
    # A supplied ALP source pair must satisfy this matrix condition before the
    # EOM-compatible maps can be used; no independent UV coefficient is fixed.
    q0,q1,q2 = (sp.Matrix(sp.symbols(prefix+'_0:4')) for prefix in ('p1','p2','p3'))
    incoming = -(q0+q1+q2)
    three_p_epsilon = sp.expand(sp.Matrix.hstack(incoming,q0,q1,q2).det(method='berkowitz'))
    assert sp.expand(three_p_epsilon)==0
    assert hp.AP.subs(dict(zip(b.ea,b.qa)),simultaneous=True)==sp.zeros(4)
    assert hp.BP.subs(dict(zip(b.ev,b.qb)),simultaneous=True)==sp.zeros(4)
    for spec in specs:
        if spec['kind']=='VVP':
            assert sp.expand(kernel(spec,A,B,U,Y)-kernel(spec,B,A,U,Y).xreplace({s:t,t:s}))==0
    # The on-shell prescription is an exact polynomial decomposition, with no
    # propagator expansion or refit of lower-order source coefficients.
    ds,dt,dh = sp.symbols('delta_s delta_t delta_h')
    s0,t0,h0 = sp.symbols('s_0 t_0 h_0')
    cs,ct,ch = sp.symbols('k_s k_t k_h')
    polynomial = cs*s+ct*t+ch*h
    expanded = sp.expand(polynomial.subs({s:s0+ds,t:t0+dt,h:h0+dh}))
    hard = expanded.subs({ds:0,dt:0,dh:0})
    soft = sp.expand(expanded-hard)
    assert sp.expand(expanded-hard-soft)==0
    normalized = {}
    W0,W1,W2 = (diagonal_weight(prefix) for prefix in ('W0','W1','W2'))
    for spec in specs:
        norm = degree_normalization(spec)
        if spec['kind']=='VJP':
            vg = -2*e*norm*kernel(spec,A,Q,U,Y)/M1
            gg = -2*sp.sqrt(2)*norm*(kernel(spec,W1,Q,U,Y)+
                   kernel(spec,W2,Q,U,Y).xreplace({s:t,t:s}))/1
            vv = 0
        else:
            vv = norm*kernel(spec,A,B,U,Y)/(M1*M2)
            vg = sp.sqrt(2)*e*norm*kernel(spec,A,W0,U,Y)/M1
            gg = 2*norm*kernel(spec,W1,W2,U,Y)
        normalized[spec['id']] = dict(aVV=str(sp.expand(vv)),aVgamma=str(sp.expand(vg)),
                                      aGammaGamma_over_e2=str(sp.expand(gg)),
                                      soft_parent_cost=spec['level']+spec['traces']-1,
                                      hard_promoted_cost=(0 if spec['mode'] in ('Ydd','fmdd') else spec['level'])+spec['traces']-1)
    # Fixed-input function-field ranks: chi, on-shell masses and W(0) are
    # coefficient-field parameters, not additional measurement directions.
    # U,Z and labelled neutral external generators are the variable sources.
    # Reparameterize only the arbitrary source-family coordinates by diagonal
    # entries. This invertible change removes sqrt(2) from the rank arithmetic;
    # actual meson generators and normalized physical gradients above are intact.
    AU = sp.diag(*sp.symbols('A_u A_d A_s'))
    BU = sp.diag(*sp.symbols('B_u B_d B_s'))
    UU = sp.diag(*sp.symbols('U_u U_d U_s'))
    response_d,response_s,mass_ratio = sp.symbols('r_d r_s mass_ratio')
    WW = sp.diag(1,response_d,response_s)
    uu,yy = canonical_source(UU,Y,True)
    snapshot = json.loads((ROOT/'data/vector_sector/alda_2026nrd_arxiv_v1.json').read_text())
    fixed_chi = {symbol:2*sp.Rational(str(value)) for symbol,value in
                 zip(CHI.diagonal(),snapshot['B0mq_gev2'])}
    channel_maps = {}
    variables = tuple(sorted(UU.free_symbols|AU.free_symbols|BU.free_symbols|{w,Z},key=str))
    for channel in ('aVV','aVgamma','aGammaGamma'):
        values = []
        for spec in specs:
            if channel=='aVV':
                value = kernel(spec,AU,BU,uu,yy) if spec['kind']=='VVP' else 0
                limits = {s:1,t:mass_ratio,h:(w-1-mass_ratio)/2}
            elif channel=='aVgamma':
                value = kernel(spec,AU,Q,uu,yy) if spec['kind']=='VJP' else kernel(spec,AU,WW,uu,yy)
                limits = {s:1,t:0,h:(w-1)/2}
            else:
                value = kernel(spec,WW,Q,uu,yy) if spec['kind']=='VJP' else kernel(spec,WW,WW,uu,yy)
                limits = {s:0,t:0,h:w/2}
            values.append(sp.expand(sp.sympify(value).subs(limits).subs(fixed_chi)))
        matrix,result = exact_map(values,ids,variables)
        assert (matrix*R.T).applyfunc(sp.expand) == sp.zeros(matrix.rows,R.rows)
        result['nested'] = {}
        for scheme in ('all_derivatives_soft','hard_vector_derivatives_promoted'):
            result['nested'][scheme] = {}
            for ceiling in (0,1,2):
                selected = [k for k,spec in enumerate(specs) if
                            ((0 if scheme=='hard_vector_derivatives_promoted' and spec['mode'] in ('Ydd','fmdd')
                              else spec['level'])+spec['traces']-1)<=ceiling]
                if len(selected)==len(ids):
                    filtered = {key:result[key] for key in ('rank','zero_images','combinations')}
                else:
                    _,filtered = exact_map([values[k] for k in selected],[ids[k] for k in selected],variables)
                result['nested'][scheme][str(ceiling)] = dict(raw_parents=len(selected),**filtered)
        result['scope'] = 'Alda:2026nrd central quark spurion fixed; generic source-family rank after invertible mass/response coordinate scaling; not a dataset or fixed ALP model rank'
        result['rank_coordinates'] = ('m_a2/M_1^2, Z_a/M_1^2, M_2^2/M_1^2' if channel=='aVV' else
            'm_a2/M_1^2, Z_a/M_1^2, W_d/W_u, W_s/W_u' if channel=='aVgamma' else
            'm_a2, Z_a, W_d/W_u, W_s/W_u')
        result['generic_domain'] = 'M_1 nonzero when present; W_u(0) nonzero for photon channels; special rank-drop loci require separate specialization'
        channel_maps[channel] = result
        (cache_status_path := Path('/private/tmp/alpaca-power-counted-channel-status.json')).write_text(
            json.dumps({key:dict(rank=value['rank'],nested={scheme:{k:entry['rank'] for k,entry in levels.items()}
                for scheme,levels in value['nested'].items()}) for key,value in channel_maps.items()},indent=2)+'\n')
    paths = ['scripts/derive_odd_power_counted_channels.py','scripts/project_odd_higher_two_boson.py',
             'scripts/derive_odd_bounded_hilbert_series.py','data/vector_sector/odd_vertex_enumeration.json',
             'data/vector_sector/odd_covariant_block_preparation.json',
             'data/vector_sector/odd_leading_action_count.json',
             'data/vector_sector/odd_three_flavor_action_quotient.json',
             'data/vector_sector/alda_2026nrd_arxiv_v1.json']
    payload = dict(scope='neutral cubic resonance kernels only; response-symbol polynomial ranks, not complete ALP predictors or data ranks',
                   invariants=dict(s='q_A^2',t='q_B^2',h='q_A dot q_B',m_a2='s+t+2h'),
                   source_convention='u_mu = i p_mu U a; i chi_minus = Y a; Z_TC = Z_a a',
                   EOM_source_specialization='Y = -2*m_a2*U + 2*Z_a*I',
                   maps=maps,relation_residuals=residuals,primitive_images=images,
                   neutral_on_shell_channel_maps=channel_maps,baseline_normalization_checks=baseline,
                   channel_rank_chi_gev2=[str(fixed_chi[symbol]) for symbol in CHI.diagonal()],
                   normalized_parent_kernels=normalized,
                   scalar_source_substitutions=dict(SM_U=sp.latex(sm_U),SM_Y=sp.latex(sm_Y),
                                                    ALP_U=sp.latex(alp_U),ALP_Y=sp.latex(alp_Y)),
                   ALP_source_EOM_defect=sp.latex(alp_source_defect),
                   target_checks=dict(a3P_single_odd_vacuum_epsilon=sp.latex(three_p_epsilon),
                                      photon_Ward_tensors_zero=True,VVP_Bose_exchange=True),
                   hard_soft_split=dict(expanded=str(expanded),hard=str(hard),soft=str(soft)),
                   deferred=['full fixed-response physical channel ranks','charged fminus kernels',
                             'anomaly/contact and LSZ basis-equivalence completion','coefficient calibration and numerical convergence'],
                   sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths})
    path = ROOT/'data/vector_sector/odd_power_counted_channels.json'
    path.write_text(json.dumps(payload,indent=2)+'\n')
    transcript = ROOT/'data/vector_sector/odd_power_counted_channels.md'
    with transcript.open('w') as handle,contextlib.redirect_stdout(handle):
        print('# Joint power-counted neutral cubic kernels\n')
        show('Scalar derivative source matrix',U)
        show('LO EOM-compatible mass source',-2*w*U+2*Z*sp.eye(3))
        show('SM derivative/mass sources',sp.BlockMatrix([[sm_U,sm_Y]]))
        show('ALP derivative/mass sources',sp.BlockMatrix([[alp_U,alp_Y]]))
        show('ALP source EOM defect: must vanish for the compatible maps',alp_source_defect)
        show('Single odd vacuum a to three scalar pseudoscalars epsilon',three_p_epsilon)
        show('Higher-derivative polynomial under a pole-reference split',expanded)
        show('Hard reference term',hard)
        show('Soft residual term',soft)
        for label,result in maps.items():
            show(label+' joint primitive rank',sp.Integer(result['rank']))
            for scheme,nested in result['nested'].items():
                show(label+' '+scheme+' nested ranks',sp.Matrix([nested[str(k)]['rank'] for k in (0,1,2)]))
        for channel,result in channel_maps.items():
            show(channel+' generic fixed-input cubic rank',sp.Integer(result['rank']))
            for scheme,nested in result['nested'].items():
                show(channel+' '+scheme+' nested ranks',sp.Matrix([nested[str(k)]['rank'] for k in (0,1,2)]))
        print('\nOriginal VJP/VVP normalized scalar kernels recovered exactly.\n')
        for label,values in residuals.items():
            print('\n'+label+' nonvanishing action-relation projections:\n')
            for i,value in values.items():show('Relation '+i,sp.sympify(value))
        print('\nAll normalized parent gradients and filtered coefficient maps are in the companion JSON.\n')
    print(path.relative_to(ROOT));print(transcript.relative_to(ROOT))


if __name__ == '__main__':
    main()
