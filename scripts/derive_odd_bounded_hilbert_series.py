"""Character checks and finite graded Hilbert series for the bounded odd action.

Run: .venv/bin/python scripts/derive_odd_bounded_hilbert_series.py
Natural units. V is an unconstrained canonical tensor nonet, chiral weight 0.
The quotient uses saved exact covariant relations: it is NOT an independent
Molien-Weyl proof of their completeness. Character integrals are independent.
Finite trace-orbit enumeration reuses the repository's tested metadata because
SymPy has no cyclic trace-word/trace-ceiling projector. All algebra uses SymPy.
"""
from pathlib import Path
import contextlib
import hashlib
import json
import sympy as sp
import enumerate_odd_breaking_vertices as b
import reduce_odd_covariant_higher as r

ROOT = Path(__file__).resolve().parents[1]
t, j, v = sp.symbols('t j v')


def show(label, expression):
    print(label + ':\n$$\n' + sp.latex(expression) + '\n$$\n')


def ct(expression, variables):
    """Exact constant term of a finite Laurent polynomial, via SymPy expansion."""
    expression = sp.expand(expression)
    for variable in variables:
        expression = expression.coeff(variable, 0)
    return sp.expand(expression)


def relation_matrix(ids, relations):
    return sp.Matrix([[sp.sympify(row.get(identifier, 0)) for identifier in ids]
                      for row in relations])


def leading_specs():
    families = b.trace_catalog(('A', 'B', 'U'))
    result = []
    def add(identifier, kind, family):
        result.append(dict(id=identifier, kind=kind, level=0,
                           traces=len(next(iter(family)))))
    for f, family in enumerate(families):
        for ell in range(6):
            add(f'VJP:F{f}:L{ell}', 'VJP', family)
    for f, family in enumerate(families):
        add(f'VJP:Y:F{f}', 'VJP', family)
    for sector, sign in [('S', 1), ('A', -1)]:
        for f, family in enumerate(b.flavor_exchange_basis(families, sign)):
            for ell in range(3):
                add(f'VVP:{sector}:F{f}:L{ell}', 'VVP', family)
    for f, family in enumerate(b.flavor_exchange_basis(families, 1)):
        add(f'VVP:Y:F{f}', 'VVP', family)
    add('VJP:fminus', 'VJP', b.trace_catalog(('A', 'B', 'X'), -1)[0])
    return result


def grade(spec):
    return (spec['kind'], (4 if spec['kind'] == 'VJP' else 2) + 2*spec['level'])


def monomial(spec):
    kind, order = grade(spec)
    return (j*v if kind == 'VJP' else v**2)*t**order


def image_dimension(selected, matrix):
    # dim image span(selected) in F/R = |selected| - dim(R intersect span(selected)).
    outside = [i for i in range(matrix.cols) if i not in selected]
    outside_rank = matrix[:, outside].rank() if outside else 0
    return len(selected) - matrix.rank() + outside_rank


def main():
    x, y = sp.symbols('x y')
    weights = (x, y, 1/(x*y))
    nonet = sum(a/c for a in weights for c in weights)
    measure = sp.prod(1-weights[k]/weights[i] for i in range(3) for k in range(i+1, 3))
    z = sp.Symbol('z')
    c_twisted_nonet = 1-z-1/z
    flavor = {}
    for n in (3, 4):
        ordinary = ct(nonet**n*measure, (x, y))
        twisted = ct(c_twisted_nonet**n*(1-z**2), (z,))
        flavor[str(n)] = dict(ordinary=int(ordinary), C_trace=int(twisted),
                             C_even=int((ordinary+twisted)/2))
    a, c = sp.symbols('a c')
    twoform = a**2+1+a**-2+c**2+1+c**-2
    vector = (a+1/a)*(c+1/c)
    haar = (1-a**2)*(1-c**2)
    tensor_odd = ct(twoform**2*vector**2*haar, (a, c))/2
    symvector = (vector**2+vector.xreplace({a:a**2, c:c**2}))/2
    symmetric_odd = ct(twoform**2*symvector*haar, (a, c))/2
    # The odd O(4) branch vanishes for each independent two-form. Thus /2
    # selects pseudoscalar contractions here. This is not a VVP Bose projection.
    assert tensor_odd == 3 and symmetric_odd == 2
    # Bose exchange is a trace insertion, not a division of the ordered count
    # by two. Trace[S(g tensor g)] = chi(g**2). On the disconnected O(4)
    # branch the square of the vector has eigenvalues z**2,z**-2,1,1.
    twoform_squared = twoform.xreplace({a:a**2, c:c**2})
    odd_twoform_squared = 2+2*(z**2+z**-2)
    odd_vector = z+1/z
    swap_X = (ct(twoform_squared*vector**2*haar, (a,c))
              -ct(odd_twoform_squared*odd_vector**2*(1-z**2), (z,)))/2
    swap_dd = (ct(twoform_squared*vector.xreplace({a:a**2,c:c**2})*haar, (a,c))
               -ct(odd_twoform_squared*(z**2+z**-2+2)*(1-z**2), (z,)))/2
    lorentz = {'leading': (2*tensor_odd, 0),
               'X': (3*tensor_odd, swap_X), 'dd': (3*tensor_odd, swap_dd)}
    bose = {name: ((dimension+trace)/2, (dimension-trace)/2)
            for name, (dimension, trace) in lorentz.items()}
    # Ordered dd is intentional. Symmetric free jets would give 2+3+2=7;
    # the two commutator slots retain curvature remainders in this seed model.
    families3 = b.trace_catalog(('A','B','U'))
    families4 = b.trace_catalog(('A','B','U','X'))
    f3 = [len(b.flavor_exchange_basis(families3, sign)) for sign in (1,-1)]
    f4 = [len(b.flavor_exchange_basis(families4, sign)) for sign in (1,-1)]
    fm = len(b.trace_catalog(('A','B','X'), -1))
    raw_character_counts = {
        'leading_VJP_derivative': len(families3)*lorentz['leading'][0],
        'leading_VVP_derivative': sum(a*b for a,b in zip(f3,bose['leading'])),
        'leading_VJP_Y': len(families3), 'leading_VVP_Y': f3[0],
        'leading_VJP_fminus_X': fm,
        'higher_VJP_X': len(families4)*lorentz['X'][0],
        'higher_VVP_X': sum(a*b for a,b in zip(f4,bose['X'])),
        'higher_VJP_Ydd': len(families3)*lorentz['dd'][0],
        'higher_VVP_Ydd': sum(a*b for a,b in zip(f3,bose['dd'])),
        'higher_VJP_fminus_Xdd': fm*lorentz['dd'][0]}
    character_raw = (j*v*t**4*sum(raw_character_counts[k] for k in
                     ('leading_VJP_derivative','leading_VJP_Y','leading_VJP_fminus_X'))
                     +v**2*t**2*(raw_character_counts['leading_VVP_derivative']+raw_character_counts['leading_VVP_Y'])
                     +j*v*t**6*sum(raw_character_counts[k] for k in
                     ('higher_VJP_X','higher_VJP_Ydd','higher_VJP_fminus_Xdd'))
                     +v**2*t**4*(raw_character_counts['higher_VVP_X']+raw_character_counts['higher_VVP_Ydd']))

    leading_path = ROOT/'data/vector_sector/odd_leading_action_count.json'
    higher_path = ROOT/'data/vector_sector/odd_three_flavor_action_quotient.json'
    leading = json.loads(leading_path.read_text())
    higher = json.loads(higher_path.read_text())
    specs = leading_specs()
    for spec in r.specifications():
        specs.append(dict(id=spec['id'], kind=spec['kind'], level=1,
                          traces=len(next(iter(spec['family']))), mode=spec['mode']))
    ids = [spec['id'] for spec in specs]
    relations = leading['relations'] + higher['all_higher_relations']
    matrix = relation_matrix(ids, relations)
    raw = sp.expand(sum(monomial(spec) for spec in specs))
    assert sp.expand(character_raw-raw) == 0
    syzygies = 0
    for relation in relations:
        supported = [spec for spec in specs if spec['id'] in relation]
        assert len({grade(spec) for spec in supported}) == 1
        syzygies += monomial(supported[0])
    assert matrix.rank() == len(relations)
    hilbert = sp.expand(raw-syzygies)
    # Reconstruct representatives by eliminating late IDs, without taking the
    # saved representative list as input. Equal ranks are basis-order invariant.
    # Choose a filtration-adapted basis: eliminate highest correction cost first.
    ordered_specs = sorted(enumerate(specs), key=lambda item:
                           (item[1]['level']+item[1]['traces']-1, item[0]))
    reverse_ids = [spec['id'] for _, spec in reversed(ordered_specs)]
    echelon, pivots = relation_matrix(reverse_ids, relations).rref()
    removed = {reverse_ids[p] for p in pivots}
    basis = [spec for spec in specs if spec['id'] not in removed]
    assert sp.expand(sum(monomial(spec) for spec in basis)-hilbert) == 0
    filtration = {}
    for ceiling in (0, 1, 2):
        selected = [i for i, spec in enumerate(specs)
                    if spec['level']+spec['traces']-1 <= ceiling]
        filtration[str(ceiling)] = dict(raw_parents=len(selected),
                                      quotient_image_dimension=image_dimension(selected, matrix))
    delta = sp.Symbol('delta')
    filtered_series = sp.expand(sum(delta**(spec['level']+spec['traces']-1) for spec in basis))
    for ceiling in (0,1,2):
        assert sum(filtered_series.coeff(delta,k) for k in range(ceiling+1)) == filtration[str(ceiling)]['quotient_image_dimension']
    topology_mixed = [i for i, rel in enumerate(relations)
                      if len({spec['traces'] for spec in specs if spec['id'] in rel}) > 1]
    output = ROOT/'data/vector_sector/odd_bounded_hilbert_series.json'
    paths = [leading_path, higher_path, Path(__file__),
             ROOT/'data/vector_sector/odd_covariant_block_preparation.json',
             ROOT/'scripts/enumerate_odd_breaking_vertices.py',
             ROOT/'scripts/reduce_odd_covariant_higher.py']
    payload = dict(method='finite graded quotient presentation; independent character checks only',
                   independent_Molien_checks=flavor,
                   character_derived_raw_counts={k:int(value) for k,value in raw_character_counts.items()},
                   Lorentz_Bose_counts={k:[int(value) for value in values] for k,values in bose.items()},
                   Lorentz_pseudoscalar_counts=dict(ordered=int(tensor_odd), symmetric=int(symmetric_odd)),
                   raw_series=str(raw), relation_series=str(sp.expand(syzygies)),
                   Hilbert_series=str(hilbert), total=int(hilbert.subs({t:1,j:1,v:1})),
                   relations_rank=matrix.rank(), trace_mixing_relation_indices=topology_mixed,
                   equal_correction_order_filtration=filtration, basis=basis,
                   associated_graded_correction_series=str(filtered_series),
                   sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    output.write_text(json.dumps(payload, indent=2)+'\n')
    transcript = ROOT/'data/vector_sector/odd_bounded_hilbert_series.md'
    with transcript.open('w') as handle, contextlib.redirect_stdout(handle):
        print('# Executed bounded Hilbert-series derivation\n')
        print('The exact relation presentation is inherited; this is not an independent Molien count of the action.\n')
        show('Nonet character', nonet)
        show('SU(3) constant-term measure', measure)
        show('C-twisted nonet character', c_twisted_nonet)
        for n, values in flavor.items():
            show('Unrestricted flavor invariants at degree '+n, sp.Matrix(list(values.values())))
        show('SO(4) two-form character', twoform)
        show('Ordered and symmetric two-vector pseudoscalars', sp.Matrix([tensor_odd, symmetric_odd]))
        show('Bose exchange traces in X and ordered dd spaces', sp.Matrix([swap_X, swap_dd]))
        for name, values in bose.items():
            show('Bose symmetric/antisymmetric Lorentz multiplicities '+name, sp.Matrix(values))
        show('Character-derived raw block multiplicities', sp.Matrix(list(raw_character_counts.values())))
        show('Raw bounded space', raw)
        show('Independent homogeneous relations', sp.expand(syzygies))
        show('Finite graded quotient Hilbert series', hilbert)
        show('Total operator count', hilbert.subs({t:1,j:1,v:1}))
        show('Filtration image dimensions through levels zero, one, two',
             sp.Matrix([filtration[str(k)]['quotient_image_dimension'] for k in (0,1,2)]))
        show('Associated graded equal-correction filtration', filtered_series)
        print('Trace-mixing relation indices:', topology_mixed)
        print('\nRepresentatives are in the companion JSON; the original normalized basis is preserved separately.\n')
    print(output.relative_to(ROOT))
    print(transcript.relative_to(ROOT))


if __name__ == '__main__':
    main()
