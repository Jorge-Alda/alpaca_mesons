"""Constant axial redefinition of the sources and the graph values of sigma -> pi+ pi- (task 4.2).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_frame_invariance.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_frame_invariance

Conventions (memory/alp/shared/murayama-radial-frame-invariance.md): partonic
source M = S + i P expanded in y = sigma/f_a and x = a/f_a,
M(x, y) = sum_{r,s} M_rs y^r x^s, coefficients of mass dimension one (GeV).
Redefinition g_R = exp(i x kappa), g_L = g_R^dagger, with a constant Hermitian
kappa, [kappa, Q] = 0; sigma is not rotated.  Starting from the connection-free
frame (C_R = C_L = 0) the new frame has
    M' = g_R M g_L^dagger,  C'_R = kappa,  C'_L = -kappa,  Theta' = Theta - x <kappa>.
Transport along x at fixed y: d/dx W = i C W, W(0) = 1; dressed source
M_hat = W_R^dagger M W_L, Theta_hat = Theta + int_0^x <(C_R - C_L)/2> dx'.
This is a reimplementation for r <= 3 of the logic of
derive_murayama_radial_frame_invariance.py; the script checks that every source
letter of the kept graphs (S_hat_30, S_hat_21, P_hat_20, ... ) lies in this
range.  Redefinitions that depend on sigma are not covered.

Checks: (1) every dressed S_hat_rs, P_hat_rs (r <= 3, s <= 3) is the same in both
frames, for noncommuting symbols; (2) the values of all kept graphs (26
one-particle-irreducible graphs and the 6 + 5 + 5 self-energies, with the D4
terms removed), polynomials in the dressed symbols with M^{-1}_{m_a^2}(C) and
lambda_{pi+-} opaque, are unchanged;
(3) negative control: inserting the raw coefficients of the new frame changes
the graph values; (4) Theta_hat_01 is unchanged, so Theta_hat_1 = 0 is frame independent.
The graph values are the termwise routing maxima of the size-rule values with
the per-term loop scales (derive_sigma_pipi_loop_traces.py); the masses,
N(m_a^2) and the Piecewise conditions do not involve the sources, so only the
trace products are transformed.
No loop is evaluated.
"""
import itertools

import sympy as sp

import sigma_pipi_loop_algebra as algebra
from sigma_pipi_loop_algebra import TraceSum

KAPPA = 'kappa'


# Noncommutative polynomials: {word (tuple of letters): coefficient}.

def poly_add(*polys):
    total = {}
    for poly in polys:
        for word, value in poly.items():
            total[word] = total.get(word, 0) + value
    return {word: sp.expand(value) for word, value in total.items() if sp.expand(value) != 0}


def poly_scale(poly, factor):
    return {word: sp.expand(factor * value) for word, value in poly.items() if sp.expand(factor * value) != 0}


def poly_mul(first, second):
    result = {}
    for word_one, value_one in first.items():
        for word_two, value_two in second.items():
            word = word_one + word_two
            result[word] = result.get(word, 0) + value_one * value_two
    return {word: sp.expand(value) for word, value in result.items() if sp.expand(value) != 0}


def poly_adjoint(poly):
    """Hermitian conjugate; all letters are Hermitian."""
    return {tuple(reversed(word)): sp.conjugate(value) for word, value in poly.items()}


# Double series {(r, s): poly}.

def series_mul(*factors, radial_order, alp_order):
    result = factors[0]
    for factor in factors[1:]:
        product = {}
        for (r1, s1), left in result.items():
            for (r2, s2), right in factor.items():
                r, s = r1 + r2, s1 + s2
                if r <= radial_order and s <= alp_order:
                    product[(r, s)] = poly_add(product.get((r, s), {}), poly_mul(left, right))
        result = product
    return result


def series_adjoint(series):
    return {key: poly_adjoint(value) for key, value in series.items()}


def exponential(sign, alp_order):
    """exp(sign i x kappa)."""
    return {(0, s): {(KAPPA,) * s: (sign * sp.I) ** s / sp.factorial(s)} for s in range(alp_order + 1)}


def transport(connection, radial_order, alp_order):
    """d/dx W = i C W, W(0, y) = 1, order by order in x at fixed y."""
    result = {(0, 0): {(): sp.Integer(1)}}
    for s in range(alp_order):
        for r in range(radial_order + 1):
            value = {}
            for (r1, s1), coefficient in connection.items():
                key = (r - r1, s - s1)
                if key in result:
                    value = poly_add(value, poly_mul(coefficient, result[key]))
            value = poly_scale(value, sp.I / (s + 1))
            if value:
                result[(r, s + 1)] = value
    return result


def original_source(radial_order, alp_order):
    return {(r, s): {(f'S{r}{s}',): sp.Integer(1), (f'P{r}{s}',): sp.I}
            for r in range(radial_order + 1) for s in range(alp_order + 1)}


def hermitian_parts(series):
    scalar = {key: poly_scale(poly_add(value, poly_adjoint(value)), sp.Rational(1, 2))
              for key, value in series.items()}
    pseudoscalar = {key: poly_scale(poly_add(value, poly_scale(poly_adjoint(value), -1)), 1 / (2 * sp.I))
                    for key, value in series.items()}
    return scalar, pseudoscalar


def frames(radial_order=3, alp_order=3):
    """Raw and dressed coefficients of the original and the redefined frame."""
    orders = dict(radial_order=radial_order, alp_order=alp_order)
    source = original_source(radial_order, alp_order)
    new_source = series_mul(exponential(+1, alp_order), source, exponential(+1, alp_order), **orders)
    right_connection = {(0, 0): {(KAPPA,): sp.Integer(1)}}
    left_connection = {(0, 0): {(KAPPA,): sp.Integer(-1)}}
    right_transport = transport(right_connection, radial_order, alp_order)
    left_transport = transport(left_connection, radial_order, alp_order)
    dressed = series_mul(series_adjoint(right_transport), new_source, left_transport, **orders)
    theta_symbol = sp.Symbol('Theta_01', real=True)
    trace_kappa = sp.Symbol('<kappa>', real=True)
    theta_new = theta_symbol - trace_kappa
    axial_connection_trace = trace_kappa  # <(C'_R - C'_L)/2> = <kappa>, integrated over x gives x <kappa>
    return {'original': hermitian_parts(source), 'raw new': hermitian_parts(new_source),
            'dressed new': hermitian_parts(dressed),
            'theta': {'original dressed': theta_symbol, 'new dressed': theta_new + axial_connection_trace,
                      'new raw': theta_new}}


def dressed_differences(data):
    """{(kind, r, s): dressed new - original} (all must vanish)."""
    differences = {}
    for index, kind in enumerate(('S', 'P')):
        original, dressed = data['original'][index], data['dressed new'][index]
        for key in original:
            differences[(kind,) + key] = poly_add(dressed.get(key, {}), poly_scale(original[key], -1))
    return differences


def substitution_map(data, which, maximum_kappa=None):
    """Letter -> polynomial for the sources of a frame ('dressed new' or 'raw new')."""
    mapping = {}
    for index, kind in enumerate(('S', 'P')):
        for (r, s), poly in data[which][index].items():
            if maximum_kappa is not None:
                poly = {word: value for word, value in poly.items() if word.count(KAPPA) <= maximum_kappa}
            mapping[f'{kind}{r}{s}'] = poly
    return mapping


def substitute(trace_sum, mapping, maximum_kappa=None):
    """Replace letters by polynomials inside every trace (exact), optionally
    dropping terms with more than maximum_kappa letters kappa."""
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        word_polys = []
        for word in key:
            poly = {(): sp.Integer(1)}
            for name in word:
                factor = mapping.get(name, {(name,): sp.Integer(1)})
                poly = poly_mul(poly, factor)
                if maximum_kappa is not None:
                    poly = {w: v for w, v in poly.items() if w.count(KAPPA) <= maximum_kappa}
            word_polys.append(poly)
        for choice in itertools.product(*(list(poly.items()) for poly in word_polys)):
            count = sum(word.count(KAPPA) for word, _ in choice)
            if maximum_kappa is not None and count > maximum_kappa:
                continue
            coefficient = value
            for _, factor in choice:
                coefficient *= factor
            result.add_term(tuple(word for word, _ in choice), coefficient)
    return result.simplified()


def sources_covered(entries, radial_order=3, alp_order=3):
    """True if every source letter of the kept T_Gamma has r <= radial_order, s <= alp_order."""
    import derive_sigma_pipi_loop_traces as traces
    return all(int(name[1]) <= radial_order and int(name[2:]) <= alp_order for name in traces.sources_in(entries))


def trace_checks(entries, data):
    """(2) and (3) for every kept T_Gamma: {label: (dressed unchanged, raw changes)}."""
    dressed_map = substitution_map(data, 'dressed new')
    raw_map = substitution_map(data, 'raw new', maximum_kappa=1)
    results = {}
    for entry in entries:
        value = entry['words']
        dressed = substitute(value, dressed_map)
        raw = substitute(value, raw_map, maximum_kappa=1)
        results[entry['label']] = ((dressed - value).is_zero(), not (raw - value).is_zero())
    return results


def main():
    data = frames(radial_order=3, alp_order=3)
    differences = dressed_differences(data)
    nonzero = {key: value for key, value in differences.items() if value}
    print(f'(1) Dressed S_hat_rs, P_hat_rs, r <= 3, s <= 3: {len(differences)} differences, nonzero: {nonzero or "none"}')
    raw_changed = sorted(key for key in differences
                         if poly_add(data['raw new'][0 if key[0] == 'S' else 1][key[1:]],
                                     poly_scale(data['original'][0 if key[0] == 'S' else 1][key[1:]], -1)))
    print('    raw coefficients changed by the redefinition:', raw_changed)
    import derive_sigma_pipi_loop_traces as traces
    entries = traces.kept_graph_traces()
    print('\nSources of the kept graphs:', traces.sources_in(entries), '; all with r <= 3, s <= 3:',
          sources_covered(entries))
    print(f'\n(2) [{len(entries)} kept graphs] T_Gamma unchanged with dressed sources; (3) changed with raw sources (first order in kappa):')
    for label, (unchanged, changed) in trace_checks(entries, data).items():
        print(f'  {label}: dressed unchanged {unchanged}; raw changes {changed}')
    theta = data['theta']
    print('\n(4) Theta_hat_01 new - original:', sp.simplify(theta['new dressed'] - theta['original dressed']),
          '; raw Theta_01 new - original:', sp.simplify(theta['new raw'] - theta['original dressed']))
    print('\nM^{-1}_{m_a^2}(C), N(m_a^2) and lambda_{pi+-} are built from S_hat_00, P_hat_01 and Theta_hat_1 only,'
          ' so they are unchanged.  Redefinitions that depend on sigma are not covered.')


if __name__ == '__main__':
    main()
