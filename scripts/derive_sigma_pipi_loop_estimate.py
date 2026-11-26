"""Assembly of the size estimate M_1^est of the one-loop sigma -> pi+ pi- amplitude and of R (tasks 4.1, 4.4).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_estimate.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_estimate

Units GeV; conventions of sigma_pipi_loop_algebra.py.

    M_1^est = sum_{Gamma 1PI} X_Gamma^max + M_0 (dZ_sigma^max + dZ_{pi+}^max + dZ_{pi-}^max)/2,
    R = |M_1^est| / |M_0|,
with M_0 = -(B_0 eps/(2 f_pi)) <S_10 {lambda_{pi+}, lambda_{pi-}}>.  X_Gamma^max
is the termwise routing maximum (D12) of the size-rule value
sum_terms S_Gamma (m_term^2)^{2-r_Gamma} T_term/(16 pi^2), with the per-term
scale m_term of D3 (largest line scale: m_sigma radial, m_a ALP-like and the
non-generator pieces of a meson line, m_P the bare generator sum); dZ^max is
the termwise routing maximum of dZ = -Pi'(m^2) taken term by term at m_term
(derive_sigma_pipi_loop_traces.py: 26 kept 1PI graphs, 6 sigma and 5 + 5 pi+-
self-energies; the D4 terms removed exactly).  Each graph value is an upper
bound per graph on the routing ambiguity of the size rule, not on the loop.
Where the largest routing depends on the masses the coefficients are
Piecewise expressions.  Everything stays exact in N(m_a^2), M^{-1}_{m_a^2}(C),
eps and the source traces; there is no expansion in eps.  R is a size ratio,
not a computed correction.  Where M_0 = 0 the external-leg terms vanish
identically and R is undefined; only M_1^est is reported.  No loop integral,
loop function or loop amplitude is evaluated.

Dimensional check: X -> c^{d_X} X with d = 1 for B_0, f_pi, the masses
(m_sigma, m_a, m_pi, m_P) and every source letter in a trace, d = 0 for eps,
N, lambda_{pi+-} and M^{-1}_{m_a^2}(C).  Every max(...) of masses is rescaled
with the exact identity max(c x, c y, ...) = c max(x, y, ...) for c > 0, and
every condition of a Piecewise coefficient is checked to be homogeneous (so
that it is unchanged by the rescaling).
"""
import sympy as sp

import sigma_pipi_loop_algebra as algebra
from sigma_pipi_loop_algebra import (B0, F_PI, M_A, M_MESON, M_PI, M_SIGMA, TraceSum)
import derive_sigma_pipi_loop_graphs as graphs
import derive_sigma_pipi_loop_traces as traces
import derive_sigma_pipi_loop_tree as tree_script


def assemble(entries=None):
    """Returns a dict with M_0, the 1PI contributions X_Gamma^max, dZ^max per
    leg (sum over its self-energies), the leg term and M_1^est (words)."""
    entries = traces.kept_graph_traces() if entries is None else entries
    tree = tree_script.tree_amplitude()
    tree_compact = tree_script.tree_amplitude_compact()
    contributions = {}
    residues = {'sigma': TraceSum(), 'pi+': TraceSum(), 'pi-': TraceSum()}
    for entry in entries:
        if entry['group'] == 'one_particle_irreducible':
            contributions[entry['label']] = entry['words']
        else:
            residues[entry['group']] = residues[entry['group']] + entry['residue']
    residues = {key: value.simplified() for key, value in residues.items()}
    leg_factor = (residues['sigma'] + residues['pi+'] + residues['pi-']).scale(sp.Rational(1, 2))
    leg_term = (tree * leg_factor).simplified()
    one_particle_irreducible = TraceSum()
    for value in contributions.values():
        one_particle_irreducible = one_particle_irreducible + value
    total = (one_particle_irreducible + leg_term).simplified()
    return {'tree': tree, 'tree_compact': tree_compact, 'contributions': contributions, 'residues': residues,
            'leg_term': leg_term, 'one_particle_irreducible': one_particle_irreducible.simplified(),
            'estimate': total, 'entries': entries}


def trace_symbol(word):
    return sp.Symbol(r'\langle ' + ''.join(algebra.letter_latex(name) for name in word) + r'\rangle')


def as_expression(trace_sum):
    """TraceSum -> sympy expression with one symbol per trace (for R)."""
    total = 0
    for key, value in trace_sum.terms.items():
        total += value * sp.Mul(*(trace_symbol(word) for word in key))
    return total


def ratio(result):
    """R = |M_1^est|/|M_0|, or None (undefined) when M_0 vanishes."""
    if result['tree'].is_zero():
        return None
    return sp.Abs(as_expression(result['estimate'])) / sp.Abs(as_expression(result['tree']))


SCALE = sp.Symbol('c', positive=True)
MASSES = (M_SIGMA, M_A, M_PI, M_MESON)


def normalized_relation(relation, scale):
    """A rescaled homogeneous condition with the common factor scale^k divided
    out; raises AssertionError if the condition is not homogeneous."""
    difference = sp.expand(relation.lhs - relation.rhs)
    for power in range(0, 5):
        reduced = sp.expand(difference / scale ** power)
        if not reduced.has(scale):
            return relation.func(sp.expand(relation.lhs / scale ** power),
                                 sp.expand(relation.rhs / scale ** power))
    raise AssertionError(f'condition not homogeneous: {relation}')


def rescale_value(value, scale=SCALE):
    """X -> c^{d_X} X for one coefficient.  Each max(...) of masses is first
    replaced by a mass symbol of dimension one, rescaled like the other masses
    (max(c x, c y, ...) = c max(x, y, ...) for c > 0) and then restored; the
    conditions of Piecewise coefficients are rescaled and normalized (they must
    be homogeneous)."""
    value = sp.sympify(value)
    maxima = sorted(value.atoms(sp.Max), key=sp.default_sort_key)
    placeholders = {atom: sp.Dummy('larger_mass', positive=True) for atom in maxima}
    value = value.subs(placeholders)
    substitution = {B0: scale * B0, F_PI: scale * F_PI}
    substitution.update({mass: scale * mass for mass in MASSES})
    substitution.update({placeholder: scale * placeholder for placeholder in placeholders.values()})
    value = value.subs(substitution, simultaneous=True)
    value = value.replace(lambda node: isinstance(node, sp.core.relational.Relational),
                          lambda node: normalized_relation(node, scale))
    return value.subs({placeholder: atom for atom, placeholder in placeholders.items()})


def rescaled(trace_sum, scale=SCALE):
    """Apply X -> c^{d_X} X to a TraceSum."""
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        result.add_term(key, rescale_value(value, scale) * scale ** algebra.source_count(key))
    return result


def scales_as(trace_sum, power):
    difference = rescaled(trace_sum) - trace_sum.scale(SCALE ** power)
    return all(graphs.expression_is_zero(value) for value in difference.terms.values())


def dimension_checks(result):
    checks = {'M_0 -> c M_0': scales_as(result['tree'], 1),
              'M_1^est -> c M_1^est': scales_as(result['estimate'], 1)}
    for label, value in result['contributions'].items():
        checks[f'1PI {label} -> c'] = scales_as(value, 1)
    for leg, value in result['residues'].items():
        checks[f'dZ_{leg} -> dZ'] = scales_as(value, 0)
    checks['R -> R (from M_0 -> c M_0 and M_1^est -> c M_1^est)'] = (checks['M_0 -> c M_0']
                                                                    and checks['M_1^est -> c M_1^est'])
    return checks


def vanishing_tree_locus(trace_sum):
    """Impose <S_10 lambda_{pi-} lambda_{pi+}> = -<S_10 lambda_{pi+} lambda_{pi-}>,
    i.e. <S_10 {lambda_{pi+}, lambda_{pi-}}> = 0."""
    plus_minus = algebra.canonical_word(('S10', 'l+', 'l-'))
    minus_plus = algebra.canonical_word(('S10', 'l-', 'l+'))
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        factor = value
        words = []
        for word in key:
            if word == minus_plus:
                factor = -factor
                words.append(plus_minus)
            else:
                words.append(word)
        result.add_term(tuple(words), factor)
    return result.simplified()


def vanishing_tree_checks(result):
    tree = vanishing_tree_locus(result['tree'])
    legs = vanishing_tree_locus(result['leg_term'])
    restricted = dict(result, tree=tree, leg_term=legs,
                      estimate=vanishing_tree_locus(result['estimate']))
    return {'M_0 vanishes': tree.is_zero(), 'leg terms vanish': legs.is_zero(),
            'R': ratio(restricted), 'M_1^est nonzero': not restricted['estimate'].is_zero()}


def graph_table(result):
    """Per kept graph: group, label, S, r, the scales m_term of its terms, the
    routing decisions and the D4 flag."""
    rows = []
    for entry in result['entries']:
        if entry['group'] == 'pi-':
            continue
        rows.append({'group': entry['group'], 'label': entry['label'], 'symmetry': entry['symmetry'],
                     'r': len(entry['blocks']), 'scales': graphs.possible_scales(entry),
                     'routing': graphs.routing_summary(entry),
                     'excluded': entry['excluded'] is not None})
    return rows


def sources_in_estimate(result):
    return sorted({name for name in result['estimate'].letters() if algebra.is_source(name)})


def main():
    result = assemble()
    print('M_0 =', tree_script.TREE_LATEX)
    print('\nM_1^est = sum over kept 1PI graphs of X^max + M_0 (dZ_sigma^max + dZ_pi+^max + dZ_pi-^max)/2')
    for leg in ('sigma', 'pi+', 'pi-'):
        print(f'\ndZ_{leg}^max: {len(result["residues"][leg].terms)} products of traces of words')
    print('\ndZ_pi+ equals dZ_pi- under lambda_pi+ <-> lambda_pi-:', len(result['residues']['pi+'].terms)
          == len(result['residues']['pi-'].terms))
    algebra.assert_no_float(result['estimate'])
    print('Number of trace products in M_1^est (words):', len(result['estimate'].terms))
    print('Sources entering M_1^est:', ', '.join(sources_in_estimate(result)))
    print('\nTable: group | label | S | r | m_term of the terms | routing decisions | D4 removed')
    for row in graph_table(result):
        print(f"  {row['group']} | {row['label']} | {row['symmetry']} | {row['r']} | {row['scales']} | "
              f"{row['routing']} | {row['excluded']}")
    print('\nDimensional checks:')
    for name, value in dimension_checks(result).items():
        print(f'  {name}: {value}')
    print('\nVanishing tree <S_10 {lambda_pi+, lambda_pi-}> = 0:', vanishing_tree_checks(result))


if __name__ == '__main__':
    main()
