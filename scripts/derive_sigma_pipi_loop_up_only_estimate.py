"""Up-only limit of the large-epsilon estimate of the one-loop sigma -> pi+ pi- amplitude:
M_1^est/M_0 -> -21 B_0 m_u eps^4/(64 pi^2 f_pi^2) (task 4.5).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_up_only_estimate.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_up_only_estimate
(about one minute: the kept graphs and their routing choices are built once).

Units: GeV for B_0, m_u, f_pi and every source letter; eps = f_pi/f_a is
dimensionless.  Conventions: those of sigma_pipi_loop_algebra.py (polar fields,
connection-free frame, <lambda_i lambda_j> = 2 delta_ij, C = B_0 P_hat_01) and
the large-eps behaviour of derive_sigma_pipi_loop_large_epsilon.py (fixed
sources and fixed quartic coupling lambda; ALP mass entry eps^2 A with
A = -2 B_0 <S_02>; m_sigma^2 = 2 lambda f_pi^2/eps^2 + eps^2 kappa_sigma with
kappa_sigma = B_0 (<S_10> - 2 <S_20>); loop scale per term, m_sigma for a
radial line, m_a for an ALP-like line and for pieces (ii)+(iii)+(iv) of a
meson line, m_P for the bare generator sum, the maximum over the lines of the
term; termwise routing maximum, chosen at large eps).

Inputs (all from derive_sigma_pipi_loop_large_epsilon.py and the graph
scripts): the kept graphs with the size-rule value of every routing, the
routing choice in the regions kappa_sigma < A and A < kappa_sigma < 2 A, the
large-eps series of every graph value and of M_0, and the leading coefficient
of M_1^est for generic fixed sources.

Steps (each printed as exact sympy text and as LaTeX):
 (a) up-only sources S_hat_{r,s} + i P_hat_{r,s} = m_u (-i)^s P_u/s! for r <= 1,
     (r, s) != (0, 0), zero for r >= 2; P_u is a symbolic projector with
     P_u^2 = P_u, S_hat_{0,0} P_u = P_u S_hat_{0,0} = m_u P_u, <P_u> = 1, and
     lambda_{pi+} lambda_{pi-} = 2 P_u, lambda_{pi-} lambda_{pi+} = 2 P_d,
     P_d^2 = P_d, P_u P_d = P_d P_u = 0, <P_d> = 1 (no matrix components).
     Assumption (up-only model, pi+ = u dbar): lambda_{pi+} = P_u lambda_{pi+} P_d
     and lambda_{pi-} = P_d lambda_{pi-} P_u, i.e. P_u lambda_{pi+} = lambda_{pi+} P_d
     = lambda_{pi+}, lambda_{pi+} P_u = P_d lambda_{pi+} = 0 (and the same with
     u <-> d for lambda_{pi-}); lambda_{pi+-}^2 = 0 (charge +-2).  A trace that
     these relations do not reduce (for instance S_hat_{0,0} next to P_d) is
     reported, and only an upper bound on the power is given.
 (b) derived inputs: A, kappa_sigma, M_0, the ratios
     <X {lambda_{pi+}, lambda_{pi-}}>/<S_hat_{1,0} {lambda_{pi+}, lambda_{pi-}}>
     and the single traces that enter the leading terms.
 (c) generic leading eps^4 terms of M_1^est/M_0 graph by graph (1PI graph:
     X_Gamma^max/M_0; self-energy: (1/2) dZ^max), their up-only values, and
     the up-only leading power of every other kept graph.
 (d) the sum and R = |M_1^est/M_0| on both sides of kappa_sigma = A, with the
     order at which max(m_sigma, m_a) is decided for the up-only sources.
 (e) agreement with the up-only limit of the full large-eps total.
No loop integral, loop function or loop amplitude is evaluated.
"""
import sympy as sp

import sigma_pipi_loop_algebra as algebra
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI, TraceSum
import derive_sigma_pipi_loop_large_epsilon as large
import derive_sigma_pipi_loop_tree as tree_script
import derive_sigma_pipi_loop_vertices as vertices

UP_QUARK_MASS = large.UP_QUARK_MASS
ALP_MASS_COEFFICIENT = large.ALP_MASS_COEFFICIENT
KAPPA_SIGMA = large.KAPPA_SIGMA
RATIO_POINTS = {'kappa_sigma < A': sp.Rational(1, 2), 'A < kappa_sigma < 2 A': sp.Rational(3, 2)}
ANTICOMMUTATOR = '{+-}'
UNIT_OF_RATIO = B0 * UP_QUARK_MASS / (sp.pi ** 2 * F_PI ** 2)


# ---------------------------------------------------------------------------
# (a) Up-only sources
# ---------------------------------------------------------------------------

def up_only_letter(name):
    """(coefficient, 'Pu') of a dressed source letter for the up-only sources, or
    None if the letter vanishes; S_hat_{0,0} stays a letter (returned as
    (1, 'S00'))."""
    kind, radial, alp = name[0], int(name[1]), int(name[2:])
    if (radial, alp) == (0, 0):
        return (1, 'S00') if kind == 'S' else None
    if radial >= 2:
        return None
    value = UP_QUARK_MASS * (-sp.I) ** alp / sp.factorial(alp)
    coefficient = sp.re(value) if kind == 'S' else sp.im(value)
    return None if coefficient == 0 else (coefficient, 'Pu')


def source_letters(entries):
    """Dressed source letters that appear in the size-rule values of the kept graphs."""
    letters = set()
    for entry in entries:
        for value in entry['routing_values']:
            letters |= {name for name in value.letters() if algebra.is_source(name)}
    return sorted(letters)


UP_ONLY_RULES = {**tree_script.UP_ONLY_RULES,
                 ('S00', 'Pu'): [(UP_QUARK_MASS, ('Pu',))], ('Pu', 'S00'): [(UP_QUARK_MASS, ('Pu',))],
                 ('Pu', 'l+'): [(1, ('l+',))], ('l+', 'Pd'): [(1, ('l+',))], ('l+', 'Pu'): [], ('Pd', 'l+'): [],
                 ('Pd', 'l-'): [(1, ('l-',))], ('l-', 'Pu'): [(1, ('l-',))], ('l-', 'Pd'): [], ('Pu', 'l-'): [],
                 ('l+', 'l+'): [], ('l-', 'l-'): []}


class UnreducedTrace(Exception):
    """A trace that the up-only relations of step (a) do not reduce to a number."""


def up_only_trace(word):
    """Number <word> for a word in P_u, P_d, S_hat_{0,0}, lambda_{pi+}, lambda_{pi-}
    with the relations of step (a)."""
    if not word:
        return algebra.NUMBER_OF_FLAVOURS
    total = 0
    for reduced, factor in tree_script._rewrite({tuple(word): 1}, UP_ONLY_RULES).items():
        rotated = tree_script._rewrite({algebra.canonical_word(reduced): 1}, UP_ONLY_RULES)
        for final, final_factor in rotated.items():
            if final not in tree_script.UP_ONLY_SINGLE_TRACES:
                raise UnreducedTrace(final)
            total += factor * final_factor * tree_script.UP_ONLY_SINGLE_TRACES[final]
    return sp.expand(total)


def drop_vanishing_letters(trace_sum):
    """Remove the products of traces that contain a source letter that vanishes
    for the up-only sources (exact)."""
    return trace_sum.filtered(lambda key: all(up_only_letter(name) is not None for word in key
                                              for name in word if algebra.is_source(name)))


def up_only_value(trace_sum):
    """Value of a TraceSum (words in source letters and lambda_{pi+-}; coefficients
    in A, kappa_sigma, B_0, f_pi, ...) for the up-only sources, with
    A = kappa_sigma = B_0 m_u."""
    total = 0
    for key, coefficient in vertices.up_only_substitution(trace_sum).terms.items():
        product = coefficient
        for word in key:
            product *= up_only_trace(word)
        total += product
    return sp.expand(sp.sympify(total).subs({ALP_MASS_COEFFICIENT: B0 * UP_QUARK_MASS,
                                             KAPPA_SIGMA: B0 * UP_QUARK_MASS}))


# ---------------------------------------------------------------------------
# Exact expressions with one symbol per (anticommutator) trace
# ---------------------------------------------------------------------------

def plain_word(word):
    return ' '.join('{l+,l-}' if name == ANTICOMMUTATOR else name for name in word)


def latex_word(word):
    return large.word_latex(word)


def trace_expression(trace_sum):
    """(sympy expression, LaTeX names) with one symbol <W> per trace; a pair of
    traces with lambda_{pi+} lambda_{pi-} and lambda_{pi-} lambda_{pi+} and equal
    coefficients is one symbol <W {l+,l-}>."""
    total, names = 0, {}
    for value, key, _ in large.anticommutator_terms(trace_sum):
        product = value
        for word in key:
            symbol = sp.Symbol('<' + plain_word(word) + '>')
            names[symbol] = r'\langle ' + latex_word(word) + r'\rangle'
            product *= symbol
        total += product
    return total, names


def anticommutator_trace_sum(word):
    """TraceSum of <W {lambda_{pi+}, lambda_{pi-}}> for a word W."""
    word = tuple(word)
    return TraceSum.single(word + ('l+', 'l-')) + TraceSum.single(word + ('l-', 'l+'))


def latex(expression, names=None):
    return sp.latex(expression, symbol_names={**large.LATEX_NAMES_LARGE_EPSILON, **(names or {})})


# ---------------------------------------------------------------------------
# Graphs, regions and totals
# ---------------------------------------------------------------------------

def region_at(entries, ratio):
    """Routing choices and branch of max(m_sigma, m_a) at r = kappa_sigma/A = ratio
    (the same choices as in the region of derive_sigma_pipi_loop_large_epsilon.py
    that contains `ratio`)."""
    return {'test': ratio, 'branch': large.branch_of(ratio),
            'choices': {entry['name']: large.routing_choices(entry['routing_values'], ratio, entry['name'])
                        for entry in entries}}


def contribution_to_ratio(coefficient, tree_coefficient):
    """Leading coefficient of the contribution of one graph to M_1^est/M_0 (an
    expression in trace symbols) from its leading coefficient in M_1^est."""
    numerator, names = trace_expression(coefficient)
    denominator, tree_names = trace_expression(tree_coefficient)
    return sp.expand(sp.cancel(numerator / denominator)), {**names, **tree_names}


def graph_contributions(total, tree_coefficient):
    """{graph: {'generic', 'names', 'up-only'}} for the graphs that attain the
    leading power of the generic total; 'up-only' is the up-only value of the
    contribution to M_1^est/M_0."""
    tree_value = up_only_value(tree_coefficient)
    result = {}
    for name, coefficient in total['by graph'].items():
        generic, names = contribution_to_ratio(coefficient, tree_coefficient)
        result[name] = {'generic': generic, 'names': names,
                        'up-only': sp.factor(up_only_value(coefficient) / tree_value)}
    return result


def up_only_leading_power(entry, region):
    """Up-only leading power of the contribution of one kept graph to
    M_1^est/M_0 (1PI: power of X - 1; self-energy: power of dZ), the order of
    the series used and a note; None for an exact zero (the graph value
    vanishes once every S_hat_{r>=2,s}, P_hat_{r>=2,s} is set to zero)."""
    value = large.chosen_value(entry['routing_values'], region['choices'][entry['name']][0])
    restricted = drop_vanishing_letters(value)
    restricted = drop_vanishing_letters(large.substitute_mass_traces(restricted))
    if restricted.is_zero():
        letters = sorted({name for key in value.terms for word in key for name in word
                          if algebra.is_source(name) and up_only_letter(name) is None})
        listed = letters[0] if len(letters) == 1 else ', '.join(letters[:-1]) + ' or ' + letters[-1]
        return {'power': None, 'note': f'zero (exact: every product of traces contains {listed})'}
    shift = -1 if entry['group'] == 'one_particle_irreducible' else 0
    precision = None
    for order in large.ORDERS:
        series = drop_vanishing_letters(large.context(order).series(restricted, region['branch']))
        precision = series.precision
        exponents = sorted({exponent for inner in series.terms.values() for exponent in inner})
        for exponent in exponents:
            try:
                value = up_only_value(series.coefficient(exponent))
            except UnreducedTrace as error:
                return {'power': None, 'bound': -exponent + shift,
                        'note': f'at most eps^{-exponent + shift} (trace {error} not reduced)'}
            if value != 0:
                return {'power': -exponent + shift, 'coefficient': value, 'order': order,
                        'note': f'eps^{-exponent + shift}'}
    return {'power': None, 'bound': -precision - 1 + shift,
            'note': f'zero up to eps^{-precision + shift} (series), below eps^{-precision + shift}'}


def all_results():
    """Every quantity printed by main()."""
    entries, _ = large.graph_entries()
    tree, tree_lead = large.tree_leading()
    tree_coefficient = tree_lead['coefficient']
    letters = source_letters(entries)
    regions = {label: region_at(entries, ratio) for label, ratio in RATIO_POINTS.items()}
    totals = {label: {sources: large.total_estimate(entries, tree, region, sources)
                      for sources in large.SOURCE_SETS}
              for label, region in regions.items()}
    contributions = {label: graph_contributions(totals[label]['generic'], tree_coefficient) for label in regions}
    powers = {label: {entry['name']: up_only_leading_power(entry, region) for entry in entries}
              for label, region in regions.items()}
    pole = large.pole_data_table()
    return {'entries': entries, 'tree': tree, 'tree coefficient': tree_coefficient, 'letters': letters,
            'regions': regions, 'totals': totals, 'contributions': contributions, 'powers': powers,
            'mass difference': pole['m_a^2 - m_sigma^2 (up-only)'],
            'mass difference generic': pole['m_a^2 - m_sigma^2']}


# ---------------------------------------------------------------------------
# Derived inputs (b)
# ---------------------------------------------------------------------------

def derived_inputs(results):
    """A, kappa_sigma, M_0, anticommutator ratios and single traces of step (b)."""
    s02, s10, s20 = TraceSum.single(('S02',)), TraceSum.single(('S10',)), TraceSum.single(('S20',))
    tree_value = up_only_value(results['tree'])
    reference = up_only_value(anticommutator_trace_sum(('S10',)))
    words, singles = set(), set()
    for label in RATIO_POINTS:
        for coefficient in results['totals'][label]['generic']['by graph'].values():
            for _, key, _ in large.anticommutator_terms(coefficient):
                for word in key:
                    if word and word[-1] == ANTICOMMUTATOR:
                        words.add(word[:-1])
                    elif 'l+' not in word and 'l-' not in word:
                        singles.add(word)
    ratios = {word: sp.factor(up_only_value(anticommutator_trace_sum(word)) / reference)
              for word in sorted(words) if word != ('S10',)}
    traces = {word: up_only_value(TraceSum.single(word)) for word in sorted(singles)}
    return {'A = -2 B_0 <S02>': up_only_value(s02.scale(-2 * B0)),
            'kappa_sigma = B_0 (<S10> - 2 <S20>)': up_only_value((s10 - s20.scale(2)).scale(B0)),
            'M_0': sp.factor(tree_value),
            '<S10 {l+,l-}>': reference, 'ratios': ratios, 'single traces': traces}


def total_value(contributions):
    return sp.factor(sum(value['up-only'] for value in contributions.values()))


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def print_item(text, expression, names=None):
    print(f'  {text}: {sp.sstr(expression)}')
    print(f'    LaTeX: {latex(expression, names)}')


BRANCH_TEXT = {large.BRANCHES[0]: 'max(m_sigma, m_a) = m_a', large.BRANCHES[1]: 'max(m_sigma, m_a) = m_sigma'}


def main():
    results = all_results()
    first, second = RATIO_POINTS
    print('(a) Up-only sources: S_hat_{r,s} + i P_hat_{r,s} = m_u (-i)^s P_u/s! for r <= 1, (r,s) != (0,0); '
          'zero for r >= 2.')
    print('    Relations: P_u^2 = P_u, S_hat_{0,0} P_u = P_u S_hat_{0,0} = m_u P_u, <P_u> = 1; '
          'lambda_{pi+} lambda_{pi-} = 2 P_u, lambda_{pi-} lambda_{pi+} = 2 P_d, P_d^2 = P_d, '
          'P_u P_d = P_d P_u = 0, <P_d> = 1; lambda_{pi+} = P_u lambda_{pi+} P_d, '
          'lambda_{pi-} = P_d lambda_{pi-} P_u, lambda_{pi+}^2 = lambda_{pi-}^2 = 0.')
    print('    Source letters of the kept graphs:', ', '.join(results['letters']))
    vanishing = []
    for name in results['letters']:
        image = up_only_letter(name)
        if image is None:
            vanishing.append(name)
        elif image[1] == 'S00':
            print(f'  {name}: stays a letter, with S00 P_u = P_u S00 = m_u P_u')
            print(r'    LaTeX: \widehat S_{0,0}P_u=P_u\widehat S_{0,0}=m_uP_u')
        else:
            print(f'  {name} = {sp.sstr(image[0])}*P_u')
            print(f'    LaTeX: {algebra.letter_latex(name)}={latex(image[0])}\\,P_u')
    print(f'  vanishing letters: {", ".join(vanishing)} = 0 (and P00 = P10 = 0 by the assumptions of the model)')
    print('    LaTeX: ' + ',\\ '.join(algebra.letter_latex(name) for name in vanishing) + '=0')

    inputs = derived_inputs(results)
    print('\n(b) Derived inputs (up-only):')
    print_item('A = -2 B_0 <S02>', inputs['A = -2 B_0 <S02>'])
    print_item('kappa_sigma = B_0 (<S10> - 2 <S20>)', inputs['kappa_sigma = B_0 (<S10> - 2 <S20>)'])
    print_item('M_0 = -(B_0 eps/(2 f_pi)) <S10 {l+,l-}>', inputs['M_0'])
    print_item('<S10 {l+,l-}>', inputs['<S10 {l+,l-}>'])
    for word, value in inputs['ratios'].items():
        print_item(f'<{plain_word(word)} {{l+,l-}}>/<S10 {{l+,l-}}>', value)
    for word, value in inputs['single traces'].items():
        print_item(f'<{plain_word(word)}>', value)

    print('\n(c) Contributions to M_1^est/M_0 at eps^4 (1PI graph: X^max/M_0; self-energy: (1/2) dZ^max).')
    total = results['totals'][first]['generic']
    print(f'  Region {first} ({BRANCH_TEXT[results["regions"][first]["branch"]]}); generic fixed sources: '
          f'M_1^est/M_0 ~ eps^{total["R power"]}, attained by {", ".join(total["attaining"])}')
    contributions = results['contributions'][first]
    for name, value in contributions.items():
        print_item(f'graph {name}, generic coefficient of eps^4', value['generic'], value['names'])
        print_item(f'graph {name}, up-only coefficient of eps^4', value['up-only'])
    survivors = [name for name, value in contributions.items() if value['up-only'] != 0]
    print('  graphs with a nonzero up-only eps^4 term:', ', '.join(survivors))
    print('  up-only leading power of the contribution to M_1^est/M_0 of every other kept graph:')
    for name, power in results['powers'][first].items():
        if name not in survivors:
            print(f'    {name}: {power["note"]}')
    total = results['totals'][second]['generic']
    print(f'  Region {second} ({BRANCH_TEXT[results["regions"][second]["branch"]]}); generic fixed sources: '
          f'M_1^est/M_0 ~ eps^{total["R power"]}, attained by {", ".join(total["attaining"])}')
    for name, value in results['contributions'][second].items():
        if sp.expand(value['generic'] - contributions[name]['generic']) != 0:
            print_item(f'graph {name}, generic coefficient of eps^4 (differs from region {first})',
                       value['generic'], value['names'])
            print_item(f'graph {name}, up-only coefficient of eps^4', value['up-only'])
    other = [name for name, value in results['contributions'][second].items() if value['up-only'] != 0]
    print('  graphs with a nonzero up-only eps^4 term:', ', '.join(other))
    same = all(results['powers'][second][name]['note'] == power['note']
               for name, power in results['powers'][first].items())
    print(f'  up-only leading powers of every other kept graph equal to those of region {first}: {same}')

    print('\n(d) Sum and R for the up-only sources:')
    for label, contributions in results['contributions'].items():
        value = total_value(contributions)
        print_item(f'region {label}: M_1^est/M_0 -> eps^4 times', value)
        print_item(f'region {label}: R = |M_1^est/M_0| -> eps^4 times (B_0 m_u > 0)', sp.Abs(value))
    mass = results['mass difference']
    print_item('m_a^2 - m_sigma^2 for generic sources, coefficient of eps^2',
               large.as_trace_expression(results['mass difference generic']['coefficient']))
    print_item(f'm_a^2 - m_sigma^2 for the up-only sources, with u = 1/eps, up to O(u^{mass["precision"] + 1})',
               mass['up-only series'])
    print(f'    The up-only sources give A = kappa_sigma (the boundary between the two regions): the eps^2 terms '
          f'cancel and max(m_sigma, m_a) is decided at order eps^{mass["power"]}, where '
          f'm_a^2 - m_sigma^2 -> {sp.sstr(mass["coefficient"].terms[()])} > 0, so max(m_sigma, m_a) = m_a.')

    print('\n(e) Up-only limit of the full large-eps total (sources r <= 1):')
    for label in RATIO_POINTS:
        value = large.up_only_ratio(results['totals'][label]['r <= 1'])
        print_item(f'region {label}: M_1^est/M_0 -> eps^4 times', value)
        print(f'    equal to the sum of (d): {sp.simplify(value - total_value(results["contributions"][label])) == 0}')


if __name__ == '__main__':
    main()
