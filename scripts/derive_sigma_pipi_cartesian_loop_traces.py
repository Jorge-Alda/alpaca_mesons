"""Cartesian frame: compact traces of the kept one-loop graphs of sigma -> pi+ pi-, with per-term loop scales and the termwise routing maximum (task 3.2).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_traces.py
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_traces.py --markdown
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_traces

OpenSpec change estimate-sigma-pipi-one-loop-murayama-cartesian.  Units GeV;
conventions of sigma_pipi_cartesian_loop_algebra.py.  The graphs, line pieces,
per-term scales (D3), D4 exclusion and size rule are those of
derive_sigma_pipi_cartesian_loop_graphs.py; everything is derived from the
Cartesian vertices of stage A.

For every kept graph and every routing the value is built in two ways:
  construction 1: vertices as words, completeness on words;
  construction 2: vertices as compact trees (sums over permutations,
                  commutators, anticommutators), completeness applied to the
                  trees before any expansion (D5);
and the two are compared as exact identities of cyclic words.

Routing (D12; assumption): for each routing rho the size-rule value is
sum_t C_t^rho t over the independent products of traces of words t (the
canonical generic-source basis, lambda already eliminated through m_sigma^2);
the graph enters with sum_t C_t^{rho*_t} t, rho*_t the routing of largest
|C_t^rho|, sign kept.  Ties of equal magnitude and opposite sign take the sign
of the sum over routings; if that sum is zero, the first routing in the
canonical order is used and the product is marked.  Mass-dependent choices are
Piecewise expressions in homogeneous inequalities between squared masses; the
comparison uses m_sigma >= 2 m_pi (open decay), so m_sigma^2 - m_pi^2 > 0.
The termwise maximum is computed with the graph-agnostic helper
termwise_maximum of derive_sigma_pipi_loop_graphs.py (it acts on a list of
sums of products of traces and contains no graph, vertex or letter).

No loop integral, scalar loop function, loop coefficient or loop amplitude is
computed.  The result is an order of magnitude per graph class and can miss
cancellations between graphs with the same epsilon scaling.
"""
import functools
import itertools
import sys

import sympy as sp

import sigma_pipi_cartesian_loop_algebra as cartesian
import derive_sigma_pipi_cartesian_loop_graphs as graphs
import derive_sigma_pipi_cartesian_loop_tree as cartesian_tree
from derive_sigma_pipi_loop_graphs import termwise_maximum    # graph-agnostic helper
from sigma_pipi_loop_algebra import M_PI, M_SIGMA, TraceSum, CompactSum

GROUP_NAMES = {'one_particle_irreducible': '1PI', 'sigma': 'sigma self-energy', 'pi': 'pi+- self-energy'}


# ---------------------------------------------------------------------------
# Graph values per routing and termwise maximum
# ---------------------------------------------------------------------------

def compact_by_scale(row, line, orientation):
    """{m_term: CompactSum} of one routing (construction 2)."""
    return graphs.routing_value(row, row['group'], line, orientation, row['selection'], construction='trees',
                                by_scale=True)


def routing_segments(maximum):
    """For the mass-dependent products: {scale: (conditions, routing indices)} if all products
    with that scale share the same branches, else None."""
    structure = {}
    for key, (kind, info) in maximum['choices'].items():
        if kind != 'mass-dependent':
            continue
        segments = tuple((segment['lower'], segment['upper'], segment['lower closed'], segment['upper closed'],
                          segment['index']) for segment in info['segments'])
        entry = (info['gamma'], segments)
        if structure.setdefault(info['scale'], entry) != entry:
            return None
    result = {}
    for scale, (gamma, segments) in structure.items():
        combination = M_SIGMA ** 2 + gamma * M_PI ** 2
        conditions = []
        for lower, upper, lower_closed, upper_closed, _ in segments:
            parts = []
            if lower != 0:
                parts.append((sp.GreaterThan if lower_closed else sp.StrictGreaterThan)(combination, lower * scale ** 2))
            if upper != sp.oo:
                parts.append((sp.LessThan if upper_closed else sp.StrictLessThan)(combination, upper * scale ** 2))
            conditions.append(sp.And(*parts))
        result[scale] = (tuple(conditions), tuple(segment[4] for segment in segments))
    return result


def displayed_form(entry):
    """Compact form of the termwise maximum: the compact value of a base routing plus, for each
    scale whose products depend on the routing, the indicator of each branch times the compact
    difference between the routing of that branch and the base routing.
    Returns {'base': index, 'base value': CompactSum, 'corrections': [(scale, condition, CompactSum)]}."""
    segments = entry['segments']
    if not segments:
        base = entry['common routing']
        return {'base': base, 'base value': entry['compact'][base], 'corrections': [],
                'base by scale': entry['compact by scale'][base]}
    first = {indices[0] for _, indices in segments.values()}
    base = first.pop() if len(first) == 1 else 0
    base_by_scale = entry['compact by scale'][base]
    corrections = []
    for scale, (conditions, indices) in segments.items():
        for condition, index in zip(conditions, indices):
            if index == base:
                continue
            other = entry['compact by scale'][index]
            difference = (other.get(scale, CompactSum()) + base_by_scale.get(scale, CompactSum()).scale(-1))
            corrections.append((scale, condition, difference.simplified()))
    base_value = CompactSum()
    for value in base_by_scale.values():
        base_value = base_value + value
    return {'base': base, 'base value': base_value.simplified(), 'corrections': corrections,
            'base by scale': base_by_scale}


def displayed_equals_maximum(entry):
    """The displayed form reproduces the termwise maximum exactly, branch by branch: for every
    combination of branches (one per scale), the base routing plus the corrections of the
    selected branches equals the maximum with every Piecewise replaced by that branch."""
    display = entry['display']
    segments = entry['segments'] or {}
    scales = list(segments)
    for combination in itertools.product(*(range(len(segments[scale][0])) for scale in scales)):
        selected = dict(zip(scales, combination))
        value = display['base value'].to_trace_sum()
        for scale, condition, difference in display['corrections']:
            if segments[scale][0][selected[scale]] == condition:
                value = value + difference.to_trace_sum()
        target = TraceSum()
        for key, coefficient in entry['maximum']['value'].terms.items():
            target.add_term(key, select_branch(coefficient, selected))
        if not (value - target).is_zero():
            return False
    return True


def branch_scale(condition, scales):
    """The loop scale whose branches a Piecewise condition describes (the largest atom present)."""
    present = [scale for scale in scales if sp.sympify(condition).has(scale)]
    outer = [scale for scale in present if not any(other != scale and other.has(scale) for other in present)]
    assert len(outer) == 1, condition
    return outer[0]


def select_branch(expression, selected):
    """Replace every Piecewise by its branch number selected[scale] (sympy may rewrite the
    conditions, e.g. the last one into True, but keeps the order of the branches)."""
    def choose(*arguments):
        scale = branch_scale(arguments[0].cond, list(selected))
        return arguments[selected[scale]].expr
    return sp.sympify(expression).replace(lambda node: isinstance(node, sp.Piecewise),
                                          lambda node: choose(*node.args))


def spread(entry):
    """Routing spread (diagnostic): for every product, the range of C_t^rho/C_t^{rho*_t} over the
    routings and over each branch (limits at the ends of the branch).  Returns (minimum, maximum,
    number of products with a zero routing coefficient somewhere)."""
    lowest, highest, zero = sp.Integer(1), sp.Integer(1), 0
    w = sp.Symbol('w', positive=True)
    for key, (kind, info) in entry['maximum']['choices'].items():
        if kind == 'independent':
            continue
        if kind == 'numerical':
            ratios = entry['maximum']['spread'][key]
            lowest, highest = min(lowest, *ratios), max(highest, *ratios)
            zero += any(ratio == 0 for ratio in ratios)
            continue
        functions = [offset + slope * w for offset, slope in zip(info['offsets'], info['slopes'])]
        found_zero = False
        for segment in info['segments']:
            chosen = functions[segment['index']]
            for function in functions:
                ratio = sp.cancel(function / chosen)
                for end in (segment['lower'], segment['upper']):
                    value = sp.limit(ratio, w, end)
                    lowest, highest = min(lowest, value), max(highest, value)
                if any(0 < root < sp.oo and segment['lower'] <= root <= segment['upper']
                       for root in sp.solve(function, w)) or function == 0:
                    found_zero = True
        zero += found_zero
    return lowest, highest, zero


@functools.lru_cache(maxsize=None)
def _entries():
    entries = []
    for group in graphs.GROUPS:
        for row in graphs.kept(graphs.table(group)):
            words = row['routing values']
            compact_scales = [compact_by_scale(row, line, orientation) for line, orientation in graphs.routings(row)]
            compact = []
            for by_scale in compact_scales:
                total = CompactSum()
                for value in by_scale.values():
                    total = total + value
                compact.append(total.simplified())
            maximum = termwise_maximum(words)
            common = next((index for index, value in enumerate(words) if (value - maximum['value']).is_zero()), None)
            entry = {**row, 'words': words, 'compact': compact, 'compact by scale': compact_scales,
                     'maximum': maximum, 'common routing': common,
                     'constructions agree': all((one.to_trace_sum() - two).is_zero()
                                                for one, two in zip(compact, words)),
                     'segments': routing_segments(maximum)}
            entry['display'] = displayed_form(entry)
            entry['spread'] = spread(entry)
            entries.append(entry)
    return tuple(entries)


def entries():
    """Kept graphs with their values per routing (words and compact), the termwise maximum
    ('maximum'['value'], a TraceSum with Piecewise coefficients where the routing choice depends
    on the masses), the displayed compact form and the routing spread.  Cached."""
    return list(_entries())


def entry_by_label(label):
    return next(entry for entry in entries() if entry['label'] == label)


# ---------------------------------------------------------------------------
# Inputs for stage C
# ---------------------------------------------------------------------------

def estimate_inputs():
    """{'one particle irreducible': {label: X_Gamma^max}, 'delta Z sigma': sum of dZ^max,
    'delta Z pi': sum of dZ^max (equal for pi+ and pi-), 'tree': M_0} as TraceSums of words."""
    result = {'one particle irreducible': {}, 'delta Z sigma': TraceSum(), 'delta Z pi': TraceSum(),
              'tree': cartesian_tree.tree_amplitude()}
    for entry in entries():
        value = entry['maximum']['value']
        if entry['group'] == 'one_particle_irreducible':
            result['one particle irreducible'][entry['label']] = value
        elif entry['group'] == 'sigma':
            result['delta Z sigma'] = result['delta Z sigma'] + value
        else:
            result['delta Z pi'] = result['delta Z pi'] + value
    return result


def pi_minus_residue_corrections():
    """dZ of the pi- leg (in = pi-, out = pi+) per kept pi graph, in routing 1."""
    values = {}
    for row in graphs.kept(graphs.table('pi')):
        total = TraceSum()
        legs = graphs.setup('pi-')
        for choice, scale, value in graphs.term_values(row, legs, 0, 1):
            if row['selection'] == 'without D4' and graphs.is_all_generator_sum(choice):
                continue
            total = total + graphs.residue_correction(row, value, scale, M_PI)
        values[row['label']] = total.simplified()
    return values


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def label_latex(label):
    if label.startswith('Csigma'):
        return r'C\sigma' + label[6:]
    if label.startswith('Cpi'):
        return r'C\pi' + label[3:]
    return label


def scale_latex(scale):
    return cartesian.coefficient_latex(scale)


def display_coefficient(value):
    return sp.factor(sp.expand(value))


def compact_lines(compact_sum):
    terms = CompactSum()
    terms.terms = {key: display_coefficient(value) for key, value in compact_sum.terms.items()}
    return cartesian.compact_latex(terms, separator=' \\\\\n&+')


def condition_latex(condition):
    return cartesian.coefficient_latex(condition).replace(r'\wedge', r',\ ')


def entry_markdown(entry):
    name = r'X_{' + label_latex(entry['label']) + r'}^{\max}' if entry['group'] == 'one_particle_irreducible' \
        else r'\delta Z_{' + label_latex(entry['label']) + r'}^{\max}'
    display = entry['display']
    parts = [compact_lines(display['base by scale'][scale])
             for scale in sorted(display['base by scale'], key=sp.default_sort_key)]
    text = '$$\n\\begin{aligned}\n' + name + '={}&' + ' \\\\\n&+'.join(parts)
    for scale, condition, difference in display['corrections']:
        text += (' \\\\\n&+\\theta\\big(' + condition_latex(condition) + '\\big)\\Big[' +
                 compact_lines(difference) + '\\Big]')
    return text + '\n\\end{aligned}\n$$'


def no_component_notation(text):
    """Output guard: no generator letters, matrix indices or unit-matrix basis E_ij."""
    forbidden = ('E_', '_{ij}', 'lambda_i', 'lambda_j', r'\lambda_{i', r'\lambda_{j', 'G1', 'G2', 'G3', "G'")
    return not any(token in text for token in forbidden)


def entry_heading(entry):
    """One line per graph: lines, S_Gamma, r_Gamma, m_term, D4 and routing summary."""
    lines = ', '.join({'S': '$s$', 'A': 'ALP-like', 'M': 'meson'}[line] for line in entry['lines'])
    scales = ','.join(scale_latex(scale) for scale in entry['scales'])
    text = (f"**${label_latex(entry['label'])}$.** Lines {lines}; $S_\\Gamma={sp.latex(entry['symmetry'])}$, "
            f"$r_\\Gamma={len(entry['blocks'])}$, $m_{{\\rm term}}\\in\\{{{scales}\\}}$")
    if entry['D4 status'] == 'excluded (D4)':
        text += '; D4 term excluded'
    counts = entry['maximum']['counts']
    if entry['common routing'] is not None:
        text += '; routing independent'
    else:
        text += (f"; {counts['mass-dependent']} of {sum(counts.values())} products depend on the routing "
                 f"(mass-dependent choice, routing {entry['display']['base'] + 1} below the boundary)")
    return text + '.'


def markdown_blocks(group):
    """Result blocks (heading and compact termwise maximum) of the kept graphs of one group."""
    blocks = []
    for entry in entries():
        if entry['group'] != group:
            continue
        text = entry_heading(entry) + '\n' + entry_markdown(entry)
        assert no_component_notation(text)
        blocks.append(text)
    return '\n\n'.join(blocks)


def main():
    if '--markdown' in sys.argv:
        for group in graphs.GROUPS:
            print(markdown_blocks(group) + '\n')
        return
    all_entries = entries()
    print('Construction 1 = construction 2 in every routing of every kept graph:',
          all(entry['constructions agree'] for entry in all_entries))
    for entry in all_entries:
        counts = entry['maximum']['counts']
        summary = (f"{entry['label']} ({GROUP_NAMES[entry['group']]}): products {sum(counts.values())}, "
                   f"routing-dependent {counts['numerical'] + counts['mass-dependent']} "
                   f"(mass-dependent {counts['mass-dependent']}, marked at a boundary "
                   f"{len(entry['maximum']['marked boundary'])}, marked on an open region "
                   f"{len(entry['maximum']['marked'])}); common routing: "
                   f"{None if entry['common routing'] is None else entry['common routing'] + 1}; "
                   f"spread {entry['spread']}")
        print(summary)
        if entry['segments']:
            for scale, (conditions, indices) in entry['segments'].items():
                print(f'    scale {scale}: ' + '; '.join(f'routing {index + 1} for {condition}'
                                                         for condition, index in zip(conditions, indices)))
    inputs = estimate_inputs()
    print('dZ_pi+ = dZ_pi- per graph (routing 1):',
          all((pi_minus_residue_corrections()[entry['label']] - entry['words'][0]).is_zero()
              for entry in all_entries if entry['group'] == 'pi'))
    print('Sum of dZ_pi^max:', cartesian.trace_sum_latex(inputs['delta Z pi']))


if __name__ == '__main__':
    main()
