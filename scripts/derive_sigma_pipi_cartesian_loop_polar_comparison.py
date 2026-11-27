"""Diagnostic comparison of the Cartesian and the polar size estimates of the one-loop sigma -> pi+ pi- amplitude (task 4.4).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_polar_comparison.py
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_polar_comparison.py --markdown
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_polar_comparison
(about ten minutes: both estimates are rebuilt from their scripts).

OpenSpec change estimate-sigma-pipi-one-loop-murayama-cartesian.  Units GeV;
eps = f_pi/f_a.  Conventions of sigma_pipi_cartesian_loop_algebra.py (Cartesian
letters S^C_{r,s}, P^C_{r,s}) and sigma_pipi_loop_algebra.py (polar letters
S_{r,s}, P_{r,s}).

This is a diagnostic made after both derivations were complete.  Both results
are order-of-magnitude estimates that can miss cancellations between graphs
with the same epsilon scaling, and differently in each parametrization.
Agreement is not expected and neither estimate is authoritative.  Only the
exact one-loop total, with tadpoles and residues included, is independent of
the parametrization.

Every number is read from the existing accessors; no polar loop data are
recomputed by new means:
  Cartesian: derive_sigma_pipi_cartesian_loop_graphs (graphs, pieces, D4 terms),
    _traces (kept graphs), _large_epsilon_graphs (regions, routing choices,
    leading terms, cancellation statement), _large_epsilon_vertices (m_a^2 root),
    _up_only_estimate, _limits (vacuum condition), _tree (A, kappa_sigma,
    m_sigma^2), _tree_equivalence (map between the coefficients).
  Polar: derive_sigma_pipi_loop_large_epsilon (all_results, survival_table,
    chosen_value, leading_trace_sum, substitute_mass_traces, alp_inverse_propagator,
    alp_pole_series) and derive_sigma_pipi_loop_up_only_estimate (all_results).

Method.
 (a) Map between the coefficients: Taylor expansion of x = sigma/f_a, y = a/f_a
     along 1 + (s + i p)/f_a = (1 + sigma/f_a) e^{i a/f_a}, through third order.
     A, kappa_sigma, C and m_sigma^2 of the Cartesian frame, written in polar
     coefficients and rewritten with the polar definitions, equal the polar ones;
     the m_a^2 series agree.  Up-only sources: the polar series S + i P =
     m_u (-i)^s P_u/s! (r <= 1) maps to S^C_10 = m_u P_u, P^C_01 = -m_u P_u and
     zero otherwise (orders 1..UP_ONLY_ORDER).  D4 terms: Cartesian C3, C15, Cpi3
     against polar 3, 15, pi2, routing by routing, after the map.
 (b) Classes: (topology, multiset of internal line types) with a meson line
     typed by the pieces that contribute to it, '(i)' (generator sum),
     '(ii)-(iv)' (ALP-sector pieces) or both.  For each class: number of graphs
     and the leading term of the sum of its members (1PI: X/M_0 as numerator over
     the tree trace; self-energies: dZ; for the polar pi legs (dZ_pi+ + dZ_pi-)/2).
 (c) Up-only: Cartesian and polar leading terms, read from the up-only scripts.
 (d) Within each class, the graphs with the largest individual power are given
     the statement of cancellation_status of the Cartesian large-eps script
     (applied in the same way to both parametrizations).
 (e) Total and per topology, per region of r = kappa_sigma/A; differences are
     reported, not attributed.

No loop integral, scalar loop function, loop coefficient or loop amplitude is
evaluated.  No numerical matrix is used.
"""
import functools
import sys

import sympy as sp

import derive_sigma_pipi_cartesian_loop_graphs as cartesian_graphs
import derive_sigma_pipi_cartesian_loop_traces as cartesian_traces
import derive_sigma_pipi_cartesian_loop_large_epsilon_graphs as cartesian_large
import derive_sigma_pipi_cartesian_loop_large_epsilon_vertices as cartesian_vertices
import derive_sigma_pipi_cartesian_loop_up_only_estimate as cartesian_up_only
import derive_sigma_pipi_cartesian_loop_limits as cartesian_limits
import derive_sigma_pipi_cartesian_loop_tree as cartesian_tree
import derive_sigma_pipi_cartesian_loop_tree_equivalence as tree_equivalence
import derive_sigma_pipi_loop_large_epsilon as polar_large
import derive_sigma_pipi_loop_up_only_estimate as polar_up_only
import sigma_pipi_loop_algebra as polar_algebra
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI, QUARTIC, TraceSum
from sigma_pipi_cartesian_loop_algebra import ALP_MASS_COEFFICIENT, KAPPA_SIGMA, UP_QUARK_MASS

RELATION_ORDER = 3            # total order of the coefficient map (the leading numerators need r + s <= 3)
UP_ONLY_ORDER = 4             # orders checked in the up-only correspondence
POLAR_TREE_TRACE = TraceSum.single(('S10', 'l+', 'l-')) + TraceSum.single(('S10', 'l-', 'l+'))
CARTESIAN_TREE_LATEX = cartesian_large.TREE_TRACE_LATEX
POLAR_TREE_LATEX = '\\langle\\widehat S_{1,0}\\{\\lambda_{\\pi^+},\\lambda_{\\pi^-}\\}\\rangle'
D4_PAIRS = (('C3', '3'), ('C15', '15'), ('Cpi3', 'pi2'))
REGION_PAIRS = (('r < 1', (0,)), ('r > 1', (1, 2, 3)))       # Cartesian region, polar region indices
TOPOLOGIES = ('tadpole', 'bubble', 'triangle', 's/sigma self-energy', 'pi self-energy')
TOPOLOGY_LATEX = {'tadpole': '1PI tadpole', 'bubble': '1PI bubble', 'triangle': '1PI triangle',
                  's/sigma self-energy': '$s$/$\\sigma$ self-energy', 'pi self-energy': '$\\pi^\\pm$ self-energy'}
LINE_LATEX = {'radial': 'radial', 'ALP-like': 'ALP-like', 'meson (i)': 'meson (i)',
              'meson (ii)-(iv)': 'meson (ii)–(iv)', 'meson (i), (ii)-(iv)': 'meson (i) and (ii)–(iv)'}


# ---------------------------------------------------------------------------
# (a) Coefficient map, parameters, regions, up-only sources, D4 terms
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def coefficient_relation(order=RELATION_ORDER):
    """{(R, Q): S^C_{R,Q} as a linear combination of polar symbols S_{r,q}}, R + Q <= order.  The script checks
    that P^C_{R,Q} is the same combination of the P_{r,q}."""
    result = {}
    for key, (scalar, pseudoscalar) in tree_equivalence.cartesian_coefficients_from_polar(order).items():
        renamed = scalar.subs({symbol: tree_equivalence.polar_symbol('P', *map(int, symbol.name[1:]))
                               for symbol in scalar.free_symbols})
        assert sp.expand(renamed - pseudoscalar) == 0, key
        result[key] = scalar
    return result


@functools.lru_cache(maxsize=None)
def letter_map():
    return tree_equivalence.letter_mapping(order=RELATION_ORDER)


def to_polar_letters(trace_sum):
    """A Cartesian TraceSum written in polar letters through the coefficient map."""
    return tree_equivalence.substitute_letters(trace_sum, letter_map())


def polar_rewriting(trace_sum):
    """The polar exact rewriting <S_02> = -A/(2 B_0), <S_10> = kappa_sigma/B_0 + 2 <S_20>."""
    return polar_large.substitute_mass_traces(trace_sum)


def polar_alp_mass_coefficient():
    """Polar A as a TraceSum, read from the polar D_a(s) = s - eps^2 A + 2 eps^2 <C M_s^{-1}(C)>."""
    total, _ = polar_large.alp_inverse_propagator()
    return TraceSum({(('S02',),): -sp.sympify(total.terms[(('S02',),)]) / EPSILON ** 2})


@functools.lru_cache(maxsize=None)
def parameter_checks():
    """Exact checks that A, kappa_sigma, C, the operator L and m_sigma^2 are the same in both frames, and that the
    large-eps m_a^2 series agree (Cartesian through u^4, polar through t = u^2)."""
    mapping = letter_map()
    cartesian_alp = to_polar_letters(cartesian_tree.alp_mass_coefficient())
    cartesian_kappa = to_polar_letters(cartesian_tree.kappa_sigma())
    radial = to_polar_letters(cartesian_tree.radial_mass_squared())          # f_a = f_pi/eps in the tree script
    u = cartesian_vertices.INVERSE_EPSILON
    cartesian_root = sp.expand(cartesian_vertices.alp_pole_root().subs(
        dict(zip(cartesian_vertices.TRACE_C_POWER, polar_large.TRACE_C_POWER))))
    polar_root = sp.expand(polar_large.alp_pole_series(order=2)['root'].subs(polar_large.INVERSE_EPSILON_SQUARED,
                                                                            u ** 2))
    difference = sp.expand(cartesian_root - polar_root)
    root_agrees = all(term.as_coeff_exponent(u)[1] >= 4 for term in sp.Add.make_args(difference)) or difference == 0
    return {
        'A: Cartesian B_0(<S^C_10> - 2<S^C_02>) equals polar -2 B_0 <S_02>':
            (cartesian_alp - polar_alp_mass_coefficient()).is_zero(),
        'A: equals the parameter A after the polar rewriting':
            (polar_rewriting(cartesian_alp) - TraceSum.scalar(ALP_MASS_COEFFICIENT)).is_zero(),
        'kappa_sigma: Cartesian B_0(<S^C_10> - 2<S^C_20>) equals the parameter kappa_sigma after the polar rewriting':
            (polar_rewriting(cartesian_kappa) - TraceSum.scalar(KAPPA_SIGMA)).is_zero(),
        'C and L: S^C_00 = S_00 and P^C_01 = P_01':
            mapping['SC00'] == [(1, 'S00')] and mapping['PC01'] == [(1, 'P01')],
        'm_sigma^2 = 2 lambda f_pi^2/eps^2 + eps^2 kappa_sigma in both':
            (polar_rewriting(radial) - TraceSum.scalar(2 * QUARTIC * F_PI ** 2 / EPSILON ** 2
                                                       + EPSILON ** 2 * KAPPA_SIGMA)).is_zero(),
        'm_a^2 series in A and c_k agree (through eps^-2)': root_agrees,
    }


@functools.lru_cache(maxsize=None)
def region_boundaries():
    """(Cartesian boundaries of the regions of r, Cartesian comparison roots, polar boundaries)."""
    merged, _ = cartesian_large.regions()
    cartesian = sorted({interval['lower'] for interval in merged} | {interval['upper'] for interval in merged},
                       key=lambda value: (value == sp.oo, value))
    polar = polar_results()['regions']
    polar_edges = sorted({region['lower'] for region in polar} | {region['upper'] for region in polar},
                         key=lambda value: (value == sp.oo, value))
    return cartesian, cartesian_large.comparison_roots(), polar_edges


def polar_up_only_source(radial, alp):
    """Polar up-only coefficient (S_{r,s}, P_{r,s}) in units of m_u P_u: S + i P = (-i)^s/s! for r <= 1,
    (r, s) != (0, 0); zero for r >= 2.  S_00 stays a letter."""
    if radial > 1:
        return 0, 0
    value = sp.expand((-sp.I) ** alp / sp.factorial(alp))
    return sp.re(value), sp.im(value)


@functools.lru_cache(maxsize=None)
def up_only_correspondence(order=UP_ONLY_ORDER):
    """{(R, Q): (S^C_{R,Q}, P^C_{R,Q}) in units of m_u P_u} for 1 <= R + Q <= order, from the polar up-only
    series; and whether S^C_00 is the polar letter S_00 alone."""
    relation = tree_equivalence.cartesian_coefficients_from_polar(order)
    substitution = {}
    for radial in range(order + 1):
        for alp in range(order + 1 - radial):
            if (radial, alp) == (0, 0):
                continue
            scalar, pseudoscalar = polar_up_only_source(radial, alp)
            substitution[tree_equivalence.polar_symbol('S', radial, alp)] = scalar
            substitution[tree_equivalence.polar_symbol('P', radial, alp)] = pseudoscalar
    substitution[tree_equivalence.polar_symbol('P', 0, 0)] = 0
    values = {key: (sp.expand(s.subs(substitution)), sp.expand(p.subs(substitution)))
              for key, (s, p) in relation.items() if key != (0, 0)}
    constant = relation[(0, 0)]
    constant_is_letter = (sp.expand(constant[0] - tree_equivalence.polar_symbol('S', 0, 0)) == 0
                          and sp.expand(constant[1] - tree_equivalence.polar_symbol('P', 0, 0)) == 0)
    return values, constant_is_letter


def up_only_correspondence_holds():
    values, constant_is_letter = up_only_correspondence()
    expected = {(1, 0): (1, 0), (0, 1): (0, -1)}
    return constant_is_letter and all(value == expected.get(key, (0, 0)) for key, value in values.items())


def canonical_list(values):
    return sorted((tuple(sorted((key, sp.srepr(sp.expand(coefficient)))
                                for key, coefficient in value.simplified().terms.items()))
                   for value in values))


@functools.lru_cache(maxsize=None)
def d4_comparison():
    """{(Cartesian label, polar name): {'same routing by routing', 'same set of routing values'}}: the
    Cartesian D4 terms in polar letters against the polar excluded terms."""
    polar = {entry['name']: entry for entry in polar_results()['entries']}
    result = {}
    for cartesian_label, polar_name in D4_PAIRS:
        row = cartesian_graphs.row_by_label(cartesian_label)
        cartesian_values = [to_polar_letters(value) for value in row['D4 term']]
        polar_values = polar[polar_name]['excluded_routings']
        same_order = len(cartesian_values) == len(polar_values) and all(
            (one - two).is_zero() for one, two in zip(cartesian_values, polar_values))
        result[(cartesian_label, polar_name)] = {
            'same routing by routing': same_order,
            'same set of routing values': canonical_list(cartesian_values) == canonical_list(polar_values)}
    return result


# ---------------------------------------------------------------------------
# Polar data (read once)
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def polar_results():
    """polar_large.all_results without the survival table (computed separately without powers)."""
    return polar_large.all_results(with_survival=False)


@functools.lru_cache(maxsize=None)
def polar_survival():
    return {row['name']: row for row in polar_large.survival_table(polar_results()['entries'], with_powers=False)}


@functools.lru_cache(maxsize=None)
def polar_up_only_results():
    return polar_up_only.all_results()


# ---------------------------------------------------------------------------
# (b) Classes of graphs
# ---------------------------------------------------------------------------

def meson_line_label(pieces):
    """'meson (i)', 'meson (ii)-(iv)' or 'meson (i), (ii)-(iv)' from the list of contributing pieces."""
    names = {piece.split()[0] for piece in pieces}
    has_generator_sum = '(i)' in names
    has_alp_sector = bool(names & {'(ii)', '(iii)', '(iv)'})
    if has_generator_sum and has_alp_sector:
        return 'meson (i), (ii)-(iv)'
    return 'meson (i)' if has_generator_sum else 'meson (ii)-(iv)'


def class_key(group, blocks, lines, pieces):
    """(topology, sorted line types); pieces: {line index: contributing pieces} for the meson lines."""
    if group == 'one_particle_irreducible':
        topology = {1: 'tadpole', 2: 'bubble', 3: 'triangle'}[len(blocks)]
    elif group == 'sigma':
        topology = 's/sigma self-energy'
    else:
        topology = 'pi self-energy'
    labels = []
    for index, line in enumerate(lines):
        if line == 'S':
            labels.append('radial')
        elif line == 'A':
            labels.append('ALP-like')
        else:
            labels.append(meson_line_label(pieces[index]))
    return topology, tuple(sorted(labels))


@functools.lru_cache(maxsize=None)
def cartesian_classes():
    """{class key: [Cartesian labels]} of the kept graphs."""
    result = {}
    for entry in cartesian_traces.entries():
        pieces = cartesian_graphs.contributing_pieces(entry)
        key = class_key(entry['group'], entry['blocks'], entry['lines'], pieces)
        result.setdefault(key, []).append(entry['label'])
    return result


def polar_group(entry):
    return 'pi' if entry['group'] in ('pi+', 'pi-') else entry['group']


@functools.lru_cache(maxsize=None)
def polar_classes():
    """{class key: [polar names]} of the kept graphs (pi+ leg; the pi- graphs have the same classes)."""
    survival = polar_survival()
    result = {}
    for entry in polar_results()['entries']:
        if entry['group'] == 'pi-':
            continue
        pieces = {}
        if entry['name'] in survival:
            pieces = {line['line']: line['survivors'] for line in survival[entry['name']]['lines']}
        key = class_key(polar_group(entry), entry['blocks'], entry['lines'], pieces)
        result.setdefault(key, []).append(entry['name'])
    return result


def class_order(key):
    topology, lines = key
    return TOPOLOGIES.index(topology), len(lines), lines


@functools.lru_cache(maxsize=None)
def all_classes():
    """(classes in both, Cartesian only, polar only), each sorted."""
    cartesian, polar = cartesian_classes(), polar_classes()
    both = sorted(set(cartesian) & set(polar), key=class_order)
    return both, sorted(set(cartesian) - set(polar), key=class_order), sorted(set(polar) - set(cartesian),
                                                                               key=class_order)


def polar_partner_names(name):
    """Polar entries that make up the share of a graph: the pi+ graph and its pi- partner."""
    return [name, 'pi-' + name[2:]] if name.startswith('pi') else [name]


# Leading terms of a sum of graphs ------------------------------------------------

def cartesian_sum_leading(labels, region_name):
    """(power of eps of the share of M_1^est/M_0, displayed TraceSum) of the sum of Cartesian graphs of one group:
    1PI: (power, numerator over the tree trace) of X/M_0; self-energies: (power, dZ)."""
    region, choices = cartesian_large.region_data(region_name)
    entries = [cartesian_traces.entry_by_label(label) for label in labels]
    groups = {entry['group'] for entry in entries}
    assert len(groups) == 1
    group = groups.pop()
    total = TraceSum()
    for entry in entries:
        total = total + cartesian_large.selected_value(entry, choices)
    total = total.simplified()
    if total.is_zero():
        return None
    power, value, _ = cartesian_large.leading_of(total, region)
    if group == 'one_particle_irreducible':
        return cartesian_large.contribution(group, power, value)
    return power, value


def polar_sum_leading(names, region_index):
    """Same for polar graphs: 1PI: (power - 1, X coefficient times -2 f_pi/B_0), the numerator over
    <S_10 {lambda_pi+, lambda_pi-}>; sigma: (power, dZ); pi: (power, (dZ_pi+ + dZ_pi-)/2)."""
    results = polar_results()
    region = results['regions'][region_index]
    entries = {entry['name']: entry for entry in results['entries']}
    selected = [entries[partner] for name in names for partner in polar_partner_names(name)]
    groups = {polar_group(entry) for entry in selected}
    assert len(groups) == 1
    group = groups.pop()
    total = TraceSum()
    for entry in selected:
        total = total + polar_large.chosen_value(entry['routing_values'], region['choices'][entry['name']][0])
    if group == 'pi':
        total = total.scale(sp.Rational(1, 2))
    total = total.simplified()
    if total.is_zero():
        return None
    leading = polar_large.leading_trace_sum(total, region['branch'])
    if leading is None:
        return None
    if group == 'one_particle_irreducible':
        return leading['power'] - 1, leading['coefficient'].scale(-2 * F_PI / B0).simplified()
    return leading['power'], leading['coefficient']


@functools.lru_cache(maxsize=None)
def class_leading(key, region_name='r < 1', polar_region_index=0):
    """{'cartesian': (power, TraceSum) or None, 'polar': ...} of the sum of the graphs of a class."""
    cartesian = cartesian_classes().get(key, [])
    polar = polar_classes().get(key, [])
    return {'cartesian': cartesian_sum_leading(tuple(cartesian), region_name) if cartesian else None,
            'polar': polar_sum_leading(tuple(polar), polar_region_index) if polar else None}


# Individual graphs and the statement on cancellations ---------------------------------

def cartesian_graph_leading(label, region_name):
    """(power of the share, displayed TraceSum): 1PI numerator of X/M_0; self-energies dZ."""
    data = cartesian_large.graph_results(region_name)[label]
    if data['group'] == 'one_particle_irreducible':
        return data['contribution']
    return data['own']


def polar_graph_leading(name, region_index):
    """Same for a polar graph, from the polar graph table of the region (pi graphs: pi+ leg)."""
    rows = {(row['name'], row['sources']): row for row in polar_results()['graphs'][region_index]}
    row = rows[(name, 'generic')]
    if row['group'] == 'one_particle_irreducible':
        if row['contribution'] is None:
            return None
        return row['contribution']['power'] - 1, row['contribution']['coefficient'].scale(-2 * F_PI / B0).simplified()
    if row['dZ'] is None:
        return None
    return row['dZ']['power'], row['dZ']['coefficient']


def cancellation_statement(leading_terms, sample_points):
    """(statement, graphs at the largest individual power): 'single graph' or the statement of
    cancellation_status of the Cartesian large-eps script on the leading numerators of those graphs."""
    present = {name: value for name, value in leading_terms.items() if value is not None}
    if not present:
        return 'exactly zero', []
    top = max(value[0] for value in present.values())
    at_top = {name: value[1] for name, value in present.items() if value[0] == top}
    if len(at_top) == 1:
        return 'single graph', sorted(at_top)
    status, _ = cartesian_large.cancellation_status(at_top, sample_points)
    return status, sorted(at_top)


@functools.lru_cache(maxsize=None)
def class_statements(key):
    """Statements on cancellations within a class, r < 1 (Cartesian) and region 1 (polar)."""
    points = tuple(cartesian_large.region_sample_points('r < 1'))
    cartesian = {label: cartesian_graph_leading(label, 'r < 1') for label in cartesian_classes().get(key, [])}
    polar = {name: polar_graph_leading(name, 0) for name in polar_classes().get(key, [])}
    return {'cartesian': cancellation_statement(cartesian, points) if cartesian else None,
            'polar': cancellation_statement(polar, points) if polar else None}


def same_leading(first, second):
    if first is None or second is None:
        return first is None and second is None
    return first[0] == second[0] and (first[1] - second[1]).simplified().is_zero()


@functools.lru_cache(maxsize=None)
def region_changes():
    """Classes whose leading class-sum term differs between r < 1 and r > 1 (Cartesian), and between polar region 1
    and any of the polar regions 2-4; and whether every leading power is the same in all regions."""
    both, _, _ = all_classes()
    cartesian_changes, polar_changes, powers_same = [], [], True
    for key in both:
        below = class_leading(key)
        above = class_leading(key, 'r > 1', 1)
        if not same_leading(below['cartesian'], above['cartesian']):
            cartesian_changes.append(key)
        powers_same &= below['cartesian'][0] == above['cartesian'][0]
        for index in (1, 2, 3):
            other = class_leading(key, 'r > 1', index)['polar']
            powers_same &= other[0] == below['polar'][0]
            if not same_leading(below['polar'], other):
                polar_changes.append(key)
                break
    return cartesian_changes, polar_changes, powers_same


# ---------------------------------------------------------------------------
# (e) Topologies and total
# ---------------------------------------------------------------------------

def topology_members(topology):
    both, cartesian_only, polar_only = all_classes()
    cartesian = [label for key in both + cartesian_only if key[0] == topology for label in cartesian_classes()[key]]
    polar = [name for key in both + polar_only if key[0] == topology for name in polar_classes()[key]]
    return tuple(cartesian), tuple(polar)


@functools.lru_cache(maxsize=None)
def topology_leading(topology, region_name, polar_region_index):
    """Leading term of the sum over a topology, and the graphs whose individual power reaches it."""
    cartesian, polar = topology_members(topology)
    cartesian_value = cartesian_sum_leading(cartesian, region_name)
    polar_value = polar_sum_leading(polar, polar_region_index)
    cartesian_graphs_at = sorted(label for label in cartesian
                                 if (cartesian_graph_leading(label, region_name) or (None,))[0] == cartesian_value[0])
    polar_graphs_at = sorted(name for name in polar
                             if (polar_graph_leading(name, polar_region_index) or (None,))[0] == polar_value[0])
    return {'cartesian': cartesian_value, 'polar': polar_value, 'cartesian graphs': cartesian_graphs_at,
            'polar graphs': polar_graphs_at}


def polar_total(region_index):
    """(power of M_1^est/M_0, numerator over <S_10 {lambda_pi+, lambda_pi-}>, graphs) of the polar total."""
    total = polar_results()['total'][region_index]['generic']
    tree = polar_results()['tree']
    assert tree['power'] == 1
    assert (tree['coefficient'] - POLAR_TREE_TRACE.scale(-B0 / (2 * F_PI))).is_zero()
    return total['R power'], total['coefficient'].scale(-2 * F_PI / B0).simplified(), list(total['attaining'])


@functools.lru_cache(maxsize=None)
def total_comparison(region_name, polar_region_index):
    """Leading terms of M_1^est/M_0 in both, and the difference of the numerators in polar letters (Cartesian
    numerator mapped to polar letters minus polar numerator, both after the polar rewriting of <S_10>, <S_02>)."""
    power, numerator, labels = cartesian_large.total_leading(region_name)
    polar_power, polar_numerator, polar_names = polar_total(polar_region_index)
    difference = None
    if power == polar_power:
        difference = (polar_rewriting(to_polar_letters(numerator)) - polar_rewriting(polar_numerator)).simplified()
    return {'cartesian': (power, numerator, labels), 'polar': (polar_power, polar_numerator, polar_names),
            'difference': difference}


# ---------------------------------------------------------------------------
# (c) Up-only sources
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def up_only_totals():
    """{'cartesian': (power, coefficient, graphs), 'polar': (power, coefficient, graphs)}; the polar values are those
    of the region kappa_sigma < A of the polar up-only script (the other side gives the same value)."""
    cartesian = cartesian_up_only.total_leading('r = 1')
    results = polar_up_only_results()
    first, second = polar_up_only.RATIO_POINTS
    values = [polar_up_only.total_value(results['contributions'][label]) for label in (first, second)]
    assert sp.expand(values[0] - values[1]) == 0
    power = results['totals'][first]['generic']['R power']
    graphs = [name for name, value in results['contributions'][first].items() if value['up-only'] != 0]
    return {'cartesian': cartesian, 'polar': (power, values[0], graphs)}


def polar_up_only_share(name):
    """(power, coefficient) of the up-only share of M_1^est/M_0 of one polar graph, from the polar up-only
    script: 1PI X/M_0, sigma dZ/2, pi (dZ_pi+ + dZ_pi-)/2; None if exactly zero, ('bound', power) if only a bound."""
    results = polar_up_only_results()
    first = next(iter(polar_up_only.RATIO_POINTS))
    tree = polar_up_only.derived_inputs(results)['M_0'] / EPSILON
    terms = []
    for partner in polar_partner_names(name):
        data = results['powers'][first][partner]
        if data['power'] is None:
            if 'bound' in data:
                return 'bound', data['bound']
            continue
        coefficient = data['coefficient']
        if partner.startswith('sigma'):
            coefficient = coefficient / 2
        elif partner.startswith('pi'):
            coefficient = coefficient / 2
        else:
            coefficient = coefficient / tree
        terms.append((data['power'], sp.factor(coefficient)))
    if not terms:
        return None
    top = max(power for power, _ in terms)
    value = sp.factor(sum(coefficient for power, coefficient in terms if power == top))
    return (top, value) if value != 0 else ('cancels', top)


def cartesian_up_only_share(label):
    return cartesian_up_only.graph_leading('r = 1')[label]


@functools.lru_cache(maxsize=None)
def up_only_classes():
    """{class key: {'cartesian': [(label, share)], 'polar': [(name, share)]}} for the classes with a nonzero
    up-only share in either parametrization."""
    both, _, _ = all_classes()
    result = {}
    for key in both:
        cartesian = [(label, cartesian_up_only_share(label)) for label in cartesian_classes()[key]]
        polar = [(name, polar_up_only_share(name)) for name in polar_classes()[key]]
        if any(value is not None for _, value in cartesian + polar):
            result[key] = {'cartesian': cartesian, 'polar': polar}
    return result


# ---------------------------------------------------------------------------
# Further facts read from the scripts
# ---------------------------------------------------------------------------

def routing_dependent_graphs():
    """(Cartesian labels, polar names) of the kept graphs whose value depends on the routing."""
    cartesian = sorted(cartesian_large.routing_dependent_products())
    polar = sorted(entry['name'] for entry in polar_results()['entries']
                   if entry['group'] != 'pi-' and polar_large.routing_dependent_keys(entry['routing_values']))
    return cartesian, polar


def vacuum_condition():
    """(power, coefficient) of the one-loop relative change of M_0 between the two vacuum conditions (Cartesian
    limits script), and the up-only largest power of the Cartesian D4 terms relative to M_0."""
    return cartesian_limits.leading_relative_change(), cartesian_limits.up_only_d4_powers()


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def polar_letter_latex(name):
    return polar_algebra.letter_latex(name)


def combination_latex(expression):
    """LaTeX of a linear combination of polar symbols S_{r,q}."""
    parts = []
    for symbol in sorted(expression.free_symbols, key=lambda s: s.name):
        coefficient = expression.coeff(symbol)
        name = f'\\widehat S_{{{symbol.name[1]},{symbol.name[2]}}}'
        if coefficient == 1:
            text = name
        elif coefficient == -1:
            text = '-' + name
        else:
            text = sp.latex(coefficient) + name
        parts.append(text)
    return '+'.join(parts).replace('+-', '-') if parts else '0'


def markdown_relation():
    lines = []
    for (radial, alp), expression in sorted(coefficient_relation().items(), key=lambda item: (sum(item[0]), item[0])):
        lines.append(f'- $\\widehat S^C_{{{radial},{alp}}}={combination_latex(expression)}$')
    return '\n'.join(lines)


def markdown_parameters():
    return '\n'.join(f'- {name}: {value}' for name, value in parameter_checks().items())


def markdown_regions():
    cartesian, roots, polar = region_boundaries()
    return (f"- Cartesian: boundaries of the regions of $r$: {', '.join(f'${sp.latex(v)}$' for v in cartesian)}; "
            f"comparison roots that are not boundaries: {', '.join(f'${sp.latex(v)}$' for v in roots)}.\n"
            f"- Polar: boundaries of the regions of $r$: {', '.join(f'${sp.latex(v)}$' for v in polar)}.")


def markdown_up_only_sources():
    values, constant_is_letter = up_only_correspondence()
    nonzero = [f'(\\widehat S^C_{{{r},{q}}},\\widehat P^C_{{{r},{q}}})=({sp.latex(s)},{sp.latex(p)})\\,m_uP_u'
               for (r, q), (s, p) in sorted(values.items()) if (s, p) != (0, 0)]
    zero = len([1 for value in values.values() if value == (0, 0)])
    return (f'- From the polar series, for $1\\le R+Q\\le{UP_ONLY_ORDER}$: ${", ".join(nonzero)}$; the other {zero} '
            f'pairs vanish. $\\widehat S^C_{{0,0}}=\\widehat S_{{0,0}}$: {constant_is_letter}. '
            f'Correspondence holds: {up_only_correspondence_holds()}.')


def markdown_d4():
    lines = ['| Cartesian | polar | equal routing by routing | equal as sets of routing values |', '|---|---|---|---|']
    for (cartesian_label, polar_name), data in d4_comparison().items():
        lines.append(f'| ${cartesian_traces.label_latex(cartesian_label)}$ | {polar_name_text(polar_name)} | '
                     f"{yes_no(data['same routing by routing'])} | {yes_no(data['same set of routing values'])} |")
    return '\n'.join(lines)


def yes_no(value):
    return 'yes' if value else 'no'


def polar_name_text(name):
    if name.startswith('sigma'):
        return 'σ' + name[5:]
    if name.startswith('pi'):
        return 'π' + name[2:]
    return name


def class_text(key):
    topology, lines = key
    return TOPOLOGY_LATEX[topology] + ': ' + ', '.join(LINE_LATEX[line] for line in lines)


def power_cell(value):
    return '0 (exact)' if value is None else f'$\\epsilon^{{{value[0]}}}$'


def statement_cell(statement):
    return '–' if statement is None else statement[0]


def markdown_classes():
    both, cartesian_only, polar_only = all_classes()
    lines = ['| class | Cartesian graphs | polar graphs | number (C / polar) | leading power of the class share, '
             '$r<1$ (C / polar) | same-power graphs (C / polar) |', '|---|---|---|---|---|---|']
    for key in both:
        leading = class_leading(key)
        statements = class_statements(key)
        cartesian_labels = cartesian_classes()[key]
        polar_names = polar_classes()[key]
        lines.append(f'| {class_text(key)} | '
                     + ', '.join(f'${cartesian_traces.label_latex(label)}$' for label in cartesian_labels) + ' | '
                     + ', '.join(polar_name_text(name) for name in polar_names)
                     + f" | {len(cartesian_labels)} / {len(polar_names)} | {power_cell(leading['cartesian'])} / "
                     f"{power_cell(leading['polar'])} | {statement_cell(statements['cartesian'])} / "
                     f"{statement_cell(statements['polar'])} |")
    lines.append('')
    lines.append('Classes in the Cartesian parametrization only: '
                 + (', '.join(class_text(key) for key in cartesian_only) or 'none') + '. Classes in the polar '
                 'parametrization only: ' + (', '.join(class_text(key) for key in polar_only) or 'none') + '.')
    return '\n'.join(lines)


def cartesian_latex(group, value):
    if value is None:
        return '0 (exact)'
    power, trace_sum = value
    if group == 'one_particle_irreducible':
        return (f'$X/M_0\\to\\epsilon^{{{power}}}\\big[' + cartesian_large.numerator_latex(trace_sum)
                + '\\big]/' + CARTESIAN_TREE_LATEX + '$')
    return f'$\\delta Z\\to\\epsilon^{{{power}}}\\big[' + cartesian_large.numerator_latex(trace_sum) + '\\big]$'


def polar_latex(group, value):
    if value is None:
        return '0 (exact)'
    power, trace_sum = value
    body = polar_large.trace_sum_latex(trace_sum, limit=200)[1:-1]
    if group == 'one_particle_irreducible':
        return f'$X/M_0\\to\\epsilon^{{{power}}}\\big[' + body + '\\big]/' + POLAR_TREE_LATEX + '$'
    name = '(\\delta Z_{\\pi^+}+\\delta Z_{\\pi^-})/2' if group == 'pi' else '\\delta Z'
    return f'${name}\\to\\epsilon^{{{power}}}\\big[' + body + '\\big]$'


def group_of_topology(topology):
    if topology == 's/sigma self-energy':
        return 'sigma'
    if topology == 'pi self-energy':
        return 'pi'
    return 'one_particle_irreducible'


def markdown_class_coefficients():
    both, _, _ = all_classes()
    blocks = []
    for key in both:
        group = group_of_topology(key[0])
        leading = class_leading(key)
        blocks.append(f'- {class_text(key)}\n  - Cartesian: {cartesian_latex(group, leading["cartesian"])}\n'
                      f'  - polar: {polar_latex(group, leading["polar"])}')
    return '\n'.join(blocks)


def markdown_region_changes():
    cartesian, polar, powers_same = region_changes()
    return ('- Leading powers of every class the same in all regions (both parametrizations): '
            f'{yes_no(powers_same)}.\n- Classes whose leading coefficient changes between $r<1$ and $r>1$, '
            'Cartesian: ' + ('; '.join(class_text(key) for key in cartesian) or 'none') + '.\n'
            '- Classes whose leading coefficient changes between polar region 1 and one of the regions 2–4: '
            + ('; '.join(class_text(key) for key in polar) or 'none') + '.')


def markdown_topologies():
    lines = ['| topology | region (C / polar) | leading power (C / polar) | graphs at that power, C | '
             'graphs at that power, polar |', '|---|---|---|---|---|']
    polar_regions = polar_results()['regions']
    for topology in TOPOLOGIES:
        for region_name, indices in REGION_PAIRS:
            for index in indices:
                data = topology_leading(topology, region_name, index)
                lines.append(f"| {TOPOLOGY_LATEX[topology]} | ${region_name}$ / "
                             f"${polar_large.region_text(polar_regions[index])}$ | {power_cell(data['cartesian'])} / "
                             f"{power_cell(data['polar'])} | "
                             + ', '.join(f'${cartesian_traces.label_latex(label)}$' for label in data['cartesian graphs'])
                             + ' | ' + ', '.join(polar_name_text(name) for name in data['polar graphs']) + ' |')
    return '\n'.join(lines)


def markdown_topology_coefficients():
    blocks = []
    for topology in TOPOLOGIES:
        group = group_of_topology(topology)
        data = topology_leading(topology, 'r < 1', 0)
        blocks.append(f'- {TOPOLOGY_LATEX[topology]} ($r<1$ / polar region 1)\n'
                      f"  - Cartesian: {cartesian_latex(group, data['cartesian'])}\n"
                      f"  - polar: {polar_latex(group, data['polar'])}")
    return '\n'.join(blocks)


def markdown_totals():
    polar_regions = polar_results()['regions']
    lines = ['| region (C / polar) | leading power of $M_1^{\\rm est}/M_0$ and $R$ (C / polar) | graphs, C | '
             'graphs, polar | numerators equal after the map |', '|---|---|---|---|---|']
    for region_name, indices in REGION_PAIRS:
        for index in indices:
            data = total_comparison(region_name, index)
            equal = data['difference'] is not None and data['difference'].is_zero()
            lines.append(f"| ${region_name}$ / ${polar_large.region_text(polar_regions[index])}$ | "
                         f"$\\epsilon^{{{data['cartesian'][0]}}}$ / $\\epsilon^{{{data['polar'][0]}}}$ | "
                         + ', '.join(f'${cartesian_traces.label_latex(label)}$' for label in data['cartesian'][2])
                         + ' | ' + ', '.join(polar_name_text(name) for name in data['polar'][2])
                         + f' | {yes_no(equal)} |')
    return '\n'.join(lines)


def markdown_total_coefficients():
    data = total_comparison('r < 1', 0)
    power, numerator, _ = data['cartesian']
    polar_power, polar_numerator, _ = data['polar']
    difference = data['difference']
    text = ('Region $r<1$ (polar region 1). Cartesian:\n$$\n\\frac{M_1^{\\rm est}}{M_0}\\to\\frac{\\epsilon^{'
            f'{power}' + '}}{' + CARTESIAN_TREE_LATEX + '}\\Big[' + cartesian_large.numerator_latex(numerator)
            + '\\Big]\n$$\nPolar:\n$$\n\\frac{M_1^{\\rm est}}{M_0}\\to\\frac{\\epsilon^{' + f'{polar_power}' + '}}{'
            + POLAR_TREE_LATEX + '}\\Big[' + polar_large.trace_sum_latex(polar_numerator, limit=200)[1:-1]
            + '\\Big]\n$$')
    if difference is not None:
        text += ('\nCartesian numerator written in polar coefficients, minus the polar numerator (both with '
                 '$\\langle\\widehat S_{1,0}\\rangle=\\kappa_\\sigma/B_0+2\\langle\\widehat S_{2,0}\\rangle$, '
                 '$\\langle\\widehat S_{0,2}\\rangle=-A/(2B_0)$):\n$$\n'
                 + (polar_large.trace_sum_latex(difference, limit=200)[1:-1] if not difference.is_zero() else '0')
                 + '\n$$')
    return text


def share_text(value):
    if value is None:
        return '0 (exact)'
    if value[0] == 'bound':
        return f'at most $\\epsilon^{{{value[1]}}}$'
    if value[0] == 'cancels':
        return f'leading terms cancel at $\\epsilon^{{{value[1]}}}$'
    return f'$\\epsilon^{{{value[0]}}}\\,\\big({sp.latex(value[1], symbol_names=polar_large.LATEX_NAMES_LARGE_EPSILON)}\\big)$'


def markdown_up_only():
    totals = up_only_totals()
    names = polar_large.LATEX_NAMES_LARGE_EPSILON
    c_power, c_value, c_graphs = totals['cartesian']
    p_power, p_value, p_graphs = totals['polar']
    lines = [f"- Cartesian: $M_1^{{\\rm est}}/M_0\\to{sp.latex(c_value, symbol_names=names)}\\,\\epsilon^{{{c_power}}}$, "
             f"from " + ', '.join(f'${cartesian_traces.label_latex(label)}$' for label in c_graphs) + '.',
             f"- Polar: $M_1^{{\\rm est}}/M_0\\to{sp.latex(p_value, symbol_names=names)}\\,\\epsilon^{{{p_power}}}$, "
             f"from " + ', '.join(polar_name_text(name) for name in p_graphs) + '.',
             '',
             '| class | Cartesian graphs: up-only share of $M_1^{\\rm est}/M_0$ | polar graphs: up-only share of '
             '$M_1^{\\rm est}/M_0$ |', '|---|---|---|']
    for key, data in up_only_classes().items():
        lines.append(f'| {class_text(key)} | '
                     + '; '.join(f'${cartesian_traces.label_latex(label)}$: {share_text(value)}'
                                 for label, value in data['cartesian']) + ' | '
                     + '; '.join(f'{polar_name_text(name)}: {share_text(value)}' for name, value in data['polar'])
                     + ' |')
    return '\n'.join(lines)


def markdown_further():
    cartesian, polar = routing_dependent_graphs()
    (power, coefficient), d4_up_only = vacuum_condition()
    names = polar_large.LATEX_NAMES_LARGE_EPSILON
    return (f'- Graphs whose value depends on the routing: Cartesian '
            + ', '.join(f'${cartesian_traces.label_latex(label)}$' for label in cartesian)
            + '; polar ' + ', '.join(polar_name_text(name) for name in polar) + '.\n'
            f'- Relative change of $M_0$ between the two vacuum conditions at one loop: '
            f'$\\epsilon^{{{power}}}\\,\\big({sp.latex(coefficient, symbol_names=names)}\\big)$.\n'
            '- Largest power of the excluded D4 terms relative to $M_0$ for up-only sources (Cartesian): '
            + ', '.join(f'${cartesian_traces.label_latex(label)}$: '
                        + ('0 (exact)' if value is None else f'$\\epsilon^{{{value}}}$')
                        for label, value in d4_up_only.items()) + '.')


def markdown_blocks():
    """Every block of the record, in order."""
    return [markdown_relation(), markdown_parameters(), markdown_regions(), markdown_up_only_sources(),
            markdown_d4(), markdown_classes(), markdown_class_coefficients(), markdown_region_changes(),
            markdown_topologies(), markdown_topology_coefficients(), markdown_totals(),
            markdown_total_coefficients(), markdown_up_only(), markdown_further()]


def main():
    if '--markdown' in sys.argv:
        for block in markdown_blocks():
            print(block + '\n')
        return
    print('Map checked to third order:', tree_equivalence.map_check(RELATION_ORDER))
    for key, value in sorted(coefficient_relation().items()):
        print(f'  S^C_{key} =', value)
    print('Parameters:', parameter_checks())
    print('Regions (Cartesian boundaries, comparison roots, polar boundaries):', region_boundaries())
    print('Up-only correspondence:', up_only_correspondence_holds(), up_only_correspondence())
    print('D4 terms:', d4_comparison())
    both, cartesian_only, polar_only = all_classes()
    print('Classes in both:', len(both), '; Cartesian only:', cartesian_only, '; polar only:', polar_only)
    for key in both:
        leading = class_leading(key)
        print(f'  {key}: C {cartesian_classes()[key]} P {polar_classes()[key]}; powers '
              f"{leading['cartesian'] and leading['cartesian'][0]} / {leading['polar'] and leading['polar'][0]}; "
              f"statements {class_statements(key)['cartesian'][0]} / {class_statements(key)['polar'][0]}")
    print('Region changes:', region_changes())
    for region_name, indices in REGION_PAIRS:
        for index in indices:
            data = total_comparison(region_name, index)
            print(f'Total {region_name} / polar region {index + 1}: powers {data["cartesian"][0]} / '
                  f'{data["polar"][0]}; graphs {data["cartesian"][2]} / {data["polar"][2]}; '
                  f'difference terms: {None if data["difference"] is None else len(data["difference"].terms)}')
    print('Up-only totals:', up_only_totals())
    print('Routing-dependent graphs:', routing_dependent_graphs())
    print('Vacuum condition:', vacuum_condition())


if __name__ == '__main__':
    main()
