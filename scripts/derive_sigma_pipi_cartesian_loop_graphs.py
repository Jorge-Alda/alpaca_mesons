"""Cartesian frame: one-loop graphs of sigma -> pi+ pi- and of the external self-energies (task 3.1).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_graphs.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_graphs

OpenSpec change estimate-sigma-pipi-one-loop-murayama-cartesian.  Units GeV;
conventions of sigma_pipi_cartesian_loop_algebra.py (Cartesian fields
phi = f_a + s + i p, eps = f_pi/f_a, all momenta incoming, d_mu -> -i k_mu,
vertex = n-th derivative of the Lagrangian, M = V for a contact graph).
Everything is built from the Cartesian vertex rules of stage A; no polar
vertex, graph or trace and no field map is used.

Graphs.  A one-loop graph is a cycle of V vertices; vertex j carries a
nonempty set of external legs and line j joins vertex j to vertex j+1.  A
vertex with no external leg would be a two-point vertex: the whole quadratic
Lagrangian is in the propagators, so it does not exist.  Hence for the three
external legs s, pi+, pi- only V = 1 (tadpole on a five-leg vertex), V = 2
(bubble, vertices with 3 and 4 legs) and V = 3 (triangle) occur; a box (V = 4)
and a 1PI bubble with a two-leg and a five-leg vertex are absent (asserted).
S_Gamma = 1/|automorphisms| (rotations and reflections of the cycle that keep
the external legs and line types), r_Gamma = V.

Internal lines (from the Cartesian quadratic Lagrangian, record
cartesian-sigma-pipi-loop-frame-and-tree.md):
  S  radial line: the field s at both ends (no meson or p component, by
     P^C_10 = 0 and <S^C_11> = 0); scale m_sigma.
  A  ALP-like line: the pole vector v_a = (2 eps M^{-1}_{m_a^2}(C)/sqrt(N), 1/sqrt(N))
     on (mesons + p) at both ends; scale m_a.
  M  meson line: residue (1_9 + 1) - v_a v_a^T on (mesons + p), split exactly
     into four pieces:
       (i)   generator sum 1_9 (U(3) generator at both ends), scale m_P;
       (ii)  p entry 1 - 1/N (p at both ends), scale m_a;
       (iii) mixed entries -2 eps M^{-1}_{m_a^2}(C)/N (p at one end,
             M^{-1}_{m_a^2}(C) at the other, both orientations), scale m_a;
       (iv)  -4 eps^2 M^{-1}_{m_a^2}(C) x M^{-1}_{m_a^2}(C)/N, scale m_a.
     m_P is one symbolic mass, the largest pole mass of the nine meson-like
     states.

Size rule (D1, D2, D3; assumptions): a term is one choice of piece on every
meson line; its scale m_term is the largest scale of its lines; k.k and k.p_j
(p_j the incoming momentum of an external leg, used as basis vector) are set
to m_term^2; the term enters with S_Gamma (m_term^2)^{2-r_Gamma}/(16 pi^2).
External products are exact: p_sigma^2 = m_sigma^2, p_pi^2 = m_pi^2,
p_1.p_2 = (m_sigma^2 - 2 m_pi^2)/2, p_sigma.p_i = -m_sigma^2/2.  For a
self-energy (legs p and -p, p^2 symbolic) the residue correction
dZ = -Pi'(m^2) is used term by term:
  dZ_term = -S_Gamma [theta(r >= 2) (m_term^2)^{1-r} T + (m_term^2)^{2-r} dT/dp^2]/(16 pi^2)
at p^2 = m_ext^2.  No loop integral, scalar loop function, loop coefficient
or loop amplitude is computed.

Generator sums (D5) are removed by U(3) completeness,
sum_i <A l_i B l_i> = 2<A><B> and sum_i <A l_i><l_i B> = 2<AB>, then the
assumptions (P^C_00 = P^C_10 = 0, <S^C_01> = <S^C_11> = 0) and the charge rule
are applied.  D11 is not imposed as matrix identities.

D4: in a graph whose lines are all meson lines, the term in which every line
takes piece (i) is excluded if, in every routing, it contains no ALP-sector
letter (M^{-1}_{m_a^2}(C), N, P^C_01, S^C_{r,q>=1}, P^C_{r,q>=1}), neither
lambda nor m_phi^2 (checked with lambda kept), and its ratio to M_0 (its dZ
for a pi leg) depends on eps only through p_sigma^2 = m_sigma^2.
"""
import functools
import itertools
import sys
from dataclasses import dataclass

import sympy as sp

import sigma_pipi_cartesian_loop_algebra as cartesian
import sigma_pipi_loop_algebra as generic
from sigma_pipi_loop_algebra import (EPSILON, M_A, M_MESON, M_PI, M_SIGMA, MOMENTUM_SQUARED, NORMALIZATION,
                                     QUARTIC, TraceSum, CompactSum, make_dot)
from sigma_pipi_cartesian_loop_algebra import MINV_LETTER, RADIAL_MASS_PARAMETER

LINE_TYPES = ('S', 'A', 'M')
GROUPS = ('one_particle_irreducible', 'sigma', 'pi')
EXTERNAL_MASS = {'one_particle_irreducible': None, 'sigma': M_SIGMA, 'pi': M_PI}


# ---------------------------------------------------------------------------
# Line pieces
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LinePiece:
    """One piece of an internal line.

    name: description; scale: loop scale (D3); terms: tuple of
    (coefficient, end at the start of the line, end at the end of the line),
    with ends 'radial', 'generator', 'p', 'MinvC' or 'ALP-like'."""
    name: str
    scale: object
    terms: tuple


RADIAL_LINE = LinePiece('radial line', M_SIGMA, ((1, 'radial', 'radial'),))
ALP_LIKE_LINE = LinePiece('ALP-like line', M_A, ((1, 'ALP-like', 'ALP-like'),))
GENERATOR_SUM = LinePiece('(i) generator sum', M_MESON, ((1, 'generator', 'generator'),))
P_ENTRY = LinePiece('(ii) p entry', M_A, ((1 - 1 / NORMALIZATION, 'p', 'p'),))
MIXED_ENTRIES = LinePiece('(iii) mixed entries', M_A, ((-2 * EPSILON / NORMALIZATION, 'p', 'MinvC'),
                                                       (-2 * EPSILON / NORMALIZATION, 'MinvC', 'p')))
MESON_BLOCK = LinePiece('(iv) meson block', M_A, ((-4 * EPSILON ** 2 / NORMALIZATION, 'MinvC', 'MinvC'),))
MESON_LINE_PIECES = (GENERATOR_SUM, P_ENTRY, MIXED_ENTRIES, MESON_BLOCK)
ALP_LIKE_REMAINDER = (P_ENTRY, MIXED_ENTRIES, MESON_BLOCK)        # (ii)+(iii)+(iv) = (0 + 1) - v_a v_a^T
# Checks only: the bare p line, and the residue (1_9 + 1) - v_a v_a^T written with three pieces.
BARE_P_LINE = LinePiece('bare p line', M_A, ((1, 'p', 'p'),))
UNIT_P_ENTRY = LinePiece('p entry of 1_9 + 1', M_A, ((1, 'p', 'p'),))
MINUS_ALP_LIKE_POLE = LinePiece('minus v_a v_a^T', M_A, ((-1, 'ALP-like', 'ALP-like'),))
THREE_PIECE_RESIDUE = (GENERATOR_SUM, UNIT_P_ENTRY, MINUS_ALP_LIKE_POLE)


def pieces_of(line_type):
    """The pieces of a line type (meson line: the four-piece split)."""
    return {'S': (RADIAL_LINE,), 'A': (ALP_LIKE_LINE,), 'M': MESON_LINE_PIECES, 'F': (BARE_P_LINE,)}[line_type]


def end_leg(end, generator, momentum):
    """Vertex leg of one end of a line piece, with incoming momentum `momentum`."""
    if end == 'radial':
        return cartesian.radial_leg(momentum)
    if end == 'generator':
        return cartesian.meson_leg(generator, momentum)
    if end == 'p':
        return cartesian.alp_field_leg(momentum)
    if end == 'MinvC':
        return cartesian.meson_leg(MINV_LETTER, momentum)
    if end == 'ALP-like':
        return cartesian.alp_like_leg(momentum)
    raise ValueError(end)


# ---------------------------------------------------------------------------
# Enumeration
# ---------------------------------------------------------------------------

def cycle_images(blocks, lines):
    """The 2V images of a cycle under rotations and reflections.  Line j joins
    vertex j to vertex j+1; a reflection reverses the order of the vertices."""
    count = len(blocks)
    images = []
    for shift in range(count):
        images.append((tuple(blocks[(shift + j) % count] for j in range(count)),
                       tuple(lines[(shift + j) % count] for j in range(count))))
        images.append((tuple(blocks[(shift - j) % count] for j in range(count)),
                       tuple(lines[(shift - j - 1) % count] for j in range(count))))
    return images


def leg_assignments(external_names, vertex_count, allow_empty=False):
    """Ways to distribute the external legs over the cycle vertices (ordered)."""
    result = []
    for assignment in itertools.product(range(vertex_count), repeat=len(external_names)):
        if not allow_empty and set(assignment) != set(range(vertex_count)):
            continue
        result.append(tuple(tuple(sorted(name for name, where in zip(external_names, assignment) if where == v))
                            for v in range(vertex_count)))
    return result


def enumerate_graphs(external_names, max_vertices=4, allow_two_point_vertices=False):
    """Distinct one-loop graphs (cycles) with the given external legs.
    Returns (graphs, {V: number of ordered leg distributions})."""
    found = {}
    distributions = {}
    for vertex_count in range(1, max_vertices + 1):
        assignments = leg_assignments(external_names, vertex_count, allow_two_point_vertices)
        distributions[vertex_count] = len(assignments)
        for blocks in assignments:
            for lines in itertools.product(LINE_TYPES, repeat=vertex_count):
                images = cycle_images(blocks, lines)
                canonical = min(images)
                if canonical in found:
                    continue
                automorphisms = sum(1 for image in images if image == (blocks, lines))
                found[canonical] = {'blocks': canonical[0], 'lines': canonical[1],
                                    'symmetry': sp.Rational(1, automorphisms)}
    order = {1: 0, 2: 1, 3: 2}
    graphs = sorted(found.values(), key=lambda graph: (order.get(len(graph['blocks']), 9), graph['blocks'],
                                                       [LINE_TYPES.index(t) for t in graph['lines']]))
    return graphs, distributions


def topology(graph):
    return {1: 'tadpole', 2: 'bubble', 3: 'triangle'}.get(len(graph['blocks']), 'box or larger')


def vertex_leg_counts(graph):
    return tuple(len(block) + 2 for block in graph['blocks'])


# ---------------------------------------------------------------------------
# Kinematics and routing
# ---------------------------------------------------------------------------

def external_dot(first, second):
    """Exact external products, all momenta incoming (p_sigma + p_1 + p_2 = 0)."""
    pair = frozenset((first, second))
    table = {frozenset(('ps',)): M_SIGMA ** 2, frozenset(('p1',)): M_PI ** 2, frozenset(('p2',)): M_PI ** 2,
             frozenset(('p1', 'p2')): (M_SIGMA ** 2 - 2 * M_PI ** 2) / 2,
             frozenset(('ps', 'p1')): -M_SIGMA ** 2 / 2, frozenset(('ps', 'p2')): -M_SIGMA ** 2 / 2,
             frozenset(('p',)): MOMENTUM_SQUARED}
    return table[pair]


def setup(group):
    """External legs {name: leg} of a group: the 1PI three-point function or a
    self-energy with incoming p and -p (for 'pi': in = pi+, out = pi-)."""
    if group == 'one_particle_irreducible':
        return {'s': cartesian.radial_leg({'ps': 1}), 'pi+': cartesian.meson_leg('l+', {'p1': 1}),
                'pi-': cartesian.meson_leg('l-', {'p2': 1})}
    if group == 'sigma':
        return {'in': cartesian.radial_leg({'p': 1}), 'out': cartesian.radial_leg({'p': -1})}
    if group == 'pi':
        return {'in': cartesian.meson_leg('l+', {'p': 1}), 'out': cartesian.meson_leg('l-', {'p': -1})}
    if group == 'pi-':
        return {'in': cartesian.meson_leg('l-', {'p': 1}), 'out': cartesian.meson_leg('l+', {'p': -1})}
    raise ValueError(group)


def external_names(group):
    return tuple(setup(group))


def routings(graph):
    """Canonical order: routing 2j+1 (2j+2) has the loop momentum on line j along (against) the line."""
    return [(line, orientation) for line in range(len(graph['lines'])) for orientation in (1, -1)]


def add_momenta(*momenta):
    total = {}
    for momentum in momenta:
        for label, value in momentum.items():
            total[label] = total.get(label, 0) + value
    return {label: value for label, value in total.items() if value != 0}


def line_momenta(graph, legs, line, orientation):
    """Momentum of every line (flowing from vertex j to vertex j+1): the loop
    momentum k on `line`, the others fixed by conservation at each vertex."""
    count = len(graph['blocks'])
    momenta = [None] * count
    momenta[line] = {'k': orientation}
    for step in range(1, count):
        j = (line + step) % count
        momenta[j] = add_momenta(momenta[(j - 1) % count], *(legs[name].momentum for name in graph['blocks'][j]))
    return momenta


def loop_dot(scale):
    """k.k -> scale^2 and k.p_j -> scale^2 (D2); external products exact."""
    def basis(first, second):
        if 'k' in (first, second):
            return scale ** 2
        return external_dot(first, second)
    return make_dot(basis)


# ---------------------------------------------------------------------------
# Values of one term
# ---------------------------------------------------------------------------

def assumptions_and_charge(value):
    if isinstance(value, TraceSum):
        return cartesian.apply_charge_rule(cartesian.apply_assumptions(value))
    return cartesian.compact_assumptions(value)


def term_scale(choice):
    """m_term: largest scale among the lines of the term (D3)."""
    return sp.Max(*{piece.scale for piece in choice})


def choice_value(graph, legs, momenta, choice, construction='words', eliminate_lambda=True):
    """T of one term (one piece per line) in one routing, without the size factor.
    Generator sums removed by completeness, then assumptions and charge rule."""
    count = len(graph['blocks'])
    dot = loop_dot(term_scale(choice))
    generators = [f'G{j + 1}' for j in range(count)]
    empty = TraceSum if construction == 'words' else CompactSum
    total = empty()
    for terms in itertools.product(*(piece.terms for piece in choice)):
        product = empty.scalar(sp.Mul(*(term[0] for term in terms)))
        for v in range(count):
            incoming, outgoing = (v - 1) % count, v
            vertex_legs = [legs[name] for name in graph['blocks'][v]]
            vertex_legs.append(end_leg(terms[incoming][2], generators[incoming], momenta[incoming]))
            vertex_legs.append(end_leg(terms[outgoing][1], generators[outgoing],
                                       {label: -value for label, value in momenta[outgoing].items()}))
            if construction == 'words':
                factor = cartesian.apply_assumptions(cartesian.vertex(vertex_legs, dot, eliminate_lambda))
                factor = cartesian.apply_charge_rule(factor)
            else:
                factor = cartesian.compact_assumptions(
                    cartesian.compact_vertex(vertex_legs, dot, eliminate_lambda))
            product = product * factor
            if not product.terms:
                break
        for j, term in enumerate(terms):
            if 'generator' in (term[1], term[2]):
                product = (generic.completeness(product, generators[j]) if construction == 'words'
                           else generic.compact_completeness(product, generators[j]))
        total = total + product
    return assumptions_and_charge(total).simplified()


def term_values(graph, legs, line, orientation, piece_lists=None, construction='words', eliminate_lambda=True):
    """[(choice, scale, T)] for every choice of one piece per line in one routing (zero terms dropped)."""
    piece_lists = piece_lists or [pieces_of(line_type) for line_type in graph['lines']]
    momenta = line_momenta(graph, legs, line, orientation)
    result = []
    for choice in itertools.product(*piece_lists):
        value = choice_value(graph, legs, momenta, choice, construction, eliminate_lambda)
        if value.terms:
            result.append((choice, term_scale(choice), value))
    return result


def is_all_generator_sum(choice):
    return all(piece is GENERATOR_SUM for piece in choice)


def size_factor(graph, scale):
    """S_Gamma (m^2)^{2-r}/(16 pi^2)."""
    return graph['symmetry'] * (scale ** 2) ** (2 - len(graph['blocks'])) / (16 * sp.pi ** 2)


def residue_correction(graph, value, scale, external_mass):
    """dZ of one term: -S [theta(r>=2) (m^2)^{1-r} T + (m^2)^{2-r} dT/dp^2]/(16 pi^2) at p^2 = m_ext^2."""
    r = len(graph['blocks'])

    def coefficient(t):
        moment = (scale ** 2) ** (1 - r) * t if r >= 2 else 0
        total = -graph['symmetry'] * (moment + (scale ** 2) ** (2 - r) * sp.diff(t, MOMENTUM_SQUARED)) / (16 * sp.pi ** 2)
        return sp.expand(total.subs(MOMENTUM_SQUARED, external_mass ** 2))
    return value.map_coefficients(coefficient).simplified()


def sized(graph, group, scale, value):
    """Size-rule value of one term: S (m^2)^{2-r} T/(16 pi^2) (1PI) or dZ (self-energy)."""
    if group == 'one_particle_irreducible':
        return value.scale(size_factor(graph, scale)).simplified()
    return residue_correction(graph, value, scale, EXTERNAL_MASS[group])


def routing_value(graph, group, line, orientation, selection='all', piece_lists=None, construction='words',
                  eliminate_lambda=True, by_scale=False):
    """Size-rule value of one routing: sum over the terms of S (m_term^2)^{2-r} T_term/(16 pi^2)
    (or dZ_term).  selection: 'all', 'D4' (only the all-generator-sum term) or 'without D4'.
    With by_scale=True: {m_term: value} (pieces with the same scale summed)."""
    empty = TraceSum if construction == 'words' else CompactSum
    totals = {}
    for choice, scale, value in term_values(graph, setup(group), line, orientation, piece_lists, construction,
                                            eliminate_lambda):
        if selection == 'D4' and not is_all_generator_sum(choice):
            continue
        if selection == 'without D4' and is_all_generator_sum(choice):
            continue
        totals[scale] = totals.get(scale, empty()) + sized(graph, group, scale, value)
    totals = {scale: value.simplified() for scale, value in totals.items()}
    if by_scale:
        return {scale: value for scale, value in totals.items() if value.terms}
    total = empty()
    for value in totals.values():
        total = total + value
    return total.simplified()


def routing_values(graph, group, selection='all', **options):
    """Size-rule value in each routing, canonical order."""
    return [routing_value(graph, group, line, orientation, selection, **options)
            for line, orientation in routings(graph)]


# ---------------------------------------------------------------------------
# D4 criterion
# ---------------------------------------------------------------------------

def is_alp_sector_letter(name):
    """M^{-1}_{m_a^2}(C), P^C_01 and every S^C_{r,q}, P^C_{r,q} with q >= 1."""
    if name == MINV_LETTER:
        return True
    return cartesian.is_cartesian_source(name) and cartesian.source_indices(name)[2] >= 1


def satisfies_d4_criterion(value, group):
    """No ALP-sector letter, no N, no lambda, no m_phi^2, and the ratio to
    M_0 = eps x (eps-free trace) depends on eps only through p_sigma^2 = m_sigma^2
    (m_sigma is a symbol independent of eps here); for a self-energy dZ itself."""
    if any(is_alp_sector_letter(name) for name in value.letters()):
        return False
    divisor = EPSILON if group == 'one_particle_irreducible' else 1
    for coefficient in value.simplified().terms.values():
        coefficient = sp.sympify(coefficient)
        if coefficient.has(NORMALIZATION) or coefficient.has(QUARTIC) or coefficient.has(RADIAL_MASS_PARAMETER):
            return False
        if sp.diff(sp.cancel(coefficient / divisor), EPSILON) != 0:
            return False
    return True


# ---------------------------------------------------------------------------
# Reasons for identically vanishing graphs
# ---------------------------------------------------------------------------

def end_kind(end):
    """Short name of a line end: s, p, generator, M^{-1}_{m_a^2}(C) or the pole vector v_a."""
    return {'radial': 's', 'generator': 'generator', 'MinvC': 'MinvC', 'p': 'p', 'ALP-like': 'v_a'}[end]


def vanishing_reasons(graph, group):
    """Why each nonzero-looking term of a vanishing graph is zero: the field-end
    pattern of a vertex that vanishes after the assumptions and the charge rule,
    or 'cancels after completeness' if every vertex is nonzero."""
    legs = setup(group)
    reasons = set()
    momenta = line_momenta(graph, legs, 0, 1)
    count = len(graph['blocks'])
    pieces = [pieces_of(line_type) for line_type in graph['lines']]
    for choice in itertools.product(*pieces):
        dot = loop_dot(term_scale(choice))
        for terms in itertools.product(*(piece.terms for piece in choice)):
            zero_vertex = None
            for v in range(count):
                incoming, outgoing = (v - 1) % count, v
                vertex_legs = [legs[name] for name in graph['blocks'][v]]
                vertex_legs.append(end_leg(terms[incoming][2], f'G{incoming + 1}', momenta[incoming]))
                vertex_legs.append(end_leg(terms[outgoing][1], f'G{outgoing + 1}',
                                           {label: -value for label, value in momenta[outgoing].items()}))
                factor = cartesian.apply_charge_rule(cartesian.apply_assumptions(
                    cartesian.vertex(vertex_legs, dot))).simplified()
                if not factor.terms:
                    external = [{'s': 's', 'in': 'in', 'out': 'out'}.get(name, name) for name in graph['blocks'][v]]
                    zero_vertex = (tuple(external) + (end_kind(terms[incoming][2]), end_kind(terms[outgoing][1])))
                    break
            if zero_vertex is not None:
                reasons.add('vertex ' + ','.join(zero_vertex) + ' vanishes')
            else:
                reasons.add('product of nonzero vertices cancels after completeness, charge rule and trace identities')
    return sorted(reasons)


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classify(graph, group):
    """Status of one graph, its D4 term and the size-rule value in every routing."""
    row = {**graph, 'group': group, 'selection': 'all', 'D4 term': None, 'D4 status': 'not meson-only'}
    raw = [TraceSum() for _ in routings(graph)]
    for index, (line, orientation) in enumerate(routings(graph)):
        for _, _, value in term_values(graph, setup(group), line, orientation):
            raw[index] = raw[index] + value
    if all(value.simplified().is_zero() for value in raw):
        return {**row, 'status': 'vanishes identically', 'routing values': [],
                'reasons': vanishing_reasons(graph, group)}
    if set(graph['lines']) == {'M'}:
        d4 = routing_values(graph, group, 'D4')
        if all(value.is_zero() for value in d4):
            row['D4 status'] = 'all-generator-sum term is zero'
        else:
            d4_lambda_kept = routing_values(graph, group, 'D4', eliminate_lambda=False)
            criterion = all(satisfies_d4_criterion(value, group) for value in d4_lambda_kept) and all(
                (one - two).is_zero() for one, two in zip(d4, d4_lambda_kept))
            if criterion:
                row.update(selection='without D4', **{'D4 term': d4, 'D4 status': 'excluded (D4)'})
            else:
                row['D4 status'] = 'all-generator-sum term kept (criterion fails)'
    values = routing_values(graph, group, row['selection'])
    row['routing values'] = values
    if all(value.is_zero() for value in values):
        status = 'zero size-rule value' if group == 'one_particle_irreducible' else 'mass shift only (dZ = 0)'
    else:
        status = 'kept' if row['selection'] == 'all' else 'kept (D4 term excluded)'
    row['status'] = status
    row['scales'] = contributing_scales(graph, group, row['selection'])
    return row


def contributing_scales(graph, group, selection):
    """Scales m_term of the terms whose size-rule value is nonzero in some routing."""
    scales = set()
    legs = setup(group)
    for line, orientation in routings(graph):
        for choice, scale, value in term_values(graph, legs, line, orientation):
            if selection == 'without D4' and is_all_generator_sum(choice):
                continue
            if sized(graph, group, scale, value).terms:
                scales.add(scale)
    return sorted(scales, key=sp.default_sort_key)


def contributing_pieces(row):
    """For each meson line of a kept graph: the pieces whose terms give a nonzero
    size-rule value in some routing (D4 term removed where excluded)."""
    graph, group = row, row['group']
    legs = setup(group)
    found = {j: set() for j, line_type in enumerate(graph['lines']) if line_type == 'M'}
    for line, orientation in routings(graph):
        for choice, scale, value in term_values(graph, legs, line, orientation):
            if row['selection'] == 'without D4' and is_all_generator_sum(choice):
                continue
            if sized(graph, group, scale, value).terms:
                for j in found:
                    found[j].add(choice[j].name.split()[0])
    return {j: sorted(names) for j, names in found.items()}


def label_rows(rows, prefix):
    for number, row in enumerate(rows, start=1):
        row['label'] = f'{prefix}{number}'
    return rows


LABEL_PREFIX = {'one_particle_irreducible': 'C', 'sigma': 'Csigma', 'pi': 'Cpi'}


@functools.lru_cache(maxsize=None)
def _table(group):
    graphs, distributions = enumerate_graphs(external_names(group))
    rows = [classify(graph, group) for graph in graphs]
    return tuple(label_rows(rows, LABEL_PREFIX[group])), distributions


def table(group):
    """Classified graphs of a group ('one_particle_irreducible', 'sigma', 'pi'), labelled
    C1, ...; Csigma1, ...; Cpi1, ... (cached; do not modify the rows)."""
    return list(_table(group)[0])


def summary():
    data = {}
    for group in GROUPS:
        data[group] = table(group)
        data[group + ' distributions'] = _table(group)[1]
    return data


def row_by_label(label):
    for group in GROUPS:
        for row in table(group):
            if row['label'] == label:
                return row
    raise KeyError(label)


def kept(rows):
    return [row for row in rows if row['status'] in ('kept', 'kept (D4 term excluded)')]


def class_counts(rows):
    """{topology: {'graphs', 'kept', 'vanishing', 'zero size-rule value or mass shift', 'D4 excluded'}}."""
    counts = {}
    for row in rows:
        entry = counts.setdefault(topology(row), {'graphs': 0, 'kept': 0, 'vanishing identically': 0,
                                                  'zero size-rule value or mass shift only': 0,
                                                  'D4 term excluded': 0})
        entry['graphs'] += 1
        if row['status'] in ('kept', 'kept (D4 term excluded)'):
            entry['kept'] += 1
        elif row['status'] == 'vanishes identically':
            entry['vanishing identically'] += 1
        else:
            entry['zero size-rule value or mass shift only'] += 1
        if row['D4 status'] == 'excluded (D4)':
            entry['D4 term excluded'] += 1
    return counts


# ---------------------------------------------------------------------------
# Absence of a box and of a 1PI (2,5) graph
# ---------------------------------------------------------------------------

def no_box_and_no_two_point_vertex():
    """(box count, number of (2,5) configurations found among the allowed graphs,
    number of configurations with a two-leg vertex if such vertices were allowed).
    A two-leg vertex is the quadratic Lagrangian, which is entirely in the propagators."""
    graphs, distributions = enumerate_graphs(external_names('one_particle_irreducible'), max_vertices=4)
    boxes = sum(1 for graph in graphs if len(graph['blocks']) >= 4)
    two_five = sum(1 for graph in graphs if sorted(vertex_leg_counts(graph)) == [2, 5])
    with_two_point, _ = enumerate_graphs(external_names('one_particle_irreducible'), max_vertices=2,
                                         allow_two_point_vertices=True)
    excluded = sum(1 for graph in with_two_point if 2 in vertex_leg_counts(graph))
    return {'box graphs': boxes, 'distributions with four vertices': distributions[4],
            '(2,5) graphs': two_five, 'configurations with a two-leg vertex (not vertices)': excluded}


# ---------------------------------------------------------------------------
# Quartic V_phi vertices (lambda) in the self-energies
# ---------------------------------------------------------------------------

def quartic_radial_potential_vertices():
    """[(group, label, vertex index)] of every four-leg vertex of every graph at which a
    quartic V_phi vertex (-6 lambda, -2 lambda xi xi, -6 lambda xi^4) is nonzero for some
    line ends.  V_phi has no meson, so every pseudoscalar leg must enter through its p
    component; the charged legs have none.  Checked with radial_potential_vertex (lambda kept)
    on the external legs and generic line ends carrying a p component."""
    found = []
    for group in GROUPS:
        legs = setup(group)
        for row in table(group):
            for v, block in enumerate(row['blocks']):
                if len(block) != 2:
                    continue
                for ends in itertools.product(('radial', 'p'), repeat=2):
                    vertex_legs = [legs[name] for name in block] + [end_leg(end, None, {}) for end in ends]
                    if cartesian.radial_potential_vertex(vertex_legs, eliminate_lambda=False).simplified().terms:
                        found.append((group, row['label'], v))
                        break
    return found


def lambda_part(row, group):
    """The lambda-dependent part of sum_terms T in every routing (lambda kept in the vertices, unsized)."""
    legs = setup(group)
    parts = []
    for line, orientation in routings(row):
        total = TraceSum()
        for _, _, value in term_values(row, legs, line, orientation, eliminate_lambda=False):
            total = total + value
        parts.append(total.map_coefficients(
            lambda value: sp.expand(value) - sp.expand(value).subs(QUARTIC, 0)).simplified())
    return parts


def kept_graphs_with_lambda():
    """Labels of the kept graphs whose T depends on lambda when lambda is kept in the vertices
    (through the cubic V_phi vertices, since the quartic ones only reach the s tadpoles)."""
    return [row['label'] for group in GROUPS for row in kept(table(group))
            if any(part.terms for part in lambda_part(row, group))]


# ---------------------------------------------------------------------------
# Four-meson kinetic vertex
# ---------------------------------------------------------------------------

def can_contain_kinetic_vertex(graph, group):
    """True if some vertex has four legs that can all carry a meson matrix: its external
    legs are pseudoscalar and both adjacent lines are ALP-like or meson lines.  The
    four-meson kinetic vertex is the only momentum-dependent Cartesian vertex."""
    legs = setup(group)
    count = len(graph['blocks'])
    for v, block in enumerate(graph['blocks']):
        if len(block) == 2 and not any(legs[name].radial for name in block):
            if graph['lines'][(v - 1) % count] in 'AM' and graph['lines'][v] in 'AM':
                return True
    return False


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
LINE_NAMES = {'S': 's', 'A': 'ALP-like', 'M': 'meson'}


def describe(row):
    blocks = ' | '.join('{' + ','.join(block) + '}' for block in row['blocks'])
    lines = ','.join(LINE_NAMES[line] for line in row['lines'])
    return f"{row['label']}: {topology(row)} [{blocks}] lines ({lines}) S={row['symmetry']} r={len(row['blocks'])}"


LEG_LATEX = {'s': 's', 'pi+': r'\pi^+', 'pi-': r'\pi^-', 'in': r'\text{in}', 'out': r'\text{out}'}
LINE_LATEX = {'S': '$s$', 'A': 'ALP-like', 'M': 'meson'}
STATUS_TEXT = {'kept': 'kept', 'kept (D4 term excluded)': 'kept, D4 term excluded',
               'vanishes identically': 'vanishes identically', 'mass shift only (dZ = 0)': 'mass shift only ($\\delta Z=0$)',
               'zero size-rule value': 'zero'}


def label_latex(label):
    if label.startswith('Csigma'):
        return r'$C\sigma' + label[6:] + '$'
    if label.startswith('Cpi'):
        return r'$C\pi' + label[3:] + '$'
    return f'${label}$'


def scale_text(scale):
    return '$' + sp.latex(scale, symbol_names=cartesian.LATEX_NAMES).replace(r'\max', r'\max') + '$'


def reason_text(reasons):
    """Vanishing reasons in words: the vertex (external legs; line ends) that vanishes by the charge rule."""
    names = {'pi+': r'\pi^+', 'pi-': r'\pi^-', 'in': r'\pi^+', 'out': r'\pi^-', 's': 's', 'p': 'p',
             'v_a': 'v_a', 'MinvC': r'\mathbb M^{-1}_{m_a^2}(\mathcal C)', 'generator': 'generator'}
    texts = []
    for reason in reasons:
        ends = reason[len('vertex '):-len(' vanishes')].split(',')
        texts.append('$V(' + ','.join(names[end] for end in ends) + ')=0$')
    return ', '.join(texts)


def markdown_table(group):
    """Markdown rows of the graph table of one group (script output for the graph-list record)."""
    lines = ['| graph | vertices (external legs) | lines | $S_\\Gamma$ | $r_\\Gamma$ | status | $m_{\\rm term}$ '
             '| pieces per meson line | routing |', '|---|---|---|---|---|---|---|---|---|']
    for row in table(group):
        vertices = ', '.join('$\\{' + ','.join(LEG_LATEX[name] for name in block) + '\\}$' for block in row['blocks'])
        line_text = ', '.join(LINE_LATEX[line] for line in row['lines'])
        status = STATUS_TEXT[row['status']]
        if row['status'] == 'vanishes identically':
            scales, pieces, routing = '–', '–', reason_text(row['reasons'])
        elif row['status'] in ('kept', 'kept (D4 term excluded)'):
            scales = ', '.join(scale_text(scale) for scale in row['scales'])
            found = contributing_pieces(row)
            pieces = '; '.join(f'line {j}: ' + ','.join(names) for j, names in found.items()) or '–'
            values = row['routing values']
            dependent = not all((value - values[0]).is_zero() for value in values)
            routing = 'depends on the routing' if dependent else 'independent'
        else:
            scales, pieces, routing = '–', '–', 'independent'
        lines.append(f"| {label_latex(row['label'])} | {vertices} | {line_text} | ${sp.latex(row['symmetry'])}$ | "
                     f"{len(row['blocks'])} | {status} | {scales} | {pieces} | {routing} |")
    return '\n'.join(lines)


def d4_markdown():
    """Compact form of every excluded D4 term, per routing where it depends on the routing."""
    blocks = []
    for group in GROUPS:
        for row in table(group):
            if row['D4 status'] != 'excluded (D4)':
                continue
            values = [routing_value(row, group, line, orientation, 'D4', construction='trees')
                      for line, orientation in routings(row)]
            independent = all((value.to_trace_sum() - row['D4 term'][0]).is_zero() for value in values)
            shown = values[:1] if independent else values
            for index, value in enumerate(shown):
                value = value.map_coefficients(lambda c: sp.factor(sp.expand(c)))
                routing_text = '' if independent else rf'\big|_{{\rho={index + 1}}}'
                name = 'D4' + routing_text
                blocks.append(f"{label_latex(row['label'])}" + (' (all routings)' if independent else '')
                              + ':\n$$\n\\begin{aligned}\n' + name + '={}&'
                              + cartesian.compact_latex(value, separator=' \\\\\n&+') + '\n\\end{aligned}\n$$')
    return '\n\n'.join(blocks)


def main():
    if '--markdown' in sys.argv:
        for group in GROUPS:
            print(markdown_table(group) + '\n')
        print(d4_markdown())
        return
    print('Absence of a box and of a 1PI (2,5) graph:', no_box_and_no_two_point_vertex())
    data = summary()
    for group in GROUPS:
        rows = data[group]
        print(f'\n== {group}: {len(rows)} graphs; leg distributions per V: {data[group + " distributions"]}')
        for topology_name, counts in class_counts(rows).items():
            print(f'   {topology_name}: {counts}')
        for row in rows:
            text = describe(row) + f" | {row['status']} | D4: {row['D4 status']}"
            if row['status'] == 'vanishes identically':
                text += ' | ' + '; '.join(row['reasons'])
            else:
                text += ' | m_term: ' + ', '.join(str(scale) for scale in row['scales'])
                text += f" | kinetic vertex possible: {can_contain_kinetic_vertex(row, group)}"
                if 'M' in row['lines'] and row in kept(rows):
                    text += f' | pieces per meson line: {contributing_pieces(row)}'
            print('  ' + text)
        for row in rows:
            if row['D4 status'] == 'excluded (D4)':
                print(f"  D4 term of {row['label']} (routing 1): {cartesian.trace_sum_latex(row['D4 term'][0])}")
    print('\nFour-leg vertices at which a quartic V_phi vertex can be nonzero:', quartic_radial_potential_vertices())
    print('Kept graphs whose T depends on lambda (cubic V_phi vertices):', kept_graphs_with_lambda())
    for row in data['sigma']:
        if topology(row) == 'tadpole':
            parts = lambda_part(row, 'sigma')
            has_p2 = any(sp.sympify(value).has(MOMENTUM_SQUARED) for part in parts for value in part.terms.values())
            print(f"  {row['label']}: lambda part of T nonzero: {any(part.terms for part in parts)}, "
                  f"depends on p^2: {has_p2}, status: {row['status']}")


if __name__ == '__main__':
    main()
