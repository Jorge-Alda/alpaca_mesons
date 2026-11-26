"""One-loop graphs of sigma -> pi+ pi-: enumeration, symmetry factors, exact split of the meson lines, per-term loop scales and termwise routing maximum (tasks 2.1, 4.4).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_graphs.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_graphs

Units GeV; conventions of sigma_pipi_loop_algebra.py.  A one-loop graph is a
cycle of V vertices; vertex v_i carries a nonempty set of external legs and
the internal lines e_i (from v_i to v_{i+1}) close the cycle.  Line types:
  S  radial line (sigma at both ends),
  A  ALP-like line (pole vector v_a at both ends),
  M  meson line, residue (1_9 + 1) - v_a v_a^T split exactly into three pieces:
       'generator sum'   the bare U(3) generator sum 1_9 at both ends,
       'ALP field'       the ALP-field entry 1 at both ends,
       'ALP-like pole'   minus v_a at both ends;
  F  bare ALP field at both ends (only for the field-basis check).
S_Gamma = 1/|Aut|, with Aut the elements of the dihedral group of the cycle
(acting on half-edges) that keep the external legs and line types.
r_Gamma = V = number of internal lines.

Selection (record murayama-sigma-pipi-loop-graph-selection.md):
 * every graph that does not vanish identically is kept; no graph and no term
   is dropped by a series in eps or 1/eps;
 * the only excluded piece (D4) is the term of a meson-only graph in which
   every meson line takes the generator sum and the only eps is the single
   sigma end of the tree; the restated criterion is checked: no letter of the
   ALP sector (M^{-1}_{m_a^2}(C), N, P_hat_01, S_hat_{r,s>=1}, P_hat_{r,s>=1})
   and the ratio to M_0 depends on eps only through p_sigma^2 = m_sigma^2;
 * D11 (assumption): CP-invariant sources with sigma CP-even, a and Phi CP-odd,
   and flavour-diagonal sources: the one-particle-reducible mixing graphs of
   the external legs vanish; they are not enumerated.
Loop scale (D3): every term (one piece per line) has the scale m_term = the
largest scale among its lines: m_sigma (radial line), m_a (ALP-like line, and
the ALP-field entry and ALP-like pole of a meson line, i.e. the pieces
(ii)+(iii)+(iv) = (0 + 1) - v_a v_a^T), m_P (bare generator sum; m_P is the
largest pole mass of the nine meson-like states, always symbolic).  Products
with the loop momentum, k.k and k.p, are set to m_term^2; the size factor is
S_Gamma (m_term^2)^{2-r_Gamma}/(16 pi^2) per term, and dZ = -Pi'(m^2) is taken
term by term at m_term.
Routing (D12): for each of the 2 r_Gamma routings (loop momentum on each
internal line, both orientations) the size-rule value is sum_t C_t^rho t; the
graph enters with sum_t C_t^{rho*_t} t, rho*_t the routing of largest
|C_t^rho| (termwise_maximum).  It is an upper bound per graph on the routing
ambiguity of the size rule, not on the loop.  Ties of equal magnitude and
opposite sign take the sign of the sum of the coefficients over all routings;
if this sum vanishes, the first routing in the canonical order is used and the
product is marked (a conventional choice at a degenerate mass relation).
The behaviour of the kept graphs at large eps is derived in
derive_sigma_pipi_loop_large_epsilon.py.  No loop integral is evaluated.
"""
import itertools

import sympy as sp

import sigma_pipi_loop_algebra as algebra
from sigma_pipi_loop_algebra import (EPSILON, M_A, M_SIGMA, MOMENTUM_SQUARED, TraceSum, CompactSum,
                                     alp_like_leg, bare_angular_leg, make_dot, meson_leg, negate,
                                     momentum_sum, radial_leg)
import derive_sigma_pipi_loop_tree as tree_script

LINE_TYPES = ('S', 'A', 'M')
GENERATOR_SUM, ALP_FIELD, ALP_LIKE_POLE = 'generator sum', 'ALP field', 'ALP-like pole'
RADIAL_LINE, ALP_LIKE_LINE = 'radial', 'ALP-like line'
M_P = algebra.M_MESON
# Loop scale of each line piece (D3).  The ALP-field entry and the ALP-like
# pole of a meson line together are pieces (ii)+(iii)+(iv) = (0 + 1) - v_a v_a^T,
# the remainder of the ALP-like propagator: scale m_a.
PIECE_SCALE = {RADIAL_LINE: M_SIGMA, ALP_LIKE_LINE: M_A, GENERATOR_SUM: M_P, ALP_FIELD: M_A, ALP_LIKE_POLE: M_A}


# ---------------------------------------------------------------------------
# Enumeration
# ---------------------------------------------------------------------------

def dihedral_images(blocks, lines):
    """All 2V images (rotations and reflections) of a cycle configuration."""
    count = len(blocks)
    images = []
    for shift in range(count):
        for reflect in (False, True):
            if reflect:
                new_blocks = tuple(blocks[(shift - i) % count] for i in range(count))
                new_lines = tuple(lines[(shift - i - 1) % count] for i in range(count))
            else:
                new_blocks = tuple(blocks[(shift + i) % count] for i in range(count))
                new_lines = tuple(lines[(shift + i) % count] for i in range(count))
            images.append((new_blocks, new_lines))
    return images


def surjections(externals, vertices):
    """Assignments of the external legs to the cycle vertices, each vertex
    receiving at least one leg (a two-point vertex does not exist)."""
    result = []
    for assignment in itertools.product(range(vertices), repeat=len(externals)):
        if set(assignment) == set(range(vertices)):
            result.append(tuple(tuple(sorted(e for e, a in zip(externals, assignment) if a == v))
                                for v in range(vertices)))
    return result


def enumerate_cycle_graphs(externals, max_vertices=5):
    """Distinct one-loop 1PI graphs with the given external legs.
    Returns (graphs, surjection count per V)."""
    counts = {}
    seen = {}
    for vertices in range(1, max_vertices + 1):
        assignments = surjections(externals, vertices)
        counts[vertices] = len(assignments)
        for blocks in assignments:
            for lines in itertools.product(LINE_TYPES, repeat=vertices):
                images = dihedral_images(blocks, lines)
                canonical = min(images)
                if canonical in seen:
                    continue
                automorphisms = sum(1 for image in images if image == (blocks, lines))
                seen[canonical] = {'blocks': canonical[0], 'lines': canonical[1],
                                   'symmetry': sp.Rational(1, automorphisms)}
    return list(seen.values()), counts


def topology(graph):
    return {1: 'tadpole', 2: 'bubble', 3: 'triangle'}.get(len(graph['blocks']), 'box or larger')


def line_content(graph):
    return ''.join(sorted(graph['lines'], key=LINE_TYPES.index))


def is_meson_only(graph):
    return set(graph['lines']) == {'M'}


def label(graph):
    blocks = ['{' + ','.join(block) + '}' for block in graph['blocks']]
    lines = graph['lines']
    parts = [blocks[i] + '-' + lines[i] + '-' for i in range(len(blocks))]
    return topology(graph) + ' ' + ''.join(parts)


# ---------------------------------------------------------------------------
# Graph values (construction 1: words; construction 2: compact trees)
# ---------------------------------------------------------------------------

def line_pieces(line_type, generator):
    """[(sign, leg factory(momentum), generator or None, piece name)] for one line type."""
    if line_type == 'S':
        return [(1, radial_leg, None, RADIAL_LINE)]
    if line_type == 'A':
        return [(1, alp_like_leg, None, ALP_LIKE_LINE)]
    if line_type == 'F':      # bare ALP field at both ends (used by the field-basis checks)
        return [(1, bare_angular_leg, None, ALP_FIELD)]
    return [(1, lambda momentum: meson_leg(generator, momentum), generator, GENERATOR_SUM),
            (1, bare_angular_leg, None, ALP_FIELD),
            (-1, alp_like_leg, None, ALP_LIKE_POLE)]


def all_routings(graph):
    """Loop momentum on each internal line, in both orientations."""
    return [(start, orientation) for start in range(len(graph['lines'])) for orientation in (1, -1)]


def routing(graph, external_momenta, start=0, orientation=1):
    """Loop momentum k on line `start`, oriented along the cycle (orientation=1)
    or against it (orientation=-1); the other lines carry +-k plus external
    momenta, fixed by momentum conservation."""
    blocks = graph['blocks']
    count = len(blocks)
    momenta = [None] * count
    momenta[start] = {'k': orientation}
    for step in range(1, count):
        index = (start + step) % count
        preceding = (index - 1) % count
        incoming = momentum_sum(*(external_momenta[name] for name in blocks[index]))
        momenta[index] = momentum_sum(momenta[preceding], incoming)
    return momenta


def loop_dot(external_basis_dot, mass):
    """k.k -> m^2, k.p -> m^2 for every external basis momentum p; external
    products exact."""
    def basis(first, second):
        if first == 'k' or second == 'k':
            return mass ** 2
        return external_basis_dot(first, second)
    return make_dot(basis)


def is_bare_choice(choice):
    """True for the piece in which every line is a meson line taking the generator sum."""
    return all(piece[3] == GENERATOR_SUM for piece in choice)


def end_leg(piece, end, momentum):
    """Leg of one end of a line piece.  piece[1] is either one leg factory (the
    same leg at both ends) or a pair (factory at the start of the line, factory
    at its end), used for pieces that differ at the two ends of the line."""
    factory = piece[1]
    if isinstance(factory, tuple):
        factory = factory[0] if end == 'start' else factory[1]
    return factory(momentum)


def piece_scale(piece):
    """Loop scale of one line piece (D3): m_sigma (radial line), m_a (ALP-like
    line and every piece of a meson line other than the bare generator sum),
    m_P (bare U(3) generator sum).  A fifth entry of the piece tuple, if
    present, overrides the table (used by the four-piece split)."""
    if len(piece) > 4:
        return piece[4]
    return PIECE_SCALE[piece[3]]


def term_scale(choice):
    """m_term = largest scale among the lines of one term (one piece per line)."""
    return sp.Max(*{piece_scale(piece) for piece in choice})


def graph_value(graph, external_legs, external_basis_dot, construction='words',
                start=0, orientation=1, pieces='all', piece_function=line_pieces, choice_filter=None,
                by_scale=False):
    """T_Gamma for one routing: product of the vertices with completeness for
    the generator sums, the assumptions and the charge rule applied.
    Every term (one piece per line) uses its own scale m_term (term_scale) in
    k.k, k.p -> m_term^2.  With by_scale=True the result is {m_term: value};
    otherwise the sum over the scales.
    pieces = 'all' (every term of the line residues), 'bare' (only the term in
    which every line takes the generator sum) or 'rest' (all but that term).
    piece_function(line type, generator) gives the pieces of a line (default:
    line_pieces); choice_filter(choice), if given, keeps only the choices of
    one piece per line for which it is true.
    The symmetry factor and the loop size are not included."""
    blocks, lines = graph['blocks'], graph['lines']
    count = len(blocks)
    external_momenta = {name: leg.momentum for name, leg in external_legs.items()}
    momenta = routing(graph, external_momenta, start, orientation)
    generators = [f'G{i + 1}' for i in range(count)]
    piece_lists = [piece_function(lines[i], generators[i]) for i in range(count)]
    empty = TraceSum if construction == 'words' else CompactSum
    totals = {}
    for choice in itertools.product(*piece_lists):
        bare = is_bare_choice(choice)
        if (pieces == 'bare' and not bare) or (pieces == 'rest' and bare):
            continue
        if choice_filter is not None and not choice_filter(choice):
            continue
        scale = term_scale(choice)
        dot = loop_dot(external_basis_dot, scale)
        sign = 1
        for piece in choice:
            sign *= piece[0]
        product = empty.scalar(sign)
        for v in range(count):
            incoming_line = (v - 1) % count
            outgoing_line = v
            legs = [external_legs[name] for name in blocks[v]]
            legs.append(end_leg(choice[incoming_line], 'end', momenta[incoming_line]))
            legs.append(end_leg(choice[outgoing_line], 'start', negate(momenta[outgoing_line])))
            if construction == 'words':
                factor = algebra.drop_vanishing_letters(algebra.vertex(legs, dot))
            else:
                factor = algebra.compact_assumptions(algebra.compact_vertex(legs, dot))
            product = product * factor
        used = [piece[2] for piece in choice if piece[2] is not None]
        for generator in used:
            if construction == 'words':
                product = algebra.completeness(product, generator)
            else:
                product = algebra.compact_completeness(product, generator)
        totals[scale] = totals.get(scale, empty()) + product
    if construction == 'words':
        totals = {scale: algebra.apply_charge_rule(algebra.apply_assumptions(value)).simplified()
                  for scale, value in totals.items()}
    else:
        totals = {scale: algebra.compact_assumptions(value).simplified() for scale, value in totals.items()}
    totals = {scale: value for scale, value in totals.items() if value.terms}
    if by_scale:
        return totals
    total = empty()
    for value in totals.values():
        total = total + value
    return total.simplified()


def size_factor(graph, mass):
    """S_Gamma (m^2)^{2-r_Gamma}/(16 pi^2)."""
    return graph['symmetry'] * (mass ** 2) ** (2 - len(graph['blocks'])) / (16 * sp.pi ** 2)


def size_rule_value(graph, external_legs, external_basis_dot, group, construction='words', start=0,
                    orientation=1, pieces='all', **options):
    """Size-rule value of one routing, every term at its own scale m_term:
    sum_terms S_Gamma (m_term^2)^{2-r}/(16 pi^2) T_term for a 1PI graph, and
    dZ = sum_terms of the derivative rule (residue_correction) at m_term for a
    self-energy.  options (piece_function, choice_filter) go to graph_value."""
    values = graph_value(graph, external_legs, external_basis_dot, construction, start=start,
                         orientation=orientation, pieces=pieces, by_scale=True, **options)
    total = TraceSum() if construction == 'words' else CompactSum()
    for scale, value in values.items():
        if group == 'one_particle_irreducible':
            total = total + value.scale(size_factor(graph, scale))
        else:
            total = total + residue_correction(graph, value, group, scale)
    return total.simplified()


def contributing_scales(graph, external_legs, external_basis_dot, group, pieces='all'):
    """Scales m_term whose terms give a nonzero size-rule value (X or dZ) in some routing."""
    scales = set()
    for start, orientation in all_routings(graph):
        values = graph_value(graph, external_legs, external_basis_dot, start=start, orientation=orientation,
                             pieces=pieces, by_scale=True)
        for scale, value in values.items():
            if group == 'one_particle_irreducible' or not residue_correction(graph, value, group, scale).is_zero():
                scales.add(scale)
    return sorted(scales, key=sp.default_sort_key)


def routing_values(graph, external_legs, external_basis_dot, group, construction='words', pieces='all', **options):
    """Size-rule value of every routing, in the canonical order of all_routings."""
    return [size_rule_value(graph, external_legs, external_basis_dot, group, construction, start=start,
                            orientation=orientation, pieces=pieces, **options)
            for start, orientation in all_routings(graph)]


# ---------------------------------------------------------------------------
# Termwise routing maximum (D12)
# ---------------------------------------------------------------------------
# For every independent product of traces t (a key of the words TraceSum, the
# canonical generic-source basis) the routing rho*_t with the largest |C_t^rho|
# is used, and the graph enters as sum_t C_t^{rho*_t} t (sign kept).  This is
# an upper bound per graph on the routing ambiguity of the size rule,
# |sum_t C_t^rho t| <= sum_t |C_t^{rho*_t}| |t| for every rho; it is not a
# bound on the loop.  Ties of equal magnitude and opposite sign take the sign
# of the sum of the coefficients over all routings; if this sum vanishes, the
# first routing in the canonical order is used and the product is marked (a
# conventional choice at a degenerate mass relation).  Where rho*_t depends on the
# masses, the coefficient is a Piecewise expression whose conditions are
# homogeneous inequalities between squared masses.
BRANCH_VARIABLE = sp.Symbol('w', positive=True)   # dimensionless ratio of squared masses


def is_number(value):
    return not sp.sympify(value).free_symbols


def tie_choice(values, candidates):
    """Index among `candidates` (indices of values of equal magnitude) chosen by
    the tie rule; returns (index, marked)."""
    distinct = {sp.nsimplify(values[i]) for i in candidates}
    if len(distinct) == 1:
        return candidates[0], False
    routing_sum = sum(values)
    if routing_sum == 0:
        return candidates[0], True
    sign = 1 if routing_sum > 0 else -1
    return next(i for i in candidates if values[i] * sign > 0), False


def largest_numerical(values):
    """(index, marked) of the largest magnitude among exact numbers."""
    magnitudes = [abs(value) for value in values]
    largest = max(magnitudes)
    candidates = [i for i, magnitude in enumerate(magnitudes) if magnitude == largest]
    return tie_choice(values, candidates)


def largest_on_half_line(offsets, slopes):
    """Largest |offsets[i] + slopes[i] w| for w > 0 (exact rationals).
    Break points are the positive roots of v_i = +-v_j and of the sum over
    routings; between them the choice is constant.  Returns merged segments
    [{'lower', 'upper', 'lower closed', 'upper closed', 'index', 'marked'}]
    covering w > 0; a break point is attached to a neighbouring segment when
    the chosen value there coincides with it.  'marked' is True if the tie rule
    had to use the first routing (sum over routings zero) somewhere in the
    segment."""
    count = len(offsets)
    w = BRANCH_VARIABLE
    functions = [offsets[i] + slopes[i] * w for i in range(count)]
    routing_sum = sum(functions)
    roots = set()
    for i in range(count):
        for j in range(i + 1, count):
            for sign in (1, -1):
                difference = sp.expand(functions[i] - sign * functions[j])
                if difference.has(w):
                    roots |= set(sp.solve(difference, w))
    if sp.expand(routing_sum).has(w):
        roots |= set(sp.solve(routing_sum, w))
    points = sorted(root for root in roots if root > 0)
    edges = [sp.Integer(0)] + points + [sp.oo]
    segments = []
    for lower, upper in zip(edges[:-1], edges[1:]):
        test = lower + 1 if upper == sp.oo else (lower + upper) / 2
        index, marked = largest_numerical([f.subs(w, test) for f in functions])
        segments.append({'lower': lower, 'upper': upper, 'lower closed': False, 'upper closed': False,
                         'index': index, 'marked': marked})
    merged = [segments[0]]
    for point, right in zip(points, segments[1:]):
        index, marked = largest_numerical([f.subs(w, point) for f in functions])
        value = functions[index].subs(w, point)
        left = merged[-1]
        if functions[left['index']].subs(w, point) == value:
            left['upper closed'] = True
            left['marked'] = left['marked'] or marked
        elif functions[right['index']].subs(w, point) == value:
            right['lower closed'] = True
            right['marked'] = right['marked'] or marked
        else:
            merged.append({'lower': point, 'upper': point, 'lower closed': True, 'upper closed': True,
                           'index': index, 'marked': marked})
            left = merged[-1]
        if left['index'] == right['index'] and (left['upper closed'] or right['lower closed']):
            left.update(upper=right['upper'], **{'upper closed': right['upper closed']},
                        marked=left['marked'] or right['marked'])
        else:
            merged.append(right)
    return merged


def scale_atoms(expression):
    """Loop-scale atoms of an expression: Max(...) objects and m_a, m_P outside them."""
    expression = sp.sympify(expression)
    maxima = expression.atoms(sp.Max)
    stripped = expression.subs({atom: sp.Dummy() for atom in maxima})
    return maxima | (stripped.free_symbols & {M_A, M_P})


def affine_reduction(coefficients):
    """Write routing coefficients that depend on the masses as
    C_rho = g (a_rho + b_rho w) times a positive factor, with
    w = (m_sigma^2 + gamma m_pi^2)/s^2, s the single loop scale.
    Returns (offsets, slopes, scale, gamma).  Raises ValueError otherwise."""
    common = sp.gcd_list([c for c in coefficients if c != 0])
    cofactors = [sp.cancel(c / common) for c in coefficients]
    scales = set()
    for cofactor in cofactors:
        scales |= scale_atoms(cofactor)
    if len(scales) != 1:
        raise ValueError(f'routing coefficients with {len(scales)} loop scales')
    scale = scales.pop()
    placeholder = sp.Dummy('scale', positive=True)
    cofactors = [sp.cancel(c.subs(scale, placeholder).subs(placeholder, 1)) for c in cofactors]
    denominator = sp.lcm_list([sp.fraction(c)[1] for c in cofactors])
    assert all(factor.is_positive for factor in sp.Mul.make_args(denominator)), denominator
    numerators = [sp.expand(c * denominator) for c in cofactors]
    radial_squared, pion_squared = sp.symbols('X Y', positive=True)
    offsets, slopes, gradients = [], [], []
    for numerator in numerators:
        polynomial = sp.Poly(numerator.subs({M_SIGMA ** 2: radial_squared, tree_script.M_PI ** 2: pion_squared}),
                             radial_squared, pion_squared)
        assert polynomial.total_degree() <= 1 and not numerator.free_symbols - {M_SIGMA, tree_script.M_PI}, numerator
        offsets.append(polynomial.coeff_monomial(1))
        gradients.append((polynomial.coeff_monomial(radial_squared), polynomial.coeff_monomial(pion_squared)))
    reference = next(gradient for gradient in gradients if gradient != (0, 0))
    assert reference[0] != 0
    gamma = reference[1] / reference[0]
    for gx, gy in gradients:
        if gx * reference[1] - gy * reference[0] != 0:
            raise ValueError('routing coefficients depend on two combinations of the masses')
        slopes.append(gx)
    # w > 0 for every physical point: m_sigma >= 2 m_pi gives m_sigma^2 + gamma m_pi^2 >= (4 + gamma) m_pi^2
    assert gamma >= -4, gamma
    return offsets, slopes, scale, gamma


def segment_condition(segment, scale, gamma):
    """Condition of a segment of w = (m_sigma^2 + gamma m_pi^2)/scale^2 as
    homogeneous relations between squared masses."""
    combination = M_SIGMA ** 2 + gamma * tree_script.M_PI ** 2
    lower, upper = segment['lower'], segment['upper']
    if lower == upper:
        return sp.Eq(combination, lower * scale ** 2)
    parts = []
    if lower != 0:
        relation = sp.GreaterThan if segment['lower closed'] else sp.StrictGreaterThan
        parts.append(relation(combination, lower * scale ** 2))
    if upper != sp.oo:
        relation = sp.LessThan if segment['upper closed'] else sp.StrictLessThan
        parts.append(relation(combination, upper * scale ** 2))
    return sp.And(*parts) if parts else sp.true


def termwise_maximum(per_routing):
    """Termwise routing maximum of a list of TraceSums (one per routing, in the
    canonical order).  Returns a dict with
      'value'     TraceSum sum_t C_t^{rho*_t} t (Piecewise where rho* depends on masses),
      'choices'   {t: description}, with the chosen routing(s),
      'counts'    numbers of independent, numerical and mass-dependent products,
      'marked'    numerical products decided by the first routing (sum over
                  routings zero),
      'marked boundary'  mass-dependent products in which the first routing
                  decides at a boundary point between two branches,
      'spread'    exact ratios C_t^rho/C_t^{rho*} for the numerical products."""
    keys = set()
    for value in per_routing:
        keys |= set(value.simplified().terms)
    result = TraceSum()
    choices, marked, marked_boundary, spread = {}, [], [], {}
    counts = {'independent': 0, 'numerical': 0, 'mass-dependent': 0}
    for key in sorted(keys, key=repr):
        coefficients = [sp.expand(value.terms.get(key, 0)) for value in per_routing]
        if all(sp.expand(c - coefficients[0]) == 0 for c in coefficients):
            counts['independent'] += 1
            result.add_term(key, coefficients[0])
            choices[key] = ('independent', None)
            continue
        reference = next(c for c in coefficients if c != 0)
        ratios = [sp.cancel(c / reference) for c in coefficients]
        if all(is_number(ratio) for ratio in ratios):
            counts['numerical'] += 1
            index, tie = largest_numerical(ratios)
            result.add_term(key, coefficients[index])
            choices[key] = ('numerical', index)
            spread[key] = [ratio / ratios[index] for ratio in ratios]
            if tie:
                marked.append(key)
            continue
        counts['mass-dependent'] += 1
        offsets, slopes, scale, gamma = affine_reduction(coefficients)
        segments = largest_on_half_line(offsets, slopes)
        branches = [(coefficients[segment['index']], segment_condition(segment, scale, gamma))
                    for segment in segments]
        result.add_term(key, sp.Piecewise(*branches) if len(branches) > 1 else branches[0][0])
        choices[key] = ('mass-dependent', {'scale': scale, 'gamma': gamma, 'segments': segments,
                                           'offsets': offsets, 'slopes': slopes})
        if any(segment['marked'] for segment in segments):
            marked_boundary.append(key)
    return {'value': result, 'choices': choices, 'counts': counts, 'marked': marked,
            'marked boundary': marked_boundary, 'spread': spread}


def common_routing(per_routing, maximum):
    """Index of the first routing whose value equals the termwise maximum for
    every product (None if the maximum is not attained by a single routing)."""
    for index, value in enumerate(per_routing):
        if (value - maximum).is_zero():
            return index
    return None


def base_routing(per_routing, maximum):
    """Routing that attains the termwise maximum for the largest number of
    products of traces (first in canonical order on ties)."""
    counts = []
    for value in per_routing:
        difference = (maximum - value).simplified()
        counts.append(sum(1 for key in maximum.terms if key not in difference.terms
                          or expression_is_zero(difference.terms[key])))
    return counts.index(max(counts))


def piecewise_branches(expression):
    """Expressions of the branches of a coefficient (itself if no Piecewise)."""
    expression = sp.piecewise_fold(sp.sympify(expression))
    if isinstance(expression, sp.Piecewise):
        return [piece for piece, _ in expression.args]
    return [expression]


def expression_is_zero(expression):
    """Exact zero test, branch by branch for Piecewise expressions."""
    return all(sp.expand(piece) == 0 for piece in piecewise_branches(sp.expand(expression)))


# ---------------------------------------------------------------------------
# External legs and residues
# ---------------------------------------------------------------------------

def three_point_setup():
    return tree_script.external_legs(), tree_script.external_basis_dot


def two_point_setup(kind):
    """External legs of a self-energy with incoming p and -p."""
    p = {'p': 1}
    if kind == 'sigma':
        legs = {'in': radial_leg(p), 'out': radial_leg(negate(p))}
    elif kind == 'pi+':
        legs = {'in': meson_leg('l+', p), 'out': meson_leg('l-', negate(p))}
    elif kind == 'pi-':
        legs = {'in': meson_leg('l-', p), 'out': meson_leg('l+', negate(p))}
    else:
        raise ValueError(kind)
    return legs, (lambda first, second: MOMENTUM_SQUARED)


EXTERNAL_MASS = {'sigma': M_SIGMA, 'pi+': tree_script.M_PI, 'pi-': tree_script.M_PI}


def residue_correction(graph, value, kind, mass):
    """dZ estimate of the terms of one self-energy graph with loop scale `mass`:
        dZ = -S [theta(r>=2) (m^2)^{2-(r+1)} T + (m^2)^{2-r} dT/dp^2] / (16 pi^2)
    at p^2 = m_ext^2 (prescription dZ = -Pi'(m^2)): the first term is the moment
    with one more propagator (a denominator differentiated), the second the
    explicit p^2 dependence of the numerator.  For tadpoles only the second
    term exists.  `value` must have been computed with k.k, k.p -> mass^2."""
    r = len(graph['blocks'])
    external = EXTERNAL_MASS[kind]
    def coefficient(t):
        derivative = sp.diff(t, MOMENTUM_SQUARED)
        moment = (mass ** 2) ** (1 - r) * t if r >= 2 else 0
        total = -graph['symmetry'] * (moment + (mass ** 2) ** (2 - r) * derivative) / (16 * sp.pi ** 2)
        return sp.expand(total.subs(MOMENTUM_SQUARED, external ** 2))
    return value.map_coefficients(coefficient).simplified()


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def ratio_to_tree_is_epsilon_independent(trace_sum, external):
    """For a 1PI value: T/M_0 is eps-independent iff T = eps x (eps-free),
    because M_0 = eps x (eps-free trace).  For a residue correction dZ (which
    multiplies M_0): dZ is eps-free."""
    divisor = EPSILON if external == 'one_particle_irreducible' else 1
    return all(sp.diff(sp.cancel(value / divisor), EPSILON) == 0
               for value in trace_sum.simplified().terms.values())


ALP_SECTOR_LETTERS = ('Minv', 'P01')


def is_alp_sector_letter(name):
    """Letters of the ALP sector: M^{-1}_{m_a^2}(C), P_hat_{0,1} and every
    S_hat_{r,s}, P_hat_{r,s} with s >= 1."""
    return name in ALP_SECTOR_LETTERS or (algebra.is_source(name) and int(name[2:]) >= 1)


def satisfies_restated_d4(trace_sum, external):
    """Restated D4 criterion (D4 addendum): no letter of the ALP sector, no
    N(m_a^2), and the ratio to M_0 (dZ itself for a self-energy) depends on eps
    only through the kinematics p_sigma^2 = m_sigma^2 (m_sigma is held as a
    symbol here)."""
    if any(is_alp_sector_letter(name) for name in trace_sum.letters()):
        return False
    if any(sp.sympify(value).has(algebra.NORMALIZATION) for value in trace_sum.terms.values()):
        return False
    return ratio_to_tree_is_epsilon_independent(trace_sum, external)


def classify(graph, legs, basis_dot, group):
    """Status, excluded piece and termwise routing maximum of one graph."""
    unsized = [graph_value(graph, legs, basis_dot, start=start, orientation=sign)
               for start, sign in all_routings(graph)]
    row = {**graph, 'group': group, 'pieces': 'all', 'excluded': None, 'excluded_routings': None,
           'excluded_status': 'none'}
    if all(value.is_zero() for value in unsized):
        return {**row, 'status': 'vanishes identically', 'value': TraceSum(),
                'routing_values': [], 'maximum': None}
    if is_meson_only(graph):
        bare_routings = routing_values(graph, legs, basis_dot, group, pieces='bare')
        bare = termwise_maximum(bare_routings)['value']
        if all(value.is_zero() for value in bare_routings):
            row['excluded_status'] = 'zero'
        elif all(satisfies_restated_d4(value, group) for value in bare_routings):
            row.update(pieces='rest', excluded=bare, excluded_routings=bare_routings,
                       excluded_status='excluded (D4)')
        else:
            row['excluded_status'] = 'kept (more than one sigma end)'
    per_routing = routing_values(graph, legs, basis_dot, group, pieces=row['pieces'])
    maximum = termwise_maximum(per_routing)
    value = maximum['value']
    row.update(routing_values=per_routing, maximum=maximum,
               scales=contributing_scales(graph, legs, basis_dot, group, row['pieces']))
    if value.is_zero():
        status = 'kept, no contribution (exact zero)' if group == 'one_particle_irreducible' \
            else 'no residue correction (exact: dZ = 0)'
        return {**row, 'status': status, 'value': value}
    status = 'kept' if row['pieces'] == 'all' else 'kept (D4 term excluded)'
    return {**row, 'status': status, 'value': value}


def one_particle_irreducible_table():
    legs, basis_dot = three_point_setup()
    graphs, counts = enumerate_cycle_graphs(('sigma', 'pi+', 'pi-'))
    rows = [{**classify(graph, legs, basis_dot, 'one_particle_irreducible'), 'label': label(graph)}
            for graph in graphs]
    return rows, counts


def self_energy_table(kind):
    legs, basis_dot = two_point_setup(kind)
    graphs, counts = enumerate_cycle_graphs(('in', 'out'))
    rows = [{**classify(graph, legs, basis_dot, kind), 'label': f'{kind} self-energy ' + label(graph),
             'kind': kind} for graph in graphs]
    return rows, counts


def summary():
    rows, counts = one_particle_irreducible_table()
    sigma_rows, two_point_counts = self_energy_table('sigma')
    plus_rows, _ = self_energy_table('pi+')
    minus_rows, _ = self_energy_table('pi-')
    return {'one_particle_irreducible': rows, 'surjections_three': counts,
            'surjections_two': two_point_counts, 'sigma': sigma_rows, 'pi+': plus_rows,
            'pi-': minus_rows}


def kept(rows):
    """Graphs that contribute to the estimate (nonzero termwise routing maximum)."""
    return [row for row in rows if row['status'] in ('kept', 'kept (D4 term excluded)')]


ALLOWED_STATUSES = ('kept', 'kept (D4 term excluded)', 'vanishes identically',
                    'no residue correction (exact: dZ = 0)', 'kept, no contribution (exact zero)')


def no_series_drop(data):
    """Every graph is kept or vanishes exactly; none is dropped by a series."""
    return all(row['status'] in ALLOWED_STATUSES
               for group in ('one_particle_irreducible', 'sigma', 'pi+', 'pi-') for row in data[group])


def possible_scales(row, piece_function=line_pieces):
    """Distinct m_term of the terms of a graph that contribute (stored by
    classify for kept graphs); for other graphs, those of all piece choices
    (D4 term excluded where removed)."""
    if row.get('scales'):
        return row['scales']
    piece_lists = [piece_function(line, f'G{i + 1}') for i, line in enumerate(row['lines'])]
    scales = set()
    for choice in itertools.product(*piece_lists):
        if row.get('pieces') == 'rest' and is_bare_choice(choice):
            continue
        scales.add(term_scale(choice))
    return sorted(scales, key=sp.default_sort_key)


def routing_summary(row):
    """Counts of the routing decisions of a kept graph."""
    maximum = row['maximum']
    return {**maximum['counts'], 'marked': len(maximum['marked']),
            'marked boundary': len(maximum['marked boundary']),
            'common routing': common_routing(row['routing_values'], row['value'])}


def main():
    data = summary()
    counts_three, counts_two = data['surjections_three'], data['surjections_two']
    print('Surjections of {sigma, pi+, pi-} onto V cycle vertices:', counts_three)
    print('Surjections of the two legs of a self-energy onto V vertices:', counts_two)
    assert counts_three[4] == counts_three[5] == 0, 'a box would need four vertices with an external leg each'
    assert counts_two[3] == 0
    print('  No box: four vertices with one external leg each cannot be built from three legs.')
    print('  No (2,5) 1PI bubble: it would need a two-point vertex with no external leg; none exists.')
    assert no_series_drop(data)
    print('  Every graph is kept or vanishes exactly; no graph is dropped by a series in eps or 1/eps.')
    rows = data['one_particle_irreducible']
    print(f'\nOne-loop 1PI graphs: {len(rows)}')
    for kind in ('tadpole', 'bubble', 'triangle'):
        subset = [row for row in rows if topology(row) == kind]
        print(f'  {kind}: {len(subset)} graphs, kept {len(kept(subset))}, '
              f'vanishing {sum(1 for row in subset if row["status"] == "vanishes identically")}, '
              f'with an excluded D4 term {sum(1 for row in kept(subset) if row["excluded"] is not None)}')
    header = ('label | S_Gamma | r | lines | m_term of the terms | status | all-generator-sum term | '
              'routing: independent/numerical/mass-dependent products, marked, common routing')
    print('\n' + header)
    all_groups = ('one_particle_irreducible', 'sigma', 'pi+', 'pi-')
    for group in all_groups:
        if group != 'one_particle_irreducible':
            print(f'\nSelf-energy graphs of {group}:')
        for row in data[group]:
            routing_text = '-'
            if row in kept([row]):
                summary_row = routing_summary(row)
                routing_text = (f"{summary_row['independent']}/{summary_row['numerical']}/"
                                f"{summary_row['mass-dependent']}, {summary_row['marked']}, "
                                f"{summary_row['common routing']}")
            print(f'  {row["label"]} | {row["symmetry"]} | {len(row["blocks"])} | {line_content(row)} | '
                  f'{possible_scales(row)} | {row["status"]} | {row["excluded_status"]} | '
                  f'{routing_text}')
    mass_dependent = sum(routing_summary(row)['mass-dependent'] for group in all_groups for row in kept(data[group]))
    marked = sum(routing_summary(row)['marked'] for group in all_groups for row in kept(data[group]))
    boundary = sum(routing_summary(row)['marked boundary'] for group in all_groups for row in kept(data[group]))
    print('\nProducts of traces whose largest routing depends on the masses:', mass_dependent,
          '; decided by the first routing (sum over routings zero): on open regions', marked,
          ', at isolated boundary points', boundary)
    total = sum(len(kept(data[group])) for group in all_groups)
    print('Kept graphs in total (1PI + sigma + pi+ + pi- self-energies):', total)


if __name__ == '__main__':
    main()
