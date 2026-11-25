"""Implicit loop-to-tree crossing with the radial mode, in phi-valued invariants.

Tasks 4.1 and 4.2 of the OpenSpec change phi-dependent-loop-tree-crossing-murayama.

Run:
    .venv/bin/python scripts/derive_murayama_radial_crossing.py
Tests:
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_murayama_radial_crossing

No loop integral, scalar loop function, loop coefficient or loop amplitude is
computed here.  Every loop size enters only through the size prescription of
memory/alp/shared/murayama-radial-loop-size-prescription.md:

    estimate = (epsilon / sqrt(N(m_a^2)))
               * sum over distinct scalar moments of
                 m_moment^2 / (16 pi^2 f_pi^2)
                 * | sum over the loop classes of that moment of
                     epsilon^p_k N(m_a^2)^(-n_k) c_k |,

with p_k the power of epsilon relative to the tree and n_k the number of
internal ALP-like pole normalizations of class k (counting record
memory/alp/shared/murayama-radial-loop-class-counting.md, whose class numbers
are used as dictionary keys below).  c_k is the reduced invariant loop
contraction of class k, as defined in
memory/alp/shared/alp-loop-crossing-exact-tree-condition.md: the
class-dependent contraction after removing the external factor, the internal
factors and the loop scale.  It depends on epsilon and on the dressed
phi-valued invariants S_hat_rs, P_hat_rs at the fixed Dalitz point, and is
kept as an undefined function: it is not evaluated.  The naturalness
assumption of memory/alp/shared/murayama-radial-naturalness-assumption.md is
what makes this magnitude an estimate of the next-order correction.

Conventions and units (Alda:2026nrd, radial-exchange record
memory/alp/three-pseudoscalar/alp-three-pseudoscalar-radial-exchange.md):
  * epsilon = f_pi/f_a (dimensionless, positive); f_pi, masses in GeV;
    squared masses and Dalitz invariants in GeV^2; amplitudes dimensionless.
  * Polar fields, dressed sources M_hat(a, sigma) = S_hat + i P_hat with
    S_hat = sum S_hat_rs (sigma/f_a)^r (a/f_a)^s (first index counts sigma).
  * Fixed independent Dalitz invariants s_12 = (p_1+p_2)^2 and
    s_23 = (p_2+p_3)^2, as in the ALP-only derivation; s_13 follows from
    the exact pole masses.
  * Fixed pole masses m_a and m_sigma; the bare ALP mass and the quartic
    coupling lambda are eliminated by the pole equation and the radial-mass
    relation (labelled choice: m_sigma, not lambda, is held fixed).

Assumption (scope fixed by the user): the dressed topological invariant
vanishes, Theta_hat = 0, so the tadpole source is C = B_0 P_hat_01 (the
Murayama model and its generalizations).  The ALP-only record allows a
nonzero Theta_hat_1; it is compared here at Theta_hat_1 = 0, where its C is
the same B_0 P_hat_01.

Scope of the numerator (decided by the user, not results):
  * Loops with meson lines only are excluded.  They carry no power of
    f_pi/f_a beyond the ordinary chiral correction.
  * Every loop class with at least one internal ALP-like or radial line is
    kept in full: no truncation in epsilon.  This includes the loop with two
    internal ALP-like propagators and the epsilon dependence of N(m_a^2) and
    of the reduced contractions.
  * The two tadpoles of the ALP-like pole at m_a (direct tadpole and the
    meson one-point source from a closed ALP-like line) share one scalar
    moment and are combined before the magnitude (prescription rule 4).
  * The radial tadpole and the two radial self-energies of ALP-like lines
    (below) all close one sigma line at one vertex: they share the sigma
    tadpole moment at m_sigma and are combined before the magnitude
    (prescription rule 4).  Every other class is a distinct moment
    (labelled choice).

Radial self-energy of ALP-like lines (added to the counting).  The polar
kinetic term (1/2)(1 + sigma/f_a)^2 (d a)^2 contains the vertex
(1/(2 f_a^2)) sigma^2 (d a)^2, with the two derivatives on the ALP legs.  A
closed sigma line on it, attached to an ALP-like line, is a self-energy of
that line:
  * on the external ALP-like leg (class 'R5a'): epsilon^2 relative to the
    tree, from 1/f_a^2 = epsilon^2/f_pi^2; no N(m_a^2)^(-1) from the loop
    line, because sigma has no meson admixture; scale m_sigma.  The
    a-component of the leg, 1/N(m_a^2) for the two ends of the vertex, is a
    factor of the leg, not of an internal ALP-like line, and is kept inside
    c_R5a (labelled choice, so that n_k keeps counting internal lines only).
  * on an internal ALP-like line of a tree exchange (class 'R5b'): the same
    epsilon^2 and scale, relative to the tree exchange that contains the
    line.  The suppression of that exchange relative to the leading tree
    term is a function of epsilon kept inside c_R5b (labelled choice).
  * The sigma tadpole on the sigma-sigma-P-P vertex of the potential (and
    on the non-derivative sigma-sigma-a-a entry S_hat_22) is the radial
    tadpole class R2; the new classes contain only the derivative vertex,
    so nothing is counted twice.
Assumption (same convention as the ALP-only record, where the external-leg
self-energies with a closed ALP-like line belong to class 2): the shift of
the ALP mass from this self-energy is absorbed into the physical pole mass
m_a, which is held fixed.  Only the residue part (derivative of the
self-energy with respect to p^2 at p^2 = m_a^2) and the momentum
dependence away from the pole enter the estimate.  The self-energy is
linear in p^2 (two derivatives on the ALP legs, the sigma tadpole does not
depend on p), so the part beyond the residue vanishes identically
(radial_alp_self_energy_parts).  The sigma tadpole itself is not evaluated.

ALP-only limit: RADIAL_COUPLING_SWITCH (printed t_sigma, value 1 in the
model) multiplies every coupling of the radial mode sigma: sigma-meson-meson,
sigma-ALP-meson, and the sigma couplings of the loops,
including the derivative sigma-sigma-a-a vertex.  Setting it to 0 turns
off the radial exchange in the tree and every loop with a radial line, at
finite m_sigma.  It is a switch on couplings, not a physical parameter.  The
limit m_sigma -> infinity is not used: the radial tadpole weight grows with
m_sigma^2 at fixed contraction.
"""
import sympy as sp

from derive_alp_loop_crossing_criterion import derive_crossing
from derive_murayama_radial_frame_invariance import (
    F_A as FRAME_F_A, F_PI as FRAME_F_PI, POLE_I, diagonal_frames,
    difference, estimate_inputs, is_zero)

EPSILON = sp.Symbol(r'\epsilon', positive=True)
F_PI = sp.Symbol(r'f_\pi', positive=True)
B_ZERO = sp.Symbol('B_0', positive=True)
MAXIMUM_RATIO = sp.Symbol(r'R_{\max}', positive=True)
ALP_POLE_MASS_SQUARED = sp.Symbol('m_a^2', positive=True)
RADIAL_MASS_SQUARED = sp.Symbol(r'm_\sigma^2', positive=True)
QUARTIC_COUPLING = sp.Symbol(r'\lambda', positive=True)
PAIR_INVARIANT_12 = sp.Symbol('s_{12}', positive=True)
PAIR_INVARIANT_23 = sp.Symbol('s_{23}', positive=True)
RADIAL_COUPLING_SWITCH = sp.Symbol('t_sigma', nonnegative=True)
PARENT_OVERLAP_TRACE = sp.Symbol(
    r'\langle[\mathbb{M}_{m_a^2}^{-1}(\mathcal{C})]^2\rangle', nonnegative=True)
DAUGHTERS = (1, 2, 3)

# Reference masses squared of the prescription record (one per scalar moment).
REFERENCE_MASS_SQUARED = {
    'alp_tadpole': ALP_POLE_MASS_SQUARED,
    'alp_bubble_triangle': sp.Symbol('m_{(3)}^2', positive=True),
    'two_alp_lines': sp.Symbol('m_{(4)}^2', positive=True),
    'meson_loop_alp_insertion': sp.Symbol('m_P^2', positive=True),
    'radial_tadpole': RADIAL_MASS_SQUARED,
    'radial_bubble_triangle': sp.Symbol('m_{(R3)}^2', positive=True),
    'radial_with_alp_line': sp.Symbol('m_{(R4a)}^2', positive=True),
    'two_radial_lines': sp.Symbol('m_{(R4b)}^2', positive=True),
}

# Counting record: class -> (scalar moment, power of epsilon relative to the
# tree, number of internal ALP-like normalizations 1/N(m_a^2), contains a
# radial line, physical description).
LOOP_CLASSES = {
    '2': ('alp_tadpole', 2, 1, False, 'tadpole of the ALP-like pole'),
    '6': ('alp_tadpole', 2, 1, False,
          'meson one-point source from a closed ALP-like line (vacuum shift)'),
    '3': ('alp_bubble_triangle', 2, 1, False,
          'bubble or triangle with one internal ALP-like propagator'),
    '4': ('two_alp_lines', 4, 2, False, 'loop with two internal ALP-like propagators'),
    '5': ('meson_loop_alp_insertion', 2, 1, False,
          'meson loop with an ALP-meson mixing insertion on an internal meson line'),
    'R2': ('radial_tadpole', 2, 0, True, 'tadpole of the radial mode'),
    'R3': ('radial_bubble_triangle', 2, 0, True,
           'bubble or triangle with one internal radial propagator'),
    'R4a': ('radial_with_alp_line', 4, 1, True,
            'loop with one internal radial and one internal ALP-like propagator'),
    'R4b': ('two_radial_lines', 4, 0, True, 'loop with two internal radial propagators'),
    'R5a': ('radial_tadpole', 2, 0, True,
            'radial self-energy of the external ALP-like leg (closed sigma line on the '
            'sigma-sigma-a-a kinetic vertex)'),
    'R5b': ('radial_tadpole', 2, 0, True,
            'radial self-energy of an internal ALP-like line of a tree exchange (closed '
            'sigma line on the sigma-sigma-a-a kinetic vertex)'),
}
# Classes with a closed sigma line on the derivative sigma-sigma-a-a vertex.
RADIAL_ALP_SELF_ENERGY_CLASSES = ('R5a', 'R5b')
ALP_ONLY_CLASSES = tuple(name for name, data in LOOP_CLASSES.items() if not data[3])
# Not part of the numerator (user decision): ordinary chiral correction.
EXCLUDED_LOOP_CLASSES = {'1': 'loop with meson lines only (ordinary chiral correction)'}


def reduced_contraction(loop_class):
    """Undefined reduced invariant contraction c_k of one loop class.

    Classes without a radial line depend on epsilon only.  A class with a
    radial line contains at least one radial coupling; one power of the
    switch is written explicitly in loop_size_estimate and the remaining
    dependence on the switch is kept inside c_k(epsilon, t_sigma)."""
    function = sp.Function(f'reduced_contraction_class_{loop_class}', real=True)
    if LOOP_CLASSES.get(loop_class, (None, None, None, False))[3]:
        return function(EPSILON, RADIAL_COUPLING_SWITCH)
    return function(EPSILON)


def alp_pole_normalization():
    """N(m_a^2) = 1 + 2 epsilon^2 <[M^{-1}_{m_a^2}(C)]^2>, fixed m_a, C = B_0 P_hat_01."""
    return 1 + 2 * EPSILON**2 * PARENT_OVERLAP_TRACE


def loop_factor(moment):
    """m_moment^2 / (16 pi^2 f_pi^2), dimensionless."""
    return REFERENCE_MASS_SQUARED[moment] / (16 * sp.pi**2 * F_PI**2)


def pull_positive_factors(expression):
    """Write each |sum| as a common factor times |remaining sum|, so that the
    positive factors epsilon^p and N^-n leave the absolute value."""
    return expression.replace(
        sp.Abs, lambda argument: sp.Abs(sp.factor_terms(sp.together(argument))))


def assemble_estimate(contraction_of_class, classes, epsilon, normalization,
                      loop_factor_of_moment, switch=1, absolute_value=sp.Abs):
    """Prescription of the size record, for given contractions of the classes.

    The same assembly is used for the symbolic estimate (task 4.1) and for the
    explicit diagonal-sector contractions of the frame check (task 4.2)."""
    combined = {}
    for loop_class in classes:
        moment, power, alp_lines, radial, _ = LOOP_CLASSES[loop_class]
        combined[moment] = combined.get(moment, 0) + (
            (switch if radial else 1) * epsilon**power * normalization**(-alp_lines)
            * contraction_of_class[loop_class])
    external_factor = epsilon / sp.sqrt(normalization)
    return external_factor * sum(loop_factor_of_moment(moment) * absolute_value(term)
                                 for moment, term in combined.items())


def loop_size_estimate(classes=tuple(LOOP_CLASSES)):
    """Estimated magnitude of the next-order correction, full in epsilon."""
    return pull_positive_factors(assemble_estimate(
        {name: reduced_contraction(name) for name in classes}, classes, EPSILON,
        alp_pole_normalization(), loop_factor, switch=RADIAL_COUPLING_SWITCH))


def radial_alp_self_energy_parts(momentum_squared, alp_component_squared, pole_mass_squared):
    """Structure of the radial self-energy of an ALP-like line.

    Per unit of (vertex coefficient 1/(2 f_a^2), written as epsilon^2/f_pi^2
    in the estimate) times (sigma tadpole), neither of which is evaluated:
    the two derivatives of sigma^2 (d a)^2 act on the two ALP ends and give
    the squared momentum p^2 of the line; each end carries the a-component
    of the pole, hence alp_component_squared.  O(1) combinatorial factors are
    not fixed here; they belong to the undefined reduced contraction.

    Returns a dict with the self-energy structure Sigma(p^2), the mass shift
    Sigma(m^2) (absorbed into the physical pole mass), the residue part
    dSigma/dp^2 at p^2 = m^2 (enters the estimate), and the remainder
    Sigma(p^2) - Sigma(m^2) - (p^2 - m^2) dSigma/dp^2(m^2), which is the
    momentum dependence beyond the residue.  On an internal line,
    (Sigma(q^2) - Sigma(m^2)) / (q^2 - m^2) cancels one propagator and equals
    the residue part plus remainder/(q^2 - m^2)."""
    self_energy = alp_component_squared * momentum_squared
    mass_shift = self_energy.subs(momentum_squared, pole_mass_squared)
    residue = sp.diff(self_energy, momentum_squared).subs(momentum_squared, pole_mass_squared)
    remainder = sp.expand(self_energy - mass_shift - (momentum_squared - pole_mass_squared) * residue)
    return {'self-energy': self_energy, 'mass shift (absorbed into m_a)': mass_shift,
            'residue part': residue, 'remainder beyond the residue': remainder}


# ---------------------------------------------------------------------------
# Exact tree with radial exchange (radial-exchange record, boxed result).
# ---------------------------------------------------------------------------

def daughter_pole_mass_squared(daughter):
    return sp.Function(f'm_{daughter}^2', positive=True)(EPSILON)


def leg_overlap_trace(leg):
    """<[M^{-1}_{m_x^2}(P_hat_01)]^2>; for the parent B_0^2 times it equals the C trace."""
    if leg == 'a':
        return PARENT_OVERLAP_TRACE / B_ZERO**2
    return sp.Function(
        rf'\langle[\mathbb{{M}}_{{m_{leg}^2}}^{{-1}}(\widehat P_{{01}})]^2\rangle')(EPSILON)


def leg_pole_mass_squared(leg):
    return ALP_POLE_MASS_SQUARED if leg == 'a' else daughter_pole_mass_squared(leg)


def leg_normalization(leg):
    return 1 + 2 * B_ZERO**2 * EPSILON**2 * leg_overlap_trace(leg)


def anticommutator_trace(first, second):
    """<S_hat_10 {M^{-1}_{m_x^2}(P_hat_01), M^{-1}_{m_y^2}(P_hat_01)}>.

    Every pair used contains a daughter, whose pole mass depends on epsilon."""
    label = (rf'\langle\widehat S_{{10}}\{{\mathbb M^{{-1}}_{{m_{first}^2}}(\widehat P_{{01}}),'
             rf'\mathbb M^{{-1}}_{{m_{second}^2}}(\widehat P_{{01}})\}}\rangle')
    return sp.Function(label)(EPSILON)


def single_trace(leg):
    """<P_hat_11 M^{-1}_{m_x^2}(P_hat_01)>."""
    name = rf'\langle\widehat P_{{11}}\mathbb M^{{-1}}_{{m_{leg}^2}}(\widehat P_{{01}})\rangle'
    return sp.Symbol(name, real=True) if leg == 'a' else sp.Function(name)(EPSILON)


def radial_coupling(first, second, momentum_squared):
    """t_sigma g_xy(q^2): radial-exchange record coupling, with f_pi/f_a = epsilon,
    multiplied by the switch on the radial couplings."""
    trace_s12 = sp.Symbol(r'\langle\widehat S_{12}\rangle', real=True)
    numerator = (
        2 * B_ZERO**3 * EPSILON**2 * anticommutator_trace(first, second)
        - 2 * B_ZERO**2 * EPSILON**2 * (single_trace(first) + single_trace(second))
        - 2 * B_ZERO * EPSILON**2 * trace_s12
        + momentum_squared - leg_pole_mass_squared(first) - leg_pole_mass_squared(second))
    f_a = F_PI / EPSILON
    return RADIAL_COUPLING_SWITCH * numerator / (
        f_a * sp.sqrt(leg_normalization(first) * leg_normalization(second)))


def pair_invariants():
    """s_12, s_23 fixed; s_13 from the exact pole masses."""
    mass_sum = ALP_POLE_MASS_SQUARED + sum(daughter_pole_mass_squared(n) for n in DAUGHTERS)
    return {(1, 2): PAIR_INVARIANT_12, (2, 3): PAIR_INVARIANT_23,
            (1, 3): mass_sum - PAIR_INVARIANT_12 - PAIR_INVARIANT_23}


def tree_without_radial_exchange():
    """Exact pole-projected tree of memory/alp/shared/alp-finite-mixing-three-meson-tree.md."""
    return sp.Function(r'\mathcal{M}^{(0)}_{\mathrm{no}\,\sigma}')(
        EPSILON, PAIR_INVARIANT_12, PAIR_INVARIANT_23)


def radial_exchange():
    """sum over channels (r|uv) of g_ar(s_uv) g_uv(s_uv) / (m_sigma^2 - s_uv)."""
    invariants = pair_invariants()
    total = 0
    for spectator in DAUGHTERS:
        pair = tuple(n for n in DAUGHTERS if n != spectator)
        momentum_squared = invariants[pair]
        total += (radial_coupling('a', spectator, momentum_squared)
                  * radial_coupling(*pair, momentum_squared)
                  / (RADIAL_MASS_SQUARED - momentum_squared))
    return total


def exact_tree():
    return tree_without_radial_exchange() + radial_exchange()


def quartic_coupling_at_fixed_radial_mass():
    """lambda from m_sigma^2 = 2 lambda f_a^2 + B_0 epsilon^2 (<S_10> - 2 <S_20>)."""
    trace_s10 = sp.Symbol(r'\langle\widehat S_{10}\rangle', real=True)
    trace_s20 = sp.Symbol(r'\langle\widehat S_{20}\rangle', real=True)
    relation = sp.Eq(RADIAL_MASS_SQUARED, 2 * QUARTIC_COUPLING * F_PI**2 / EPSILON**2
                     + B_ZERO * EPSILON**2 * (trace_s10 - 2 * trace_s20))
    return relation, sp.solve(relation, QUARTIC_COUPLING)[0]


# ---------------------------------------------------------------------------
# Task 4.1: implicit crossing and its ALP-only limit.
# ---------------------------------------------------------------------------

def crossing_equation(classes=tuple(LOOP_CLASSES)):
    """F = estimate - R_max |M^(0)| = 0, and the ratio R."""
    numerator = loop_size_estimate(classes)
    denominator = sp.Abs(exact_tree())
    return numerator - MAXIMUM_RATIO * denominator, numerator / denominator


def leading_numerator_truncation(numerator):
    """Truncation of the ALP-only record: epsilon^3 lim_{epsilon'->0} numerator/epsilon'^3.

    Positive powers of epsilon are moved out of each magnitude; the remaining
    argument is assumed continuous at epsilon = 0 (smooth branch of the
    ALP-only record), so the limit is taken by substitution."""
    return EPSILON**3 * sp.expand(pull_positive_factors(numerator) / EPSILON**3).subs(EPSILON, 0)


def alp_only_limit(expression):
    """Radial couplings switched off, m_sigma finite."""
    return expression.subs(RADIAL_COUPLING_SWITCH, 0)


def alp_only_record_equation():
    """Accepted equation of alp-loop-crossing-exact-tree-condition.md at Theta_hat_1 = 0:
    epsilon^3 lim |M^(1)|/epsilon^3 - R_max |M^(0)|, with the leading ALP-only classes."""
    leading = (loop_factor('alp_tadpole') * sp.Abs(
                   reduced_contraction('2').subs(EPSILON, 0) + reduced_contraction('6').subs(EPSILON, 0))
               + loop_factor('alp_bubble_triangle') * sp.Abs(reduced_contraction('3').subs(EPSILON, 0))
               + loop_factor('meson_loop_alp_insertion') * sp.Abs(reduced_contraction('5').subs(EPSILON, 0)))
    return EPSILON**3 * leading - MAXIMUM_RATIO * sp.Abs(tree_without_radial_exchange())


def alp_only_difference_terms():
    """Terms of (ALP-only limit of the full equation) - (record equation).

    Returns a list of (loop, origin, expression).  For a moment present in the
    record, with combined contraction c(epsilon) and n = 1 internal ALP-like line,
        epsilon^3 L |c(epsilon)| N^(-3/2) - epsilon^3 L |c(0)|
          = epsilon^3 L (N^(-3/2) - 1) |c(epsilon)|          [pole normalizations]
          + epsilon^3 L (|c(epsilon)| - |c(0)|)              [epsilon dependence of c],
    an exact split.  The loop with two internal ALP-like propagators is absent
    from the record and enters as a whole."""
    normalization = alp_pole_normalization()
    groups = {}
    for name in ALP_ONLY_CLASSES:
        groups.setdefault(LOOP_CLASSES[name][0], []).append(name)
    terms = []
    for moment, names in groups.items():
        _, power, alp_lines, _, _ = LOOP_CLASSES[names[0]]
        loop = ' plus '.join(LOOP_CLASSES[name][4] for name in names)
        contraction = sum(reduced_contraction(name) for name in names)
        prefactor = EPSILON**(power + 1) * loop_factor(moment)
        full_magnitude = prefactor * normalization**(-alp_lines - sp.Rational(1, 2)) * sp.Abs(contraction)
        if power > 2:
            terms.append((loop, 'whole contribution, absent from the record', full_magnitude))
            continue
        terms.append((loop, 'pole normalizations N(m_a^2) of the external and internal ALP-like lines',
                      prefactor * (normalization**(-alp_lines - sp.Rational(1, 2)) - 1)
                      * sp.Abs(contraction)))
        terms.append((loop, 'dependence of the reduced contraction on epsilon at finite mixing',
                      prefactor * (sp.Abs(contraction) - sp.Abs(contraction.subs(EPSILON, 0)))))
    return terms


def leading_relative_power(term):
    """Lowest power of epsilon of a difference term, relative to the leading tree
    term (which is linear in epsilon).

    Used only to label the term.  Each reduced contraction is replaced by
    c_0 + c_2 epsilon^2 with positive symbols: smooth branch with a nonzero
    leading coefficient and dependence on epsilon^2, as in the ALP-only
    record.  Positive signs only fix the branch of |.|; they do not affect the
    power unless c_2 = 0, in which case the power is higher."""
    replacements = {}
    for function in term.atoms(sp.Function):
        if function.func.__name__.startswith('reduced_contraction_class_'):
            label = function.func.__name__.split('_')[-1]
            leading = sp.Symbol(f'c_{label}_0', positive=True)
            correction = sp.Symbol(f'c_{label}_2', positive=True) * EPSILON**2
            replacements[function] = leading + correction.subs(EPSILON, function.args[0])
    smooth = term.xreplace(replacements)
    expansion = sp.expand(sp.series(smooth, EPSILON, 0, 10).removeO())
    return min(sp.Poly(expansion, EPSILON).monoms())[0] - 1


def alp_only_check():
    """Residual after the record truncation, and residual of the term split.

    Returns (truncation residual, full difference, split residual, terms)."""
    equation, _ = crossing_equation()
    limit_equation = alp_only_limit(equation)
    record = alp_only_record_equation()
    estimate = alp_only_limit(loop_size_estimate())
    truncated_equation = limit_equation - estimate + leading_numerator_truncation(estimate)
    truncation_residual = sp.simplify(truncated_equation - record)
    full_difference = pull_positive_factors(limit_equation - record)
    terms = alp_only_difference_terms()
    split_residual = sp.simplify(full_difference - sum(term for _, _, term in terms))
    return truncation_residual, full_difference, split_residual, terms


def single_line_ratio_comparison():
    """Only the tadpole of the ALP-like pole, radial couplings off, no truncation,
    against the ratio of derive_alp_loop_crossing_criterion.py with its
    coefficient identified."""
    _, ratio = crossing_equation(('2',))
    ratio = alp_only_limit(ratio)
    (_, _, _, _, _, criterion_ratio, criterion_epsilon, _, _, criterion_overlap,
     *_) = derive_crossing()
    reduced_tree = sp.sqrt(alp_pole_normalization()) / EPSILON * tree_without_radial_exchange()
    coefficient = loop_factor('alp_tadpole') * sp.Abs(reduced_contraction('2')) / sp.Abs(reduced_tree)
    criterion_coefficient = sp.Function('dimensionless_loop_coefficient')(criterion_epsilon**2)
    mapped = criterion_ratio.subs(criterion_coefficient, coefficient).subs(
        {criterion_epsilon: EPSILON, criterion_overlap: PARENT_OVERLAP_TRACE})
    return sp.simplify(ratio - mapped)


# ---------------------------------------------------------------------------
# Task 4.2: frame invariance of R on the estimate's contractions, Theta_hat = 0.
# ---------------------------------------------------------------------------

def ratio_from_estimate_inputs(inputs):
    """R built with assemble_estimate from the diagonal-sector entries of the
    frame-invariance script (task 3.4), radial couplings on.

    Each contraction is a combination, with unspecified frame-independent
    coefficients, of the invariant entries that the counting record assigns to
    the class; the loops with two internal lines use products of the entries
    of the one-line loops.  The tree without radial exchange is a
    frame-independent symbol (its invariance is the ALP-only result); the
    radial exchange uses g_ij and m_sigma^2.  POLE_I plays the role of m_a^2.

    The radial self-energies of ALP-like lines use the residue part of
    radial_alp_self_energy_parts with p^2 = m_i^2 and a-component squared
    1/N(m_i^2); for the internal line it multiplies the ALP-like tree
    exchange, a frame-independent symbol for the same reason as the tree
    without radial exchange."""
    epsilon = FRAME_F_PI / FRAME_F_A
    coefficient = sp.symbols('unspecified_coefficient_1:17', real=True)
    alp_momentum_squared = sp.Symbol('p_a^2', positive=True)
    self_energy_residue = radial_alp_self_energy_parts(
        alp_momentum_squared, 1 / inputs['N(m_i^2)'], POLE_I)['residue part']
    alp_like_tree_exchange = sp.Symbol('alp_like_tree_exchange', real=True)
    alp_pole_term = inputs['ALP pole term tr(C M^{-1}(C)) at m_i^2']
    trace_s02 = inputs['tr S_hat_02']
    radial_alp_meson = inputs['sigma-a-P contraction tr(P_hat_11 Omega_i)']
    radial_meson_meson = FRAME_F_A * inputs['sigma-PP coupling g_ij']
    contractions = {
        '2': coefficient[0] * alp_pole_term + coefficient[1] * trace_s02,
        '6': coefficient[2] * alp_pole_term,
        '3': coefficient[3] * alp_pole_term,
        '4': coefficient[4] * alp_pole_term**2,
        '5': coefficient[5] * trace_s02,
        'R2': (coefficient[6] * inputs['R2: tr(S_hat_20 {Omega_i, Omega_j})']
               + coefficient[7] * inputs['R2: tr(P_hat_21 (Omega_i + Omega_j))']
               + coefficient[8] * inputs['R2: tr S_hat_22']
               + coefficient[9] * inputs['R2: tr S_hat_20 (sigma^2 a^0 vertex)']),
        'R3': coefficient[10] * radial_alp_meson + coefficient[11] * radial_meson_meson,
        'R4a': coefficient[12] * radial_alp_meson * alp_pole_term,
        'R4b': coefficient[13] * radial_meson_meson**2,
        'R5a': coefficient[14] * self_energy_residue,
        'R5b': coefficient[15] * self_energy_residue * alp_like_tree_exchange,
    }
    scales = {'alp_tadpole': POLE_I, 'radial_tadpole': inputs['m_sigma^2']}

    def frame_loop_factor(moment):
        scale = scales.get(moment, sp.Symbol(f'reference_mass_squared_{moment}', positive=True))
        return scale / (16 * sp.pi**2 * FRAME_F_PI**2)

    numerator = assemble_estimate(
        contractions, tuple(LOOP_CLASSES), epsilon, inputs['N(m_i^2)'], frame_loop_factor,
        absolute_value=lambda term: sp.Abs(term, evaluate=False))
    tree_no_radial = sp.Symbol('tree_without_radial_exchange', real=True)
    pair_invariant = sp.Symbol('s_pair', positive=True)
    tree = tree_no_radial + 3 * inputs['sigma-PP coupling g_ij']**2 / (inputs['m_sigma^2'] - pair_invariant)
    return numerator / sp.Abs(tree, evaluate=False)


def frames_without_topological_invariant():
    """Diagonal-sector frames of task 3.4 with Theta_hat_01 = 0 in the original
    frame.  The dressed coefficient stays 0 after the redefinition."""
    frames = diagonal_frames(alp_order=2, radial_order=2)
    theta_one = sp.Symbol(r'\Theta_{01}', real=True)
    frames['theta'] = {name: sp.sympify(value).subs(theta_one, 0)
                       for name, value in frames['theta'].items()}
    return frames


def frame_invariance_check():
    """R(dressed transformed frame) - R(original), both from dressed invariants."""
    frames = frames_without_topological_invariant()
    ratios = {name: ratio_from_estimate_inputs(estimate_inputs(frames[name], frames['theta'][name]))
              for name in ('original', 'dressed_after')}
    dressed_difference = difference(ratios['dressed_after'], ratios['original'])
    return dressed_difference, is_zero(dressed_difference), frames['theta']


def main():
    print('Assumption: Theta_hat = 0, C = B_0 P_hat_01.')
    print('Excluded from the numerator:', '; '.join(EXCLUDED_LOOP_CLASSES.values()))
    print('Included in full (no truncation in epsilon):')
    for name, (moment, power, alp_lines, _, description) in LOOP_CLASSES.items():
        print(f'    {description}: epsilon^{power} relative to the tree, '
              f'N(m_a^2)^-{alp_lines} from internal ALP-like lines, scale {REFERENCE_MASS_SQUARED[moment]}')
    print('Radial tadpole and radial self-energies of ALP-like lines share the sigma tadpole moment')
    print('at m_sigma and are combined inside one magnitude (prescription rule 4).')
    momentum_squared = sp.Symbol('p_a^2', positive=True)
    parts = radial_alp_self_energy_parts(momentum_squared, 1 / alp_pole_normalization(),
                                         ALP_POLE_MASS_SQUARED)
    print('Radial self-energy of an ALP-like line, from (1/(2 f_a^2)) sigma^2 (d a)^2 with a closed')
    print('sigma line, per unit of vertex coefficient and sigma tadpole (neither evaluated):')
    for name, value in parts.items():
        print(f'    {name}: {sp.latex(value)}')
    print('Assumption (ALP-only convention): the mass shift is absorbed into the physical m_a;')
    print('only the residue part and the momentum dependence beyond it enter (the latter is zero).')
    print('\nTask 4.1: implicit crossing  estimate = R_max |M_no_sigma + radial exchange|')
    print('Estimate (numerator):')
    print(sp.latex(loop_size_estimate()), '\n')
    print('Radial exchange, channel (3|12): g_a3(s_12) g_12(s_12) / (m_sigma^2 - s_12) with')
    print('t_sigma g_12(s_12) =', sp.latex(radial_coupling(1, 2, PAIR_INVARIANT_12)), '\n')
    relation, quartic = quartic_coupling_at_fixed_radial_mass()
    print('Radial mass held fixed; quartic coupling eliminated:')
    print(sp.latex(sp.Eq(QUARTIC_COUPLING, quartic)), '\n')

    truncation_residual, full_difference, split_residual, terms = alp_only_check()
    print('ALP-only limit: radial couplings switched off (t_sigma = 0), m_sigma finite.')
    print('Full ALP-only equation minus the ALP-only record equation, term by term')
    print('(powers of epsilon relative to the leading tree term, which is linear in epsilon):')
    for loop, origin, term in terms:
        print(f'  * {loop} -- {origin}: leading power epsilon^{leading_relative_power(term)}, '
              f'removed by the record truncation: {leading_numerator_truncation(term) == 0}')
        print('      ', sp.latex(term))
    print('Difference minus the sum of these terms:', split_residual)
    print('After the record truncation of the numerator, minus the record equation:',
          truncation_residual, '\n')
    print('Tadpole of the ALP-like pole alone versus derive_alp_loop_crossing_criterion ratio:',
          single_line_ratio_comparison(), '\n')
    print('Large m_sigma at fixed contractions: radial exchange term ->',
          sp.limit(radial_exchange(), RADIAL_MASS_SQUARED, sp.oo),
          '; radial tadpole weight ->', sp.limit(loop_factor('radial_tadpole'), RADIAL_MASS_SQUARED, sp.oo))
    print('(the radial tadpole does not decouple: the ALP-only limit is t_sigma = 0, not m_sigma -> oo)\n')

    dressed_difference, dressed_zero, theta = frame_invariance_check()
    print('Task 4.2 (Theta_hat = 0): Theta_01 in the frames:', theta)
    print('R(dressed transformed frame) - R(original frame) =', dressed_difference,
          '; zero:', dressed_zero)


if __name__ == '__main__':
    main()
