"""Cartesian-frame vertices and flavour-trace tools for the size estimate of the one-loop sigma -> pi+ pi- amplitude.

This module is a library; it prints nothing.  It is used by the scripts
derive_sigma_pipi_cartesian_loop_*.py (OpenSpec change
estimate-sigma-pipi-one-loop-murayama-cartesian).

Conventions (Alda:2026nrd, DiLuzio:2026poh):

* Natural units, GeV.  epsilon = f_pi / f_a.  N_f = 3, <lambda_i lambda_j> = 2 delta_ij,
  lambda_0 = sqrt(2/3) 1.
* Cartesian fields phi = f_a + s + i p.  Lagrangian
      L = f_pi^2/4 <u_mu u^mu> + f_pi^2/4 <chi_+> - m_1^2/12 <Phi>^2
          + 1/2 (d s)^2 + 1/2 (d p)^2 - V_phi(s, p),
      chi = 2 B_0 M_hat(s, p),
      M_hat(s, p) = sum_{r,q} (S^C_{r,q} + i P^C_{r,q}) (s/f_a)^r (p/f_a)^q,
      V_phi = m_phi^2/2 |phi|^2 + lambda/4 |phi|^4,  |phi|^2 = (f_a + s)^2 + p^2,
  with U = exp(i Phi/f_pi).  S^C_{r,q}, P^C_{r,q} are Hermitian 3x3 matrices of
  mass dimension one, generic and noncommuting.  (In the memory records the
  second index is written s, as in the OpenSpec change.)
* A vertex V(x_1, ..., x_n) is the n-th derivative of L with respect to the
  fields of the legs; d_mu -> -i k_mu with all momenta incoming, so M = V for
  a contact graph.  A leg is a radial leg (field s) or a pseudoscalar leg
  carrying a meson matrix lambda_x and an ALP-field component xi_x (the
  coefficient of the Cartesian ALP field p).

Letters (strings inside flavour traces):
  'SC{r}{q}', 'PC{r}{q}'  Cartesian sources S^C_{r,q}, P^C_{r,q}
                          (never confused with the polar letters 'S{r}{s}', 'P{r}{s}');
  'l+', 'l-'              lambda_{pi+}, lambda_{pi-};
  'MinvC'                 M^{-1}_{m_a^2}(C) with C = B_0 P^C_{0,1} (Cartesian pole data);
  'Pu', 'Pd'              up and down projectors (up-only limit);
  'G...'                  internal U(3) generators summed by completeness.

The generic trace tools (TraceSum, CompactSum, canonical words, completeness,
charge rule, trees) are imported from sigma_pipi_loop_algebra.py.  The vertex
functions of that module belong to the polar frame and are not used here:
every Cartesian vertex is written anew below.

The module never evaluates a loop integral.
"""
import itertools
import math
from dataclasses import dataclass, field

import sympy as sp

import sigma_pipi_loop_algebra as generic
from sigma_pipi_loop_algebra import (B0, EPSILON, F_PI, M_SIGMA, QUARTIC, TraceSum, CompactSum,
                                     canonical_word, letter, symbolic_dot, make_dot, tree_key)

SINGLET_MASS_SQUARED = sp.Symbol('m_1^2', positive=True)       # anomalous singlet mass m_1^2
RADIAL_MASS_PARAMETER = sp.Symbol('m_phi^2', real=True)          # m_phi^2 of V_phi
SQUARED_MOMENTUM = sp.Symbol('s_k', real=True)                   # k^2 of a two-point function
UP_QUARK_MASS = sp.Symbol('m_u', positive=True)
ALP_MASS_COEFFICIENT = sp.Symbol('A', positive=True)             # A = B_0(<S^C_10> - 2<S^C_02>)
KAPPA_SIGMA = sp.Symbol('kappa_sigma', positive=True)            # kappa_sigma = B_0(<S^C_10> - 2<S^C_20>)
NORMALIZATION = generic.NORMALIZATION                            # N(m_a^2)

MINV_LETTER = 'MinvC'


# ---------------------------------------------------------------------------
# Letters
# ---------------------------------------------------------------------------

def source_letter(kind, radial, alp):
    """Letter of S^C_{radial,alp} (kind 'S') or P^C_{radial,alp} (kind 'P')."""
    return f'{kind}C{radial}{alp}'


def is_cartesian_source(name):
    return len(name) >= 4 and name[0] in 'SP' and name[1] == 'C' and name[2:].isdigit()


def source_indices(name):
    """(kind, r, q) of a Cartesian source letter."""
    return name[0], int(name[2]), int(name[3:])


def source_count(key):
    return sum(1 for word in key for name in word if is_cartesian_source(name))


# Assumptions 1 and 2 of the Cartesian frame (derived in
# derive_sigma_pipi_cartesian_loop_tree.py): P^C_00 = P^C_10 = 0 as matrix
# identities, <S^C_01> = <S^C_11> = 0 as trace identities.
VANISHING_LETTERS = ('PC00', 'PC10')
VANISHING_SINGLE_TRACES = (('SC01',), ('SC11',))


def apply_assumptions(trace_sum):
    """Remove the terms that vanish by P^C_00 = P^C_10 = 0, <S^C_01> = <S^C_11> = 0."""
    return trace_sum.filtered(
        lambda key: not any(name in VANISHING_LETTERS for word in key for name in word)
        and not any(word in VANISHING_SINGLE_TRACES for word in key))


def apply_charge_rule(trace_sum):
    """A trace of ad_Q eigen-matrices with nonzero total charge vanishes (generic rule)."""
    return generic.apply_charge_rule(trace_sum)


def compact_assumptions(compact_sum):
    """Assumptions 1 and 2 and the charge rule on compact traces."""
    def keep(key):
        for tree in key:
            names = generic.tree_letters(tree)
            if any(name in VANISHING_LETTERS for name in names):
                return False
            if tree in (('L', 'SC01'), ('L', 'SC11')):
                return False
            charges = [generic.letter_charge(name) for name in names]
            if None not in charges and sum(charges) != 0:
                return False
        return True
    return compact_sum.filtered(keep)


# ---------------------------------------------------------------------------
# Legs
# ---------------------------------------------------------------------------

@dataclass
class CartesianLeg:
    """A vertex leg in the Cartesian frame.

    radial: a leg of the radial field s (no meson matrix, no ALP-field component).
    meson: tuple of (coefficient, letter), the meson matrix of the leg.
    alp_field: coefficient of the Cartesian ALP field p in the leg.
    momentum: {basis label: coefficient}, incoming."""
    radial: bool = False
    meson: tuple = ()
    alp_field: object = 0
    momentum: dict = field(default_factory=dict)


def radial_leg(momentum):
    return CartesianLeg(radial=True, momentum=dict(momentum))


def meson_leg(name, momentum, coefficient=1, alp_field=0):
    return CartesianLeg(meson=((coefficient, name),), alp_field=alp_field, momentum=dict(momentum))


def alp_field_leg(momentum, coefficient=sp.Integer(1)):
    """A leg of the Cartesian ALP field p alone."""
    return CartesianLeg(alp_field=coefficient, momentum=dict(momentum))


def alp_like_leg(momentum):
    """The ALP-like pole vector v_a = (2 eps M^{-1}_{m_a^2}(C)/sqrt(N), 1/sqrt(N)) on (mesons + p)."""
    return CartesianLeg(meson=((2 * EPSILON / sp.sqrt(NORMALIZATION), MINV_LETTER),),
                        alp_field=1 / sp.sqrt(NORMALIZATION), momentum=dict(momentum))


def generic_legs(number_of_legs, radial):
    """radial legs first (momenta q1, ...), then generic pseudoscalar legs
    (x_j, xi_j) carrying both a meson matrix and an ALP-field component."""
    legs = []
    for index in range(number_of_legs):
        label = f'q{index + 1}'
        if index < radial:
            legs.append(radial_leg({label: 1}))
        else:
            legs.append(meson_leg(f'x{index + 1}', {label: 1},
                                  alp_field=sp.Symbol(f'xi{index + 1}', real=True)))
    return legs


# ---------------------------------------------------------------------------
# Construction 2: closed vertex rules
# ---------------------------------------------------------------------------

def _leg_roles(pseudoscalar_legs):
    options = []
    for leg in pseudoscalar_legs:
        roles = []
        if leg.alp_field != 0:
            roles.append(('alp', leg.alp_field, None))
        for coefficient, name in leg.meson:
            roles.append(('meson', coefficient, name))
        options.append(roles)
    return options


def source_vertex_coefficient(number_of_legs, radial, alp, meson):
    """B_0 f_pi^{2-n} r! q! eps^{r+q} (-1)^{m(m-1)/2}/m!  (r radial ends, q ALP-field ends, m meson ends)."""
    return (B0 * F_PI ** (2 - number_of_legs) * math.factorial(radial) * math.factorial(alp)
            * EPSILON ** (radial + alp) * sp.Integer(-1) ** (meson * (meson - 1) // 2) / math.factorial(meson))


def source_vertex(legs):
    """Source part (f_pi^2/4)<chi_+> of the vertex, word level:
    sum over the choices (R, A, P) of B_0 f_pi^{2-n} r! q! eps^{r+q} (-1)^{m(m-1)/2}/m!
    prod_{A} xi  sum_{perms} <X^C_{r,q} lambda_{pi_1} ... lambda_{pi_m}>, X = S (m even), P (m odd)."""
    n = len(legs)
    radial = sum(1 for leg in legs if leg.radial)
    pseudoscalar = [leg for leg in legs if not leg.radial]
    result = TraceSum()
    for choice in itertools.product(*_leg_roles(pseudoscalar)):
        alp = sum(1 for role in choice if role[0] == 'alp')
        names = [role[2] for role in choice if role[0] == 'meson']
        m = len(names)
        coefficient = source_vertex_coefficient(n, radial, alp, m)
        for role in choice:
            coefficient *= role[1]
        source = source_letter('S' if m % 2 == 0 else 'P', radial, alp)
        for order in itertools.permutations(names):
            result.add_term(((source,) + order,), coefficient)
    return result


def radial_lambda_times_fa():
    """lambda f_a = eps [m_sigma^2 - eps^2 B_0(<S^C_10> - 2<S^C_20>)]/(2 f_pi), from
    m_sigma^2 = 2 lambda f_a^2 + eps^2 B_0(<S^C_10> - 2<S^C_20>) (Cartesian quadratic Lagrangian)."""
    return (TraceSum.scalar(EPSILON * M_SIGMA ** 2 / (2 * F_PI))
            + TraceSum.single(('SC10',), -EPSILON ** 3 * B0 / (2 * F_PI))
            + TraceSum.single(('SC20',), EPSILON ** 3 * B0 / F_PI))


# -d^n V_phi / ds^r dp^q at the vacuum for n >= 3, in units of lambda f_a (cubic) and lambda (quartic):
#   cubic:   -lambda f_a (s^3 + s p^2)       -> (r,q) = (3,0): -6 lambda f_a, (1,2): -2 lambda f_a
#   quartic: -(lambda/4)(s^2 + p^2)^2         -> (4,0): -6 lambda, (2,2): -2 lambda, (0,4): -6 lambda
RADIAL_POTENTIAL_CUBIC = {(3, 0): -6, (1, 2): -2}
RADIAL_POTENTIAL_QUARTIC = {(4, 0): -6, (2, 2): -2, (0, 4): -6}


def radial_potential_vertex(legs, eliminate_lambda=True):
    """V_phi part.  V_phi does not contain the mesons, so every pseudoscalar leg
    enters through its ALP-field component: q = number of pseudoscalar legs and
    the vertex is prod xi times the table above.  No quintic term."""
    n = len(legs)
    radial = sum(1 for leg in legs if leg.radial)
    pseudoscalar = [leg for leg in legs if not leg.radial]
    alp = len(pseudoscalar)
    product = sp.Integer(1)
    for leg in pseudoscalar:
        product *= leg.alp_field
    if product == 0:
        return TraceSum()
    lambda_fa = (radial_lambda_times_fa() if eliminate_lambda
                 else TraceSum.scalar(QUARTIC * F_PI / EPSILON))
    if n == 3 and (radial, alp) in RADIAL_POTENTIAL_CUBIC:
        return lambda_fa.scale(RADIAL_POTENTIAL_CUBIC[(radial, alp)] * product)
    if n == 4 and (radial, alp) in RADIAL_POTENTIAL_QUARTIC:
        # lambda = (lambda f_a) eps / f_pi
        return lambda_fa.scale(RADIAL_POTENTIAL_QUARTIC[(radial, alp)] * product * EPSILON / F_PI)
    return TraceSum()


def meson_kinetic_vertex(legs, dot):
    """Four-meson part of f_pi^2/4 <u_mu u^mu>:
    1/(48 f_pi^2) sum_{perms} (-k_{pi_2}.k_{pi_4}) <[l_{pi_1}, l_{pi_2}][l_{pi_3}, l_{pi_4}]>.
    No cubic or quintic term; s and p do not enter (f_pi constant)."""
    if len(legs) != 4 or any(leg.radial for leg in legs):
        return TraceSum()
    result = TraceSum()
    for order in itertools.permutations(range(4)):
        momentum_factor = -dot(legs[order[1]].momentum, legs[order[3]].momentum) / (48 * F_PI ** 2)
        for choice in itertools.product(*(legs[i].meson for i in order)):
            coefficient = momentum_factor
            for value, _ in choice:
                coefficient *= value
            first, second, third, fourth = (name for _, name in choice)
            for word, sign in (((first, second, third, fourth), 1), ((first, second, fourth, third), -1),
                               ((second, first, third, fourth), -1), ((second, first, fourth, third), 1)):
                result.add_term((word,), sign * coefficient)
    return result


def vertex(legs, dot, eliminate_lambda=True):
    """Full Cartesian vertex rule (word level): source + V_phi + four-meson kinetic.
    The kinetic term of (s, p) gives no vertex."""
    return (source_vertex(legs) + radial_potential_vertex(legs, eliminate_lambda)
            + meson_kinetic_vertex(legs, dot))


def compact_vertex(legs, dot, eliminate_lambda=True):
    """The same rule written with trees: <X^C_{r,q} sum_perm(lambda ...)> and
    <[l, l][l, l]>; V_phi terms are scalars or single traces of source letters."""
    n = len(legs)
    radial = sum(1 for leg in legs if leg.radial)
    pseudoscalar = [leg for leg in legs if not leg.radial]
    result = CompactSum()
    for choice in itertools.product(*_leg_roles(pseudoscalar)):
        alp = sum(1 for role in choice if role[0] == 'alp')
        names = [role[2] for role in choice if role[0] == 'meson']
        m = len(names)
        coefficient = source_vertex_coefficient(n, radial, alp, m)
        for role in choice:
            coefficient *= role[1]
        source = letter(source_letter('S' if m % 2 == 0 else 'P', radial, alp))
        tree = ('P', source, ('S', *(letter(name) for name in names))) if names else source
        result.add_term((tree,), coefficient)
    if n == 4 and radial == 0:
        for order in itertools.permutations(range(4)):
            momentum_factor = -dot(legs[order[1]].momentum, legs[order[3]].momentum) / (48 * F_PI ** 2)
            for choice in itertools.product(*(legs[i].meson for i in order)):
                coefficient = momentum_factor
                for value, _ in choice:
                    coefficient *= value
                first, second, third, fourth = (letter(name) for _, name in choice)
                result.add_term((('P', ('C', first, second), ('C', third, fourth)),), coefficient)
    for key, value in radial_potential_vertex(legs, eliminate_lambda).terms.items():
        result.add_term(tuple(('L', word[0]) for word in key), value)
    return result


# ---------------------------------------------------------------------------
# Construction 1: direct noncommutative expansion of the Cartesian Lagrangian
# ---------------------------------------------------------------------------

def _multilinear_trace_sum(expression, fields, markers=None):
    """Multilinear coefficient of prod(fields) of a traced noncommutative
    expression, as a TraceSum.  markers: momentum symbols kappa_j; exactly two
    per term are replaced by the dot product of the corresponding momenta."""
    result = TraceSum()
    for term in sp.Add.make_args(sp.expand(expression)):
        commuting, ordered = term.args_cnc()
        commuting = sp.Mul(*commuting)
        powers = commuting.as_powers_dict()
        if any(powers.get(one_field, 0) != 1 for one_field in fields):
            continue
        coefficient = commuting.subs({one_field: 1 for one_field in fields})
        if markers:
            labels = [index for index, marker in enumerate(markers)
                      for _ in range(int(powers.get(marker, 0)))]
            if len(labels) != 2:
                continue
            coefficient = coefficient.subs({marker: 1 for marker in markers})
            coefficient *= symbolic_dot({f'q{labels[0] + 1}': 1}, {f'q{labels[1] + 1}': 1})
        word = []
        for factor in ordered:
            base, exponent = factor.as_base_exp()
            word += [base.name] * int(exponent)
        result.add_term((tuple(word),), coefficient)
    return result


def _commuting_multilinear(expression, fields, markers=None):
    total = 0
    for term in sp.Add.make_args(sp.expand(expression)):
        powers = term.as_powers_dict()
        if any(powers.get(one_field, 0) != 1 for one_field in fields):
            continue
        value = term.subs({one_field: 1 for one_field in fields})
        if markers:
            labels = [index for index, marker in enumerate(markers)
                      for _ in range(int(powers.get(marker, 0)))]
            if len(labels) != 2:
                continue
            value = value.subs({marker: 1 for marker in markers})
            value *= symbolic_dot({f'q{labels[0] + 1}': 1}, {f'q{labels[1] + 1}': 1})
        total += value
    return sp.expand(total)


def direct_vertex(number_of_legs, radial, include_lambda=True, include_quadratic=False,
                  up_only=False):
    """Construction 1 for the pattern (radial s legs first, then generic
    pseudoscalar legs (x_j, xi_j)).

    The Lagrangian is expanded with s = sum_{radial j} phi_j,
    p = sum_{pseudoscalar j} xi_j phi_j, Phi = sum_{pseudoscalar j} phi_j x_j and
    d_mu phi_j -> -i kappa_j phi_j; the coefficient of phi_1 ... phi_n is returned
    as a TraceSum.  The source part is taken in the exponential form
    (f_pi^2/4) 2 B_0 <M_hat U^dagger + U M_hat^dagger>, U = exp(i Phi/f_pi),
    with M_hat = sum (S^C_{r,q} + i P^C_{r,q})(s/f_a)^r (p/f_a)^q (or, for
    up_only=True, M_hat = m_u P_u phi^*/f_a with phi^*/f_a = 1 + (s - i p)/f_a).
    The kinetic terms 1/2 (d s)^2 + 1/2 (d p)^2, f_pi^2/4 <d U d U^dagger>, the
    anomalous term -m_1^2/12 <Phi>^2 and -V_phi (with m_phi^2 only if
    include_quadratic) are expanded as well."""
    n = number_of_legs
    fields = sp.symbols(f'phi1:{n + 1}')
    markers = sp.symbols(f'kappa1:{n + 1}')
    matrices = [sp.Symbol(f'x{j + 1}', commutative=False) for j in range(n)]
    alp_components = [sp.Symbol(f'xi{j + 1}', real=True) for j in range(n)]
    pseudoscalar = range(radial, n)
    radial_field = sum(fields[:radial], sp.Integer(0))
    alp_field = sum((alp_components[j] * fields[j] for j in pseudoscalar), sp.Integer(0))
    meson = sum((fields[j] * matrices[j] for j in pseudoscalar), sp.Integer(0))
    f_a = F_PI / EPSILON
    result = TraceSum()

    # Source part, exponential form.  Only the order-m term of U contributes for m = n - r - q.
    def meson_power(power):
        return sp.expand(meson ** power) if power else sp.Integer(1)

    pseudoscalar_count = n - radial
    if up_only:
        projector = sp.Symbol('Pu', commutative=False)
        terms = [((0, 0), UP_QUARK_MASS * projector, UP_QUARK_MASS * projector),
                 ((1, 0), UP_QUARK_MASS * projector, UP_QUARK_MASS * projector),
                 ((0, 1), -sp.I * UP_QUARK_MASS * projector, sp.I * UP_QUARK_MASS * projector)]
    else:
        terms = []
        for q in range(pseudoscalar_count + 1):
            scalar_part = sp.Symbol(source_letter('S', radial, q), commutative=False)
            pseudoscalar_part = sp.Symbol(source_letter('P', radial, q), commutative=False)
            terms.append(((radial, q), scalar_part + sp.I * pseudoscalar_part,
                          scalar_part - sp.I * pseudoscalar_part))
    for (r, q), source, source_dagger in terms:
        m = n - r - q
        if m < 0:
            continue
        field_factor = (radial_field / f_a) ** r * (alp_field / f_a) ** q
        power = meson_power(m)
        u_term = (sp.I / F_PI) ** m / math.factorial(m) * power          # order-m term of U
        u_dagger_term = (-sp.I / F_PI) ** m / math.factorial(m) * power  # order-m term of U^dagger
        expression = F_PI ** 2 / 4 * 2 * B0 * field_factor * (source * u_dagger_term + u_term * source_dagger)
        result = result + _multilinear_trace_sum(expression, fields)

    # Meson kinetic part f_pi^2/4 <d U d U^dagger>.
    if radial == 0 and n % 2 == 0:
        derivative_meson = sum((-sp.I * markers[j] * fields[j] * matrices[j] for j in pseudoscalar),
                               sp.Integer(0))

        def derivative_of_power(power):
            return sum((meson ** k * derivative_meson * meson ** (power - 1 - k) for k in range(power)),
                       sp.Integer(0))
        kinetic = 0
        for k in range(1, n):
            kinetic += ((sp.I / F_PI) ** k / math.factorial(k) * derivative_of_power(k)
                        * (-sp.I / F_PI) ** (n - k) / math.factorial(n - k) * derivative_of_power(n - k))
        result = result + _multilinear_trace_sum(F_PI ** 2 / 4 * kinetic, fields, markers)

    # Anomalous singlet term -m_1^2/12 <Phi>^2 (independent of s and p).
    if n == 2 and radial == 0:
        # <Phi>^2 = (phi_1 <x_1> + phi_2 <x_2>)^2 -> coefficient of phi_1 phi_2: 2 <x_1><x_2>
        result = result + TraceSum({((matrices[0].name,), (matrices[1].name,)): -SINGLET_MASS_SQUARED / 6})

    # Kinetic term of (s, p) and -V_phi (commuting; the trace of the empty word is not used).
    derivative_radial = sum((-sp.I * markers[j] * fields[j] for j in range(radial)), sp.Integer(0))
    derivative_alp = sum((-sp.I * markers[j] * alp_components[j] * fields[j] for j in pseudoscalar),
                         sp.Integer(0))
    scalar = _commuting_multilinear(sp.Rational(1, 2) * (derivative_radial ** 2 + derivative_alp ** 2),
                                    fields, markers)
    modulus_squared = (f_a + radial_field) ** 2 + alp_field ** 2
    potential = sp.Integer(0)
    if include_lambda:
        potential += QUARTIC / 4 * modulus_squared ** 2
    if include_quadratic:
        potential += RADIAL_MASS_PARAMETER / 2 * modulus_squared
    scalar += _commuting_multilinear(-potential, fields)
    if scalar != 0:
        result = result + TraceSum.scalar(scalar)
    return result


# ---------------------------------------------------------------------------
# Up-only limit (D-A4): M_hat = m_u P_u phi^*/f_a
# ---------------------------------------------------------------------------
# S^C_00 = S^C_10 = m_u P_u, P^C_01 = -m_u P_u, every other S^C, P^C zero.


def up_only_letter(name):
    """(coefficient, replacement letter) of a Cartesian source letter in the up-only
    limit, or None if the letter vanishes."""
    if name in ('SC00', 'SC10'):
        return UP_QUARK_MASS, 'Pu'
    if name == 'PC01':
        return -UP_QUARK_MASS, 'Pu'
    return None


def _reduce_projector_word(word):
    """Use P_u P_u = P_u on a cyclic word; return the reduced canonical word."""
    word = list(canonical_word(word))
    changed = True
    while changed and len(word) > 1:
        changed = False
        for index in range(len(word)):
            if word[index] == 'Pu' and word[(index + 1) % len(word)] == 'Pu':
                del word[index]
                changed = True
                break
    return tuple(word)


def up_only_substitution(trace_sum):
    """Up-only limit of a TraceSum: source letters replaced by +-m_u P_u or zero,
    then P_u P_u = P_u and <P_u> = 1.  Other letters are kept."""
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        factor = value
        words = []
        vanishes = False
        for word in key:
            new_word = []
            for name in word:
                if is_cartesian_source(name):
                    replacement = up_only_letter(name)
                    if replacement is None:
                        vanishes = True
                        break
                    factor *= replacement[0]
                    new_word.append(replacement[1])
                else:
                    new_word.append(name)
            if vanishes:
                break
            reduced = _reduce_projector_word(new_word)
            if reduced == ('Pu',):
                continue                      # <P_u> = 1
            words.append(reduced)
        if not vanishes:
            result.add_term(tuple(words), factor)
    return result.simplified()


# ---------------------------------------------------------------------------
# LaTeX output (compact matrix notation only)
# ---------------------------------------------------------------------------

def letter_latex(name):
    if name in ('l+', 'l-'):
        return r'\lambda_{\pi^' + name[1] + '}'
    if name == MINV_LETTER:
        return r'\mathbb M^{-1}_{m_a^2}(\mathcal C)'
    if is_cartesian_source(name):
        kind, r, q = source_indices(name)
        return r'\widehat ' + kind + '^C_{' + f'{r},{q}' + '}'
    if name == 'Pu':
        return 'P_u'
    if name == 'Pd':
        return 'P_d'
    return name


def tree_latex(tree):
    kind = tree[0]
    if kind == 'L':
        return letter_latex(tree[1])
    if kind == 'U':
        return r'\mathbf 1'
    if kind == 'C':
        return '[' + tree_latex(tree[1]) + ',' + tree_latex(tree[2]) + ']'
    if kind == 'A':
        return r'\{' + tree_latex(tree[1]) + ',' + tree_latex(tree[2]) + r'\}'
    if kind == 'S':
        return r'\sum_{\rm perm}\big(' + r'\,'.join(tree_latex(child) for child in tree[1:]) + r'\big)'
    if kind == 'P':
        return r'\,'.join('(' + tree_latex(child) + ')' if child[0] == '+' else tree_latex(child)
                          for child in tree[1:])
    if kind == '+':
        text = ''
        for coefficient, child in tree[1]:
            if coefficient == 1:
                text += ('+' if text else '') + tree_latex(child)
            elif coefficient == -1:
                text += '-' + tree_latex(child)
            else:
                coefficient_text = sp.latex(coefficient)
                sign = '' if (not text or coefficient_text.startswith('-')) else '+'
                text += sign + coefficient_text + r'\,' + tree_latex(child)
        return text
    raise ValueError(tree)


LATEX_NAMES = {**generic.LATEX_NAMES, SINGLET_MASS_SQUARED: 'm_1^2', UP_QUARK_MASS: 'm_u',
               ALP_MASS_COEFFICIENT: 'A', KAPPA_SIGMA: r'\kappa_\sigma', B0: 'B_0'}


def coefficient_latex(value):
    return sp.latex(value, symbol_names=LATEX_NAMES)


def compact_latex(compact_sum, separator='\n'):
    lines = []
    for key, value in sorted(compact_sum.terms.items(), key=lambda item: tuple(map(tree_key, item[0]))):
        traces = ''.join(r'\langle ' + tree_latex(tree) + r'\rangle' for tree in key)
        lines.append(r'\left(' + coefficient_latex(value) + r'\right)' + traces)
    return separator.join(lines) if lines else '0'


def trace_sum_latex(trace_sum):
    """LaTeX of a TraceSum with words written as products of letters (no indices)."""
    parts = []
    for key, value in sorted(trace_sum.simplified().terms.items()):
        traces = ''.join(r'\langle ' + r'\,'.join(letter_latex(name) for name in word) + r'\rangle'
                         for word in key)
        parts.append(r'\left(' + coefficient_latex(value) + r'\right)' + traces)
    return ' + '.join(parts) if parts else '0'


__all__ = ['CartesianLeg', 'radial_leg', 'meson_leg', 'alp_field_leg', 'alp_like_leg', 'generic_legs',
           'source_vertex', 'radial_potential_vertex', 'meson_kinetic_vertex', 'vertex', 'compact_vertex',
           'direct_vertex', 'radial_lambda_times_fa', 'apply_assumptions', 'apply_charge_rule',
           'compact_assumptions', 'up_only_substitution', 'make_dot', 'symbolic_dot']
