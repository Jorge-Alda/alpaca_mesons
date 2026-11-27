"""Up-only specialization of the implicit loop-to-tree crossing with the radial mode.

Tasks 5.1 and 5.2 of the OpenSpec change phi-dependent-loop-tree-crossing-murayama.

Run:
    .venv/bin/python scripts/derive_murayama_radial_up_only_specialization.py
Tests:
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_murayama_radial_up_only_specialization

No loop integral, scalar loop function, loop coefficient or loop amplitude is
computed here.  The loop-size estimate of the general crossing
(scripts/derive_murayama_radial_crossing.py) is specialized by substitution:
its reduced contractions c_k stay undefined functions, and only the invariant
entries from which they are built (frame-invariance script, task 3.4) are
given their up-only values.

Model (memory/alp/meson-mixing/finite-mixing/spectrum-and-poles/quadratic-mixing/only-up-finite-mixing-inputs.md and
memory/alp/radial-mode/couplings-and-decays/finite-mixing/murayama-sigma-two-pseudoscalar-decays.md):
  * (n_u, n_d, n_s) = (1, 0, 0), f_a = f_phi, theta = 0, CP-conserving vacuum.
  * m_u = kappa_u f_a / sqrt(2) (DiLuzio:2026poh), hence
        f_a = sqrt(2) m_u / kappa_u,   epsilon = f_pi / f_a = f_pi kappa_u / (sqrt(2) m_u).
  * Partonic source -qbar_R M q_L - h.c. with, in the polar representative,
        M_hat = diag[m_u (1 + sigma/f_a) exp(-i a/f_a), m_d, m_s],
    expanded as S_hat = sum S_hat_rs (sigma/f_a)^r (a/f_a)^s (first index
    counts sigma), likewise P_hat.  Coefficients have mass dimension one.
  * Theta_hat = 0, so C = B_0 P_hat_01 (scope fixed by the user).
  * M_s(X) = B_0 {S_hat_00, X} + (m_1^2/3) <X> 1 - s X on neutral flavour
    matrices; m_1^2 is the singlet mass parameter of the only-up inputs record.

Units: masses, B_0, f_pi, f_a in GeV; squared masses, s and Dalitz
invariants in GeV^2; kappa_u, lambda and epsilon dimensionless.

kappa_u and the quartic coupling lambda stay symbols throughout; nothing
is substituted numerically.

Radial self-energy of ALP-like lines (closed sigma line on the derivative
vertex (1/(2 f_a^2)) sigma^2 (d a)^2 of the polar kinetic term):
  * on the external ALP-like leg (class R5a of the general script): kept,
    epsilon^2 relative to the tree, scale m_sigma, combined with the radial
    tadpole in one magnitude.  Its residue part is a_component^2 = 1/N(m_a^2)
    = 1/F'(m_a^2); the mass shift m_a^2/F'(m_a^2) is absorbed into the
    physical m_a (assumption, ALP-only convention).
  * on an internal ALP-like line of a tree exchange (class R5b): absent in
    the up-only model, because every cubic vertex vanishes and the tree has
    no ALP-like exchange (murayama-cubic-vertex-parity-rule.md).  Its
    reduced contraction is set to zero by specialize (a structural zero, not
    an evaluation of a loop).
"""
import sympy as sp

import derive_murayama_radial_crossing as general
import derive_murayama_radial_frame_invariance as frame
from derive_murayama_radial_crossing import (
    ALP_POLE_MASS_SQUARED, B_ZERO, DAUGHTERS, EPSILON, F_PI, MAXIMUM_RATIO,
    PAIR_INVARIANT_12, PAIR_INVARIANT_23, PARENT_OVERLAP_TRACE, QUARTIC_COUPLING,
    RADIAL_COUPLING_SWITCH, RADIAL_MASS_SQUARED)

UP_QUARK_MASS = sp.Symbol('m_u', positive=True)
DOWN_QUARK_MASS = sp.Symbol('m_d', positive=True)
STRANGE_QUARK_MASS = sp.Symbol('m_s', positive=True)
SINGLET_MASS_SQUARED = sp.Symbol('m_1^2', positive=True)
KAPPA_U = sp.Symbol(r'\kappa_u', positive=True)
SQUARED_MOMENTUM = sp.Symbol('s', real=True)      # argument of M_s, GeV^2
FLAVOURS = ('u', 'd', 's')
QUARK_MASSES = {'u': UP_QUARK_MASS, 'd': DOWN_QUARK_MASS, 's': STRANGE_QUARK_MASS}
LEGS = ('a',) + DAUGHTERS

# Derivative F'(s) of the secular function, evaluated at a pole mass.  It is
# the pole normalization N(m^2) (checked in normalization_is_secular_derivative).
SECULAR_DERIVATIVE = sp.Function("F'", positive=True)


# ---------------------------------------------------------------------------
# Model parameters in terms of kappa_u.
# ---------------------------------------------------------------------------

def alp_decay_constant():
    """f_a = sqrt(2) m_u / kappa_u (GeV)."""
    return sp.sqrt(2) * UP_QUARK_MASS / KAPPA_U


def up_only_epsilon():
    """epsilon = f_pi / f_a = f_pi kappa_u / (sqrt(2) m_u)."""
    return F_PI / alp_decay_constant()


def up_only_source_coefficients(radial_order=2, alp_order=2):
    """Up-quark entries (S_hat_rs, P_hat_rs) of m_u (1 + y) exp(-i x).

    y = sigma/f_a, x = a/f_a.  The d and s entries are constant (S_hat_00
    only); every coefficient with r + s > 0 is proportional to e_1."""
    x, y = sp.symbols('x y', real=True)
    source = sp.expand(sp.series(UP_QUARK_MASS * (1 + y) * sp.exp(-sp.I * x), x, 0,
                                 alp_order + 1).removeO())
    coefficients = {}
    for r in range(radial_order + 1):
        for s in range(alp_order + 1):
            value = sp.expand(source.coeff(y, r).coeff(x, s))
            coefficients[(r, s)] = (sp.re(value), sp.im(value))
    return coefficients


def flavour_entries(kind, r, s):
    """Diagonal flavour entries of S_hat_rs (kind 'S') or P_hat_rs (kind 'P')."""
    up_entry = up_only_source_coefficients()[(r, s)][0 if kind == 'S' else 1]
    entries = {'u': up_entry, 'd': 0, 's': 0}
    if kind == 'S' and (r, s) == (0, 0):
        entries.update(d=DOWN_QUARK_MASS, s=STRANGE_QUARK_MASS)
    return entries


def flavour_trace(*factors):
    return sum(sp.Mul(*[factor[q] for factor in factors]) for q in FLAVOURS)


# ---------------------------------------------------------------------------
# Pole operator, secular function and pole normalization.
# ---------------------------------------------------------------------------

def pole_inverse(source, squared_momentum=SQUARED_MOMENTUM):
    """Diagonal solution X of M_s(X) = source (diagonal source entries)."""
    unknowns = sp.symbols('X_u X_d X_s')
    equations = [(2 * B_ZERO * QUARK_MASSES[q] - squared_momentum) * unknowns[k]
                 + SINGLET_MASS_SQUARED / 3 * sum(unknowns) - source[q]
                 for k, q in enumerate(FLAVOURS)]
    solution = sp.solve(equations, unknowns, dict=True)[0]
    return {q: sp.factor(solution[unknowns[k]]) for k, q in enumerate(FLAVOURS)}


def secular_function(squared_momentum=SQUARED_MOMENTUM, epsilon=EPSILON):
    """F(s) = s + 2 B_0 epsilon^2 [<S_hat_02> + B_0 <P_hat_01 M_s^{-1}(P_hat_01)>].

    The ALP-like and meson poles of the neutral block are the zeros of F
    (memory/alp/radial-mode/couplings-and-decays/finite-mixing/murayama-sigma-two-pseudoscalar-decays.md)."""
    inverse = pole_inverse(flavour_entries('P', 0, 1), squared_momentum)
    return (squared_momentum + 2 * B_ZERO * epsilon**2 * (
        flavour_trace(flavour_entries('S', 0, 2))
        + B_ZERO * flavour_trace(flavour_entries('P', 0, 1), inverse)))


def pole_normalization(squared_momentum=SQUARED_MOMENTUM, epsilon=EPSILON):
    """N(s) = 1 + 2 B_0^2 epsilon^2 <[M_s^{-1}(P_hat_01)]^2>."""
    inverse = pole_inverse(flavour_entries('P', 0, 1), squared_momentum)
    return 1 + 2 * B_ZERO**2 * epsilon**2 * flavour_trace(inverse, inverse)


def normalization_is_secular_derivative():
    """N(s) - dF/ds; zero identically."""
    return sp.simplify(pole_normalization() - sp.diff(secular_function(), SQUARED_MOMENTUM))


def up_entry_at_pole(pole_mass_squared, epsilon=EPSILON):
    """(M^{-1}_{m^2}(P_hat_01))_uu at a pole: (m^2/(B_0 m_u epsilon^2) - 1)/(2 B_0)."""
    return (pole_mass_squared / (B_ZERO * UP_QUARK_MASS * epsilon**2) - 1) / (2 * B_ZERO)


def up_entry_pole_identity():
    """(2 B_0 X_u(s) + 1) B_0 m_u epsilon^2 - s + F(s); zero identically.

    At F(m^2) = 0 it gives up_entry_at_pole(m^2)."""
    up_entry = pole_inverse(flavour_entries('P', 0, 1))['u']
    return sp.simplify((2 * B_ZERO * up_entry + 1) * B_ZERO * UP_QUARK_MASS * EPSILON**2
                       - SQUARED_MOMENTUM + secular_function())


# ---------------------------------------------------------------------------
# Task 5.1: up-only tree and the specialized crossing.
# ---------------------------------------------------------------------------

def contact_amplitude():
    """Four-point contact V^(4)_{a123} of the up-only tree.

    Kept as an undefined function: its forms are in
    memory/alp/three-pseudoscalar/finite-mixing/murayama-unmixed-basis-quartic-vertex.md, and
    it vanishes for three distinct neutral daughters
    (murayama-distinct-neutral-contact-cancellation.md)."""
    return sp.Function(r'\mathcal{V}^{(4)}_{a123}')(EPSILON, PAIR_INVARIANT_12, PAIR_INVARIANT_23)


def leg_normalization_up_only(leg):
    """N(m_leg^2) = F'(m_leg^2)."""
    return SECULAR_DERIVATIVE(general.leg_pole_mass_squared(leg))


def up_only_trace_substitution():
    """Up-only values of the traces in the general radial coupling g_xy(q^2).

    S_hat_10 = m_u e_1, P_hat_11 = -m_u e_1, S_hat_12 = -m_u/2 e_1 and
    M^{-1}(P_hat_01) is diagonal, so only its up entry enters.  The up entry
    is fixed by the pole equation (up_entry_at_pole).  The overlap traces of
    the normalizations are replaced so that N(m^2) = F'(m^2)."""
    up_entry = {leg: up_entry_at_pole(general.leg_pole_mass_squared(leg)) for leg in LEGS}
    radial_scalar = flavour_entries('S', 1, 0)
    radial_pseudoscalar = flavour_entries('P', 1, 1)
    substitution = {sp.Symbol(r'\langle\widehat S_{12}\rangle', real=True):
                    flavour_trace(flavour_entries('S', 1, 2))}
    for leg in LEGS:
        substitution[general.single_trace(leg)] = radial_pseudoscalar['u'] * up_entry[leg]
        if leg != 'a':
            substitution[general.leg_overlap_trace(leg)] = (
                (leg_normalization_up_only(leg) - 1) / (2 * B_ZERO**2 * EPSILON**2))
        for other in LEGS:
            if other != leg:
                substitution[general.anticommutator_trace(leg, other)] = (
                    2 * radial_scalar['u'] * up_entry[leg] * up_entry[other])
    substitution[PARENT_OVERLAP_TRACE] = (leg_normalization_up_only('a') - 1) / (2 * EPSILON**2)
    return substitution


def closed_form_radial_coupling(first, second, momentum_squared):
    """t_sigma g_xy(q^2) of the only-up three-pseudoscalar record:
    [q^2 - m_x^2 - m_y^2 + m_x^2 m_y^2/(B_0 m_u epsilon^2)] / (f_a sqrt(N_x N_y))."""
    first_mass, second_mass = (general.leg_pole_mass_squared(leg) for leg in (first, second))
    bracket = (momentum_squared - first_mass - second_mass
               + first_mass * second_mass / (B_ZERO * UP_QUARK_MASS * EPSILON**2))
    return RADIAL_COUPLING_SWITCH * bracket * EPSILON / F_PI / sp.sqrt(
        leg_normalization_up_only(first) * leg_normalization_up_only(second))


def up_only_tree():
    """V^(4)_{a123} + sum_(r|uv) g_ar(s_uv) g_uv(s_uv) / (m_sigma^2 - s_uv).

    No meson or ALP exchange: all cubic vertices vanish in the up-only model
    (memory/alp/three-pseudoscalar/murayama-cubic-vertex-parity-rule.md)."""
    invariants = general.pair_invariants()
    exchange = 0
    for spectator in DAUGHTERS:
        pair = tuple(n for n in DAUGHTERS if n != spectator)
        momentum_squared = invariants[pair]
        exchange += (closed_form_radial_coupling('a', spectator, momentum_squared)
                     * closed_form_radial_coupling(*pair, momentum_squared)
                     / (RADIAL_MASS_SQUARED - momentum_squared))
    return contact_amplitude() + exchange


def structural_zero_contractions():
    """Reduced contractions that vanish in the up-only model because the
    diagram does not exist: the radial self-energy of an internal ALP-like
    line needs an ALP-like tree exchange, and the up-only tree has none."""
    return {general.reduced_contraction('R5b'): 0}


def specialize(expression):
    """General expression -> up-only, still in epsilon and m_sigma^2."""
    return expression.subs(up_only_trace_substitution()).subs(
        general.tree_without_radial_exchange(), contact_amplitude()).subs(
        structural_zero_contractions())


def tree_substitution_residual():
    """specialize(general exact tree) - up_only_tree(); zero."""
    return sp.simplify(specialize(general.exact_tree()) - up_only_tree())


def radial_coupling_residuals():
    """specialize(general g_xy(q^2)) - closed form, for every pair of legs used."""
    momentum_squared = sp.Symbol('q^2', real=True)
    pairs = [('a', n) for n in DAUGHTERS] + [(1, 2), (1, 3), (2, 3)]
    return {pair: sp.simplify(specialize(general.radial_coupling(*pair, momentum_squared))
                              - closed_form_radial_coupling(*pair, momentum_squared))
            for pair in pairs}


def up_only_crossing():
    """(F, R) of the up-only crossing, F = estimate - R_max |tree| = 0.

    epsilon -> f_pi kappa_u/(sqrt(2) m_u), m_sigma^2 -> radial_mass_squared(),
    N(m^2) -> F'(m^2).  The reduced contractions c_k remain undefined
    functions (of kappa_u through epsilon); m_a^2 is the ALP-like root of
    F(s) = 0 (task 5.2) and stays the symbol m_a^2."""
    equation, ratio = general.crossing_equation()
    replacements = {EPSILON: up_only_epsilon(), RADIAL_MASS_SQUARED: radial_mass_squared()}
    return (specialize(equation).subs(replacements), specialize(ratio).subs(replacements))


def up_only_crossing_compact():
    """Same equation written with the closed-form tree and epsilon kept as a
    symbol, for printing; epsilon = f_pi kappa_u/(sqrt(2) m_u) is understood."""
    equation = specialize(general.loop_size_estimate()) - MAXIMUM_RATIO * sp.Abs(up_only_tree())
    return equation.subs(RADIAL_MASS_SQUARED, radial_mass_squared())


def secular_polynomial_is_irreducible():
    """True if the secular polynomial has no factor over the rational
    functions of the parameters (so no root separates from the others)."""
    _, factors = sp.factor_list(secular_polynomial().as_expr())
    return len(factors) == 1 and factors[0][1] == 1


def up_only_contraction_entries():
    """Up-only values of the invariant entries from which the estimate's
    contractions are built (estimate_inputs of the frame-invariance script).

    Returns (entries, identities): entries in the symbols of that script
    (f_a, f_pi, m_0^2 playing the role of m_1^2, m_i^2, m_j^2, lambda), and the
    residuals of the identities that fix them at the poles.  Each residual is
    zero identically in m_i^2, m_j^2; F(m^2) is the secular function in the
    same symbols, so at the poles F = 0 gives up_only_contraction_values_at_poles."""
    frames = frame.diagonal_frames(alp_order=2, radial_order=2)
    n_u, n_d, n_s = frames['charges']
    charges = {n_u: 1, n_d: 0, n_s: 0}
    coefficients = {q: tuple({key: sp.sympify(value).subs(charges) for key, value in part.items()}
                             for part in frames['original'][q]) for q in frame.FLAVOURS}
    entries = frame.estimate_inputs(coefficients, 0)
    epsilon = frame.F_PI / frame.F_A
    m_u = frames['masses'][0]
    scale = frame.B0 * m_u * epsilon**2
    scalar_zero = {q: coefficients[q][0][(0, 0)] for q in frame.FLAVOURS}
    pseudoscalar_one = {q: coefficients[q][1][(0, 1)] for q in frame.FLAVOURS}

    def secular(pole):
        inverse = frame.pole_inverse(scalar_zero, pseudoscalar_one, pole)
        return pole + 2 * frame.B0 * epsilon**2 * (
            entries['tr S_hat_02']
            + frame.B0 * sum(pseudoscalar_one[q] * inverse[q] for q in frame.FLAVOURS))

    def normalization(pole):
        inverse = frame.pole_inverse(scalar_zero, pseudoscalar_one, pole)
        return 1 + 2 * epsilon**2 * sum((frame.B0 * value)**2 for value in inverse.values())

    pole_i, pole_j = frame.POLE_I, frame.POLE_J
    secular_i, secular_j = secular(pole_i), secular(pole_j)
    coupling_bracket = entries['sigma-PP coupling g_ij'] * frame.F_A * sp.sqrt(
        normalization(pole_i) * normalization(pole_j))
    identities = {
        'tr S_hat_02 + m_u/2': entries['tr S_hat_02'] + m_u / 2,
        'm_sigma^2 - 2 lambda f_a^2 - B_0 m_u epsilon^2':
            entries['m_sigma^2'] - 2 * frame.QUARTIC * frame.F_A**2 - scale,
        'tr(C M^-1(C)) - (F(m_i^2) - m_i^2 + B_0 m_u epsilon^2)/(2 epsilon^2)':
            entries['ALP pole term tr(C M^{-1}(C)) at m_i^2']
            - (secular_i - pole_i + scale) / (2 * epsilon**2),
        'tr(P_hat_11 Omega_i) + m_u (m_i^2 - F(m_i^2) - B_0 m_u epsilon^2)/(2 B_0 m_u epsilon^2)':
            entries['sigma-a-P contraction tr(P_hat_11 Omega_i)']
            + m_u * (pole_i - secular_i - scale) / (2 * scale),
        'f_a sqrt(N_i N_j) g_ij - (m_i^2 - F_i)(m_j^2 - F_j)/(B_0 m_u epsilon^2) - m_sigma^2 + m_i^2 + m_j^2':
            coupling_bracket - (pole_i - secular_i) * (pole_j - secular_j) / scale
            - entries['m_sigma^2'] + pole_i + pole_j,
    }
    # Radial self-energy of the external ALP-like leg: its residue part is
    # 1/N(m_i^2), and N(m_i^2) is the derivative of the secular function.
    identities["N(m_i^2) - F'(m_i^2) (residue part of the radial self-energy of the ALP leg is 1/N)"] = (
        entries['N(m_i^2)'] - sp.diff(secular_i, pole_i))
    for name in ('R2: tr(S_hat_20 {Omega_i, Omega_j})', 'R2: tr(P_hat_21 (Omega_i + Omega_j))',
                 'R2: tr S_hat_22', 'R2: tr S_hat_20 (sigma^2 a^0 vertex)'):
        identities[name] = entries[name]
    return entries, {name: sp.simplify(value) for name, value in identities.items()}


def up_only_radial_alp_self_energy_parts():
    """Parts of the radial self-energy of the external ALP-like leg in the
    up-only model, from general.radial_alp_self_energy_parts with
    a-component squared 1/F'(m_a^2) and p_a^2 = m_a^2 on shell."""
    momentum_squared = sp.Symbol('p_a^2', positive=True)
    return general.radial_alp_self_energy_parts(
        momentum_squared, 1 / SECULAR_DERIVATIVE(ALP_POLE_MASS_SQUARED), ALP_POLE_MASS_SQUARED)


def up_only_contraction_values_at_poles():
    """The same entries at poles (F = 0), in the symbols of the crossing script."""
    epsilon, m_u = EPSILON, UP_QUARK_MASS
    scale = B_ZERO * m_u * epsilon**2

    def alp_pole_term(mass_squared):
        return (scale - mass_squared) / (2 * epsilon**2)

    def radial_alp_meson(mass_squared):
        return -m_u * (mass_squared - scale) / (2 * scale)

    return {
        'tr(C M^-1_{m_a^2}(C))': alp_pole_term(ALP_POLE_MASS_SQUARED),
        'tr S_hat_02': -m_u / 2,
        'tr(P_hat_11 B_0 M^-1_{m_a^2}(P_hat_01))': radial_alp_meson(ALP_POLE_MASS_SQUARED),
        'm_sigma^2': radial_mass_squared(),
        'R2 entries built from S_hat_20, P_hat_21, S_hat_22': 0,
        'R5a residue part 1/N(m_a^2) (radial self-energy of the external ALP leg)':
            1 / SECULAR_DERIVATIVE(ALP_POLE_MASS_SQUARED),
        'R5a mass shift m_a^2/N(m_a^2) (absorbed into the physical m_a)':
            ALP_POLE_MASS_SQUARED / SECULAR_DERIVATIVE(ALP_POLE_MASS_SQUARED),
        'R5b (radial self-energy of an internal ALP-like line)': 0,
    }


# ---------------------------------------------------------------------------
# Task 5.2: m_sigma and m_a in kappa_u and lambda.
# ---------------------------------------------------------------------------

def radial_mass_squared():
    """m_sigma^2 = 2 lambda f_a^2 + B_0 epsilon^2 (<S_hat_10> - 2 <S_hat_20>), up-only.

    Returned as B_0 f_pi^2 kappa_u^2/(2 m_u) + 4 lambda m_u^2/kappa_u^2."""
    relation, _ = general.quartic_coupling_at_fixed_radial_mass()
    traces = {sp.Symbol(r'\langle\widehat S_{10}\rangle', real=True): flavour_trace(flavour_entries('S', 1, 0)),
              sp.Symbol(r'\langle\widehat S_{20}\rangle', real=True): flavour_trace(flavour_entries('S', 2, 0))}
    return sp.expand(relation.rhs.subs(traces).subs(EPSILON, up_only_epsilon()))


def radial_mass_terms():
    """(first term, lambda term, ratio lambda term / first term)."""
    mass_squared = radial_mass_squared()
    lambda_term = sp.expand(QUARTIC_COUPLING * sp.diff(mass_squared, QUARTIC_COUPLING))
    first_term = sp.expand(mass_squared - lambda_term)
    return first_term, lambda_term, sp.simplify(lambda_term / first_term)


def expected_radial_mass_ratio():
    return 8 * QUARTIC_COUPLING * UP_QUARK_MASS**3 / (B_ZERO * F_PI**2 * KAPPA_U**4)


def secular_polynomial():
    """Numerator of F(s) at epsilon = f_pi kappa_u/(sqrt(2) m_u), as a polynomial in s.

    Its roots are the ALP-like pole m_a^2 and the neutral meson poles."""
    numerator, _ = sp.fraction(sp.together(secular_function(epsilon=up_only_epsilon())))
    return sp.Poly(sp.expand(numerator), SQUARED_MOMENTUM)


def alp_like_root_small_kappa_limit():
    """Limit only (not the result): m_a^2 = epsilon^2 B_0 m_u [1 - 2 B_0 m_u w(0)] + O(epsilon^4),
    with w(0) = <e_1 M_0^{-1}(e_1)>, from F(epsilon^2 t) = 0 at leading order in epsilon."""
    rescaled = sp.Symbol('t', real=True)
    leading = sp.limit(secular_function(rescaled * EPSILON**2) / EPSILON**2, EPSILON, 0)
    return sp.factor(EPSILON**2 * sp.solve(leading, rescaled)[0])


def main():
    print('Up-only model: f_a = sqrt(2) m_u / kappa_u, epsilon =', up_only_epsilon())
    print('Source coefficients (up entry; first index sigma, second a):')
    for key, (scalar, pseudoscalar) in up_only_source_coefficients().items():
        print(f'    S_hat_{key[0]}{key[1]} = {scalar},  P_hat_{key[0]}{key[1]} = {pseudoscalar}')
    print('\nN(s) - dF/ds =', normalization_is_secular_derivative())
    print('(2 B_0 X_u(s) + 1) B_0 m_u epsilon^2 - s + F(s) =', up_entry_pole_identity())

    print('\nTask 5.1')
    print('General radial couplings with up-only traces minus closed form:',
          radial_coupling_residuals())
    print('General exact tree with up-only traces minus up-only tree:', tree_substitution_residual())
    print('Up-only tree (contact V4 kept; zero for three distinct neutral daughters):')
    print(sp.latex(up_only_tree()), '\n')
    entries, identities = up_only_contraction_entries()
    print('Invariant entries of the contractions, up-only identities (residuals):')
    for name, value in identities.items():
        print(f'    {name}: {value}')
    print('Values at the poles:')
    for name, value in up_only_contraction_values_at_poles().items():
        print(f'    {name} = {sp.latex(value)}')
    print('Radial self-energy of the external ALP leg (sigma^2 (d a)^2 / (2 f_a^2), closed sigma line):')
    for name, value in up_only_radial_alp_self_energy_parts().items():
        print(f'    {name}: {sp.latex(value)}')
    print('    The mass shift is absorbed into the physical m_a (assumption, ALP-only convention).')
    print('    The radial self-energy of an internal ALP-like line is absent: no ALP-like tree exchange.')
    print('\nUp-only crossing estimate - R_max |tree| = 0, with epsilon = f_pi kappa_u/(sqrt(2) m_u),')
    print("N(m^2) = F'(m^2) and m_sigma^2 in kappa_u, lambda:")
    print(sp.latex(up_only_crossing_compact()), '\n')

    print('Task 5.2')
    first, lambda_term, ratio = radial_mass_terms()
    print('m_sigma^2 =', sp.latex(radial_mass_squared()))
    print('lambda term / first term =', sp.latex(ratio),
          '; equals 8 lambda m_u^3/(B_0 f_pi^2 kappa_u^4):',
          sp.simplify(ratio - expected_radial_mass_ratio()) == 0)
    polynomial = secular_polynomial()
    print('Secular polynomial in s: degree', polynomial.degree(),
          '; depends on lambda:', polynomial.as_expr().has(QUARTIC_COUPLING))
    for power, coefficient in zip(range(polynomial.degree(), -1, -1), polynomial.all_coeffs()):
        print(f'    s^{power}:', sp.latex(sp.factor(coefficient)))
    print('Irreducible over the rational functions of the parameters:',
          secular_polynomial_is_irreducible())
    print('m_a^2 is the ALP-like root of this quartic. It is kept implicit: a radical')
    print('(Ferrari) form exists for degree 4 but is not used here.')
    print('Small-kappa_u limit of the ALP-like root (limit only):',
          sp.latex(alp_like_root_small_kappa_limit()))


if __name__ == '__main__':
    main()
