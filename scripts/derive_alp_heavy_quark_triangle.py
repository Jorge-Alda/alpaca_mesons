"""Derive the one-loop heavy-quark contribution to the doubly virtual ALP two-photon form factor.

Run .venv/bin/python scripts/derive_alp_heavy_quark_triangle.py [--test].

Units and conventions: natural units, masses and momenta in GeV, the form
factor F_{a gamma* gamma*} in 1/GeV. Metric (+---), epsilon^{0123}=+1
(epsilon_{0123}=-1), gamma5 = i gamma0 gamma1 gamma2 gamma3, so that
tr(gamma0 gamma1 gamma2 gamma3 gamma5) = -4i. The form factor is defined by
i M = i e^2 R F_{a gamma* gamma*}(q1^2, q2^2), with
R = epsilon_{mu nu rho sigma} eps1^mu eps2^nu q1^rho q2^sigma, an incoming
ALP of momentum p = q1 + q2 (p^2 = m_a^2) and outgoing photons.

Lagrangian (Alda:2026nrd, eq. lagr_partons, restricted to one heavy quark h
of mass m_h, charge Q_h, N_c colors, diagonal couplings):
  L = (d_mu a / f_a) hbar gamma^mu (c_hh^V + c_hh^A gamma5) h
      + i (a / f_a) 2 m_h c_hh^P hbar gamma5 h
      - e Q_h A_mu hbar gamma^mu h.
Feynman rules: quark propagator i (kslash + m_h) / (k^2 - m_h^2 + i0), photon
vertex -i e Q_h gamma^mu, vertex of the pseudoscalar coupling
-(2 m_h c_hh^P / f_a) gamma5, vertex of the derivative coupling with an
incoming ALP (c_hh^A / f_a) pslash gamma5, factor -1 for the closed quark loop.

Inputs:
- memory/sources/Alda:2026nrd/main_arXiv.tex, eq. lagr_partons: couplings.
- memory/alp/lepton-pair/alp-high-energy-tff-conventions.md: tensor R, sign
  and dimension of F, momentum labels.
- memory/rcht/odd-parity/odd-alp-source-lsz-assembly.md and
  memory/alp/lepton-pair/alp-doubly-virtual-tff.md: the local photon term
  F = C_{a gamma gamma} / (4 pi^2 f_a), C = c_gamma + 2 N_c tr[Q^2 c^A] + ...,
  used here only as the reference for the overall sign.
- memory/alp/lepton-pair/alp-heavy-flavor-thresholds.md: local expansion
  F = 2 i Q_h^2 p_sigma <0|hbar gamma^sigma gamma5 h|a(p)> / ((q^2 - m_h^2) m_a^2).

Approximations and assumptions (all explicit):
- one quark loop, leading order in alpha_s and in 1/f_a, free heavy quark;
- the vector coupling c_hh^V does not contribute (Part 6): its trace has no
  gamma5 and therefore no Levi-Civita tensor, and the two photon orderings
  cancel (charge conjugation, Furry theorem). The cancellation uses the change
  of loop momentum k -> -k - p in a linearly divergent integral; it is assumed
  that the regularization allows it (dimensional regularization, or the
  Pauli-Villars difference, whose integrand falls as k^{-5});
- the derivative coupling is regulated with one Pauli-Villars quark of mass
  M, taken to infinity. The momentum shifts used in the axial Ward identity
  are allowed only for the regulated difference;
- Part 5: the vacuum-to-ALP matrix element of the axial current is
  ultraviolet divergent. Integrals of the form int d^4k k^sigma/(k^2-m_h^2)
  are set to zero (true in dimensional regularization), and the
  logarithmically divergent scalar integral is cut off at Euclidean loop
  momentum Lambda. The constant that accompanies ln(Lambda^2) depends on this
  choice; only the terms that are not polynomial in m_h^2 are independent of it;
- expansions in m_a^2 are Taylor expansions, valid for m_a^2 < 4 m_h^2;
- Part 7 (m_a > 2 m_h): the quark propagator i(kslash+m_h)/(k^2-m_h^2+i0)
  gives the Feynman-parameter denominator with m_h^2 -> m_h^2 - i0
  (equivalently p^2 -> p^2 + i0). The relation
  1/(D - i0) = principal value of 1/D + i pi delta(D) is used as an input and
  is tested by quadrature at complex mass. With this prescription the
  imaginary part of the loop integral is >= 0.

Notation (approved by the user): with tau_a = 4 m_h^2/m_a^2,
tau_1 = 4 m_h^2/q_1^2, tau_2 = 4 m_h^2/q_2^2 (a real photon, q_i^2 = 0, means
tau_i -> infinity) the dimensionless loop function is
  B_4(tau_a, tau_1, tau_2) = 1 - int_0^1 dx int_0^{1-x} dy 2 m_h^2 /
      [m_h^2 - x y m_a^2 - x (1-x-y) q_2^2 - y (1-x-y) q_1^2 - i0],
that is, one minus the "loop integral" of Parts 1-7, and
  F_h = (2 N_c Q_h^2/(4 pi^2 f_a)) [c_hh^P + (c_hh^A - c_hh^P) B_4].
The name B_3 is not used for it because the literature already has
B_3(tau_1, tau_2) = 1 + tau_1 tau_2/(tau_1 - tau_2) [f^2(tau_1) - f^2(tau_2)]
(Bauer:2021mvw Eq. (2.82), Alda:2024cxn Eq. (A1)) for one real photon, and
B_1(tau) = 1 - tau f^2(tau) (Bauer:2021mvw Eq. (2.71)) for two real photons;
see memory/extracts/Bauer:2021mvw/alp-photon-loop-functions-b1-b3.md and
memory/extracts/Alda:2024cxn/alp-photon-loop-function-b3-offshell.md. Part 8
shows B_4(tau_a, tau_1, infinity) = B_3(tau_a, tau_1) and
B_4(tau_a, infinity, infinity) = B_1(tau_a), the symmetry of B_4 under any
exchange of its three arguments, and its limits. The literature gives
f(tau) for tau >= 1 and 0 < tau < 1 only; the form for tau < 0 (spacelike
momentum) used here, f = (i/2) ln[(sqrt(1-tau)+1)/(sqrt(1-tau)-1)], is derived
in Part 8 from the closed form of Part 1 and is not taken from the literature.

The closed forms are proven by differentiation with respect to the momentum
invariant and comparison at zero; the antiderivatives used are checked by
differentiation. Numerical quadrature (mpmath, fixed precision) is used only
in the unit tests; the two-dimensional quadratures at complex mass of Part 7
take most of the test time (about one minute in total for Part 7 and a
similar time for the complex-mass points of Part 8).
"""

from functools import lru_cache
import sys
import unittest

import mpmath
import sympy as sp
from sympy.physics.matrices import mgamma


def show(label, expression):
    """Print expressions through the standard SymPy LaTeX renderer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


METRIC = sp.diag(1, -1, -1, -1)
GAMMA = tuple(mgamma(index) for index in range(4))
GAMMA5 = mgamma(5)
IDENTITY = sp.eye(4)

quark_mass = sp.Symbol("m_h", positive=True)
regulator_mass = sp.Symbol("M", positive=True)
decay_constant = sp.Symbol("f_a", positive=True)
color_number = sp.Symbol("N_c", positive=True)
quark_charge = sp.Symbol("Q_h", real=True)
electric_charge = sp.Symbol("e", positive=True)
axial_coupling = sp.Symbol("c_{hh}^{A}", real=True)
pseudoscalar_coupling = sp.Symbol("c_{hh}^{P}", real=True)
photon_coupling = sp.Symbol("c_{\\gamma}", real=True)
rotation_parameter = sp.Symbol("kappa", real=True)
first_virtuality = sp.Symbol("q_1^2", real=True)
second_virtuality = sp.Symbol("q_2^2", real=True)
alp_momentum_squared = sp.Symbol("p^2", real=True)
alp_mass = sp.Symbol("m_a", positive=True)
spacelike_magnitude = sp.Symbol(r"\left|q^{2}\right|", positive=True)
cutoff = sp.Symbol("Lambda", positive=True)
parameter_x, parameter_y = sp.symbols("x y", nonnegative=True)
parameter_u = sp.Symbol("u", nonnegative=True)
# Arguments of B_4: 4 m_h^2/m_a^2, 4 m_h^2/q_1^2, 4 m_h^2/q_2^2 (Part 8).
tau_alp = sp.Symbol("tau_a", real=True)
tau_first_photon = sp.Symbol("tau_1", real=True)
tau_second_photon = sp.Symbol("tau_2", real=True)
# One argument 4 m_h^2/q^2 at spacelike momentum, q^2 < 0.
tau_spacelike = sp.Symbol("tau", negative=True)

# Common factor 2 N_c Q_h^2 / (4 pi^2 f_a) of the reference local term.
REFERENCE_FACTOR = 2 * color_number * quark_charge**2 / (4 * sp.pi**2 * decay_constant)


def lower(vector):
    """Lower the index of a contravariant four-vector."""
    return METRIC * vector


def dot(first, second):
    """Minkowski product of two contravariant four-vectors."""
    return sp.expand((first.T * METRIC * second)[0])


def slash(vector):
    """Contract a contravariant four-vector with the Dirac matrices."""
    lowered = lower(vector)
    return sum((lowered[index] * GAMMA[index] for index in range(4)), sp.zeros(4))


def four_vector(name):
    """Real symbolic contravariant four-vector."""
    return sp.Matrix(sp.symbols(f"{name}_0 {name}_1 {name}_2 {name}_3", real=True))


def feynman_denominator(mass, virtuality_one, virtuality_two, momentum_squared,
                        x=parameter_x, y=parameter_y):
    """Denominator m^2 - x y p^2 - x (1-x-y) q2^2 - y (1-x-y) q1^2 of the triangle."""
    return (mass**2 - x * y * momentum_squared - x * (1 - x - y) * virtuality_two
            - y * (1 - x - y) * virtuality_one)


def loop_integrand(mass, virtuality_one, virtuality_two, momentum_squared,
                   x=parameter_x, y=parameter_y):
    """Integrand 2 m^2 / denominator of the dimensionless loop integral."""
    return 2 * mass**2 / feynman_denominator(mass, virtuality_one, virtuality_two,
                                             momentum_squared, x, y)


def loop_integral(mass=quark_mass, virtuality_one=first_virtuality,
                  virtuality_two=second_virtuality, momentum_squared=alp_momentum_squared):
    """Unevaluated dimensionless loop integral over 0 < x < 1, 0 < y < 1 - x."""
    return sp.Integral(loop_integrand(mass, virtuality_one, virtuality_two, momentum_squared),
                       (parameter_y, 0, 1 - parameter_x), (parameter_x, 0, 1))


def loop_integral_numeric(mass, virtuality_one, virtuality_two, momentum_squared):
    """Fixed-precision quadrature of the loop integral; used only by the tests."""
    mpmath.mp.dps = 25

    def inner(x_value):
        return mpmath.quad(lambda y_value: 2 * mass**2 / (
            mass**2 - x_value * y_value * momentum_squared
            - x_value * (1 - x_value - y_value) * virtuality_two
            - y_value * (1 - x_value - y_value) * virtuality_one), [0, 1 - x_value])

    return mpmath.quad(inner, [0, 1])


@lru_cache(maxsize=1)
def derive_pseudoscalar_triangle():
    """Part 1: triangle with the pseudoscalar coupling, both photon orderings.

    Momentum routing of the first ordering: the quark enters the ALP vertex
    with momentum k and leaves with k + p, emits photon 1 (index mu) and then
    photon 2 (index nu). The closed loop is
      -N_c int d^4k/(2 pi)^4 tr[V_a S(k) V_2 S(k + q2) V_1 S(k + p)].
    The second ordering exchanges the two photons.
    """
    loop_momentum = four_vector("k")
    first_photon = four_vector("q1")
    second_photon = four_vector("q2")
    first_polarization = four_vector("eps1")
    second_polarization = four_vector("eps2")
    alp_momentum = first_photon + second_photon

    mass_matrix = quark_mass * IDENTITY
    start = GAMMA5 * (slash(loop_momentum) + mass_matrix)
    after_second = slash(loop_momentum) + slash(second_photon) + mass_matrix
    after_first = slash(loop_momentum) + slash(first_photon) + mass_matrix
    end = slash(loop_momentum) + slash(alp_momentum) + mass_matrix
    numerator_first_ordering = tuple(tuple(sp.expand(sp.trace(
        start * GAMMA[second] * after_second * GAMMA[first] * end))
        for second in range(4)) for first in range(4))
    numerator_second_ordering = tuple(tuple(sp.expand(sp.trace(
        start * GAMMA[first] * after_first * GAMMA[second] * end))
        for second in range(4)) for first in range(4))
    first_lower, second_lower = lower(first_photon), lower(second_photon)
    # 4 i m epsilon^{mu nu rho sigma} q1_rho q2_sigma, independent of k.
    numerator_expected = tuple(tuple(sp.expand(4 * sp.I * quark_mass * sum(
        sp.LeviCivita(first, second, third, fourth) * first_lower[third] * second_lower[fourth]
        for third in range(4) for fourth in range(4)))
        for second in range(4)) for first in range(4))

    # R = epsilon_{mu nu rho sigma} eps1^mu eps2^nu q1^rho q2^sigma with
    # epsilon_{0123} = -1, compared with the contraction of the numerator.
    photon_tensor = sp.expand(-sum(
        sp.LeviCivita(first, second, third, fourth) * first_polarization[first]
        * second_polarization[second] * first_photon[third] * second_photon[fourth]
        for first in range(4) for second in range(4)
        for third in range(4) for fourth in range(4)))
    polarization_one_lower = lower(first_polarization)
    polarization_two_lower = lower(second_polarization)
    contracted_numerator = sp.expand(sum(
        polarization_one_lower[first] * polarization_two_lower[second]
        * numerator_expected[first][second] for first in range(4) for second in range(4)))
    numerator_coefficient = sp.cancel(contracted_numerator / photon_tensor)

    # Feynman parameters: x on the propagator of momentum k (between the
    # photon-2 vertex and the ALP vertex), y on the propagator of momentum
    # k + p (between the ALP vertex and the photon-1 vertex), 1 - x - y on the
    # propagator of momentum k + q2 (between the two photon vertices).
    invariants = (dot(first_photon, first_photon), dot(second_photon, second_photon),
                  dot(alp_momentum, alp_momentum))
    weight_between_photons = 1 - parameter_x - parameter_y
    combined_first_ordering = (
        parameter_x * (dot(loop_momentum, loop_momentum) - quark_mass**2)
        + parameter_y * (dot(loop_momentum + alp_momentum, loop_momentum + alp_momentum)
                         - quark_mass**2)
        + weight_between_photons * (dot(loop_momentum + second_photon,
                                        loop_momentum + second_photon) - quark_mass**2))
    shifted_first = loop_momentum + parameter_y * alp_momentum + weight_between_photons * second_photon
    denominator_first_ordering = feynman_denominator(quark_mass, invariants[0], invariants[1],
                                                     invariants[2])
    combination_check_first = sp.expand(combined_first_ordering - dot(shifted_first, shifted_first)
                                        + denominator_first_ordering)
    combined_second_ordering = (
        parameter_x * (dot(loop_momentum, loop_momentum) - quark_mass**2)
        + parameter_y * (dot(loop_momentum + alp_momentum, loop_momentum + alp_momentum)
                         - quark_mass**2)
        + weight_between_photons * (dot(loop_momentum + first_photon,
                                        loop_momentum + first_photon) - quark_mass**2))
    shifted_second = loop_momentum + parameter_y * alp_momentum + weight_between_photons * first_photon
    denominator_second_ordering = feynman_denominator(quark_mass, invariants[1], invariants[0],
                                                      invariants[2])
    combination_check_second = sp.expand(combined_second_ordering - dot(shifted_second, shifted_second)
                                         + denominator_second_ordering)

    # 1/(A B C) = 2 int_0^1 dx int_0^{1-x} dy [x A + y B + (1-x-y) C]^{-3}.
    first_factor, second_factor, third_factor = sp.symbols("A B C", positive=True)
    combined = (parameter_x * first_factor + parameter_y * second_factor
                + (1 - parameter_x - parameter_y) * third_factor)
    inner_antiderivative = -1 / ((second_factor - third_factor) * combined**2)
    inner_antiderivative_check = sp.simplify(sp.diff(inner_antiderivative, parameter_y) - 2 / combined**3)
    inner_result = sp.simplify(inner_antiderivative.subs(parameter_y, 1 - parameter_x)
                               - inner_antiderivative.subs(parameter_y, 0))
    outer_integrand = sp.simplify(inner_result)
    first_line = parameter_x * first_factor + (1 - parameter_x) * second_factor
    second_line = parameter_x * first_factor + (1 - parameter_x) * third_factor
    outer_antiderivative = (1 / (second_factor - third_factor)) * (
        1 / ((first_factor - second_factor) * first_line)
        - 1 / ((first_factor - third_factor) * second_line))
    outer_antiderivative_check = sp.simplify(sp.diff(outer_antiderivative, parameter_x) - outer_integrand)
    three_factor_identity = sp.simplify(outer_antiderivative.subs(parameter_x, 1)
                                        - outer_antiderivative.subs(parameter_x, 0)
                                        - 1 / (first_factor * second_factor * third_factor))

    # Shifted loop momentum l, rotated to Euclidean components: d^4l = i d^4l_E,
    # (l^2 - Delta)^3 = -(l_E^2 + Delta)^3, d^4l_E = 2 pi^2 l_E^3 dl_E.
    radial = sp.Symbol("l_E", positive=True)
    positive_denominator = sp.Symbol("Delta", positive=True)
    radial_integral = sp.integrate(radial**3 / (radial**2 + positive_denominator)**3, (radial, 0, sp.oo))
    momentum_integral = sp.simplify(sp.I * (-1) * 2 * sp.pi**2 / (2 * sp.pi)**4 * radial_integral)

    # Assembly of i M for one ordering, in units of R times the Feynman-parameter
    # integral of 1/denominator.
    closed_loop_sign = -1
    alp_vertex = -2 * quark_mass * pseudoscalar_coupling / decay_constant
    photon_vertex = -sp.I * electric_charge * quark_charge
    propagator_numerator_factor = sp.I
    one_ordering = (closed_loop_sign * color_number * alp_vertex * photon_vertex**2
                    * propagator_numerator_factor**3 * numerator_coefficient
                    * 2 * momentum_integral * positive_denominator)
    # Both orderings give the same integral after exchanging x and y.
    amplitude_coefficient = 2 * one_ordering
    # i M = i e^2 R F and the loop integral is 2 m^2 int 1/denominator.
    form_factor_per_loop_integral = sp.simplify(
        amplitude_coefficient / (sp.I * electric_charge**2) / (2 * quark_mass**2))
    overall_sign = sp.simplify(form_factor_per_loop_integral
                               / (REFERENCE_FACTOR * pseudoscalar_coupling))
    pseudoscalar_form_factor = form_factor_per_loop_integral * loop_integral()

    # Local photon operator c_gamma (alpha/(4 pi)) (a/f_a) F_{mu nu} Ftilde^{mu nu}
    # with Ftilde^{mu nu} = (1/2) epsilon^{mu nu rho sigma} F_{rho sigma},
    # epsilon^{0123}=+1 and alpha = e^2/(4 pi): two outgoing photons give
    # i M = i c_gamma e^2/(16 pi^2 f_a) * 2 * 2 * i^2 epsilon^{mu nu rho sigma}
    #       q1_mu eps1_nu q2_rho eps2_sigma.
    dual_contraction = sp.expand(sum(
        sp.LeviCivita(first, second, third, fourth) * first_lower[first]
        * polarization_one_lower[second] * second_lower[third] * polarization_two_lower[fourth]
        for first in range(4) for second in range(4)
        for third in range(4) for fourth in range(4)))
    local_photon_form_factor = sp.simplify(
        photon_coupling * electric_charge**2 / (16 * sp.pi**2 * decay_constant)
        * 4 * sp.I**2 * sp.cancel(dual_contraction / photon_tensor) / electric_charge**2)
    return locals()


@lru_cache(maxsize=1)
def derive_loop_integral_properties():
    """Part 1 (continued): symmetry, mass limits and the real-photon closed form."""
    integrand = loop_integrand(quark_mass, first_virtuality, second_virtuality, alp_momentum_squared)
    exchanged = loop_integrand(quark_mass, second_virtuality, first_virtuality, alp_momentum_squared,
                               x=parameter_y, y=parameter_x)
    photon_exchange_difference = sp.simplify(integrand - exchanged)

    # Large mass: 2 m^2/denominator = 2/(1 - N/m^2) with N the momentum-dependent
    # part; the geometric series converges uniformly when all invariants are
    # small compared with m^2, so it can be integrated term by term.
    momentum_part = sp.expand(quark_mass**2 - feynman_denominator(
        quark_mass, first_virtuality, second_virtuality, alp_momentum_squared))

    def triangle_integral(expression):
        return sp.integrate(sp.integrate(expression, (parameter_y, 0, 1 - parameter_x)),
                            (parameter_x, 0, 1))

    heavy_mass_limit = triangle_integral(sp.Integer(2))
    heavy_mass_first_correction = sp.factor(triangle_integral(2 * momentum_part) / quark_mass**2)
    heavy_mass_second_correction = sp.factor(triangle_integral(2 * momentum_part**2) / quark_mass**4)
    geometric_remainder = sp.simplify(
        integrand - 2 - 2 * momentum_part / quark_mass**2 - 2 * momentum_part**2 / quark_mass**4
        - 2 * momentum_part**3 / (quark_mass**4 * (quark_mass**2 - momentum_part)))

    # One nonzero invariant. With the ratio a = invariant/m^2 the loop integral is
    # (2/a) g(a), g(a) = -int_0^1 dx ln(1 - a x (1-x))/x.
    ratio = sp.Symbol("a", real=True)
    # Real photons: q1^2 = q2^2 = 0, the y integration is elementary.
    real_integrand = loop_integrand(quark_mass, 0, 0, ratio * quark_mass**2)
    real_inner_antiderivative = -2 * sp.log(1 - ratio * parameter_x * parameter_y) / (ratio * parameter_x)
    real_inner_check = sp.simplify(sp.diff(real_inner_antiderivative, parameter_y) - real_integrand)
    logarithm_integrand = -sp.log(1 - ratio * parameter_x * (1 - parameter_x)) / parameter_x
    real_inner_result_check = sp.simplify(
        real_inner_antiderivative.subs(parameter_y, 1 - parameter_x)
        - real_inner_antiderivative.subs(parameter_y, 0) - 2 * logarithm_integrand / ratio)
    # dg/da: integrand (1-x)/(1 - a x(1-x)); its average with x -> 1-x is
    # 1/(2 (1 - a x(1-x))).
    derivative_integrand = sp.simplify(sp.diff(logarithm_integrand, ratio))
    symmetric_integrand = 1 / (1 - ratio * parameter_x * (1 - parameter_x))
    symmetrization_check = sp.simplify(
        (derivative_integrand + derivative_integrand.subs(parameter_x, 1 - parameter_x)) / 2
        - symmetric_integrand / 2)

    # Below threshold, 0 < a < 4: a = 4 t^2/(1 + t^2), t > 0, so that
    # tau = 4/a = (1 + t^2)/t^2 and arcsin(1/sqrt(tau)) = arctan(t).
    tangent = sp.Symbol("t", positive=True)
    ratio_below_threshold = 4 * tangent**2 / (1 + tangent**2)
    below_antiderivative = (1 + tangent**2) / (2 * tangent) * sp.atan(tangent * (2 * parameter_x - 1))
    below_antiderivative_check = sp.simplify(
        sp.diff(below_antiderivative, parameter_x)
        - symmetric_integrand.subs(ratio, ratio_below_threshold))
    below_derivative = sp.simplify((below_antiderivative.subs(parameter_x, 1)
                                    - below_antiderivative.subs(parameter_x, 0)) / 2)
    below_claim = 2 * sp.atan(tangent)**2
    below_claim_derivative_check = sp.simplify(
        sp.diff(below_claim, tangent) / sp.diff(ratio_below_threshold, tangent) - below_derivative)
    below_claim_at_zero = below_claim.subs(tangent, 0)
    real_photon_in_tangent = sp.simplify(2 * below_claim / ratio_below_threshold)
    # arcsin(1/sqrt(tau)) = arcsin(t/sqrt(1+t^2)) = arctan(t): equal derivatives, equal at t = 0.
    arcsine_of_tangent = sp.asin(tangent / sp.sqrt(1 + tangent**2))
    arcsine_derivative_check = sp.simplify(sp.diff(arcsine_of_tangent, tangent) - 1 / (1 + tangent**2))
    arcsine_at_zero = arcsine_of_tangent.subs(tangent, 0)

    # Spacelike invariant, a = -(v-1)^2/v with v > 1.
    exponential = sp.Symbol("v", positive=True)
    ratio_spacelike = -(exponential - 1)**2 / exponential
    spacelike_antiderivative = exponential / (exponential**2 - 1) * sp.log(
        (1 + (exponential - 1) * parameter_x) / (exponential - (exponential - 1) * parameter_x))
    spacelike_antiderivative_check = sp.simplify(
        sp.diff(spacelike_antiderivative, parameter_x)
        - symmetric_integrand.subs(ratio, ratio_spacelike))
    spacelike_symmetric_integral = sp.simplify(sp.expand_log(
        spacelike_antiderivative.subs(parameter_x, 1) - spacelike_antiderivative.subs(parameter_x, 0),
        force=True))
    spacelike_claim = -sp.log(exponential)**2 / 2
    spacelike_claim_derivative_check = sp.simplify(
        sp.diff(spacelike_claim, exponential) / sp.diff(ratio_spacelike, exponential)
        - spacelike_symmetric_integral / 2)
    spacelike_claim_at_zero = spacelike_claim.subs(exponential, 1)
    one_spacelike_invariant = sp.simplify(2 * spacelike_claim / ratio_spacelike)

    # Printed real-photon form with tau = 4 m_h^2/m_a^2 >= 1.
    threshold_ratio = 4 * quark_mass**2 / alp_mass**2
    real_photon_closed_form = threshold_ratio * sp.asin(1 / sp.sqrt(threshold_ratio))**2
    return locals()


@lru_cache(maxsize=1)
def derive_derivative_triangle():
    """Part 2: triangle with the derivative coupling from the regulated axial Ward identity."""
    loop_momentum = four_vector("k")
    first_photon = four_vector("q1")
    second_photon = four_vector("q2")
    alp_momentum = first_photon + second_photon
    mass_matrix = quark_mass * IDENTITY

    # pslash gamma5 = (kslash + pslash - m) gamma5 + gamma5 (kslash - m) + 2 m gamma5,
    # and (kslash - m)(kslash + m) = k^2 - m^2, so that between the propagators
    # S(k+p) pslash gamma5 S(k) = i gamma5 S(k) + S(k+p) i gamma5 + 2 m S(k+p) gamma5 S(k).
    vertex_identity = (slash(alp_momentum) * GAMMA5
                       - (slash(loop_momentum) + slash(alp_momentum) - mass_matrix) * GAMMA5
                       - GAMMA5 * (slash(loop_momentum) - mass_matrix)
                       - 2 * quark_mass * GAMMA5).applyfunc(sp.expand)
    inverse_propagator_identity = (
        (slash(loop_momentum) - mass_matrix) * (slash(loop_momentum) + mass_matrix)
        - (dot(loop_momentum, loop_momentum) - quark_mass**2) * IDENTITY).applyfunc(sp.expand)

    # The two terms without the ALP vertex have two propagators each:
    # tr[gamma^nu (kslash+q2slash+m) gamma^mu gamma5 (kslash+m)] and
    # tr[gamma^nu (kslash+q2slash+m) gamma^mu (kslash+pslash+m) gamma5].
    second_lower, first_lower = lower(second_photon), lower(first_photon)
    loop_lower = lower(loop_momentum)
    shifted_lower = lower(loop_momentum + second_photon)
    after_second = slash(loop_momentum) + slash(second_photon) + mass_matrix
    remainder_one = tuple(tuple(sp.expand(sp.trace(
        GAMMA[second] * after_second * GAMMA[first] * GAMMA5 * (slash(loop_momentum) + mass_matrix)))
        for second in range(4)) for first in range(4))
    remainder_two = tuple(tuple(sp.expand(sp.trace(
        GAMMA[second] * after_second * GAMMA[first]
        * (slash(loop_momentum) + slash(alp_momentum) + mass_matrix) * GAMMA5))
        for second in range(4)) for first in range(4))

    def epsilon_two_vectors(first, second, vector_one, vector_two):
        return sp.expand(sum(sp.LeviCivita(first, second, third, fourth)
                             * vector_one[third] * vector_two[fourth]
                             for third in range(4) for fourth in range(4)))

    # Both are independent of the mass and antisymmetric in two momenta, of which
    # only one is external: -4i epsilon^{mu nu rho sigma} k_rho q2_sigma and
    # -4i epsilon^{mu nu rho sigma} (k+q2)_rho q1_sigma.
    remainder_one_expected = tuple(tuple(
        -4 * sp.I * epsilon_two_vectors(first, second, loop_lower, second_lower)
        for second in range(4)) for first in range(4))
    remainder_two_expected = tuple(tuple(
        -4 * sp.I * epsilon_two_vectors(first, second, shifted_lower, first_lower)
        for second in range(4)) for first in range(4))
    # After the Feynman-parameter shift k = l - x q2 (first) or k + q2 = l - x q1
    # (second) the numerators are odd in l. The shift is allowed for the
    # difference of the quark and the Pauli-Villars regulator: the numerators
    # are equal and the difference of the denominators falls as (M^2-m^2)/k^6,
    # so the integrand falls as k^{-5}.
    shifted = four_vector("l")
    shift_substitution_one = {loop_momentum[index]: shifted[index] - parameter_x * second_photon[index]
                              for index in range(4)}
    remainder_one_shifted = tuple(tuple(
        sp.expand(remainder_one[first][second].subs(shift_substitution_one, simultaneous=True))
        for second in range(4)) for first in range(4))
    shift_substitution_two = {
        loop_momentum[index]: shifted[index] - second_photon[index] - parameter_x * first_photon[index]
        for index in range(4)}
    remainder_two_shifted = tuple(tuple(
        sp.expand(remainder_two[first][second].subs(shift_substitution_two, simultaneous=True))
        for second in range(4)) for first in range(4))
    reflection = {component: -component for component in shifted}

    # Regulated amplitude. A vertex (2 m / f_a) gamma5 gives the form factor
    # -(overall sign) (reference factor) (loop integral of mass m); this is the
    # pseudoscalar result with -2 m c^P replaced by 2 m.
    pseudoscalar = derive_pseudoscalar_triangle()
    properties = derive_loop_integral_properties()
    loop_integral_symbol = sp.Symbol("loop_integral_quark")
    loop_integral_regulator = sp.Symbol("loop_integral_regulator")
    unit_mass_vertex_quark = -pseudoscalar["overall_sign"] * REFERENCE_FACTOR * loop_integral_symbol
    unit_mass_vertex_regulator = -pseudoscalar["overall_sign"] * REFERENCE_FACTOR * loop_integral_regulator
    # Derivative vertex (c^A/f_a) pslash gamma5; the regulator loop is subtracted.
    regulated_axial = axial_coupling * (unit_mass_vertex_quark - unit_mass_vertex_regulator)
    axial_form_factor_symbolic = sp.factor(regulated_axial.subs(
        loop_integral_regulator, properties["heavy_mass_limit"]))
    axial_form_factor = axial_form_factor_symbolic.subs(loop_integral_symbol, loop_integral())
    axial_massless_limit = axial_form_factor_symbolic.subs(loop_integral_symbol, 0)
    reference_massless = REFERENCE_FACTOR * axial_coupling
    return locals()


@lru_cache(maxsize=1)
def derive_combined_form_factor():
    """Part 3: sum of both couplings, invariance under the ALP-dependent axial rotation, limits."""
    pseudoscalar = derive_pseudoscalar_triangle()
    derivative = derive_derivative_triangle()
    properties = derive_loop_integral_properties()
    loop_integral_symbol = derivative["loop_integral_symbol"]
    pseudoscalar_symbolic = (pseudoscalar["form_factor_per_loop_integral"] * loop_integral_symbol)
    combined_symbolic = sp.simplify(derivative["axial_form_factor_symbolic"] + pseudoscalar_symbolic)
    expected_symbolic = REFERENCE_FACTOR * (
        axial_coupling - (axial_coupling - pseudoscalar_coupling) * loop_integral_symbol)
    combined_difference = sp.simplify(combined_symbolic - expected_symbolic)
    combined_form_factor = expected_symbolic.subs(loop_integral_symbol, loop_integral(
        momentum_squared=alp_mass**2))

    # Rotation h -> exp(-i kappa a/f_a gamma5) h: c^A -> c^A + kappa,
    # c^P -> c^P + kappa, c_gamma -> c_gamma - 2 N_c Q_h^2 kappa.
    total_symbolic = photon_coupling / (4 * sp.pi**2 * decay_constant) + combined_symbolic
    rotated_total = total_symbolic.subs({
        axial_coupling: axial_coupling + rotation_parameter,
        pseudoscalar_coupling: pseudoscalar_coupling + rotation_parameter,
        photon_coupling: photon_coupling - 2 * color_number * quark_charge**2 * rotation_parameter},
        simultaneous=True)
    rotation_difference = sp.simplify(rotated_total - total_symbolic)

    decoupling_limit = sp.simplify(combined_symbolic.subs(loop_integral_symbol,
                                                          properties["heavy_mass_limit"]))
    massless_limit = sp.simplify(combined_symbolic.subs(loop_integral_symbol, 0))
    real_photon_limit = expected_symbolic.subs(loop_integral_symbol,
                                               properties["real_photon_closed_form"])
    return locals()


@lru_cache(maxsize=1)
def derive_kinematic_limits():
    """Part 4: equal spacelike virtualities, the first ALP-mass correction, one real photon, low energy."""
    properties = derive_loop_integral_properties()
    mass_ratio = sp.Symbol("r", positive=True)  # m_h^2/|q^2|, used only inside this function
    to_physical = {mass_ratio: quark_mass**2 / spacelike_magnitude}

    # (a) q1^2 = q2^2 = -|q^2|, p^2 = 0: the denominator depends on u = x + y only.
    equal_integrand = loop_integrand(quark_mass, -spacelike_magnitude, -spacelike_magnitude, 0)
    in_sum_variable = sp.simplify(equal_integrand.subs(parameter_y, parameter_u - parameter_x))
    independent_of_x = sp.simplify(sp.diff(in_sum_variable, parameter_x))
    one_dimensional_integrand = sp.simplify(sp.integrate(in_sum_variable, (parameter_x, 0, parameter_u)))
    expected_one_dimensional = (2 * quark_mass**2 * parameter_u
                                / (quark_mass**2 + spacelike_magnitude * parameter_u * (1 - parameter_u)))
    one_dimensional_check = sp.simplify(one_dimensional_integrand - expected_one_dimensional)
    # Average with u -> 1-u: m^2/(m^2 + |q^2| u (1-u)).
    symmetric_one_dimensional = quark_mass**2 / (
        quark_mass**2 + spacelike_magnitude * parameter_u * (1 - parameter_u))
    equal_symmetrization_check = sp.simplify(
        (expected_one_dimensional + expected_one_dimensional.subs(parameter_u, 1 - parameter_u)) / 2
        - symmetric_one_dimensional)
    equal_virtuality_integral = sp.Integral(symmetric_one_dimensional, (parameter_u, 0, 1))
    # Closed form from the antiderivative
    # (m^2/(|q^2| root)) ln[(root - 1 + 2u)/(root + 1 - 2u)], root = sqrt(1 + 4 m^2/|q^2|).
    root = sp.sqrt(1 + 4 * mass_ratio)
    equal_antiderivative = mass_ratio / root * sp.log((root - 1 + 2 * parameter_u)
                                                      / (root + 1 - 2 * parameter_u))
    equal_antiderivative_check = sp.simplify(
        sp.diff(equal_antiderivative, parameter_u)
        - mass_ratio / (mass_ratio + parameter_u * (1 - parameter_u)))
    equal_closed_ratio = 2 * mass_ratio / root * sp.log((root + 1) / (root - 1))
    logarithm_of_roots = sp.Symbol("logarithm_of_roots")  # placeholder used only to collect terms
    equal_endpoint_check = sp.simplify(sp.expand_log(
        equal_antiderivative.subs(parameter_u, 1) - equal_antiderivative.subs(parameter_u, 0)
        - equal_closed_ratio, force=True))
    equal_closed_form = equal_closed_ratio.subs(to_physical)
    logarithm = sp.Symbol("ln_inverse_ratio")  # stands for ln(|q^2|/m_h^2) inside the expansions

    def expand_small_ratio(expression, order):
        series = sp.series(expression, mass_ratio, 0, order).removeO()
        with_symbol = sp.expand(sp.expand_log(series, force=True).subs(sp.log(mass_ratio), -logarithm))
        # Keep powers of the ratio below the requested order, logarithms included.
        return sp.expand(sum(with_symbol.coeff(mass_ratio, power) * mass_ratio**power
                             for power in range(order)))

    physical_logarithm = sp.log(spacelike_magnitude / quark_mass**2)

    def physical(expression):
        return expression.subs(logarithm, physical_logarithm).subs(to_physical)

    equal_expansion_ratio = expand_small_ratio(equal_closed_ratio, 3)
    equal_expansion = physical(equal_expansion_ratio)
    equal_massless_limit = sp.limit(equal_closed_ratio, mass_ratio, 0, "+")
    equal_heavy_expansion = sp.series(equal_closed_ratio.subs(mass_ratio, 1 / sp.Symbol("s", positive=True)),
                                      sp.Symbol("s", positive=True), 0, 2).removeO()

    # (b) first order in p^2 = m_a^2 at equal spacelike virtualities:
    # d(integrand)/d(p^2) at p^2 = 0, reduced to u = x + y.
    first_order_integrand = sp.diff(loop_integrand(
        quark_mass, -spacelike_magnitude, -spacelike_magnitude, alp_momentum_squared),
        alp_momentum_squared).subs(alp_momentum_squared, 0)
    first_order_in_sum = sp.simplify(first_order_integrand.subs(parameter_y, parameter_u - parameter_x))
    first_order_one_dimensional = sp.factor(sp.integrate(first_order_in_sum, (parameter_x, 0, parameter_u)))
    first_order_expected = (quark_mass**2 * parameter_u**3
                            / (3 * (quark_mass**2 + spacelike_magnitude * parameter_u * (1 - parameter_u))**2))
    first_order_check = sp.simplify(first_order_one_dimensional - first_order_expected)
    alp_mass_correction_integral = alp_mass**2 * sp.Integral(first_order_expected, (parameter_u, 0, 1))
    # In units of 1/|q^2|: (r/3) u^3/(r + u(1-u))^2. Its average with u -> 1-u is
    # (r/6) [(1+3r)/(r+u(1-u))^2 - 3/(r+u(1-u))], and
    # int 1/(r+u(1-u))^2 = -d/dr int 1/(r+u(1-u)) = -d/dr (equal-virtuality result / r).
    scaled_first_order = mass_ratio * parameter_u**3 / (3 * (mass_ratio + parameter_u * (1 - parameter_u))**2)
    scaled_symmetric = mass_ratio / 6 * (
        (1 + 3 * mass_ratio) / (mass_ratio + parameter_u * (1 - parameter_u))**2
        - 3 / (mass_ratio + parameter_u * (1 - parameter_u)))
    first_order_symmetrization_check = sp.simplify(
        (scaled_first_order + scaled_first_order.subs(parameter_u, 1 - parameter_u)) / 2 - scaled_symmetric)
    scaling_check = sp.simplify(
        (first_order_expected * spacelike_magnitude).subs(
            quark_mass, sp.sqrt(mass_ratio * spacelike_magnitude)) - scaled_first_order)
    alp_mass_coefficient_ratio = sp.collect(sp.expand(sp.simplify(mass_ratio / 6 * (
        -(1 + 3 * mass_ratio) * sp.diff(equal_closed_ratio / mass_ratio, mass_ratio)
        - 3 * equal_closed_ratio / mass_ratio)
        .subs(sp.log((root + 1) / (root - 1)), logarithm_of_roots))), logarithm_of_roots).subs(
            logarithm_of_roots, sp.log((root + 1) / (root - 1)))
    alp_mass_coefficient_expansion_ratio = expand_small_ratio(alp_mass_coefficient_ratio, 2)
    alp_mass_correction_closed = alp_mass**2 / spacelike_magnitude * alp_mass_coefficient_ratio.subs(to_physical)
    alp_mass_correction_expansion = (alp_mass**2 / spacelike_magnitude
                                     * physical(alp_mass_coefficient_expansion_ratio))

    # (c) q1^2 = -|q^2|, q2^2 = 0, p^2 = 0: the x integration is elementary and
    # gives the function of one spacelike invariant of Part 1.
    single_integrand = loop_integrand(quark_mass, -spacelike_magnitude, 0, 0)
    single_inner_antiderivative = -2 * quark_mass**2 / (spacelike_magnitude * parameter_y) * sp.log(
        quark_mass**2 + spacelike_magnitude * parameter_y * (1 - parameter_x - parameter_y))
    single_inner_check = sp.simplify(sp.diff(single_inner_antiderivative, parameter_x) - single_integrand)
    ratio = properties["ratio"]
    single_inner_result_check = sp.simplify(sp.expand_log(
        single_inner_antiderivative.subs(parameter_x, 1 - parameter_y)
        - single_inner_antiderivative.subs(parameter_x, 0)
        - (2 * properties["logarithm_integrand"] / ratio).subs(
            {ratio: -spacelike_magnitude / quark_mass**2, parameter_x: parameter_y}, simultaneous=True),
        force=True))
    single_closed_ratio = mass_ratio * sp.log((root + 1) / (root - 1))**2
    exponential = properties["exponential"]
    exponential_of_ratio = (root + 1) / (root - 1)
    single_matches_one_invariant = sp.simplify(
        properties["one_spacelike_invariant"]
        - (exponential / (exponential - 1)**2) * sp.log(exponential)**2)
    ratio_matches = sp.simplify(exponential_of_ratio / (exponential_of_ratio - 1)**2 - mass_ratio)
    single_closed_form = single_closed_ratio.subs(to_physical)
    single_expansion_ratio = expand_small_ratio(single_closed_ratio, 3)
    single_expansion = physical(single_expansion_ratio)
    single_massless_limit = sp.limit(single_closed_ratio, mass_ratio, 0, "+")

    # (d) low energy: first two terms of the expansion in inverse powers of m_h^2.
    low_energy_expansion = (properties["heavy_mass_limit"] + properties["heavy_mass_first_correction"])
    low_energy_next_term = properties["heavy_mass_second_correction"]
    return locals()


@lru_cache(maxsize=1)
def derive_axial_current_matrix_element():
    """Part 5: one-loop vacuum-to-ALP matrix element of hbar gamma^sigma gamma5 h and the local expansion.

    <0|hbar gamma^sigma gamma5 h(0)|a(p)> =
      -N_c int d^4k/(2 pi)^4 tr[gamma^sigma gamma5 S(k+p) V_a S(k)],
    with V_a = (c^A/f_a) pslash gamma5 - (2 m c^P/f_a) gamma5. No photon is
    attached, so there is no anomaly term in this two-point function.
    """
    loop_momentum = four_vector("k")
    alp_momentum = four_vector("p")
    mass_matrix = quark_mass * IDENTITY
    first_denominator = dot(loop_momentum, loop_momentum) - quark_mass**2
    second_denominator = dot(loop_momentum + alp_momentum, loop_momentum + alp_momentum) - quark_mass**2
    outgoing = slash(loop_momentum) + slash(alp_momentum) + mass_matrix
    incoming = slash(loop_momentum) + mass_matrix
    derivative_numerator = tuple(sp.expand(sp.trace(
        GAMMA[index] * GAMMA5 * outgoing * slash(alp_momentum) * GAMMA5 * incoming)) for index in range(4))
    # 4 k^sigma D(k+p) - 4 (k+p)^sigma D(k) - 8 m^2 p^sigma.
    derivative_numerator_expected = tuple(sp.expand(
        4 * loop_momentum[index] * second_denominator
        - 4 * (loop_momentum[index] + alp_momentum[index]) * first_denominator
        - 8 * quark_mass**2 * alp_momentum[index]) for index in range(4))
    pseudoscalar_numerator = tuple(sp.expand(sp.trace(
        GAMMA[index] * GAMMA5 * outgoing * GAMMA5 * incoming)) for index in range(4))
    pseudoscalar_numerator_expected = tuple(sp.expand(-4 * quark_mass * alp_momentum[index])
                                            for index in range(4))

    # Scalar integral: 1/(D(k) D(k+p)) = int_0^1 dx [(k + x p)^2 - (m^2 - x(1-x)p^2)]^{-2}.
    two_point_combination_check = sp.expand(
        parameter_x * second_denominator + (1 - parameter_x) * first_denominator
        - dot(loop_momentum + parameter_x * alp_momentum, loop_momentum + parameter_x * alp_momentum)
        + quark_mass**2 - parameter_x * (1 - parameter_x) * dot(alp_momentum, alp_momentum))
    radial = sp.Symbol("l_E", positive=True)
    positive_denominator = sp.Symbol("Delta", positive=True)
    cut_radial_integral = sp.integrate(radial**3 / (radial**2 + positive_denominator)**2, (radial, 0, cutoff))
    large_cutoff_form = (sp.log(cutoff**2 / positive_denominator) - 1) / 2
    cutoff_remainder = sp.limit(sp.simplify(cut_radial_integral - large_cutoff_form), cutoff, sp.oo)
    # d^4l = i d^4l_E, (l^2 - Delta)^2 = (l_E^2 + Delta)^2.
    scalar_integral_factor = sp.I * 2 * sp.pi**2 / (2 * sp.pi)**4
    two_point_denominator = quark_mass**2 - parameter_x * (1 - parameter_x) * alp_mass**2
    scalar_integrand = 2 * large_cutoff_form.subs(positive_denominator, two_point_denominator)
    scalar_integral = sp.Integral(scalar_integrand, (parameter_x, 0, 1))
    # Taylor expansion in m_a^2/m_h^2 to first order, integrated term by term.
    mass_ratio_squared = sp.Symbol("b", positive=True)  # m_a^2/m_h^2 inside this function
    scalar_integrand_expanded = (sp.log(cutoff**2 / quark_mass**2) - 1
                                 - sp.series(sp.log(1 - parameter_x * (1 - parameter_x) * mass_ratio_squared),
                                             mass_ratio_squared, 0, 2).removeO())
    scalar_integral_expanded = sp.expand(sp.integrate(scalar_integrand_expanded, (parameter_x, 0, 1)).subs(
        mass_ratio_squared, alp_mass**2 / quark_mass**2))

    # Coefficient of p^sigma: -N_c * i^2 * [(c^A/f_a)(-8 m^2) + (-2 m c^P/f_a)(-4 m)]
    # times the scalar integral; the terms k^sigma/D are set to zero (assumption).
    closed_loop_sign = -1
    propagator_numerator_factor = sp.I
    vertex_sum = (axial_coupling / decay_constant * (-8 * quark_mass**2)
                  + (-2 * quark_mass * pseudoscalar_coupling / decay_constant) * (-4 * quark_mass))
    matrix_element_coefficient_per_integral = sp.factor(
        closed_loop_sign * color_number * propagator_numerator_factor**2 * vertex_sum
        * scalar_integral_factor / 2)
    matrix_element_coefficient = matrix_element_coefficient_per_integral * scalar_integral
    matrix_element_coefficient_expanded = matrix_element_coefficient_per_integral * scalar_integral_expanded

    # Local expansion: F = 2 i Q_h^2 p_sigma <...> / ((q^2 - m_h^2) m_a^2) with
    # <...> = p^sigma * coefficient, p^2 = m_a^2 and q^2 = -|q^2|.
    local_per_integral = sp.factor(2 * sp.I * quark_charge**2 * matrix_element_coefficient_per_integral
                                   / (-spacelike_magnitude - quark_mass**2))
    local_form_factor = local_per_integral * scalar_integral
    local_leading = sp.expand(sp.limit(local_per_integral * spacelike_magnitude, spacelike_magnitude, sp.oo)
                              / spacelike_magnitude * scalar_integral_expanded)

    # Direct triangle: mass-dependent part -(reference factor)(c^A - c^P) x loop
    # integral, with the leading terms of Part 4 (a) and (b) for |q^2| >> m_h^2.
    limits = derive_kinematic_limits()
    mass_ratio = limits["mass_ratio"]
    leading_loop_integral = (
        limits["equal_expansion_ratio"].coeff(mass_ratio, 1) * mass_ratio
        + alp_mass**2 / spacelike_magnitude * limits["alp_mass_coefficient_expansion_ratio"].coeff(mass_ratio, 0))
    leading_loop_integral = limits["physical"](leading_loop_integral)
    direct_leading = sp.expand(-REFERENCE_FACTOR * (axial_coupling - pseudoscalar_coupling)
                               * leading_loop_integral)
    remainder = sp.simplify(sp.expand_log(direct_leading - local_leading, force=True))
    remainder_expected = (-REFERENCE_FACTOR * (axial_coupling - pseudoscalar_coupling)
                          * 2 * quark_mass**2 / spacelike_magnitude
                          * (sp.log(spacelike_magnitude / cutoff**2) + 1))
    remainder_check = sp.simplify(sp.expand_log(remainder - remainder_expected, force=True))
    return locals()


def loop_integral_complex_mass_numeric(mass_squared, virtuality_one, virtuality_two, momentum_squared):
    """Quadrature of the loop integral at complex m^2 (15 digits); used only by the tests."""
    with mpmath.workdps(15):
        def inner(x_value):
            return mpmath.quad(lambda y_value: 2 * mass_squared / (
                mass_squared - x_value * y_value * momentum_squared
                - x_value * (1 - x_value - y_value) * virtuality_two
                - y_value * (1 - x_value - y_value) * virtuality_one), [0, 1 - x_value])

        return mpmath.quad(inner, [0, 1])


def extrapolate_to_zero(nodes, values):
    """Value at zero of the polynomial through (nodes, values) (Neville scheme); used only by the tests."""
    table = list(values)
    for level in range(1, len(nodes)):
        table = [(nodes[index + level] * table[index] - nodes[index] * table[index + 1])
                 / (nodes[index + level] - nodes[index]) for index in range(len(nodes) - level)]
    return table[0]


def loop_integral_real_mass_limit(virtuality_one, virtuality_two, momentum_squared,
                                  imaginary_fractions=(0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)):
    """Two-dimensional quadratures at m^2 = (1 - i eta) GeV^2 extrapolated to eta = 0.

    The loop integral is analytic in m^2 for Im m^2 < 0, so the values at
    eta > 0 lie on a Taylor series around eta = 0 as long as eta is smaller
    than the distances to the threshold and to m^2 = 0, that is
    eta < min(p^2/(4 GeV^2) - 1, 1). The polynomial through the eight default
    values reproduces the value at eta = 0 to about 1e-6 (measured, not a bound).
    Returns the extrapolated value and the list of (eta, quadrature) pairs.
    """
    nodes = [mpmath.mpf(str(fraction)) for fraction in imaginary_fractions]
    values = [loop_integral_complex_mass_numeric(1 - mpmath.mpc(0, 1) * node, virtuality_one,
                                                 virtuality_two, momentum_squared) for node in nodes]
    return extrapolate_to_zero(nodes, values), list(zip(nodes, values))


def levi_civita_coefficient(numerator, first_photon, second_photon):
    """Coefficient of epsilon^{mu nu rho sigma} q1_rho q2_sigma in a numerator tensor N^{mu nu}.

    N^{mu nu} is contracted with epsilon_{mu nu alpha beta} q1^alpha q2^beta
    (epsilon_{0123} = -1) and divided by the contraction of the structure with
    itself, -2 (q1^2 q2^2 - (q1.q2)^2). The metric tensor and products
    a^mu b^nu of the three momenta k, q1, q2 drop out of the contraction,
    because epsilon_{mu nu alpha beta} a^mu b^nu q1^alpha q2^beta vanishes
    when a and b are taken from k, q1, q2.
    """
    contraction = sp.expand(-sum(
        sp.LeviCivita(first, second, third, fourth) * numerator[first][second]
        * first_photon[third] * second_photon[fourth]
        for first in range(4) for second in range(4)
        for third in range(4) for fourth in range(4)))
    normalization = sp.expand(-2 * (dot(first_photon, first_photon) * dot(second_photon, second_photon)
                                    - dot(first_photon, second_photon)**2))
    return sp.cancel(contraction / normalization)


@lru_cache(maxsize=1)
def derive_vector_triangle():
    """Part 6: triangle with the vector coupling (d_mu a/f_a) c_hh^V hbar gamma^mu h.

    Vertex with an incoming ALP: (c_hh^V/f_a) pslash. Same momentum routing as
    in Part 1; the numerators are
      first ordering:  tr[pslash (kslash+m) gamma^nu (kslash+q2slash+m) gamma^mu (kslash+pslash+m)],
      second ordering: tr[pslash (kslash+m) gamma^mu (kslash+q1slash+m) gamma^nu (kslash+pslash+m)].
    """
    loop_momentum = four_vector("k")
    first_photon = four_vector("q1")
    second_photon = four_vector("q2")
    alp_momentum = first_photon + second_photon
    mass_matrix = quark_mass * IDENTITY

    start = (slash(alp_momentum) * (slash(loop_momentum) + mass_matrix)).applyfunc(sp.expand)
    after_second = slash(loop_momentum) + slash(second_photon) + mass_matrix
    after_first = slash(loop_momentum) + slash(first_photon) + mass_matrix
    end = slash(loop_momentum) + slash(alp_momentum) + mass_matrix
    numerator_first_ordering = tuple(tuple(sp.expand(sp.trace(
        start * GAMMA[second] * after_second * GAMMA[first] * end))
        for second in range(4)) for first in range(4))
    numerator_second_ordering = tuple(tuple(sp.expand(sp.trace(
        start * GAMMA[first] * after_first * GAMMA[second] * end))
        for second in range(4)) for first in range(4))

    # (i) No Levi-Civita tensor. Under the reflection of the space components of
    # k, q1, q2 a tensor built from the metric and the momenta obeys
    # N^{mu nu} -> s_mu s_nu N^{mu nu} (s_0 = +1, s_i = -1), while any term with
    # one Levi-Civita tensor obeys N^{mu nu} -> -s_mu s_nu N^{mu nu}.
    space_reflection = {vector[index]: -vector[index]
                        for vector in (loop_momentum, first_photon, second_photon) for index in (1, 2, 3)}

    def index_sign(index):
        return 1 if index == 0 else -1

    def reflection_odd_part(numerator):
        return tuple(tuple(sp.expand(
            (numerator[first][second] - index_sign(first) * index_sign(second)
             * numerator[first][second].subs(space_reflection, simultaneous=True)) / 2)
            for second in range(4)) for first in range(4))

    def reflection_even_part(numerator):
        return tuple(tuple(sp.expand(
            (numerator[first][second] + index_sign(first) * index_sign(second)
             * numerator[first][second].subs(space_reflection, simultaneous=True)) / 2)
            for second in range(4)) for first in range(4))

    odd_part_first_ordering = reflection_odd_part(numerator_first_ordering)
    odd_part_second_ordering = reflection_odd_part(numerator_second_ordering)
    levi_civita_coefficient_first = levi_civita_coefficient(numerator_first_ordering, first_photon, second_photon)
    levi_civita_coefficient_second = levi_civita_coefficient(numerator_second_ordering, first_photon, second_photon)
    # The same two operations applied to the pseudoscalar numerator of Part 1,
    # to show that they do detect the Levi-Civita structure.
    pseudoscalar = derive_pseudoscalar_triangle()
    pseudoscalar_even_part = reflection_even_part(pseudoscalar["numerator_first_ordering"])
    pseudoscalar_levi_civita_coefficient = levi_civita_coefficient(
        pseudoscalar["numerator_first_ordering"], first_photon, second_photon)
    numerators_are_real = not any(
        entry.has(sp.I) for numerator in (numerator_first_ordering, numerator_second_ordering)
        for row in numerator for entry in row)

    # (ii) Charge conjugation. C = i gamma^2 gamma^0 obeys
    # C gamma^mu C^{-1} = -(gamma^mu)^T and C gamma5 C^{-1} = +(gamma5)^T. The
    # trace of the transposed product, with three vector vertices, gives
    # N_first^{mu nu}(k) = -N_second^{mu nu}(-k - p).
    charge_conjugation = sp.I * GAMMA[2] * GAMMA[0]
    charge_conjugation_inverse = charge_conjugation.inv()
    vector_conjugation = tuple(
        (charge_conjugation * GAMMA[index] * charge_conjugation_inverse + GAMMA[index].T).applyfunc(sp.expand)
        for index in range(4))
    pseudoscalar_conjugation = (charge_conjugation * GAMMA5 * charge_conjugation_inverse
                                - GAMMA5.T).applyfunc(sp.expand)
    momentum_reversal = {loop_momentum[index]: -loop_momentum[index] - alp_momentum[index]
                         for index in range(4)}
    ordering_sum = tuple(tuple(sp.expand(
        numerator_first_ordering[first][second]
        + numerator_second_ordering[first][second].subs(momentum_reversal, simultaneous=True))
        for second in range(4)) for first in range(4))
    # With gamma5 at the ALP vertex the sign is opposite: the orderings add.
    pseudoscalar_ordering_difference = tuple(tuple(sp.expand(
        pseudoscalar["numerator_first_ordering"][first][second]
        - pseudoscalar["numerator_second_ordering"][first][second].subs(momentum_reversal, simultaneous=True))
        for second in range(4)) for first in range(4))

    # The three propagator denominators of the second ordering, taken at
    # k -> -k - p, are those of the first ordering.
    def propagator_denominator(momentum):
        return dot(momentum, momentum) - quark_mass**2

    denominators_first_ordering = (
        propagator_denominator(loop_momentum), propagator_denominator(loop_momentum + second_photon),
        propagator_denominator(loop_momentum + alp_momentum))
    denominators_second_ordering = (
        propagator_denominator(loop_momentum), propagator_denominator(loop_momentum + first_photon),
        propagator_denominator(loop_momentum + alp_momentum))
    denominator_difference = sp.expand(
        sp.Mul(*denominators_first_ordering)
        - sp.Mul(*denominators_second_ordering).subs(momentum_reversal, simultaneous=True))
    vector_form_factor = sp.Integer(0)
    return locals()


@lru_cache(maxsize=1)
def derive_above_threshold():
    """Part 7: the loop integral for m_a > 2 m_h, with m_h^2 -> m_h^2 - i0 in the denominator.

    Real photons: with a = m_a^2/m_h^2 > 4 the loop integral is (2/a) g(a) and,
    from Part 1, dg/da = (1/2) int_0^1 dx/(1 - a x(1-x) - i0). The parameter
    w > 1 inside this function is defined by a = (w+1)^2/w, that is
    w = (1 + sqrt(1 - 4 m_h^2/m_a^2))/(1 - sqrt(1 - 4 m_h^2/m_a^2)).
    """
    properties = derive_loop_integral_properties()
    ratio = properties["ratio"]
    symmetric_integrand = properties["symmetric_integrand"]
    logarithm_argument = sp.Symbol("w", positive=True)
    ratio_above_threshold = (logarithm_argument + 1)**2 / logarithm_argument

    # Zeros of 1 - a x(1-x) inside 0 < x < 1 and the slopes there.
    parameter_denominator = (1 - ratio * parameter_x * (1 - parameter_x)).subs(ratio, ratio_above_threshold)
    lower_zero = 1 / (logarithm_argument + 1)
    upper_zero = logarithm_argument / (logarithm_argument + 1)
    zero_checks = (sp.simplify(parameter_denominator.subs(parameter_x, lower_zero)),
                   sp.simplify(parameter_denominator.subs(parameter_x, upper_zero)))
    slope_at_lower_zero = sp.simplify(sp.diff(parameter_denominator, parameter_x).subs(parameter_x, lower_zero))
    slope_at_upper_zero = sp.simplify(sp.diff(parameter_denominator, parameter_x).subs(parameter_x, upper_zero))
    # For w > 1 the slope is negative at the lower zero and positive at the upper one.
    slope_magnitude = (logarithm_argument**2 - 1) / logarithm_argument
    slope_checks = (sp.simplify(slope_at_lower_zero + slope_magnitude),
                    sp.simplify(slope_at_upper_zero - slope_magnitude))
    # 1/(D - i0) = principal value 1/D + i pi delta(D); int delta(D) dx = sum over zeros of 1/|slope|.
    delta_function_integral = 2 / slope_magnitude
    # Principal value: antiderivative (w/(w^2-1)) ln|((w+1)x - w)/((w+1)x - 1)|; the
    # logarithmic terms cancel between the two sides of each zero when equal
    # intervals are removed around it (tested numerically).
    principal_value_antiderivative = logarithm_argument / (logarithm_argument**2 - 1) * sp.log(
        ((logarithm_argument + 1) * parameter_x - logarithm_argument)
        / ((logarithm_argument + 1) * parameter_x - 1))
    principal_value_antiderivative_check = sp.simplify(
        sp.diff(principal_value_antiderivative, parameter_x)
        - symmetric_integrand.subs(ratio, ratio_above_threshold))
    principal_value = sp.simplify(sp.expand_log(
        principal_value_antiderivative.subs(parameter_x, 1)
        - principal_value_antiderivative.subs(parameter_x, 0), force=True))
    derivative_above_threshold = (principal_value + sp.I * sp.pi * delta_function_integral) / 2

    # Claim g = -(1/2) (ln w - i pi)^2: equal derivative, equal value at the
    # threshold w = 1 (a = 4), where the form below threshold gives 2 arcsin(1)^2.
    # dg/da behaves as 1/sqrt(a - 4) at the threshold, which is integrable.
    above_claim = -(sp.log(logarithm_argument) - sp.I * sp.pi)**2 / 2
    above_claim_derivative_check = sp.simplify(
        sp.diff(above_claim, logarithm_argument) / sp.diff(ratio_above_threshold, logarithm_argument)
        - derivative_above_threshold)
    tangent = properties["tangent"]
    below_value_at_threshold = sp.limit(properties["below_claim"], tangent, sp.oo)
    above_claim_at_threshold_check = sp.simplify(above_claim.subs(logarithm_argument, 1)
                                                 - below_value_at_threshold)
    # The same function is the form for one spacelike invariant, a = -(v-1)^2/v,
    # g = -(1/2) ln(v)^2, continued to v = w exp(-i pi).
    exponential = properties["exponential"]
    spacelike_continuation_check = (
        sp.simplify(properties["ratio_spacelike"].subs(exponential, -logarithm_argument) - ratio_above_threshold),
        sp.expand(properties["spacelike_claim"].subs(sp.log(exponential), sp.log(logarithm_argument) - sp.I * sp.pi)
                  - above_claim))
    above_in_logarithm_argument = sp.simplify(2 * above_claim / ratio_above_threshold)
    massless_limit_above_threshold = sp.limit(above_in_logarithm_argument, logarithm_argument, sp.oo)

    # Printed form in m_h and m_a: 4 m_h^2/m_a^2 = 4 w/(w+1)^2 and
    # sqrt(1 - 4 m_h^2/m_a^2) = (w-1)/(w+1) (the quark velocity in the ALP rest frame).
    def closed_form_parts(threshold_ratio, quark_velocity):
        logarithm = sp.log((1 + quark_velocity) / (1 - quark_velocity))
        complete = threshold_ratio * (sp.pi / 2 + sp.I * logarithm / 2)**2
        real_part = threshold_ratio * (sp.pi**2 - logarithm**2) / 4
        imaginary_part = sp.pi * threshold_ratio * logarithm / 2
        return complete, real_part, imaginary_part

    threshold_ratio_in_argument = 4 * logarithm_argument / (logarithm_argument + 1)**2
    velocity_in_argument = (logarithm_argument - 1) / (logarithm_argument + 1)
    velocity_relation_check = sp.simplify(velocity_in_argument**2 - 1 + threshold_ratio_in_argument)
    in_argument = closed_form_parts(threshold_ratio_in_argument, velocity_in_argument)
    closed_form_check = sp.simplify(sp.expand(sp.expand_log(in_argument[0], force=True)
                                              - above_in_logarithm_argument))
    parts_check = sp.expand(in_argument[0] - in_argument[1] - sp.I * in_argument[2])
    threshold_ratio = 4 * quark_mass**2 / alp_mass**2
    above_closed_form, above_real_part, above_imaginary_part = closed_form_parts(
        threshold_ratio, sp.sqrt(1 - threshold_ratio))

    # Threshold. distance = sqrt(|4 m_h^2/m_a^2 - 1|), used only inside this function.
    distance = sp.Symbol("d", positive=True)
    arcsine_form = sp.asin(1 / sp.sqrt(1 + distance**2))
    arcsine_form_checks = (sp.simplify(sp.diff(arcsine_form, distance) + 1 / (1 + distance**2)),
                           arcsine_form.subs(distance, 0) - sp.pi / 2)
    half_logarithm = sp.log((1 + distance) / (1 - distance)) / 2
    half_logarithm_checks = (sp.simplify(sp.diff(half_logarithm - sp.atanh(distance), distance)),
                             half_logarithm.subs(distance, 0))
    # Below threshold, 4 m_h^2/m_a^2 = 1 + d^2; above threshold, 4 m_h^2/m_a^2 = 1 - d^2.
    below_near_threshold = (1 + distance**2) * (sp.pi / 2 - sp.atan(distance))**2
    above_near_threshold = (1 - distance**2) * (sp.pi / 2 + sp.I * sp.atanh(distance))**2
    # m_h^2 -> m_h^2 - i0 turns sqrt(4 m_h^2/m_a^2 - 1) into -i sqrt(1 - 4 m_h^2/m_a^2).
    continuation_check = sp.simplify(below_near_threshold.subs(distance, -sp.I * distance)
                                     - above_near_threshold)
    below_series = sp.series(below_near_threshold, distance, 0, 3).removeO()
    above_series = sp.series(above_near_threshold, distance, 0, 3).removeO()
    value_at_threshold = (below_near_threshold.subs(distance, 0), above_near_threshold.subs(distance, 0))
    # Derivative with respect to 4 m_h^2/m_a^2, multiplied by the distance.
    below_slope_times_distance = sp.limit(
        sp.diff(below_near_threshold, distance) / (2 * distance) * distance, distance, 0, "+")
    above_slope = sp.diff(above_near_threshold, distance) / (-2 * distance)
    above_slope_real = sp.simplify(sp.re(sp.expand(above_slope)))
    above_slope_imaginary = sp.simplify(sp.im(sp.expand(above_slope)))
    above_real_slope_at_threshold = sp.limit(above_slope_real, distance, 0, "+")
    above_imaginary_slope_times_distance = sp.limit(above_slope_imaginary * distance, distance, 0, "+")

    # Spacelike photons, q1^2 = q2^2 = -|q^2|, p^2 = m_a^2 > 4 m_h^2. With
    # u = x + y the denominator is m_h^2 + |q^2| u(1-u) - m_a^2 x (u - x).
    spacelike_denominator = feynman_denominator(quark_mass, -spacelike_magnitude, -spacelike_magnitude,
                                                alp_mass**2)
    photon_part = quark_mass**2 + spacelike_magnitude * parameter_u * (1 - parameter_u)
    denominator_in_sum_variable = photon_part - alp_mass**2 * parameter_x * (parameter_u - parameter_x)
    sum_variable_check = sp.expand(spacelike_denominator.subs(parameter_y, parameter_u - parameter_x)
                                   - denominator_in_sum_variable)
    # Zeros in x at fixed u: x = (u -+ sqrt(u^2 - 4 photon_part/m_a^2))/2, both in
    # 0 < x < u when the square root is real (their sum is u, their product is
    # positive). The slope there is -+ sqrt(m_a^4 u^2 - 4 m_a^2 photon_part).
    slope_squared = alp_mass**4 * parameter_u**2 - 4 * alp_mass**2 * photon_part
    zeros_in_x = tuple((parameter_u + sign * sp.sqrt(slope_squared) / alp_mass**2) / 2 for sign in (-1, 1))
    spacelike_zero_checks = tuple(sp.simplify(sp.expand(denominator_in_sum_variable.subs(parameter_x, zero)))
                                  for zero in zeros_in_x)
    spacelike_slope_checks = tuple(sp.simplify(sp.expand(
        sp.diff(denominator_in_sum_variable, parameter_x).subs(parameter_x, zero)**2 - slope_squared))
        for zero in zeros_in_x)
    # Im = pi int du int dx 2 m_h^2 delta(denominator): two zeros, each 1/|slope|.
    spacelike_imaginary_integrand = sp.pi * 2 * quark_mass**2 * 2 / sp.sqrt(slope_squared)
    combined_scale = alp_mass**2 + 4 * spacelike_magnitude
    slope_squared_factored = alp_mass**2 * (combined_scale * parameter_u**2
                                            - 4 * spacelike_magnitude * parameter_u - 4 * quark_mass**2)
    slope_squared_check = sp.expand(slope_squared - slope_squared_factored)
    lower_limit = (2 * spacelike_magnitude + 2 * sp.sqrt(
        spacelike_magnitude**2 + quark_mass**2 * combined_scale)) / combined_scale
    lower_limit_check = sp.simplify(sp.expand(slope_squared_factored.subs(parameter_u, lower_limit)))
    slope_squared_at_one = sp.factor(slope_squared_factored.subs(parameter_u, 1))
    spacelike_imaginary_integral = sp.Integral(
        4 * sp.pi * quark_mass**2 / sp.sqrt(slope_squared_factored), (parameter_u, lower_limit, 1))
    spacelike_integrand_check = sp.simplify(
        spacelike_imaginary_integrand**2 - (4 * sp.pi * quark_mass**2)**2 / slope_squared_factored)
    quadratic = combined_scale * parameter_u**2 - 4 * spacelike_magnitude * parameter_u - 4 * quark_mass**2
    spacelike_antiderivative_argument = (2 * combined_scale * parameter_u - 4 * spacelike_magnitude
                                         + 2 * sp.sqrt(combined_scale) * sp.sqrt(quadratic))
    spacelike_imaginary_antiderivative = (4 * sp.pi * quark_mass**2 / (alp_mass * sp.sqrt(combined_scale))
                                          * sp.log(spacelike_antiderivative_argument))
    # d/du ln(argument) = sqrt(m_a^2 + 4|q^2|)/sqrt(quadratic), written without quotients.
    spacelike_imaginary_antiderivative_check = sp.expand(
        sp.diff(spacelike_antiderivative_argument, parameter_u) * sp.sqrt(quadratic)
        - sp.sqrt(combined_scale) * spacelike_antiderivative_argument)
    argument_at_one = (2 * alp_mass**2 + 4 * spacelike_magnitude
                       + 2 * sp.sqrt(combined_scale) * sp.sqrt(alp_mass**2 - 4 * quark_mass**2))
    argument_at_lower_limit = 4 * sp.sqrt(spacelike_magnitude**2 + quark_mass**2 * combined_scale)
    spacelike_imaginary_closed_form = (4 * sp.pi * quark_mass**2 / (alp_mass * sp.sqrt(combined_scale))
                                       * sp.log(argument_at_one / argument_at_lower_limit))
    return locals()


def literature_function_f_numeric(tau):
    """f(tau) of Bauer:2021mvw Eq. (2.71), with the form for tau < 0 derived in Part 8; used only by the tests.

    Real argument: arcsin(1/sqrt(tau)) for tau >= 1,
    pi/2 + (i/2) ln[(1+sqrt(1-tau))/(1-sqrt(1-tau))] for 0 < tau < 1 and
    (i/2) ln[(sqrt(1-tau)+1)/(sqrt(1-tau)-1)] for tau < 0.
    Argument with a nonzero imaginary part (tau = 4 m^2/q^2 at complex m^2 in
    the lower half plane): the single expression
    (i/2) ln[(sqrt(1-tau)+1)/(sqrt(1-tau)-1)] with principal branches.
    """
    tau = mpmath.mpmathify(tau)
    imaginary_unit = mpmath.mpc(0, 1)
    if mpmath.im(tau) != 0:
        root = mpmath.sqrt(1 - tau)
        return imaginary_unit / 2 * mpmath.log((root + 1) / (root - 1))
    tau = mpmath.re(tau)
    if tau >= 1:
        return mpmath.asin(1 / mpmath.sqrt(tau))
    root = mpmath.sqrt(1 - tau)
    if tau > 0:
        return mpmath.pi / 2 + imaginary_unit / 2 * mpmath.log((1 + root) / (1 - root))
    return imaginary_unit / 2 * mpmath.log((root + 1) / (root - 1))


def literature_b1_numeric(tau):
    """B_1(tau) = 1 - tau f^2(tau), Bauer:2021mvw Eq. (2.71); used only by the tests."""
    return 1 - tau * literature_function_f_numeric(tau)**2


def literature_b3_numeric(first_tau, second_tau):
    """B_3 of Bauer:2021mvw Eq. (2.82) and Alda:2024cxn Eq. (A1) for different arguments; used only by the tests."""
    return 1 + first_tau * second_tau / (first_tau - second_tau) * (
        literature_function_f_numeric(first_tau)**2 - literature_function_f_numeric(second_tau)**2)


def b4_numeric(mass_squared, virtuality_one, virtuality_two, alp_mass_squared):
    """B_4 = 1 - loop integral by two-dimensional quadrature (15 digits); used only by the tests.

    Arguments in GeV^2. A real mass_squared requires a denominator without
    zeros; above a threshold use a complex mass_squared with negative
    imaginary part.
    """
    return 1 - loop_integral_complex_mass_numeric(mass_squared, virtuality_one, virtuality_two,
                                                  alp_mass_squared)


def one_real_photon_bracket_numeric(mass_squared, alp_momentum_squared_value, virtual_photon_momentum_squared):
    """Bracket of Bauer:2021mvw Eq. (2.73) by quadrature (15 digits); used only by the tests.

    1 - int_0^1 dx int_0^1 dy m^2/(m^2 - x(1-x)(k^2 - 2 y k.q)), with
    k^2 the squared ALP momentum and k^2 - 2 k.q = p^2 the squared momentum
    of the virtual photon (the other photon is real).
    """
    with mpmath.workdps(15):
        def inner(x_value):
            return mpmath.quad(lambda y_value: mass_squared / (
                mass_squared - x_value * (1 - x_value) * (
                    (1 - y_value) * alp_momentum_squared_value + y_value * virtual_photon_momentum_squared)),
                [0, 1])

        return 1 - mpmath.quad(inner, [0, 1])


@lru_cache(maxsize=1)
def derive_b4_properties():
    """Part 8: the function B_4(tau_a, tau_1, tau_2) = 1 - loop integral and its properties.

    tau_a = 4 m_h^2/m_a^2, tau_1 = 4 m_h^2/q_1^2, tau_2 = 4 m_h^2/q_2^2; a real
    photon corresponds to tau_i -> infinity, that is 1/tau_i = 0. The
    prescription m_h^2 -> m_h^2 - i0 of Part 7 means tau -> tau - i0 for a
    timelike invariant and tau -> tau + i0 for a spacelike one.
    """
    properties = derive_loop_integral_properties()
    derivative = derive_derivative_triangle()
    combined = derive_combined_form_factor()
    limits = derive_kinematic_limits()
    above = derive_above_threshold()
    loop_integral_symbol = derivative["loop_integral_symbol"]
    b4_function = sp.Function("B_4")
    b3_function = sp.Function("B_3")
    b1_function = sp.Function("B_1")
    f_function = sp.Function("f")

    # 1. Definition and equivalence with the form factor of Part 3.
    b4_in_loop_integral = 1 - loop_integral_symbol
    form_factor_with_b4 = REFERENCE_FACTOR * (
        pseudoscalar_coupling + (axial_coupling - pseudoscalar_coupling) * b4_in_loop_integral)
    equivalence_with_derived_sum = sp.simplify(form_factor_with_b4 - combined["combined_symbolic"])
    equivalence_with_printed_form = sp.simplify(form_factor_with_b4 - combined["expected_symbolic"])
    tau_definitions = {tau_alp: 4 * quark_mass**2 / alp_mass**2,
                       tau_first_photon: 4 * quark_mass**2 / first_virtuality,
                       tau_second_photon: 4 * quark_mass**2 / second_virtuality}
    weight_between_photons = 1 - parameter_x - parameter_y
    tau_integrand = 2 / (1 - 4 * parameter_x * parameter_y / tau_alp
                         - 4 * parameter_x * weight_between_photons / tau_second_photon
                         - 4 * parameter_y * weight_between_photons / tau_first_photon)
    # The integrand depends on the masses and momenta only through the three ratios.
    tau_integrand_check = sp.simplify(tau_integrand.subs(tau_definitions) - loop_integrand(
        quark_mass, first_virtuality, second_virtuality, alp_mass**2))
    b4_arguments = (tau_alp, tau_first_photon, tau_second_photon)
    b4_definition = sp.Eq(b4_function(*b4_arguments), 1 - sp.Integral(
        tau_integrand, (parameter_y, 0, 1 - parameter_x), (parameter_x, 0, 1)))
    b4_definition_physical = sp.Eq(b4_function(*b4_arguments), 1 - loop_integral(momentum_squared=alp_mass**2))
    form_factor_definition = REFERENCE_FACTOR * (
        pseudoscalar_coupling + (axial_coupling - pseudoscalar_coupling) * b4_function(*b4_arguments))

    # 2. Permutation symmetry. With the third Feynman parameter 1 - x - y the
    # three invariants multiply the three products of two parameters:
    # m_a^2 <-> x y, q_2^2 <-> x (1-x-y), q_1^2 <-> y (1-x-y). Each exchange of
    # two invariants is undone by exchanging two of the three parameters, a
    # linear map of the integration triangle onto itself with Jacobian of
    # modulus one.
    relabelings = {
        "tau_1 <-> tau_2": ((tau_alp, tau_second_photon, tau_first_photon), (parameter_y, parameter_x)),
        "tau_a <-> tau_1": ((tau_first_photon, tau_alp, tau_second_photon),
                            (weight_between_photons, parameter_y)),
        "tau_a <-> tau_2": ((tau_second_photon, tau_first_photon, tau_alp),
                            (parameter_x, weight_between_photons)),
    }
    triangle_corners = ((0, 0), (1, 0), (0, 1))
    permutation_checks = {}
    for name, (exchanged_arguments, new_parameters) in relabelings.items():
        exchanged_integrand = tau_integrand.subs(dict(zip(b4_arguments, exchanged_arguments)), simultaneous=True)
        relabeled_integrand = tau_integrand.subs(
            {parameter_x: new_parameters[0], parameter_y: new_parameters[1]}, simultaneous=True)
        jacobian = sp.Matrix(new_parameters).jacobian(sp.Matrix([parameter_x, parameter_y])).det()
        mapped_corners = {tuple(parameter.subs({parameter_x: corner[0], parameter_y: corner[1]})
                                for parameter in new_parameters) for corner in triangle_corners}
        permutation_checks[name] = (sp.simplify(exchanged_integrand - relabeled_integrand),
                                    sp.Abs(jacobian) - 1,
                                    0 if mapped_corners == set(triangle_corners) else 1)

    # 3. The function f(tau) in the three ranges of a real argument.
    # tau >= 1 (timelike below threshold): Part 1, loop integral = tau arcsin^2(1/sqrt(tau)).
    threshold_ratio = 4 * quark_mass**2 / alp_mass**2
    f_below_threshold = sp.asin(1 / sp.sqrt(threshold_ratio))
    b1_below_threshold_check = sp.simplify(
        threshold_ratio * f_below_threshold**2 - properties["real_photon_closed_form"])
    # tau < 0 (spacelike): Part 1 proves, for invariant/m^2 = -(v-1)^2/v with
    # v > 1, loop integral = (v/(v-1)^2) ln^2(v). Then tau = -4 v/(v-1)^2,
    # sqrt(1 - tau) = (v+1)/(v-1) > 0 and (sqrt(1-tau)+1)/(sqrt(1-tau)-1) = v,
    # so loop integral = tau f^2 requires f^2 = -ln^2(v)/4: f = (i/2) ln(v).
    exponential = properties["exponential"]
    tau_of_exponential = 4 / properties["ratio_spacelike"]
    root_of_exponential = (exponential + 1) / (exponential - 1)
    spacelike_root_check = sp.simplify(1 - tau_of_exponential - root_of_exponential**2)
    spacelike_argument_check = sp.simplify((root_of_exponential + 1) / (root_of_exponential - 1) - exponential)
    f_spacelike_in_exponential = sp.I / 2 * sp.log(exponential)
    f_spacelike_check = sp.simplify(tau_of_exponential * f_spacelike_in_exponential**2
                                    - properties["one_spacelike_invariant"])
    f_squared_spacelike_imaginary_part = sp.im(sp.expand(f_spacelike_in_exponential**2))
    root_spacelike = sp.sqrt(1 - tau_spacelike)
    f_spacelike = sp.I / 2 * sp.log((root_spacelike + 1) / (root_spacelike - 1))
    # 0 < tau < 1 (timelike above threshold): Part 7 continues v -> w exp(-i pi),
    # w = (1+sqrt(1-tau))/(1-sqrt(1-tau)) > 1, so f = (i/2)(ln w - i pi).
    logarithm_argument = above["logarithm_argument"]
    f_above_in_argument = sp.I / 2 * (sp.log(logarithm_argument) - sp.I * sp.pi)
    f_above_form_check = sp.expand(f_above_in_argument - sp.pi / 2 - sp.I / 2 * sp.log(logarithm_argument))
    f_above_check = sp.simplify(sp.expand(
        above["threshold_ratio_in_argument"] * f_above_in_argument**2 - above["above_in_logarithm_argument"]))
    root_above = sp.sqrt(1 - tau_alp)
    f_above_threshold = sp.pi / 2 + sp.I / 2 * sp.log((1 + root_above) / (1 - root_above))

    # 3 (i), (iii). One real photon, q_2^2 = 0. The substitution
    # x = (1 - y) xi, 0 < xi < 1, and the average with y -> 1 - y turn the
    # triangle integrand into the integrand of Bauer:2021mvw Eq. (2.73) with
    # k^2 = m_a^2 and p^2 = q_1^2 (their x is y here, their y is 1 - xi).
    line_fraction = sp.Symbol("xi", nonnegative=True)
    one_real_photon_integrand = loop_integrand(quark_mass, first_virtuality, 0, alp_mass**2)
    substituted_integrand = sp.simplify(
        one_real_photon_integrand.subs(parameter_x, (1 - parameter_y) * line_fraction) * (1 - parameter_y))
    interpolated_invariant = line_fraction * alp_mass**2 + (1 - line_fraction) * first_virtuality
    one_real_photon_literature_integrand = quark_mass**2 / (
        quark_mass**2 - parameter_y * (1 - parameter_y) * interpolated_invariant)
    one_real_photon_integrand_check = sp.simplify(
        (substituted_integrand + substituted_integrand.subs(parameter_y, 1 - parameter_y)) / 2
        - one_real_photon_literature_integrand)
    # The y integral at fixed xi is int_0^1 dy/(1 - a y(1-y)) = 2 dg/da at
    # a = interpolated invariant/m^2 (Part 1), and d/dxi g(a) = (dg/da)(m_a^2 - q_1^2)/m^2,
    # so the xi integral gives 2 m^2 [g(m_a^2/m^2) - g(q_1^2/m^2)]/(m_a^2 - q_1^2).
    # Parts 1 and 7 give g(a) = (a/2) x (loop integral for one invariant) = 2 f^2(4/a).
    ratio = properties["ratio"]
    y_integrand_check = sp.simplify(
        one_real_photon_literature_integrand
        - properties["symmetric_integrand"].subs({ratio: interpolated_invariant / quark_mass**2,
                                                  parameter_x: parameter_y}, simultaneous=True))
    logarithm_integral_function = sp.Function("g")
    chain_rule_check = sp.simplify(
        sp.diff(logarithm_integral_function(interpolated_invariant / quark_mass**2), line_fraction)
        - (alp_mass**2 - first_virtuality) / quark_mass**2
        * sp.diff(logarithm_integral_function(ratio), ratio).subs(ratio, interpolated_invariant / quark_mass**2))
    f_squared_alp, f_squared_first_photon = sp.symbols("f_squared_alp f_squared_first_photon")
    one_real_photon_result = 1 - 2 * quark_mass**2 / (alp_mass**2 - first_virtuality) * (
        2 * f_squared_alp - 2 * f_squared_first_photon)
    literature_b3 = 1 + tau_alp * tau_first_photon / (tau_alp - tau_first_photon) * (
        f_squared_alp - f_squared_first_photon)
    one_real_photon_b3_check = sp.simplify(one_real_photon_result - literature_b3.subs(tau_definitions))
    b3_printed = 1 + tau_alp * tau_first_photon / (tau_alp - tau_first_photon) * (
        f_function(tau_alp)**2 - f_function(tau_first_photon)**2)
    b1_printed = 1 - tau_alp * f_function(tau_alp)**2
    # Equal arguments: B_3(tau, tau) = 1 + tau^2 d(f^2)/dtau. For tau < 0:
    b3_equal_arguments_spacelike = sp.simplify(1 + tau_spacelike**2 * sp.diff(f_spacelike**2, tau_spacelike))

    # 4. Closed forms in the tau variables.
    # (a) tau_a -> infinity, tau_1 = tau_2 = tau < 0: Part 4 (a) with m_h^2/|q^2| = -tau/4.
    mass_ratio = limits["mass_ratio"]
    equal_spacelike_closed_form = 1 + tau_spacelike / (2 * root_spacelike) * sp.log(
        (root_spacelike + 1) / (root_spacelike - 1))
    equal_spacelike_check = sp.simplify(
        equal_spacelike_closed_form - (1 - limits["equal_closed_ratio"].subs(mass_ratio, -tau_spacelike / 4)))
    # By the permutation symmetry this is also B_3(tau, tau).
    equal_arguments_check = sp.simplify(b3_equal_arguments_spacelike - equal_spacelike_closed_form)
    # (b) tau_a -> infinity, tau_2 -> infinity, tau_1 = tau < 0: Part 4 (c).
    single_spacelike_closed_form = 1 - tau_spacelike * f_spacelike**2
    single_spacelike_check = sp.simplify(sp.expand(
        single_spacelike_closed_form - (1 - limits["single_closed_ratio"].subs(mass_ratio, -tau_spacelike / 4))))
    # (c) real photons: B_1(tau_a), in the two ranges of tau_a (Parts 1 and 7).
    real_photons_below_threshold = 1 - properties["real_photon_closed_form"]
    real_photons_above_threshold = 1 - above["above_closed_form"]
    real_photons_above_check = sp.simplify(
        (1 - tau_alp * f_above_threshold**2).subs(tau_alp, threshold_ratio) - real_photons_above_threshold)

    # 5. Limits in the physical variables.
    # m_h -> infinity: 1 - (expansion of Part 1), with p^2 = m_a^2.
    heavy_expansion = sp.expand((1 - properties["heavy_mass_limit"] - properties["heavy_mass_first_correction"]
                                 - properties["heavy_mass_second_correction"]).subs(
                                     alp_momentum_squared, alp_mass**2))
    invariants = (alp_mass**2, first_virtuality, second_virtuality)
    heavy_expansion_expected = (
        -sum(invariants) / (12 * quark_mass**2)
        - (sum(invariant**2 for invariant in invariants)
           + invariants[0] * invariants[1] + invariants[0] * invariants[2] + invariants[1] * invariants[2])
        / (90 * quark_mass**4))
    heavy_expansion_check = sp.simplify(heavy_expansion - heavy_expansion_expected)
    inverse_arguments = tuple(1 / argument for argument in b4_arguments)
    heavy_expansion_in_tau = (
        -sum(inverse_arguments) / 3
        - sp.Rational(8, 45) * (sum(inverse**2 for inverse in inverse_arguments)
                                + inverse_arguments[0] * inverse_arguments[1]
                                + inverse_arguments[0] * inverse_arguments[2]
                                + inverse_arguments[1] * inverse_arguments[2]))
    heavy_expansion_in_tau_check = sp.simplify(heavy_expansion_in_tau.subs(tau_definitions)
                                               - heavy_expansion_expected)
    # m_h -> 0, or |q^2| -> infinity at fixed m_h: B_4 -> 1. First correction in
    # the three cases expanded in Part 4.
    physical = limits["physical"]
    massless_limits = (1 - limits["equal_massless_limit"], 1 - limits["single_massless_limit"],
                       1 - above["massless_limit_above_threshold"])
    equal_first_correction = physical(-limits["equal_expansion_ratio"].coeff(mass_ratio, 1) * mass_ratio)
    single_first_correction = physical(-limits["single_expansion_ratio"].coeff(mass_ratio, 1) * mass_ratio)
    alp_mass_first_correction = (-alp_mass**2 / spacelike_magnitude
                                 * limits["alp_mass_coefficient_expansion_ratio"].coeff(mass_ratio, 0))

    # 6. Imaginary part. B_4 = 1 - loop integral, so Im B_4 = -Im(loop integral) <= 0.
    real_photons_real_part = 1 - above["above_real_part"]
    real_photons_imaginary_part = -above["above_imaginary_part"]
    real_photons_parts_check = sp.simplify(sp.expand(
        real_photons_above_threshold - real_photons_real_part - sp.I * real_photons_imaginary_part))
    spacelike_photons_imaginary_part = -above["spacelike_imaginary_closed_form"]
    return locals()


class PseudoscalarTriangleTests(unittest.TestCase):
    """Part 1: Dirac traces, Feynman parameters, sign and the loop integral."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_pseudoscalar_triangle()
        cls.properties = derive_loop_integral_properties()

    def test_gamma5_orientation(self):
        self.assertEqual(sp.trace(GAMMA[0] * GAMMA[1] * GAMMA[2] * GAMMA[3] * GAMMA5), -4 * sp.I)
        self.assertEqual(GAMMA5, sp.I * GAMMA[0] * GAMMA[1] * GAMMA[2] * GAMMA[3])

    def test_numerators_of_both_orderings(self):
        for first in range(4):
            for second in range(4):
                expected = self.result["numerator_expected"][first][second]
                self.assertEqual(sp.expand(self.result["numerator_first_ordering"][first][second] - expected), 0)
                self.assertEqual(sp.expand(self.result["numerator_second_ordering"][first][second] - expected), 0)
        self.assertEqual(self.result["numerator_coefficient"], 4 * sp.I * quark_mass)

    def test_feynman_parameters(self):
        self.assertEqual(self.result["combination_check_first"], 0)
        self.assertEqual(self.result["combination_check_second"], 0)
        self.assertEqual(self.result["inner_antiderivative_check"], 0)
        self.assertEqual(self.result["outer_antiderivative_check"], 0)
        self.assertEqual(self.result["three_factor_identity"], 0)
        # The second ordering is the first with x and y exchanged.
        first_ordering = feynman_denominator(quark_mass, first_virtuality, second_virtuality,
                                             alp_momentum_squared)
        second_ordering = feynman_denominator(quark_mass, second_virtuality, first_virtuality,
                                              alp_momentum_squared, x=parameter_y, y=parameter_x)
        self.assertEqual(sp.expand(first_ordering - second_ordering), 0)

    def test_momentum_integral_and_sign(self):
        self.assertEqual(self.result["radial_integral"], 1 / (4 * self.result["positive_denominator"]))
        self.assertEqual(sp.simplify(self.result["momentum_integral"]
                                     + sp.I / (32 * sp.pi**2 * self.result["positive_denominator"])), 0)
        self.assertEqual(self.result["overall_sign"], 1)
        self.assertEqual(sp.simplify(self.result["form_factor_per_loop_integral"]
                                     - REFERENCE_FACTOR * pseudoscalar_coupling), 0)

    def test_local_photon_operator_sign(self):
        self.assertEqual(sp.simplify(self.result["local_photon_form_factor"]
                                     - photon_coupling / (4 * sp.pi**2 * decay_constant)), 0)

    def test_symmetry_and_heavy_mass_expansion(self):
        properties = self.properties
        self.assertEqual(properties["photon_exchange_difference"], 0)
        self.assertEqual(properties["heavy_mass_limit"], 1)
        self.assertEqual(properties["geometric_remainder"], 0)
        self.assertEqual(sp.simplify(properties["heavy_mass_first_correction"]
                                     - (first_virtuality + second_virtuality + alp_momentum_squared)
                                     / (12 * quark_mass**2)), 0)

    def test_real_photon_closed_form(self):
        properties = self.properties
        for name in ("real_inner_check", "real_inner_result_check", "symmetrization_check",
                     "below_antiderivative_check", "below_claim_derivative_check", "below_claim_at_zero"):
            self.assertEqual(properties[name], 0, name)
        tangent = properties["tangent"]
        self.assertEqual(sp.simplify(properties["real_photon_in_tangent"]
                                     - (1 + tangent**2) * sp.atan(tangent)**2 / tangent**2), 0)
        # tau = (1+t^2)/t^2 and arcsin(1/sqrt(tau)) = arctan(t).
        self.assertEqual(properties["arcsine_derivative_check"], 0)
        self.assertEqual(properties["arcsine_at_zero"], 0)
        # Quadrature of the two-dimensional integral, m_h = 1 GeV.
        for alp_mass_value in (mpmath.mpf("0.5"), mpmath.mpf("1"), mpmath.mpf("1.9")):
            closed = properties["real_photon_closed_form"].subs({quark_mass: 1, alp_mass: sp.Float(str(alp_mass_value), 30)})
            numeric = loop_integral_numeric(1, 0, 0, alp_mass_value**2)
            self.assertAlmostEqual(float(numeric), float(closed), places=12)
        # tau = 4: 4 arcsin(1/2)^2 = pi^2/9.
        self.assertEqual(sp.simplify(properties["real_photon_closed_form"].subs({quark_mass: 1, alp_mass: 1})
                                     - sp.pi**2 / 9), 0)

    def test_one_spacelike_invariant(self):
        properties = self.properties
        for name in ("spacelike_antiderivative_check", "spacelike_claim_derivative_check",
                     "spacelike_claim_at_zero"):
            self.assertEqual(properties[name], 0, name)


class DerivativeTriangleTests(unittest.TestCase):
    """Part 2: regulated axial Ward identity."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_derivative_triangle()

    def test_vertex_identities(self):
        self.assertEqual(self.result["vertex_identity"], sp.zeros(4))
        self.assertEqual(self.result["inverse_propagator_identity"], sp.zeros(4))

    def test_terms_without_alp_vertex(self):
        result = self.result
        for first in range(4):
            for second in range(4):
                self.assertEqual(sp.expand(result["remainder_one"][first][second]
                                           - result["remainder_one_expected"][first][second]), 0)
                self.assertEqual(sp.expand(result["remainder_two"][first][second]
                                           - result["remainder_two_expected"][first][second]), 0)
                self.assertEqual(sp.diff(result["remainder_one"][first][second], quark_mass), 0)
                self.assertEqual(sp.diff(result["remainder_two"][first][second], quark_mass), 0)
                for name in ("remainder_one_shifted", "remainder_two_shifted"):
                    shifted = result[name][first][second]
                    self.assertEqual(sp.expand(shifted + shifted.subs(result["reflection"], simultaneous=True)), 0)

    def test_regulated_amplitude(self):
        result = self.result
        expected = REFERENCE_FACTOR * axial_coupling * (1 - result["loop_integral_symbol"])
        self.assertEqual(sp.simplify(result["axial_form_factor_symbolic"] - expected), 0)
        self.assertEqual(sp.simplify(result["axial_massless_limit"] - result["reference_massless"]), 0)


class CombinedFormFactorTests(unittest.TestCase):
    """Part 3: combination, rotation invariance and limits."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_combined_form_factor()

    def test_combination_and_rotation(self):
        self.assertEqual(self.result["combined_difference"], 0)
        self.assertEqual(self.result["rotation_difference"], 0)

    def test_limits(self):
        self.assertEqual(sp.simplify(self.result["decoupling_limit"]
                                     - REFERENCE_FACTOR * pseudoscalar_coupling), 0)
        self.assertEqual(sp.simplify(self.result["massless_limit"] - REFERENCE_FACTOR * axial_coupling), 0)
        threshold_ratio = 4 * quark_mass**2 / alp_mass**2
        expected = REFERENCE_FACTOR * (axial_coupling - (axial_coupling - pseudoscalar_coupling)
                                       * threshold_ratio * sp.asin(1 / sp.sqrt(threshold_ratio))**2)
        self.assertEqual(sp.simplify(self.result["real_photon_limit"] - expected), 0)


class KinematicLimitTests(unittest.TestCase):
    """Part 4: closed forms and expansions."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_kinematic_limits()

    def test_equal_virtualities(self):
        result = self.result
        for name in ("independent_of_x", "one_dimensional_check", "equal_symmetrization_check",
                     "equal_antiderivative_check", "equal_endpoint_check"):
            self.assertEqual(result[name], 0, name)
        ratio, logarithm = result["mass_ratio"], result["logarithm"]
        self.assertEqual(sp.expand(result["equal_expansion_ratio"]
                                   - (2 * ratio * logarithm + ratio**2 * (4 - 4 * logarithm))), 0)
        self.assertEqual(result["equal_massless_limit"], 0)
        heavy_variable = sp.Symbol("s", positive=True)
        self.assertEqual(sp.expand(result["equal_heavy_expansion"] - (1 - heavy_variable / 6)), 0)
        for mass_value, magnitude in ((1, 3), (1, 400), (2, 1)):
            closed = result["equal_closed_form"].subs({quark_mass: mass_value, spacelike_magnitude: magnitude})
            numeric = loop_integral_numeric(mass_value, -magnitude, -magnitude, 0)
            self.assertAlmostEqual(float(numeric), float(closed), places=12)

    def test_alp_mass_correction(self):
        result = self.result
        for name in ("first_order_check", "first_order_symmetrization_check", "scaling_check"):
            self.assertEqual(result[name], 0, name)
        ratio, logarithm = result["mass_ratio"], result["logarithm"]
        self.assertEqual(sp.expand(result["alp_mass_coefficient_expansion_ratio"]
                                   - (sp.Rational(1, 3) - ratio * (1 + logarithm) / 3)), 0)
        mpmath.mp.dps = 25
        for mass_value, magnitude in ((1, 3), (1, 50)):
            closed = result["alp_mass_coefficient_ratio"].subs(ratio, sp.Rational(mass_value**2, magnitude))
            numeric = mpmath.quad(lambda u: mass_value**2 * u**3 / (3 * (mass_value**2 + magnitude * u * (1 - u))**2),
                                  [0, 1]) * magnitude
            self.assertAlmostEqual(float(numeric), float(closed), places=12)
        # Finite difference of the two-dimensional integral in p^2, m_h = 1 GeV, |q^2| = 3 GeV^2.
        step = mpmath.mpf("1e-6")
        slope = (loop_integral_numeric(1, -3, -3, step) - loop_integral_numeric(1, -3, -3, -step)) / (2 * step)
        closed = result["alp_mass_coefficient_ratio"].subs(ratio, sp.Rational(1, 3)) / 3
        self.assertAlmostEqual(float(slope), float(closed), places=9)

    def test_one_real_photon(self):
        result = self.result
        for name in ("single_inner_check", "single_inner_result_check", "single_matches_one_invariant",
                     "ratio_matches"):
            self.assertEqual(result[name], 0, name)
        ratio, logarithm = result["mass_ratio"], result["logarithm"]
        self.assertEqual(sp.expand(result["single_expansion_ratio"]
                                   - (ratio * logarithm**2 + 4 * ratio**2 * logarithm)), 0)
        self.assertEqual(result["single_massless_limit"], 0)
        for mass_value, magnitude in ((1, 3), (1, 400)):
            closed = result["single_closed_form"].subs({quark_mass: mass_value, spacelike_magnitude: magnitude})
            self.assertAlmostEqual(float(loop_integral_numeric(mass_value, -magnitude, 0, 0)), float(closed), places=12)
            self.assertAlmostEqual(float(loop_integral_numeric(mass_value, 0, -magnitude, 0)), float(closed), places=12)

    def test_low_energy(self):
        expected = 1 + (first_virtuality + second_virtuality + alp_momentum_squared) / (12 * quark_mass**2)
        self.assertEqual(sp.simplify(self.result["low_energy_expansion"] - expected), 0)
        point = {first_virtuality: -1, second_virtuality: -2, alp_momentum_squared: sp.Rational(1, 2),
                 quark_mass: 10}
        numeric = loop_integral_numeric(10, -1, -2, mpmath.mpf("0.5"))
        self.assertAlmostEqual(float(numeric), float(expected.subs(point)), delta=1e-5)
        with_next_term = (expected + self.result["low_energy_next_term"]).subs(point)
        self.assertAlmostEqual(float(numeric), float(with_next_term), places=7)


class AxialCurrentMatrixElementTests(unittest.TestCase):
    """Part 5: two-point function and comparison with the direct triangle."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_axial_current_matrix_element()

    def test_numerators(self):
        result = self.result
        for index in range(4):
            self.assertEqual(sp.expand(result["derivative_numerator"][index]
                                       - result["derivative_numerator_expected"][index]), 0)
            self.assertEqual(sp.expand(result["pseudoscalar_numerator"][index]
                                       - result["pseudoscalar_numerator_expected"][index]), 0)

    def test_scalar_integral(self):
        result = self.result
        self.assertEqual(result["two_point_combination_check"], 0)
        self.assertEqual(result["cutoff_remainder"], 0)
        expected = sp.log(cutoff**2 / quark_mass**2) - 1 + alp_mass**2 / (6 * quark_mass**2)
        self.assertEqual(sp.simplify(sp.expand_log(result["scalar_integral_expanded"] - expected, force=True)), 0)

    def test_matrix_element_and_local_expansion(self):
        result = self.result
        expected = (-sp.I * color_number * quark_mass**2 / (2 * sp.pi**2 * decay_constant)
                    * (axial_coupling - pseudoscalar_coupling))
        self.assertEqual(sp.simplify(result["matrix_element_coefficient_per_integral"] - expected), 0)
        expected_local = (-REFERENCE_FACTOR * (axial_coupling - pseudoscalar_coupling)
                          * 2 * quark_mass**2 / (spacelike_magnitude + quark_mass**2))
        self.assertEqual(sp.simplify(result["local_per_integral"] - expected_local), 0)

    def test_comparison_with_direct_triangle(self):
        result = self.result
        self.assertEqual(result["remainder_check"], 0)
        # The remainder contains neither ln(m_h^2) nor m_a^2: it is polynomial in m_h^2.
        scaled = sp.expand_log(result["remainder"] * spacelike_magnitude / quark_mass**2, force=True)
        self.assertEqual(sp.simplify(sp.diff(scaled, quark_mass)), 0)
        self.assertEqual(sp.simplify(sp.diff(scaled, alp_mass)), 0)


class VectorCouplingTests(unittest.TestCase):
    """Part 6: the vector coupling gives no contribution."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_vector_triangle()

    def test_no_levi_civita_structure(self):
        result = self.result
        self.assertTrue(result["numerators_are_real"])
        for first in range(4):
            for second in range(4):
                self.assertEqual(result["odd_part_first_ordering"][first][second], 0)
                self.assertEqual(result["odd_part_second_ordering"][first][second], 0)
                # The pseudoscalar numerator is the opposite case: it has no even part.
                self.assertEqual(result["pseudoscalar_even_part"][first][second], 0)
        self.assertEqual(result["levi_civita_coefficient_first"], 0)
        self.assertEqual(result["levi_civita_coefficient_second"], 0)
        self.assertEqual(sp.simplify(result["pseudoscalar_levi_civita_coefficient"] - 4 * sp.I * quark_mass), 0)
        # The vector numerators themselves are not zero.
        self.assertNotEqual(result["numerator_first_ordering"][0][0], 0)
        self.assertNotEqual(result["numerator_first_ordering"][1][2], 0)

    def test_photon_orderings_cancel(self):
        result = self.result
        for index in range(4):
            self.assertEqual(result["vector_conjugation"][index], sp.zeros(4))
        self.assertEqual(result["pseudoscalar_conjugation"], sp.zeros(4))
        self.assertEqual(result["denominator_difference"], 0)
        for first in range(4):
            for second in range(4):
                self.assertEqual(result["ordering_sum"][first][second], 0)
                self.assertEqual(result["pseudoscalar_ordering_difference"][first][second], 0)
        self.assertEqual(result["vector_form_factor"], 0)


class AboveThresholdTests(unittest.TestCase):
    """Part 7: m_a > 2 m_h, imaginary part from m_h^2 -> m_h^2 - i0."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_above_threshold()

    @staticmethod
    def closed_form(threshold_ratio):
        """Closed form at complex 4 m^2/m_a^2 in the lower half plane (principal branches)."""
        velocity = mpmath.sqrt(1 - threshold_ratio)
        return threshold_ratio * (mpmath.pi / 2 + mpmath.mpc(0, 1) / 2
                                  * mpmath.log((1 + velocity) / (1 - velocity)))**2

    def test_closed_form_derivation(self):
        result = self.result
        for name in ("zero_checks", "slope_checks", "spacelike_continuation_check", "arcsine_form_checks",
                     "half_logarithm_checks"):
            self.assertEqual(tuple(result[name]), (0, 0), name)
        for name in ("principal_value_antiderivative_check", "above_claim_derivative_check",
                     "above_claim_at_threshold_check", "velocity_relation_check", "closed_form_check",
                     "parts_check", "massless_limit_above_threshold"):
            self.assertEqual(result[name], 0, name)
        argument = result["logarithm_argument"]
        self.assertEqual(sp.simplify(result["principal_value"]
                                     + 2 * argument * sp.log(argument) / (argument**2 - 1)), 0)
        self.assertEqual(sp.simplify(result["delta_function_integral"] - 2 * argument / (argument**2 - 1)), 0)
        # Principal value by quadrature with equal intervals removed around the zeros.
        mpmath.mp.dps = 25
        argument_value, excluded = mpmath.mpf(3), mpmath.mpf("1e-9")
        ratio_value = (argument_value + 1)**2 / argument_value
        zeros = (1 / (argument_value + 1), argument_value / (argument_value + 1))
        numeric = sum(mpmath.quad(lambda x: 1 / (1 - ratio_value * x * (1 - x)), interval) for interval in (
            [0, zeros[0] - excluded], [zeros[0] + excluded, zeros[1] - excluded], [zeros[1] + excluded, 1]))
        self.assertAlmostEqual(float(numeric), float(result["principal_value"].subs(argument, 3)), places=7)

    def test_real_and_imaginary_parts_by_quadrature(self):
        """ln(D - i0) = ln|D| - i pi for D < 0 in the one-dimensional integral of Part 1, 25 digits."""
        mpmath.mp.dps = 25
        result = self.result
        for threshold_ratio in ("0.05", "0.25", "0.5", "0.8", "0.99"):
            ratio_value = 4 / mpmath.mpf(threshold_ratio)
            velocity = mpmath.sqrt(1 - mpmath.mpf(threshold_ratio))
            zeros = ((1 - velocity) / 2, (1 + velocity) / 2)
            real_numeric = -2 / ratio_value * mpmath.quad(
                lambda x: mpmath.log(abs(1 - ratio_value * x * (1 - x))) / x, [0, zeros[0], zeros[1], 1])
            imaginary_numeric = 2 * mpmath.pi / ratio_value * mpmath.quad(lambda x: 1 / x, [zeros[0], zeros[1]])
            point = {quark_mass: 1, alp_mass: sp.sqrt(4 / sp.Rational(threshold_ratio))}
            self.assertAlmostEqual(float(real_numeric), float(result["above_real_part"].subs(point)), places=10)
            self.assertAlmostEqual(float(imaginary_numeric),
                                   float(result["above_imaginary_part"].subs(point)), places=12)
            self.assertGreater(float(imaginary_numeric), 0)
            closed = self.closed_form(mpmath.mpf(threshold_ratio))
            self.assertAlmostEqual(float(real_numeric), float(closed.real), places=10)
            self.assertAlmostEqual(float(imaginary_numeric), float(closed.imag), places=12)

    def test_two_dimensional_integral_at_complex_mass(self):
        """Full two-dimensional quadrature at m^2 = (1 - i eta) GeV^2, then eta -> 0."""
        # m_a^2 = 16 GeV^2 (4 m^2/m_a^2 = 1/4): every node, and the extrapolation.
        extrapolated, samples = loop_integral_real_mass_limit(0, 0, 16)
        for node, value in samples:
            self.assertLess(abs(value - self.closed_form(4 * (1 - mpmath.mpc(0, 1) * node) / 16)), 1e-8)
        closed = self.closed_form(mpmath.mpf(1) / 4)
        self.assertLess(abs(extrapolated - closed), 1e-5)
        self.assertGreater(extrapolated.imag, 0)
        # Other values of 4 m^2/m_a^2 on both sides of the threshold, eta = 0.2.
        imaginary_fraction = mpmath.mpf("0.2")
        for momentum_squared in (8, 5, mpmath.mpf("2.5")):
            mass_squared = 1 - mpmath.mpc(0, 1) * imaginary_fraction
            numeric = loop_integral_complex_mass_numeric(mass_squared, 0, 0, momentum_squared)
            self.assertLess(abs(numeric - self.closed_form(4 * mass_squared / momentum_squared)), 1e-9)

    def test_sign_of_the_imaginary_part(self):
        """m^2 - i0 gives Im > 0; the opposite prescription gives the complex conjugate."""
        with mpmath.workdps(15):
            lower = loop_integral_complex_mass_numeric(1 - mpmath.mpc(0, "0.3"), 0, 0, 16)
            upper = loop_integral_complex_mass_numeric(1 + mpmath.mpc(0, "0.3"), 0, 0, 16)
            self.assertGreater(lower.imag, 0)
            self.assertLess(abs(upper - mpmath.conj(lower)), 1e-12)
        point = {quark_mass: 1, alp_mass: 4}
        self.assertGreater(float(self.result["above_imaginary_part"].subs(point)), 0)

    def test_threshold(self):
        result = self.result
        distance = result["distance"]
        self.assertEqual(result["continuation_check"], 0)
        self.assertEqual(result["value_at_threshold"], (sp.pi**2 / 4, sp.pi**2 / 4))
        self.assertEqual(sp.expand(result["below_series"] - (
            sp.pi**2 / 4 - sp.pi * distance + (1 + sp.pi**2 / 4) * distance**2)), 0)
        self.assertEqual(sp.expand(result["above_series"] - (
            sp.pi**2 / 4 + sp.I * sp.pi * distance - (1 + sp.pi**2 / 4) * distance**2)), 0)
        self.assertEqual(result["below_slope_times_distance"], -sp.pi / 2)
        self.assertEqual(sp.simplify(result["above_real_slope_at_threshold"] - 1 - sp.pi**2 / 4), 0)
        self.assertEqual(result["above_imaginary_slope_times_distance"], -sp.pi / 2)
        # Closed forms on the two sides of 4 m_h^2/m_a^2 = 1, m_h = 1 GeV.
        mpmath.mp.dps = 25
        properties = derive_loop_integral_properties()
        offset = sp.Rational(1, 10**12)
        below = properties["real_photon_closed_form"].subs({quark_mass: 1, alp_mass: 2 / sp.sqrt(1 + offset)})
        above = result["above_closed_form"].subs({quark_mass: 1, alp_mass: 2 / sp.sqrt(1 - offset)})
        self.assertAlmostEqual(float(below), float(sp.pi**2 / 4), delta=1e-5)
        self.assertAlmostEqual(complex(above).real, float(sp.pi**2 / 4), delta=1e-5)
        self.assertAlmostEqual(complex(above).imag, 0, delta=1e-5)
        # The closed form at 4 m^2/m_a^2 - i0 > 1 is the arcsine form.
        for threshold_ratio in ("1.5", "4"):
            shifted = mpmath.mpc(threshold_ratio, "-1e-20")
            arcsine = mpmath.mpf(threshold_ratio) * mpmath.asin(1 / mpmath.sqrt(mpmath.mpf(threshold_ratio)))**2
            self.assertLess(abs(self.closed_form(shifted) - arcsine), 1e-15)

    def test_spacelike_photons_imaginary_part(self):
        result = self.result
        for name in ("sum_variable_check", "slope_squared_check", "lower_limit_check", "spacelike_integrand_check",
                     "spacelike_imaginary_antiderivative_check"):
            self.assertEqual(result[name], 0, name)
        self.assertEqual(tuple(result["spacelike_zero_checks"]), (0, 0))
        self.assertEqual(tuple(result["spacelike_slope_checks"]), (0, 0))
        self.assertEqual(sp.simplify(result["slope_squared_at_one"]
                                     - alp_mass**2 * (alp_mass**2 - 4 * quark_mass**2)), 0)
        # Endpoint values of the argument of the logarithm.
        self.assertEqual(sp.simplify(result["spacelike_antiderivative_argument"].subs(parameter_u, 1)
                                     - result["argument_at_one"]), 0)
        # At the lower limit the square root in the argument vanishes (lower_limit_check).
        combined_scale = result["combined_scale"]
        self.assertEqual(sp.simplify(2 * combined_scale * result["lower_limit"] - 4 * spacelike_magnitude
                                     - result["argument_at_lower_limit"]), 0)
        mpmath.mp.dps = 25
        for mass_value, magnitude, alp_mass_value in ((1, 3, 4), (1, sp.Rational(1, 2), 3), (1, 40, 5)):
            point = {quark_mass: mass_value, spacelike_magnitude: magnitude, alp_mass: alp_mass_value}
            integrand = sp.lambdify(parameter_u, result["spacelike_imaginary_integral"].function.subs(point),
                                    "mpmath")
            one_dimensional = mpmath.quad(integrand, [mpmath.mpf(sp.N(result["lower_limit"].subs(point), 30)), 1])
            closed = sp.N(result["spacelike_imaginary_closed_form"].subs(point), 30)
            self.assertAlmostEqual(float(one_dimensional.real), float(closed), places=10)
        # |q^2| -> 0 gives the imaginary part for real photons.
        real_photon_limit = result["spacelike_imaginary_closed_form"].subs(spacelike_magnitude, 0)
        for alp_mass_value in (3, 4, 10):
            point = {quark_mass: 1, alp_mass: alp_mass_value}
            self.assertAlmostEqual(float(real_photon_limit.subs(point)),
                                   float(result["above_imaginary_part"].subs(point)), places=12)
        # Full two-dimensional integral, m_h = 1 GeV, |q^2| = 3 GeV^2, m_a = 4 GeV.
        extrapolated, _ = loop_integral_real_mass_limit(-3, -3, 16)
        point = {quark_mass: 1, spacelike_magnitude: 3, alp_mass: 4}
        self.assertAlmostEqual(float(extrapolated.imag),
                               float(result["spacelike_imaginary_closed_form"].subs(point)), delta=1e-5)


class B4PropertyTests(unittest.TestCase):
    """Part 8: the function B_4, its symmetry, reductions to B_3 and B_1, closed forms and limits."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive_b4_properties()

    def assert_close(self, first, second, tolerance):
        self.assertLess(abs(mpmath.mpmathify(first) - mpmath.mpmathify(second)), tolerance)

    def test_equivalence_with_form_factor(self):
        result = self.result
        for name in ("equivalence_with_derived_sum", "equivalence_with_printed_form", "tau_integrand_check"):
            self.assertEqual(result[name], 0, name)
        # Numerical value of the complete form factor in both forms, m_h = 1 GeV.
        loop_integral_value = loop_integral_complex_mass_numeric(1, -2, mpmath.mpf("0.5"), mpmath.mpf("1.5"))
        b4_value = b4_numeric(1, -2, mpmath.mpf("0.5"), mpmath.mpf("1.5"))
        axial_value, pseudoscalar_value = mpmath.mpf("0.7"), mpmath.mpf("-1.3")
        self.assert_close(pseudoscalar_value + (axial_value - pseudoscalar_value) * b4_value,
                          axial_value - (axial_value - pseudoscalar_value) * loop_integral_value, 1e-13)

    def test_permutation_symmetry(self):
        self.assertEqual(set(self.result["permutation_checks"]),
                         {"tau_1 <-> tau_2", "tau_a <-> tau_1", "tau_a <-> tau_2"})
        for name, checks in self.result["permutation_checks"].items():
            self.assertEqual(tuple(checks), (0, 0, 0), name)
        # All six orderings of three different invariants (GeV^2), m_h = 1 GeV:
        # timelike below threshold and spacelike values.
        orderings = ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0))
        for invariants in ((mpmath.mpf("1.3"), mpmath.mpf("-2.7"), mpmath.mpf("0.6")),
                           (mpmath.mpf("0.5"), mpmath.mpf(-3), mpmath.mpf("-7.5")),
                           (mpmath.mpf("2.2"), mpmath.mpf("1.1"), mpmath.mpf(-40))):
            values = [b4_numeric(1, *(invariants[index] for index in ordering)) for ordering in orderings]
            for value in values[1:]:
                self.assert_close(value, values[0], 1e-11)
        # One invariant above threshold (9 GeV^2), m^2 = (1 - 0.3 i) GeV^2.
        mass_squared = 1 - mpmath.mpc(0, "0.3")
        invariants = (mpmath.mpf(-2), mpmath.mpf("0.7"), mpmath.mpf(9))
        values = [b4_numeric(mass_squared, *(invariants[index] for index in ordering))
                  for ordering in ((0, 1, 2), (2, 0, 1), (1, 2, 0))]
        self.assertGreater(abs(values[0].imag), 0.1)
        for value in values[1:]:
            self.assert_close(value, values[0], 1e-10)

    def test_function_f_in_the_three_ranges(self):
        result = self.result
        for name in ("b1_below_threshold_check", "spacelike_root_check", "spacelike_argument_check",
                     "f_spacelike_check", "f_squared_spacelike_imaginary_part", "f_above_form_check",
                     "f_above_check"):
            self.assertEqual(result[name], 0, name)
        mpmath.mp.dps = 25
        # The three real forms are boundary values of one function of
        # tau = 4 (m^2 - i0)/q^2: tau - i0 for q^2 > 0, tau + i0 for q^2 < 0.
        for tau in ("7.5", "1.2", "0.9", "0.3", "0.01"):
            shifted = mpmath.mpc(tau, "-1e-22")
            self.assert_close(literature_function_f_numeric(shifted),
                              literature_function_f_numeric(mpmath.mpf(tau)), 1e-10)
        for tau in ("-0.02", "-1", "-35"):
            shifted = mpmath.mpc(tau, "1e-22")
            real_form = literature_function_f_numeric(mpmath.mpf(tau))
            self.assert_close(literature_function_f_numeric(shifted), real_form, 1e-18)
            # f is imaginary and f^2 is real and negative.
            self.assertEqual(real_form.real, 0)
            self.assertLess((real_form**2).real, 0)
            self.assertEqual((real_form**2).imag, 0)
            # The form printed in the literature for 0 < tau < 1, used at tau < 0, has a
            # negative argument of the logarithm. With ln(-r) = ln(r) + i pi it is the
            # form derived here; with ln(-r) = ln(r) - i pi it is not (f^2 is not real).
            root = mpmath.sqrt(1 - mpmath.mpf(tau))
            magnitude_logarithm = mpmath.log((root + 1) / (root - 1))
            imaginary_unit = mpmath.mpc(0, 1)
            with_plus = mpmath.pi / 2 + imaginary_unit / 2 * (magnitude_logarithm + imaginary_unit * mpmath.pi)
            with_minus = mpmath.pi / 2 + imaginary_unit / 2 * (magnitude_logarithm - imaginary_unit * mpmath.pi)
            self.assert_close(with_plus, real_form, 1e-20)
            self.assertGreater(abs((with_minus**2).imag), 0.1)
        # The symbolic form for tau < 0 agrees with the numerical one.
        symbolic = sp.N(result["f_spacelike"].subs(tau_spacelike, -3), 25)
        self.assert_close(mpmath.mpc(*symbolic.as_real_imag()), literature_function_f_numeric(-3), 1e-20)
        # B_1 at tau < 0 against the two-dimensional integral with one spacelike invariant.
        for invariant in (mpmath.mpf(-3), mpmath.mpf(-400)):
            for placement in ((invariant, 0, 0), (0, invariant, 0), (0, 0, invariant)):
                self.assert_close(b4_numeric(1, *placement), literature_b1_numeric(4 / invariant), 1e-10)

    def test_reduction_to_b3(self):
        result = self.result
        for name in ("one_real_photon_integrand_check", "y_integrand_check", "chain_rule_check",
                     "one_real_photon_b3_check"):
            self.assertEqual(result[name], 0, name)
        # Real mass: m_a below threshold; q_1^2 spacelike or timelike below threshold. (m_a^2, q_1^2) in GeV^2.
        for alp_mass_squared, virtuality in ((1, -3), (2, "1.5"), ("0.5", -40), ("3.5", "0.2"), ("0.3", "-0.1")):
            alp_mass_squared, virtuality = mpmath.mpf(alp_mass_squared), mpmath.mpf(virtuality)
            literature = literature_b3_numeric(4 / alp_mass_squared, 4 / virtuality)
            self.assertEqual(mpmath.im(literature), 0)
            self.assert_close(b4_numeric(1, virtuality, 0, alp_mass_squared), literature, 1e-10)
            self.assert_close(b4_numeric(1, 0, virtuality, alp_mass_squared), literature, 1e-10)
            self.assert_close(one_real_photon_bracket_numeric(1, alp_mass_squared, virtuality), literature, 1e-10)
        # tau_a < 1 (m_a^2 = 16 GeV^2) at m^2 = (1 - 0.3 i) GeV^2, q_1^2 spacelike or timelike.
        mass_squared = 1 - mpmath.mpc(0, "0.3")
        for virtuality in (mpmath.mpf(-3), mpmath.mpf(6)):
            literature = literature_b3_numeric(4 * mass_squared / 16, 4 * mass_squared / virtuality)
            self.assert_close(b4_numeric(mass_squared, virtuality, 0, 16), literature, 1e-9)
        self.assert_close(one_real_photon_bracket_numeric(mass_squared, 16, -3),
                          literature_b3_numeric(4 * mass_squared / 16, 4 * mass_squared / mpmath.mpf(-3)), 1e-9)
        # The same at real mass m = 1 GeV. In the integral of Bauer:2021mvw Eq. (2.73) the
        # y integration is elementary:
        # bracket = 1 + m^2/(k^2 - p^2) int_0^1 dx [ln(m^2 - x(1-x) k^2 - i0) - ln(m^2 - x(1-x) p^2 - i0)]/(x(1-x)),
        # with ln(D - i0) = ln|D| - i pi for D < 0.
        mpmath.mp.dps = 25

        def logarithm_with_prescription(value):
            return mpmath.log(abs(value)) - (mpmath.mpc(0, 1) * mpmath.pi if value < 0 else 0)

        for alp_mass_squared, virtuality in ((16, -3), (16, 6), (5, "0.5"), (9, 30), ("2.5", 7)):
            alp_mass_squared, virtuality = mpmath.mpf(alp_mass_squared), mpmath.mpf(virtuality)

            def integrand(x_value):
                weight = x_value * (1 - x_value)
                return (logarithm_with_prescription(1 - weight * alp_mass_squared)
                        - logarithm_with_prescription(1 - weight * virtuality)) / weight

            # Break points at the zeros of the two arguments of the logarithms.
            break_points = {mpmath.mpf(0), mpmath.mpf(1)}
            for invariant in (alp_mass_squared, virtuality):
                if invariant > 4:
                    velocity = mpmath.sqrt(1 - 4 / invariant)
                    break_points |= {(1 - velocity) / 2, (1 + velocity) / 2}
            one_dimensional = 1 + mpmath.quad(integrand, sorted(break_points)) / (alp_mass_squared - virtuality)
            literature = literature_b3_numeric(4 / alp_mass_squared, 4 / virtuality)
            self.assert_close(one_dimensional, literature, 1e-12)
            self.assertNotEqual(mpmath.im(literature), 0)
        # Only the ALP above threshold: Im B_3 = Im B_4 < 0.
        self.assertLess(mpmath.im(literature_b3_numeric(mpmath.mpf(4) / 16, mpmath.mpf(4) / -3)), 0)

    def test_b3_at_equal_arguments(self):
        result = self.result
        self.assertEqual(result["equal_arguments_check"], 0)
        mpmath.mp.dps = 25
        # tau < 0: closed form, B_3 at nearly equal arguments, and the integral with m_a^2 = q_1^2 < 0.
        for tau in (mpmath.mpf("-0.8"), mpmath.mpf(-6)):
            closed = mpmath.mpf(sp.N(result["equal_spacelike_closed_form"].subs(
                tau_spacelike, sp.Rational(str(tau))), 25))
            self.assert_close(literature_b3_numeric(tau, tau * (1 + mpmath.mpf("1e-10"))), closed, 1e-9)
            self.assert_close(b4_numeric(1, 4 / tau, 0, 4 / tau), closed, 1e-10)
        # tau > 1: B_3(tau, tau) = 1 - tau f(tau)/sqrt(tau - 1).
        for tau in (mpmath.mpf("1.6"), mpmath.mpf(9)):
            closed = 1 - tau * literature_function_f_numeric(tau) / mpmath.sqrt(tau - 1)
            self.assert_close(literature_b3_numeric(tau, tau * (1 + mpmath.mpf("1e-10"))), closed, 1e-9)
            self.assert_close(b4_numeric(1, 4 / tau, 0, 4 / tau), closed, 1e-10)

    def test_reduction_to_b1(self):
        result = self.result
        self.assertEqual(result["real_photons_above_check"], 0)
        threshold_ratio = 4 * quark_mass**2 / alp_mass**2
        self.assertEqual(sp.simplify(result["real_photons_below_threshold"]
                                     - (1 - threshold_ratio * sp.asin(1 / sp.sqrt(threshold_ratio))**2)), 0)
        for alp_mass_squared in (mpmath.mpf("0.25"), mpmath.mpf(1), mpmath.mpf("3.6")):
            self.assert_close(b4_numeric(1, 0, 0, alp_mass_squared), literature_b1_numeric(4 / alp_mass_squared),
                              1e-10)
        # B_3 with a large second argument.
        self.assert_close(literature_b3_numeric(mpmath.mpf(3), mpmath.mpf("1e12")),
                          literature_b1_numeric(mpmath.mpf(3)), 1e-11)
        # Above threshold, symbolic closed form against the numerical B_1.
        mpmath.mp.dps = 25
        for alp_mass_value in (3, 4, 10):
            symbolic = sp.N(result["real_photons_above_threshold"].subs({quark_mass: 1, alp_mass: alp_mass_value}),
                            25)
            self.assert_close(mpmath.mpc(*symbolic.as_real_imag()),
                              literature_b1_numeric(mpmath.mpf(4) / alp_mass_value**2), 1e-20)

    def test_closed_forms_in_tau(self):
        result = self.result
        for name in ("equal_spacelike_check", "single_spacelike_check"):
            self.assertEqual(result[name], 0, name)
        limits = derive_kinematic_limits()
        mpmath.mp.dps = 25
        for magnitude in (3, 400, sp.Rational(1, 4)):
            tau_value = -sp.Integer(4) / magnitude
            # (a) equal spacelike virtualities, m_a = 0.
            closed = sp.N(result["equal_spacelike_closed_form"].subs(tau_spacelike, tau_value), 25)
            earlier = sp.N(1 - limits["equal_closed_form"].subs({quark_mass: 1, spacelike_magnitude: magnitude}), 25)
            self.assert_close(closed, earlier, 1e-20)
            virtuality = -mpmath.mpf(sp.N(magnitude, 25))
            self.assert_close(b4_numeric(1, virtuality, virtuality, 0), closed, 1e-10)
            # (b) one spacelike photon, one real photon, m_a = 0.
            closed = sp.N(result["single_spacelike_closed_form"].subs(tau_spacelike, tau_value), 25)
            earlier = sp.N(1 - limits["single_closed_form"].subs({quark_mass: 1, spacelike_magnitude: magnitude}), 25)
            self.assert_close(closed, earlier, 1e-20)
            self.assert_close(closed, literature_b1_numeric(mpmath.mpf(sp.N(tau_value, 25))), 1e-20)
        # (b) timelike below threshold, and (c) real photons: B_1 of the only finite argument.
        for invariant in (mpmath.mpf("0.7"), mpmath.mpf("3.2")):
            expected = literature_b1_numeric(4 / invariant)
            self.assert_close(b4_numeric(1, invariant, 0, 0), expected, 1e-10)
            self.assert_close(b4_numeric(1, 0, invariant, 0), expected, 1e-10)
            self.assert_close(b4_numeric(1, 0, 0, invariant), expected, 1e-10)

    def test_limits(self):
        result = self.result
        for name in ("heavy_expansion_check", "heavy_expansion_in_tau_check"):
            self.assertEqual(result[name], 0, name)
        self.assertEqual(tuple(result["massless_limits"]), (1, 1, 1))
        logarithm = sp.log(spacelike_magnitude / quark_mass**2)
        self.assertEqual(sp.simplify(result["equal_first_correction"]
                                     + 2 * quark_mass**2 / spacelike_magnitude * logarithm), 0)
        self.assertEqual(sp.simplify(result["single_first_correction"]
                                     + quark_mass**2 / spacelike_magnitude * logarithm**2), 0)
        self.assertEqual(sp.simplify(result["alp_mass_first_correction"]
                                     + alp_mass**2 / (3 * spacelike_magnitude)), 0)
        # m_h = 10 GeV: the two terms of the expansion against the integral.
        point = {quark_mass: 10, alp_mass: sp.sqrt(sp.Rational(1, 2)), first_virtuality: -1, second_virtuality: -2}
        numeric = b4_numeric(100, -1, -2, mpmath.mpf("0.5"))
        self.assert_close(numeric, float(result["heavy_expansion_expected"].subs(point)), 1e-7)
        first_term = -(sp.Rational(1, 2) - 1 - 2) / 1200
        self.assertGreater(abs(numeric - float(first_term)), 1e-7)
        self.assert_close(numeric, float(first_term), 1e-5)
        # Small m_h^2/|q^2| (m_h = 1 GeV, m_a^2 = 0.5 GeV^2): B_4 approaches 1, and the
        # distance to the first corrections falls faster than the corrections.
        previous_distance = None
        for magnitude in (100, 10000):
            point = {quark_mass: 1, alp_mass: sp.sqrt(sp.Rational(1, 2)), spacelike_magnitude: magnitude}
            equal = b4_numeric(1, -magnitude, -magnitude, mpmath.mpf("0.5"))
            leading = float((1 + result["equal_first_correction"] + result["alp_mass_first_correction"]).subs(point))
            self.assertLess(abs(equal - leading), 60.0 / magnitude**2 * mpmath.log(magnitude))
            single = b4_numeric(1, -magnitude, 0, 0)
            leading_single = float((1 + result["single_first_correction"]).subs(point))
            self.assertLess(abs(single - leading_single), 60.0 / magnitude**2 * mpmath.log(magnitude))
            distance = abs(1 - equal)
            if previous_distance is not None:
                self.assertLess(distance, previous_distance / 20)
            previous_distance = distance
        self.assertLess(previous_distance, 2e-3)

    def test_imaginary_part(self):
        result = self.result
        self.assertEqual(result["real_photons_parts_check"], 0)
        mpmath.mp.dps = 25
        # Real photons above threshold: closed form, sign, and the logarithmic integral of Part 7.
        for threshold_ratio in ("0.05", "0.25", "0.8", "0.99"):
            tau = mpmath.mpf(threshold_ratio)
            literature = literature_b1_numeric(tau)
            point = {quark_mass: 1, alp_mass: sp.sqrt(4 / sp.Rational(threshold_ratio))}
            self.assert_close(literature.real, float(result["real_photons_real_part"].subs(point)), 1e-12)
            self.assert_close(literature.imag, float(result["real_photons_imaginary_part"].subs(point)), 1e-12)
            self.assertLess(literature.imag, 0)
            velocity = mpmath.sqrt(1 - tau)
            imaginary_numeric = -2 * mpmath.pi * tau / 4 * mpmath.quad(
                lambda x: 1 / x, [(1 - velocity) / 2, (1 + velocity) / 2])
            self.assert_close(imaginary_numeric, literature.imag, 1e-18)
        # Below every threshold B_4 is real.
        self.assertEqual(mpmath.im(literature_b3_numeric(mpmath.mpf(2), mpmath.mpf(-3))), 0)
        # Two-dimensional integral at m^2 = (1 - 0.3 i) GeV^2 and its conjugate prescription.
        mass_squared = 1 - mpmath.mpc(0, "0.3")
        value = b4_numeric(mass_squared, 0, 0, 16)
        self.assert_close(value, literature_b1_numeric(4 * mass_squared / 16), 1e-9)
        self.assertLess(value.imag, 0)
        # Spacelike photons, m_a above threshold: sign of the closed form of Part 7.
        point = {quark_mass: 1, spacelike_magnitude: 3, alp_mass: 4}
        self.assertLess(float(result["spacelike_photons_imaginary_part"].subs(point)), 0)


def main():
    pseudoscalar = derive_pseudoscalar_triangle()
    properties = derive_loop_integral_properties()
    derivative = derive_derivative_triangle()
    combined = derive_combined_form_factor()
    limits = derive_kinematic_limits()
    matrix_element = derive_axial_current_matrix_element()
    vector = derive_vector_triangle()
    above = derive_above_threshold()
    b4 = derive_b4_properties()

    print("One heavy quark h; metric (+---), epsilon^{0123}=+1, gamma5=i gamma0 gamma1 gamma2 gamma3;")
    print("i M = i e^2 R F, R = epsilon_{mu nu rho sigma} eps1^mu eps2^nu q1^rho q2^sigma; F in 1/GeV.\n")

    print("PART 1: pseudoscalar coupling")
    show("numerator trace of either ordering, mu=1 nu=2 (equals 4 i m_h epsilon^{12 rho sigma} q1_rho q2_sigma)",
         pseudoscalar["numerator_first_ordering"][1][2])
    show("numerator contracted with the polarizations, divided by R", pseudoscalar["numerator_coefficient"])
    show("Feynman-parameter denominator (x: quark line between photon 2 and ALP, pairs with q_2^2; "
         "y: quark line between ALP and photon 1, pairs with q_1^2)",
         feynman_denominator(quark_mass, first_virtuality, second_virtuality, alp_momentum_squared))
    show("int d^4l/(2 pi)^4 (l^2 - Delta)^{-3}", pseudoscalar["momentum_integral"])
    show("dimensionless loop integral", loop_integral())
    show("overall sign relative to 2 N_c Q_h^2 c_hh^P/(4 pi^2 f_a) times the loop integral",
         pseudoscalar["overall_sign"])
    show("form factor from the pseudoscalar coupling", pseudoscalar["pseudoscalar_form_factor"])
    show("form factor of c_gamma (alpha/4 pi)(a/f_a) F Ftilde with epsilon^{0123}=+1 in the dual",
         pseudoscalar["local_photon_form_factor"])
    show("loop integral minus itself with photons exchanged (integrand, x <-> y)",
         properties["photon_exchange_difference"])
    show("loop integral for m_h -> infinity", properties["heavy_mass_limit"])
    show("first correction in 1/m_h^2", properties["heavy_mass_first_correction"])
    show("second correction in 1/m_h^4", properties["heavy_mass_second_correction"])
    show("loop integral for real photons, 4 m_h^2 >= m_a^2", properties["real_photon_closed_form"])
    print("m_h -> 0: the loop integral vanishes in the closed forms of Part 4 (equal spacelike")
    print("virtualities and one real photon). For real photons above threshold see Part 7.\n")

    print("PART 2: derivative coupling, Pauli-Villars regulator")
    show("pslash gamma5 - (kslash+pslash-m) gamma5 - gamma5 (kslash-m) - 2 m gamma5", derivative["vertex_identity"])
    show("two-propagator numerator, mu=1 nu=2 (mass independent, linear in k)", derivative["remainder_one"][1][2])
    show("form factor from the derivative coupling", derivative["axial_form_factor"])
    show("its massless limit minus the reference 2 N_c Q_h^2 c_hh^A/(4 pi^2 f_a)",
         sp.simplify(derivative["axial_massless_limit"] - derivative["reference_massless"]))

    print("PART 3: sum")
    show("heavy-quark form factor", combined["combined_form_factor"])
    show("change of c_gamma/(4 pi^2 f_a) + form factor under the axial rotation", combined["rotation_difference"])
    show("m_h -> infinity", combined["decoupling_limit"])
    show("m_h -> 0", combined["massless_limit"])
    show("real photons", combined["real_photon_limit"])

    print("PART 4: kinematic limits of the loop integral")
    show("(a) q_1^2 = q_2^2 = -|q^2|, p^2 = 0: one-dimensional integral", limits["equal_virtuality_integral"])
    show("(a) closed form", limits["equal_closed_form"])
    show("(a) expansion for |q^2| >> m_h^2", limits["equal_expansion"])
    show("(b) term linear in m_a^2: integral", limits["alp_mass_correction_integral"])
    show("(b) term linear in m_a^2: closed form", limits["alp_mass_correction_closed"])
    show("(b) term linear in m_a^2 for |q^2| >> m_h^2 (requires m_a^2 < 4 m_h^2)",
         limits["alp_mass_correction_expansion"])
    show("(c) q_1^2 = -|q^2|, q_2^2 = 0, p^2 = 0: closed form", limits["single_closed_form"])
    show("(c) expansion for |q^2| >> m_h^2", limits["single_expansion"])
    show("(d) |q_i^2|, |p^2| << m_h^2", limits["low_energy_expansion"])
    show("(d) next term", limits["low_energy_next_term"])

    print("PART 5: local expansion with the one-loop axial-current matrix element")
    show("coefficient of p^sigma in <0|hbar gamma^sigma gamma5 h|a(p)> (Euclidean cutoff Lambda)",
         matrix_element["matrix_element_coefficient"])
    show("the same to first order in m_a^2/m_h^2", sp.factor(matrix_element["matrix_element_coefficient_expanded"]))
    show("local term 2 i Q_h^2 p_sigma<...>/((q^2-m_h^2) m_a^2) at q^2 = -|q^2|", matrix_element["local_form_factor"])
    show("local term, leading power in 1/|q^2|", matrix_element["local_leading"])
    show("direct triangle, mass-dependent part at leading power in 1/|q^2|", matrix_element["direct_leading"])
    show("direct minus local", matrix_element["remainder_expected"])
    print("The local term reproduces the sign, the factor 2 m_h^2/|q^2|, the logarithm of m_h^2 and the")
    print("term m_a^2/(3|q^2|). The difference is polynomial in m_h^2; its constant depends on the cutoff")
    print("prescription of the matrix element and belongs to loop momenta of order |q^2|, not derived here.")


    print("\nPART 6: vector coupling (d_mu a/f_a) c_hh^V hbar gamma^mu h, vertex (c_hh^V/f_a) pslash")
    show("numerator trace of the first ordering, mu=1 nu=2 (real, no Levi-Civita tensor)",
         vector["numerator_first_ordering"][1][2])
    show("part of the first-ordering numerator that is odd under space reflection, mu=1 nu=2",
         vector["odd_part_first_ordering"][1][2])
    show("coefficient of epsilon^{mu nu rho sigma} q1_rho q2_sigma, first ordering",
         vector["levi_civita_coefficient_first"])
    show("coefficient of epsilon^{mu nu rho sigma} q1_rho q2_sigma, second ordering",
         vector["levi_civita_coefficient_second"])
    show("the same projection applied to the pseudoscalar numerator of Part 1",
         vector["pseudoscalar_levi_civita_coefficient"])
    show("first-ordering numerator at k plus second-ordering numerator at -k-p, mu=1 nu=2",
         vector["ordering_sum"][1][2])
    show("product of the three propagator denominators: first ordering at k minus second ordering at -k-p",
         vector["denominator_difference"])
    show("form factor from the vector coupling", vector["vector_form_factor"])
    print("The cancellation of the two orderings changes the loop momentum in a linearly divergent")
    print("integral; it holds in a regularization that allows this change (assumption stated above).\n")

    print("PART 7: m_a > 2 m_h, denominator with m_h^2 -> m_h^2 - i0 (from the propagator 1/(k^2 - m_h^2 + i0))")
    print("Sign convention: 1/(D - i0) = principal value 1/D + i pi delta(D), so the imaginary part of the")
    print("loop integral is positive; the form factor has the imaginary part")
    print("-2 N_c Q_h^2 (c_hh^A - c_hh^P)/(4 pi^2 f_a) times it.\n")
    show("loop integral for real photons, m_a > 2 m_h", above["above_closed_form"])
    show("its real part", above["above_real_part"])
    show("its imaginary part", above["above_imaginary_part"])
    show("its limit m_h -> 0", above["massless_limit_above_threshold"])
    distance = above["distance"]
    threshold_ratio = 4 * quark_mass**2 / alp_mass**2
    show("m_a slightly below 2 m_h, up to terms of order (4 m_h^2/m_a^2 - 1)^(3/2)",
         above["below_series"].subs(distance, sp.sqrt(threshold_ratio - 1)))
    show("m_a slightly above 2 m_h, up to terms of order (1 - 4 m_h^2/m_a^2)^(3/2)",
         above["above_series"].subs(distance, sp.sqrt(1 - threshold_ratio)))
    show("sqrt(4 m_h^2/m_a^2 - 1) times the derivative with respect to 4 m_h^2/m_a^2, at the threshold from below",
         above["below_slope_times_distance"])
    show("sqrt(1 - 4 m_h^2/m_a^2) times the derivative of the imaginary part with respect to 4 m_h^2/m_a^2, "
         "at the threshold from above", above["above_imaginary_slope_times_distance"])
    show("derivative of the real part with respect to 4 m_h^2/m_a^2 at the threshold, from above threshold",
         above["above_real_slope_at_threshold"])
    print("The value is continuous, pi^2/4 on both sides; the derivative diverges as the inverse square root")
    print("of the distance to the threshold (real part on the side m_a < 2 m_h, imaginary part on the side")
    print("m_a > 2 m_h).\n")
    show("imaginary part for q_1^2 = q_2^2 = -|q^2|, p^2 = m_a^2 > 4 m_h^2: one-dimensional integral",
         above["spacelike_imaginary_integral"])
    show("the same in closed form", above["spacelike_imaginary_closed_form"])

    print("PART 8: B_4(tau_a, tau_1, tau_2), tau_a = 4 m_h^2/m_a^2, tau_1 = 4 m_h^2/q_1^2, tau_2 = 4 m_h^2/q_2^2")
    print("(real photon: tau_i -> infinity; m_h^2 -> m_h^2 - i0 in the denominator). B_3 and B_1 are the")
    print("functions of Bauer:2021mvw Eqs. (2.82), (2.71) and Alda:2024cxn Eq. (A1).\n")
    b4_function, b3_function, b1_function, f_function = (
        b4["b4_function"], b4["b3_function"], b4["b1_function"], b4["f_function"])
    show("definition", b4["b4_definition_physical"])
    show("definition in the tau variables", b4["b4_definition"])
    show("heavy-quark form factor", b4["form_factor_definition"])
    show("this form minus the form factor of Part 3", b4["equivalence_with_derived_sum"])
    show("symmetry under the exchange of any two arguments", sp.Eq(
        sp.Eq(b4_function(tau_alp, tau_first_photon, tau_second_photon),
              b4_function(tau_first_photon, tau_alp, tau_second_photon), evaluate=False),
        sp.Eq(b4_function(tau_second_photon, tau_first_photon, tau_alp),
              b4_function(tau_alp, tau_second_photon, tau_first_photon), evaluate=False), evaluate=False))
    for name, checks in b4["permutation_checks"].items():
        print(f"{name}: integrand with exchanged arguments minus integrand with exchanged Feynman parameters = "
              f"{checks[0]}; |Jacobian| - 1 = {checks[1]}")
    print()
    generic_tau = sp.Symbol("tau", real=True)
    generic_root = sp.sqrt(1 - generic_tau)
    show("f(tau), tau >= 1", sp.Eq(f_function(generic_tau), sp.asin(1 / sp.sqrt(generic_tau))))
    show("f(tau), 0 < tau < 1", sp.Eq(f_function(generic_tau), sp.pi / 2 + sp.I / 2 * sp.log(
        (1 + generic_root) / (1 - generic_root))))
    show("f(tau), tau < 0 (derived here from Part 1; f^2 is real and negative)",
         sp.Eq(f_function(tau_spacelike), b4["f_spacelike"]))
    show("one real photon", sp.Eq(b4_function(tau_alp, tau_first_photon, sp.oo),
                                  sp.Eq(b3_function(tau_alp, tau_first_photon), b4["b3_printed"], evaluate=False),
                                  evaluate=False))
    show("two real photons", sp.Eq(b4_function(tau_alp, sp.oo, sp.oo),
                                   sp.Eq(b1_function(tau_alp), b4["b1_printed"], evaluate=False), evaluate=False))
    show("one real photon: triangle integrand after x = (1-y) xi and the average with y -> 1-y "
         "(integrand of Bauer:2021mvw Eq. (2.73) with k^2 = m_a^2, p^2 = q_1^2)",
         b4["one_real_photon_literature_integrand"])
    show("B_4(infinity, tau, tau) = B_3(tau, tau), tau < 0",
         sp.Eq(b4_function(sp.oo, tau_spacelike, tau_spacelike), b4["equal_spacelike_closed_form"]))
    show("B_4(infinity, tau, infinity) = B_1(tau), tau < 0",
         sp.Eq(b4_function(sp.oo, tau_spacelike, sp.oo), b4["single_spacelike_closed_form"]))
    show("B_4(tau_a, infinity, infinity) = B_1(tau_a), tau_a >= 1", b4["real_photons_below_threshold"])
    show("B_4(tau_a, infinity, infinity) = B_1(tau_a), 0 < tau_a < 1", b4["real_photons_above_threshold"])
    show("its real part", b4["real_photons_real_part"])
    show("its imaginary part (negative)", b4["real_photons_imaginary_part"])
    show("imaginary part of B_4 for q_1^2 = q_2^2 = -|q^2|, m_a > 2 m_h", b4["spacelike_photons_imaginary_part"])
    show("B_4 for m_h -> infinity, up to terms of order 1/m_h^6", b4["heavy_expansion_expected"])
    show("the same in the tau variables", b4["heavy_expansion_in_tau"])
    show("B_4 for m_h -> 0 (equal spacelike virtualities; one real photon; real photons)",
         sp.Tuple(*b4["massless_limits"]))
    show("B_4 - 1 for |q^2| >> m_h^2, q_1^2 = q_2^2 = -|q^2|, m_a = 0: first term", b4["equal_first_correction"])
    show("B_4 - 1 for |q^2| >> m_h^2, q_1^2 = -|q^2|, q_2^2 = 0, m_a = 0: first term",
         b4["single_first_correction"])
    show("B_4 - 1 for |q^2| >> m_h^2, q_1^2 = q_2^2 = -|q^2|: first term in m_a^2 (m_a < 2 m_h)",
         b4["alp_mass_first_correction"])
    print("B_4 depends only on the three ratios, so |q^2| -> infinity at fixed m_h and m_h -> 0 at fixed")
    print("|q^2| are the same limit when m_a = 0; in both B_4 -> 1 and F_h -> 2 N_c Q_h^2 c_hh^A/(4 pi^2 f_a).")



if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
