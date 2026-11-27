"""Cartesian cross-check of the weights of the radial-mode loop-size estimate.

Task 6.1 of the OpenSpec change phi-dependent-loop-tree-crossing-murayama.

Run:
    .venv/bin/python scripts/derive_murayama_radial_cartesian_weights.py
Tests:
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_murayama_radial_cartesian_weights

What is checked.  The weight given by the size prescription
(memory/alp/radial-mode/finite-mixing/loop-estimates/power-counting/murayama-radial-loop-size-prescription.md) to a loop
class is
    (epsilon/sqrt(N(m_a^2))) * epsilon^p N(m_a^2)^(-n) * m^2/(16 pi^2 f_pi^2) * |c_k|,
where c_k is built from tree-level vertices with legs on exact poles.  This
script computes those tree-level vertices in polar fields (sigma, a) and in
Cartesian fields (s, p), with the exact pole projectors of the up-only model,
and compares them.  It does not compute any loop integral, scalar loop
function, loop coefficient, combinatorial factor of a loop diagram, or loop
amplitude: only vertices and one tree-level four-point amplitude (control).

Model and conventions (memory/alp/radial-mode/couplings-and-decays/finite-mixing/murayama-sigma-two-pseudoscalar-decays.md,
memory/alp/decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-polar-cartesian-amplitude-equivalence.md):
  * Up-only model, f_a = f_phi, theta = 0, Theta_hat = 0.
    f_a phi = (f_a + sigma) exp(i a/f_a)/sqrt(2) = (f_a + s + i p)/sqrt(2).
  * Only the up-quark entry Phi_uu of the neutral meson matrix couples to a,
    sigma (or p, s).  With U = exp(i Phi/f_pi) and diagonal Phi,
        polar:     L ⊃ B_0 f_pi^2 m_u (1 + sigma/f_a) cos(a/f_a + Phi_uu/f_pi)
                       + (1/2)(1 + sigma/f_a)^2 (d a)^2 - V(f_a + sigma)
        Cartesian: L ⊃ B_0 f_pi^2 m_u [(1 + s/f_a) cos(Phi_uu/f_pi) - (p/f_a) sin(Phi_uu/f_pi)]
                       + (1/2)(d p)^2 - V(sqrt((f_a + s)^2 + p^2))
    with V(rho) = m_phi^2 rho^2/2 + lambda rho^4/4 and, from stationarity,
    m_phi^2 = B_0 m_u epsilon^2 - lambda f_a^2, epsilon = f_pi/f_a.  The
    d and s quark parts and the meson kinetic term are the same in both
    field choices and drop out of every difference.
  * Exact pole projectors (no expansion in epsilon): leg r has ALP component
    1/sqrt(N_r) and Phi_uu component 2 B_0 epsilon X_u(m_r^2)/sqrt(N_r), with
    X_u(m_r^2) = (m_r^2/(B_0 m_u epsilon^2) - 1)/(2 B_0) from the pole equation
    (derive_murayama_radial_up_only_specialization.up_entry_at_pole).  The
    quadratic Lagrangians of the two field choices coincide (checked), so the
    projectors are the same.  N_r stays a symbol.
  * Vertex function: coefficient of the product of the leg fields, summed
    over assignments; d_mu -> i p_mu with all momenta incoming, so that two
    derivatives give -p_x.p_y.  Pair invariant w_xy = (p_x + p_y)^2.

Units: GeV for masses, f_pi, f_a, B_0; GeV^2 for squared masses and momenta.
"""
import itertools

import sympy as sp

from derive_murayama_radial_crossing import B_ZERO, EPSILON, F_PI, QUARTIC_COUPLING
from derive_murayama_radial_up_only_specialization import UP_QUARK_MASS, up_entry_at_pole

F_A = F_PI / EPSILON
RADIAL_MASS_SQUARED_UP_ONLY = B_ZERO * UP_QUARK_MASS * EPSILON**2 + 2 * QUARTIC_COUPLING * F_A**2
ANGULAR_MASS_PARAMETER = B_ZERO * UP_QUARK_MASS * EPSILON**2 - QUARTIC_COUPLING * F_A**2  # m_phi^2
POLE_MASSES = sp.symbols('m_1^2 m_2^2 m_3^2 m_4^2', positive=True)
NORMALIZATIONS = sp.symbols('N_1 N_2 N_3 N_4', positive=True)
LEG_FIELDS = sp.symbols('P_1 P_2 P_3 P_4', real=True)
RADIAL_LEG_FIELDS = sp.symbols('R_1 R_2 R_3 R_4', real=True)
ANGLE, RADIAL, MESON = sp.symbols('a sigma Phi_uu', real=True)


def alp_component(leg):
    return 1 / sp.sqrt(NORMALIZATIONS[leg])


def meson_component(leg):
    """Phi_uu component of the exact pole projector of leg."""
    return 2 * B_ZERO * EPSILON * up_entry_at_pole(POLE_MASSES[leg]) / sp.sqrt(NORMALIZATIONS[leg])


def potential(radius_squared):
    return ANGULAR_MASS_PARAMETER * radius_squared / 2 + QUARTIC_COUPLING * radius_squared**2 / 4


def polar_lagrangian(alp, radial, meson):
    """Non-derivative a- and sigma-dependent part (up entry)."""
    return (B_ZERO * F_PI**2 * UP_QUARK_MASS * (1 + radial / F_A) * sp.cos(alp / F_A + meson / F_PI)
            - potential((F_A + radial)**2))


def cartesian_lagrangian(angular, radial, meson):
    """Non-derivative p- and s-dependent part (up entry)."""
    return (B_ZERO * F_PI**2 * UP_QUARK_MASS * ((1 + radial / F_A) * sp.cos(meson / F_PI)
                                               - angular / F_A * sp.sin(meson / F_PI))
            - potential((F_A + radial)**2 + angular**2))


def polar_kinetic_prefactor(radial):
    """Coefficient of (d a)^2 in polar fields: (1/2)(1 + sigma/f_a)^2."""
    return (1 + radial / F_A)**2 / 2


# ---------------------------------------------------------------------------
# Field map checks.
# ---------------------------------------------------------------------------

def field_map_residuals():
    """Polar minus Cartesian Lagrangian after s, p -> functions of sigma, a.

    Returns (non-derivative residual, kinetic residual); both zero."""
    radial_cartesian = (F_A + RADIAL) * sp.cos(ANGLE / F_A) - F_A
    angular_cartesian = (F_A + RADIAL) * sp.sin(ANGLE / F_A)
    potential_residual = sp.simplify(sp.expand_trig(
        polar_lagrangian(ANGLE, RADIAL, MESON)
        - cartesian_lagrangian(angular_cartesian, radial_cartesian, MESON)))
    radial_derivative, angle_derivative = sp.symbols('dsigma da', real=True)
    jacobian = sp.Matrix([radial_cartesian, angular_cartesian]).jacobian([RADIAL, ANGLE])
    cartesian_derivatives = jacobian * sp.Matrix([radial_derivative, angle_derivative])
    kinetic_residual = sp.simplify(
        (cartesian_derivatives[0]**2 + cartesian_derivatives[1]**2) / 2
        - radial_derivative**2 / 2 - polar_kinetic_prefactor(RADIAL) * angle_derivative**2)
    return potential_residual, kinetic_residual


def quadratic_part_residual():
    """Quadratic parts of the two Lagrangians, fields identified linearly; zero.

    Also checks that the linear (tadpole) terms of sigma and s vanish."""
    point = {ANGLE: 0, RADIAL: 0, MESON: 0}
    fields = (ANGLE, RADIAL, MESON)
    polar = polar_lagrangian(*fields)
    cartesian = cartesian_lagrangian(*fields)
    hessian_difference = sp.simplify(sp.hessian(polar, fields).subs(point)
                                     - sp.hessian(cartesian, fields).subs(point))
    tadpoles = [sp.simplify(sp.diff(lagrangian, RADIAL).subs(point)) for lagrangian in (polar, cartesian)]
    return hessian_difference, tadpoles


# ---------------------------------------------------------------------------
# Vertex functions with exact pole projectors.
# ---------------------------------------------------------------------------

def non_derivative_vertex(lagrangian, pole_legs, radial_legs):
    """Coefficient of prod(pole leg fields) prod(radial leg fields), summed over assignments."""
    alp = sum(alp_component(leg) * LEG_FIELDS[leg] for leg in pole_legs)
    meson = sum(meson_component(leg) * LEG_FIELDS[leg] for leg in pole_legs)
    radial = sum(RADIAL_LEG_FIELDS[k] for k in range(radial_legs))
    expression = lagrangian(alp, radial, meson)
    variables = [LEG_FIELDS[leg] for leg in pole_legs] + list(RADIAL_LEG_FIELDS[:radial_legs])
    zero = {variable: 0 for variable in variables}
    for variable in variables:
        expression = sp.diff(expression, variable)
    return sp.simplify(expression.subs(zero))


def polar_kinetic_vertex(first, second, radial_legs, momentum_product):
    """(d a)^2 vertex: d^n/dsigma^n of (1/2)(1+sigma/f_a)^2 times 2 a_x a_y (-p_x.p_y)."""
    prefactor = sp.diff(polar_kinetic_prefactor(RADIAL), RADIAL, radial_legs).subs(RADIAL, 0)
    return -2 * prefactor * alp_component(first) * alp_component(second) * momentum_product


def on_shell_momentum_product(first, second, pair_invariant):
    """p_x.p_y = (w_xy - m_x^2 - m_y^2)/2 for legs on their poles."""
    return (pair_invariant - POLE_MASSES[first] - POLE_MASSES[second]) / 2


def radial_two_pole_vertices(momentum_squared):
    """sigma P_1 P_2 at off-shell radial momentum q^2: (polar, Cartesian)."""
    polar = (non_derivative_vertex(polar_lagrangian, (0, 1), 1)
             + polar_kinetic_vertex(0, 1, 1, on_shell_momentum_product(0, 1, momentum_squared)))
    cartesian = non_derivative_vertex(cartesian_lagrangian, (0, 1), 1)
    return sp.simplify(polar), sp.simplify(cartesian)


def two_radial_two_pole_vertices(momentum_product):
    """sigma sigma P_1 P_2 with p_1.p_2 kept as a symbol: (polar, Cartesian)."""
    polar = (non_derivative_vertex(polar_lagrangian, (0, 1), 2)
             + polar_kinetic_vertex(0, 1, 2, momentum_product))
    cartesian = non_derivative_vertex(cartesian_lagrangian, (0, 1), 2)
    return sp.simplify(polar), sp.simplify(cartesian)


def four_pole_vertex_difference():
    """Polar minus Cartesian four-pole contact vertex (no derivative terms differ)."""
    legs = (0, 1, 2, 3)
    return sp.factor(sp.simplify(non_derivative_vertex(polar_lagrangian, legs, 0)
                                 - non_derivative_vertex(cartesian_lagrangian, legs, 0)))


def radial_self_coupling_differences():
    """sigma^n self-couplings, n = 1..4, polar minus Cartesian; all zero."""
    return {n: sp.simplify(non_derivative_vertex(polar_lagrangian, (), n)
                           - non_derivative_vertex(cartesian_lagrangian, (), n)) for n in range(1, 5)}


def odd_vertices():
    """Vertices with an odd number of pole legs, both field choices; all zero."""
    results = {}
    for pole_legs, radial_legs in (((0,), 1), ((0,), 2), ((0, 1, 2), 0), ((0, 1, 2), 1)):
        for name, lagrangian in (('polar', polar_lagrangian), ('Cartesian', cartesian_lagrangian)):
            results[(name, len(pole_legs), radial_legs)] = non_derivative_vertex(
                lagrangian, pole_legs, radial_legs)
    return results


# ---------------------------------------------------------------------------
# Control: on-shell four-pole amplitude with radial exchange.
# ---------------------------------------------------------------------------

def four_pole_amplitude_difference():
    """Polar minus Cartesian tree amplitude of four pole legs, on shell.

    M = V4 + sum over the three channels of V_sigma(x,y; w) V_sigma(z,t; w) / (m_sigma^2 - w).
    Only the parts that differ are kept (the meson-only parts are equal).
    Pole exchanges vanish in both field choices (odd_vertices).  The third
    pair invariant follows from w_12 + w_13 + w_14 = sum m_r^2."""
    w_12, w_13 = sp.symbols('w_12 w_13', real=True)
    w_14 = sum(POLE_MASSES) - w_12 - w_13
    total = four_pole_vertex_difference()
    for (pair, other), invariant in zip((((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))),
                                        (w_12, w_13, w_14)):
        polar = [non_derivative_vertex(polar_lagrangian, legs, 1)
                 + polar_kinetic_vertex(*legs, 1, on_shell_momentum_product(*legs, invariant))
                 for legs in (pair, other)]
        cartesian = [non_derivative_vertex(cartesian_lagrangian, legs, 1) for legs in (pair, other)]
        total += (polar[0] * polar[1] - cartesian[0] * cartesian[1]) / (
            RADIAL_MASS_SQUARED_UP_ONLY - invariant)
    return sp.simplify(sp.together(total))


# ---------------------------------------------------------------------------
# Weights of the estimate.
# ---------------------------------------------------------------------------

def common_factor(pole_legs, inverse_f_a_power):
    """prod over pole legs of 1/sqrt(N_r), times f_a^-k: the factor that carries
    the pole normalizations and the epsilon label of the weight."""
    return sp.Mul(*[alp_component(leg) for leg in pole_legs]) / F_A**inverse_f_a_power


def weight_comparison():
    """Per vertex: (common factor, polar bracket, Cartesian bracket, difference).

    bracket = vertex / common factor.  The common factor is the same in both
    field choices; a bracket free of N_r confirms it carries all pole
    normalizations."""
    momentum_squared = sp.Symbol('q^2', real=True)
    momentum_product = sp.Symbol('p_1.p_2', real=True)
    rows = {}
    factor = common_factor((0, 1), 1)
    polar, cartesian = radial_two_pole_vertices(momentum_squared)
    rows['sigma P P at q^2 (tree exchange, R3, R4a, R4b)'] = (
        factor, sp.factor(polar / factor), sp.factor(cartesian / factor),
        sp.factor((polar - cartesian) / factor))
    factor = common_factor((0, 1), 2)
    polar, cartesian = two_radial_two_pole_vertices(momentum_product)
    rows['sigma sigma P P (R2)'] = (
        factor, sp.factor(polar / factor), sp.factor(cartesian / factor),
        sp.factor((polar - cartesian) / factor))
    factor = common_factor((0, 1, 2, 3), 2)
    difference = four_pole_vertex_difference()
    rows['P P P P contact (classes 2, 3, 6 with ALP-like legs)'] = (
        factor, None, None, sp.factor(sp.simplify(difference / factor)))
    return rows


def radial_alp_self_energy_weight():
    """Listing only: derivative sigma sigma A A vertex of the polar kinetic term
    with both pole legs on the ALP-like pole (leg 1), momenta p and -p
    incoming, so p_x.p_y = -p_a^2, on shell p_a^2 = m_1^2.  A closed sigma line
    on it is the radial self-energy of the ALP-like leg (classes R5a, R5b of
    derive_murayama_radial_crossing.py).  Returns (common factor, polar
    bracket).  The Cartesian field choice has no derivative sigma sigma A A
    vertex; its counterpart is part of the non-derivative s s P P vertex."""
    vertex = polar_kinetic_vertex(0, 0, 2, -POLE_MASSES[0])
    factor = common_factor((0, 0), 2)
    return factor, sp.factor(vertex / factor)


def main():
    potential_residual, kinetic_residual = field_map_residuals()
    print('Field map polar <-> Cartesian: non-derivative residual', potential_residual,
          '; kinetic residual', kinetic_residual)
    hessian_difference, tadpoles = quadratic_part_residual()
    print('Quadratic parts, polar minus Cartesian:', hessian_difference.tolist(),
          '; radial tadpoles (polar, Cartesian):', tadpoles)
    print('Radial self-couplings, polar minus Cartesian:', radial_self_coupling_differences())
    print('Vertices with an odd number of pole legs:', set(odd_vertices().values()))
    print('\nWeights: common factor prod 1/sqrt(N_r) / f_a^k, and bracket = vertex / factor')
    for name, (factor, polar, cartesian, difference) in weight_comparison().items():
        print(f'  {name}:')
        print('      common factor:', sp.latex(factor))
        if polar is not None:
            print('      polar bracket:    ', sp.latex(polar))
            print('      Cartesian bracket:', sp.latex(cartesian))
        print('      polar - Cartesian:', sp.latex(difference))
    factor, bracket = radial_alp_self_energy_weight()
    print('  sigma sigma A A, derivative vertex, both legs on the ALP-like pole, p_a^2 = m_1^2')
    print('  (radial self-energy of the ALP-like leg, R5a/R5b; listed, not compared):')
    print('      common factor:', sp.latex(factor))
    print('      polar bracket:    ', sp.latex(bracket))
    print('      Cartesian: no derivative vertex; counterpart inside the non-derivative s s P P vertex')
    print('\nm_sigma^2 (up-only) =', sp.latex(RADIAL_MASS_SQUARED_UP_ONLY))
    print('Control: on-shell four-pole amplitude, polar minus Cartesian:',
          four_pole_amplitude_difference())
    print('\nThe common factors (pole normalizations and powers of 1/f_a) agree; the')
    print('brackets of the off-shell vertices differ.  The weights |c_k| of the')
    print('estimate therefore depend on the field choice; only the on-shell sum agrees.')


if __name__ == '__main__':
    main()
