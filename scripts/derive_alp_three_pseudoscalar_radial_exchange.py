#!/usr/bin/env python3
"""Exact LO tree amplitude a_phys -> P_i P_j P_k with radial-mode (sigma) exchange.

Model: DiLuzio:2026poh (dynamical Peccei-Quinn field), polar coordinates, general Hermitian
source coefficients, conventions of Alda:2026nrd and of
memory/alp/shared/radial-decay-polar-cartesian-general-sources.md:

* M(a, sigma) = S + i P, S = sum_{r,s} S_rs (sigma/f_a)^r (a/f_a)^s (same for P); the
  coefficients are Hermitian 3x3 matrices of mass dimension one.
* Phi = sum_k P_k lambda_k, <lambda_k lambda_l> = 2 delta_kl, U = exp(i Phi/f_pi), chi = 2 B0 M,
  L = (f_pi^2/4) <d U d U^dagger + chi U^dagger + chi^dagger U> - (m_1^2/12) <Phi>^2
      + (1/2)(d sigma)^2 + (1/2)(1 + sigma/f_a)^2 (d a)^2 - V_phi,
  V_phi = (m_phi^2/2)(f_a + sigma)^2 + (lambda/4)(f_a + sigma)^4.
* Cartesian fields: f_a + s + i p = (f_a + sigma) exp(i a/f_a), L_kin = (1/2)(d s)^2 + (1/2)(d p)^2.
* M_s(X) = B0 {S_00, X} + (m_1^2/3) <X> 1 - s X.
* Vertices omit the Feynman-rule factor i; all-outgoing momenta, the parent carries -p.
  The tree amplitude is M = V4 + sum_channels V3 (M^2 - s K)^{-1} V3 (Alda:2026nrd, S = 1 + i T).

Units: GeV for masses, fields, B0, f_pi, f_a and the source coefficients; GeV^2 for s_ij.

Assumptions (labelled as in the memory record):
  1. P_00 = 0, <S_01> = 0; m_phi^2 fixed by radial stationarity.
  2. P_10 = 0, <S_11> = 0 (no quadratic radial mixing).
  3. The anomalous singlet term does not depend on sigma or a.
  4. Tree level, LO U(3), no expansion in f_pi/f_a, connection-free polar representative.
  5. Regular simple poles; the trace forms need a nonzero angular component of the leg.

Content:
  (1) sympy: no quadratic term in the symmetric chiral map; no quadratic radial mixing in polar
      coordinates; off-shell polar sigma-x-y vertex; Cartesian quartic contacts at s = 0; the
      mapping of the quadratic blocks onto the finite-mixing three-pseudoscalar record.
  (2) numerics: tree amplitude from the exact Lagrangian (polar or Cartesian, exact map expanded
      to quartic order) with Hessian poles and complete contraction, and the analytic result
      (record amplitude without sigma + trace-form sigma exchange).

Run:
    .venv/bin/python scripts/derive_alp_three_pseudoscalar_radial_exchange.py
Tests:
    .venv/bin/python -m unittest scripts.test_derive_alp_three_pseudoscalar_radial_exchange
"""
import itertools
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp

_SCRIPTS = str(Path(__file__).resolve().parent)
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

from derive_radial_decay_polar_cartesian_general import (  # noqa: E402
    LAMBDA, NUMBER_OF_FIELDS, ANGULAR_COMPONENT_THRESHOLD, B_0, f_a, f_pi, quartic,
    derive_general_cubic_lagrangians, impose_stationarity, inverse_mass_operator,
    mass_operator_matrix, pole_normalization, polar_vertex_trace_form, radial_mass_squared_trace_form,
    random_hermitian, trace)
from derive_alp_finite_mixing_three_pseudoscalars import (  # noqa: E402
    PAIRINGS, aligned_source_coefficients, compact_amplitude, component_quadratic_blocks)

MAXIMUM_SOURCE_ORDER = 4  # r + s <= 4 is complete for four-point vertices

# ---------------------------------------------------------------------------
# (1) Symbolic derivations
# ---------------------------------------------------------------------------


def derive_symmetric_map_second_order():
    """log(e^{tA} e^{tB} e^{tA}) to second order in t for non-commuting A, B.

    Returns the second-order coefficient (zero: the symmetric chiral map has no quadratic term).
    """
    first, second = sp.symbols('A B', commutative=False)
    t = sp.Symbol('t')

    def exponential(x):
        return 1 + t * x + t ** 2 * x ** 2 / 2

    product = sp.expand(exponential(first) * exponential(second) * exponential(first))
    deviation = product - 1
    logarithm = sp.expand(deviation - deviation ** 2 / 2)
    first_order = logarithm.coeff(t, 1)
    second_order = logarithm.coeff(t, 2)
    return sp.expand(first_order - 2 * first - second), sp.expand(second_order)


def _source_lagrangian(radial_field, angular_field, scale, orders, traces):
    """B0 f_pi^2 <S> + B0 f_pi <P Phi> - (B0/2) <S Phi^2> with trace symbols; scale counts Phi."""
    scalar, pseudoscalar_meson, scalar_meson_pair = traces
    total = 0
    for (r, s) in orders:
        monomial = (radial_field / f_a) ** r * (angular_field / f_a) ** s
        total += (B_0 * f_pi ** 2 * scalar[(r, s)] + B_0 * f_pi * scale * pseudoscalar_meson[(r, s)]
                  - B_0 / 2 * scale ** 2 * scalar_meson_pair[(r, s)]) * monomial
    return total


def _trace_symbols(maximum_order=MAXIMUM_SOURCE_ORDER):
    orders = [(r, s) for r in range(maximum_order + 1) for s in range(maximum_order + 1) if r + s <= maximum_order]
    scalar = {o: sp.Symbol(f'tr_S_{o[0]}{o[1]}', real=True) for o in orders}
    pseudoscalar_meson = {o: sp.Symbol(f'tr_P_{o[0]}{o[1]}_Phi', real=True) for o in orders}
    scalar_meson_pair = {o: sp.Symbol(f'tr_S_{o[0]}{o[1]}_PhiPhi', real=True) for o in orders}
    return orders, (scalar, pseudoscalar_meson, scalar_meson_pair)


def derive_polar_radial_mixing():
    """Quadratic sigma-a and sigma-Phi terms of the polar Lagrangian (trace symbols).

    The polar kinetic term (1/2)(1 + sigma/f_a)^2 (d a)^2 has no quadratic sigma-a part.
    Returns (coefficient of sigma a, dict of coefficients of sigma <P_rs Phi>).
    """
    orders, traces = _trace_symbols(2)
    sigma, phase, scale, count = sp.symbols('sigma a epsilon field_count', real=True)
    lagrangian = _source_lagrangian(count * sigma, count * phase, count * scale, orders, traces)
    quadratic = sp.expand(sp.series(lagrangian, count, 0, 3).removeO()).coeff(count, 2)
    sigma_angular = sp.simplify(quadratic.coeff(scale, 0).coeff(sigma, 1).coeff(phase, 1))
    sigma_meson = {o: sp.simplify(quadratic.coeff(scale, 1).coeff(traces[1][o]).coeff(sigma, 1))
                   for o in orders}
    return sigma_angular, {o: v for o, v in sigma_meson.items() if v != 0}


def derive_off_shell_radial_vertex():
    """Polar sigma-x-y vertex for general columns (Lambda_x, xi_x), (Lambda_y, xi_y), off shell.

    Uses the polar cubic coefficients of derive_general_cubic_lagrangians.  With all-outgoing
    momenta k_x, k_y and the derivative term sigma (d a)^2/f_a, the vertex is
    V = (coefficients) - (2/f_a) xi_x xi_y k_x.k_y, and g = -V.  Returns
    (f_a g with 2 k_x.k_y written through q^2, difference to the conjectured form, difference to
    the on-shell g of the sigma -> PP record with m_sigma^2 -> q^2, parent-pair check).
    """
    data = derive_general_cubic_lagrangians()
    anticommutator_s10 = sp.Symbol('tr_S10_anticomm_Lx_Ly', real=True)  # <S_10 {Lambda_x, Lambda_y}>
    p11_x, p11_y = sp.symbols('tr_P11_Lx tr_P11_Ly', real=True)
    xi_x, xi_y = sp.symbols('xi_x xi_y', real=True)
    s12 = data['trace_scalar'][(1, 2)]
    q_squared, m_x_squared, m_y_squared, radial_mass_squared = sp.symbols('q2 m_x2 m_y2 m_sigma2', real=True)
    two_kx_ky = q_squared - m_x_squared - m_y_squared
    vertex = (data['polar_sigma_meson_pair'][(1, 0)] * anticommutator_s10
              + data['polar_sigma_a_meson'][(1, 1)] * (xi_x * p11_y + xi_y * p11_x)
              + 2 * data['polar_sigma_a_a'] * xi_x * xi_y
              - data['polar_kinetic_sigma_da_da'] * xi_x * xi_y * two_kx_ky)
    f_a_g = sp.expand(-f_a * vertex)
    conjectured = (B_0 / 2 * anticommutator_s10 - B_0 * f_pi / f_a * (xi_x * p11_y + xi_y * p11_x)
                   - 2 * B_0 * f_pi ** 2 / f_a ** 2 * xi_x * xi_y * s12 + xi_x * xi_y * two_kx_ky)
    on_shell = conjectured.subs(q_squared, radial_mass_squared)
    # Parent pair: k_a = -p, k_k = p_k; q = p - p_k = p_i + p_j, q^2 = s_ij.
    m_a_squared, m_k_squared, p_dot_pk = sp.symbols('m_a2 m_k2 p_dot_pk', real=True)
    s_ij = m_a_squared + m_k_squared - 2 * p_dot_pk
    parent_pair = sp.expand(2 * (-1) * p_dot_pk - (s_ij - m_a_squared - m_k_squared))
    return (f_a_g, sp.expand(f_a_g - conjectured),
            sp.expand(f_a_g.subs(q_squared, radial_mass_squared) - on_shell), parent_pair)


def derive_coordinate_map_quartic():
    """sigma(s, p) and a(s, p) of the exact map, expanded to total degree four."""
    radial, angular, scale = sp.symbols('s p epsilon', real=True)
    sigma_exact = sp.sqrt((f_a + radial) ** 2 + angular ** 2) - f_a
    angular_exact = f_a * sp.atan(angular / (f_a + radial))
    expansions = []
    for exact in (sigma_exact, angular_exact):
        expansions.append(sp.expand(sp.series(exact.subs({radial: scale * radial, angular: scale * angular}),
                                              scale, 0, 5).removeO().subs(scale, 1)))
    return tuple(expansions)


def derive_cartesian_quartic_contacts():
    """Pure (p, Phi) quartic terms of the Cartesian Lagrangian at s = 0, minus the polar ones.

    Polar: a^4, a^3 <P Phi>, a^2 <S Phi^2> terms of the source Lagrangian (a -> p).  Cartesian:
    substitute the exact map (s = 0) into the source Lagrangian and V_phi at stationarity.
    Only Phi up to second order is needed (Phi^3 and Phi^4 multiply at most one power of a).
    Returns dict with keys 'p^4', 'p^3 Phi', 'p^2 Phi^2' (differences), and the polar a^4 term.
    """
    orders, traces = _trace_symbols()
    scalar = traces[0]
    sigma_series, angular_series = derive_coordinate_map_quartic()
    radial, angular, scale, count = sp.symbols('s p epsilon field_count', real=True)
    mass_phi_squared = B_0 * f_pi ** 2 * scalar[(1, 0)] / f_a ** 2 - quartic * f_a ** 2
    sigma_of_p = sigma_series.subs(radial, 0).subs(angular, count * angular)
    angle_of_p = angular_series.subs(radial, 0).subs(angular, count * angular)
    cartesian = (_source_lagrangian(sigma_of_p, angle_of_p, count * scale, orders, traces)
                 - mass_phi_squared / 2 * (f_a ** 2 + (count * angular) ** 2)
                 - quartic / 4 * (f_a ** 2 + (count * angular) ** 2) ** 2)
    polar = _source_lagrangian(0, count * angular, count * scale, orders, traces)
    difference = sp.expand(sp.series(cartesian - polar, count, 0, 5).removeO()).coeff(count, 4)
    results = {
        'p^4': sp.simplify(difference.coeff(scale, 0).coeff(angular, 4)),
        'p^3 Phi': sp.simplify(difference.coeff(scale, 1).coeff(angular, 3)),
        'p^2 Phi^2': sp.simplify(difference.coeff(scale, 2).coeff(angular, 2)),
        'polar a^4': sp.simplify(sp.expand(sp.series(polar, count, 0, 5).removeO()).coeff(count, 4).coeff(angular, 4)),
    }
    quadratic = sp.expand(sp.series(cartesian - polar, count, 0, 3).removeO()).coeff(count, 2)
    results['p^2 (must vanish)'] = sp.simplify(quadratic)
    return results


def derive_record_mapping():
    """Quadratic blocks of the polar Lagrangian against those of the finite-mixing record.

    Record (Alda:2026nrd conventions): M^2_{I a} = -f_pi <C lambda_I>/f_a,
    M^2_{aa} = m_{a,0}^2 - (f_pi^2/4) <d_a^2 chi_+|_0> + (2 f_pi^2 m_0^2/(N_f f_a^2)) c_GA^2 (N_f = 3),
    d_a^n chi_+|_0 = 4 B0 S_n/f_a^n with x-derivatives S_n = n! S_0n.
    Polar: L_2 contains (B0 f_pi/f_a) a <P_01 Phi> + (B0 f_pi^2/f_a^2) a^2 <S_02>.
    Returns (C solved from the mixing entry in units of <P_01 lambda_I>, M^2_aa difference with
    m_{a,0} = c_GA = 0).
    """
    tr_p01_lambda, tr_c_lambda, tr_s02 = sp.symbols('tr_P01_lambda tr_C_lambda tr_S02', real=True)
    bare_mass_squared, gluon_coupling, anomaly_mass_squared = sp.symbols('m_a0_2 c_GA m_0_2', real=True)
    polar_mixing = -B_0 * f_pi / f_a * tr_p01_lambda
    record_mixing = -f_pi * tr_c_lambda / f_a
    coupling = sp.solve(sp.Eq(polar_mixing, record_mixing), tr_c_lambda)[0] / tr_p01_lambda
    second_derivative = 4 * B_0 * sp.factorial(2) * tr_s02 / f_a ** 2
    record_alp_entry = (bare_mass_squared - f_pi ** 2 * second_derivative / 4
                        + 2 * f_pi ** 2 * anomaly_mass_squared * gluon_coupling ** 2 / (3 * f_a ** 2))
    polar_alp_entry = -2 * B_0 * f_pi ** 2 * tr_s02 / f_a ** 2
    return sp.simplify(coupling), sp.simplify(record_alp_entry.subs({bare_mass_squared: 0, gluon_coupling: 0})
                                              - polar_alp_entry)


# ---------------------------------------------------------------------------
# (2) Numerical construction from the exact Lagrangian
# ---------------------------------------------------------------------------

_DISJOINT_PAIRS = {}


def _disjoint_pairs(size):
    if size not in _DISJOINT_PAIRS:
        _DISJOINT_PAIRS[size] = [(left, right) for left in range(size) for right in range(size) if left & right == 0]
    return _DISJOINT_PAIRS[size]


class TruncatedExpansion:
    """Polynomial in formal variables t_1..t_n with t_l^2 = 0; scalar or matrix coefficients.

    The coefficient with all bits set of a function evaluated at x = sum_l t_l v_l is its n-th
    derivative contracted with v_1 ... v_n (n <= 4 here).
    """

    def __init__(self, coefficients):
        self.coefficients = coefficients  # shape (2^n, ...) indexed by bit mask

    @property
    def size(self):
        return self.coefficients.shape[0]

    def constant_like(self, value):
        value = np.asarray(value, dtype=complex)
        coefficients = np.zeros((self.size,) + value.shape, dtype=complex)
        coefficients[0] = value
        return TruncatedExpansion(coefficients)

    def _wrap(self, other):
        return other if isinstance(other, TruncatedExpansion) else self.constant_like(other)

    def __add__(self, other):
        return TruncatedExpansion(self.coefficients + self._wrap(other).coefficients)

    __radd__ = __add__

    def __neg__(self):
        return TruncatedExpansion(-self.coefficients)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def _product(self, other, multiply):
        first, second = self.coefficients, other.coefficients
        shape = np.shape(multiply(first[0], second[0]))
        result = np.zeros((self.size,) + shape, dtype=complex)
        for left, right in _disjoint_pairs(self.size):
            result[left | right] += multiply(first[left], second[right])
        return TruncatedExpansion(result)

    def __mul__(self, other):
        if not isinstance(other, TruncatedExpansion):
            return TruncatedExpansion(self.coefficients * other)
        return self._product(other, np.multiply)

    __rmul__ = __mul__

    def __truediv__(self, number):
        return TruncatedExpansion(self.coefficients / number)

    def matmul(self, other):
        return self._product(other, np.matmul)

    def times_matrix(self, matrix):
        return TruncatedExpansion(np.einsum('m,ij->mij', self.coefficients, matrix))

    def trace(self):
        return TruncatedExpansion(np.trace(self.coefficients, axis1=1, axis2=2))

    def dagger(self):
        return TruncatedExpansion(np.conj(np.swapaxes(self.coefficients, 1, 2)))

    def apply(self, derivatives):
        """f(x0 + delta) from f(x0), f'(x0), ..., f''''(x0) (scalar expansions)."""
        delta = self - self.coefficients[0]
        result = self.constant_like(derivatives[0])
        power = self.constant_like(1.0)
        for order in range(1, len(derivatives)):
            power = power * delta
            result = result + power * (derivatives[order] / math.factorial(order))
        return result


def derivative_product(first, second, mask_dot, matrices=False):
    """Expansion of d_mu F d^mu H for plane-wave legs: (i K_1).(i K_2) = -K_1.K_2 per mask pair."""
    shape = (first.coefficients[0] @ second.coefficients[0]).shape if matrices else ()
    result = np.zeros((first.size,) + shape, dtype=complex)
    for left, right in _disjoint_pairs(first.size):
        if left == 0 or right == 0:
            continue
        product = (first.coefficients[left] @ second.coefficients[right] if matrices
                   else first.coefficients[left] * second.coefficients[right])
        result[left | right] -= mask_dot[left, right] * product
    return TruncatedExpansion(result)


def matrix_exponential(generator):
    """exp(X) for a matrix expansion with vanishing constant part (X^5 = 0)."""
    result = generator.constant_like(np.eye(3))
    power = generator.constant_like(np.eye(3))
    for order in range(1, 5):
        power = power.matmul(generator)
        result = result + power / math.factorial(order)
    return result


def field_point(directions):
    """Field expansions x^A = sum_l t_l v_l^A for the given direction vectors."""
    size = 1 << len(directions)
    coefficients = np.zeros((size, NUMBER_OF_FIELDS))
    for position, direction in enumerate(directions):
        coefficients[1 << position] = direction
    return [TruncatedExpansion(coefficients[:, index].astype(complex)) for index in range(NUMBER_OF_FIELDS)]


def source_matrices(sources, radial_polar, angular_polar, decay_constant):
    """S(sigma, a) and P(sigma, a) as matrix expansions."""
    scalar = radial_polar.constant_like(np.zeros((3, 3)))
    pseudoscalar = radial_polar.constant_like(np.zeros((3, 3)))
    radial_ratio, angular_ratio = radial_polar / decay_constant, angular_polar / decay_constant
    for (r, s) in set(sources['S']) | set(sources['P']):
        monomial = radial_polar.constant_like(1.0)
        for _ in range(r):
            monomial = monomial * radial_ratio
        for _ in range(s):
            monomial = monomial * angular_ratio
        if (r, s) in sources['S']:
            scalar = scalar + monomial.times_matrix(sources['S'][(r, s)])
        if (r, s) in sources['P']:
            pseudoscalar = pseudoscalar + monomial.times_matrix(sources['P'][(r, s)])
    return scalar, pseudoscalar


def polar_fields_from_cartesian(radial, angular, decay_constant):
    """sigma(s, p) and a(s, p) of the exact map, exact to fourth order in the formal variables."""
    rho_squared = (radial + decay_constant) * (radial + decay_constant) + angular * angular
    x0 = rho_squared.coefficients[0].real
    square_root = [x0 ** 0.5, 0.5 * x0 ** -0.5, -0.25 * x0 ** -1.5, 0.375 * x0 ** -2.5, -0.9375 * x0 ** -3.5]
    radial_polar = rho_squared.apply(square_root) - decay_constant
    inverse = (radial + decay_constant).apply([(-1) ** n * math.factorial(n) / decay_constant ** (n + 1)
                                               for n in range(5)])
    angular_polar = (angular * inverse).apply([0.0, 1.0, 0.0, -2.0, 0.0]) * decay_constant
    return radial_polar, angular_polar, rho_squared


def lagrangian(fields, parameters, coordinates, mask_dot):
    """Complete Lagrangian (kinetic and potential) on field expansions with plane-wave momenta."""
    decay_constant, pion_decay_constant = parameters['f_a'], parameters['f_pi']
    radial, angular = fields[0], fields[1]
    if coordinates == 'polar':
        radial_polar, angular_polar = radial, angular
        rho_squared = (radial + decay_constant) * (radial + decay_constant)
        kinetic = (derivative_product(radial, radial, mask_dot) * 0.5
                   + (1 + radial / decay_constant) * (1 + radial / decay_constant)
                   * derivative_product(angular, angular, mask_dot) * 0.5)
    else:
        radial_polar, angular_polar, rho_squared = polar_fields_from_cartesian(radial, angular, decay_constant)
        kinetic = (derivative_product(radial, radial, mask_dot) * 0.5
                   + derivative_product(angular, angular, mask_dot) * 0.5)
    meson_matrix = radial.constant_like(np.zeros((3, 3)))
    for k in range(9):
        meson_matrix = meson_matrix + fields[2 + k].times_matrix(LAMBDA[k])
    chiral_field = matrix_exponential(meson_matrix * (1j / pion_decay_constant))
    chiral_field_dagger = matrix_exponential(meson_matrix * (-1j / pion_decay_constant))
    kinetic = kinetic + derivative_product(chiral_field, chiral_field_dagger, mask_dot, matrices=True).trace() * (
        pion_decay_constant ** 2 / 4)
    scalar, pseudoscalar = source_matrices(parameters['sources'], radial_polar, angular_polar, decay_constant)
    chi = (scalar + pseudoscalar * 1j) * (2 * parameters['B_0'])
    chi_plus_trace = (chi.matmul(chiral_field_dagger) + chi.dagger().matmul(chiral_field)).trace()
    singlet_trace = meson_matrix.trace()
    radial_potential = rho_squared * (parameters['m_phi_squared'] / 2) + rho_squared * rho_squared * (
        parameters['lambda'] / 4)
    return (kinetic + chi_plus_trace * (pion_decay_constant ** 2 / 4)
            - singlet_trace * singlet_trace * (parameters['m_1_squared'] / 12) - radial_potential)


def vertex(parameters, coordinates, directions, momentum_coefficients, gram):
    """n-point vertex (no factor i) on directions with momenta sum_e c_e k_e; gram = k_e.k_f."""
    size = 1 << len(directions)
    mask_momenta = np.zeros((size, gram.shape[0]))
    for mask in range(size):
        for leg in range(len(directions)):
            if mask >> leg & 1:
                mask_momenta[mask] += momentum_coefficients[leg]
    mask_dot = mask_momenta @ gram @ mask_momenta.T
    return lagrangian(field_point(directions), parameters, coordinates, mask_dot).coefficients[-1].real


def quadratic_forms(parameters, coordinates):
    """(gradient, kinetic matrix K, mass matrix M^2) at the vacuum, all eleven fields."""
    identity = np.eye(NUMBER_OF_FIELDS)
    zero_gram = np.zeros((2, 2))
    unit_gram = np.array([[0.0, -1.0], [-1.0, 0.0]])  # k_1.k_2 = -1: kinetic vertex = K_ab
    gradient = np.array([vertex(parameters, coordinates, [identity[a]], [np.ones(1)], np.zeros((1, 1)))
                         for a in range(NUMBER_OF_FIELDS)])
    mass_squared = np.zeros((NUMBER_OF_FIELDS, NUMBER_OF_FIELDS))
    kinetic = np.zeros((NUMBER_OF_FIELDS, NUMBER_OF_FIELDS))
    for a in range(NUMBER_OF_FIELDS):
        for b in range(a, NUMBER_OF_FIELDS):
            directions = [identity[a], identity[b]]
            potential = vertex(parameters, coordinates, directions, np.eye(2), zero_gram)
            mass_squared[a, b] = mass_squared[b, a] = -potential
            kinetic[a, b] = kinetic[b, a] = vertex(parameters, coordinates, directions, np.eye(2), unit_gram) - potential
    return gradient, kinetic, mass_squared


def pole_vectors(mass_squared):
    """Pole masses and unit vectors of the angular-meson block (K = 1), angular component >= 0."""
    masses_squared, vectors = np.linalg.eigh(mass_squared[1:, 1:])
    signs = np.where(np.abs(vectors[0]) > ANGULAR_COMPONENT_THRESHOLD, np.sign(vectors[0]), 1.0)
    embedded = np.zeros((NUMBER_OF_FIELDS, len(masses_squared)))
    embedded[1:] = vectors * signs
    return masses_squared, embedded


def decay_gram(parent_mass_squared, daughter_masses_squared, invariants):
    """Gram matrix of all-outgoing momenta (k_0, k_1, k_2, k_3) = (-p, p_1, p_2, p_3).

    invariants maps frozenset({u, v}) to s_uv = (p_u + p_v)^2 for daughters 1, 2, 3.
    """
    gram = np.zeros((4, 4))
    gram[0, 0] = parent_mass_squared
    for r in (1, 2, 3):
        gram[r, r] = daughter_masses_squared[r - 1]
        others = frozenset({1, 2, 3} - {r})
        gram[0, r] = gram[r, 0] = -(parent_mass_squared + daughter_masses_squared[r - 1] - invariants[others]) / 2
    for u, v in itertools.combinations((1, 2, 3), 2):
        gram[u, v] = gram[v, u] = (invariants[frozenset({u, v})] - daughter_masses_squared[u - 1]
                                   - daughter_masses_squared[v - 1]) / 2
    return gram


def dalitz_invariants(parent_mass_squared, daughter_masses_squared, s12, s13):
    """Invariants with s_23 fixed by s_12 + s_13 + s_23 = m_a^2 + m_1^2 + m_2^2 + m_3^2."""
    s23 = parent_mass_squared + sum(daughter_masses_squared) - s12 - s13
    return {frozenset({1, 2}): s12, frozenset({1, 3}): s13, frozenset({2, 3}): s23}


def direct_amplitude(parameters, coordinates, mass_squared_matrix, kinetic_matrix, masses_squared, vectors,
                     legs, invariants):
    """Tree amplitude from the exact Lagrangian: contact plus the three exchanges.

    legs = (parent, daughter 1, daughter 2, daughter 3) pole indices.  The exchange uses the full
    eleven-field propagator (M^2 - s K)^{-1}.  Returns dict with 'total', 'contact', 'exchange'
    (all fields), 'radial_exchange' (internal radial field only) and 'radial_offdiagonal'.
    """
    parent_mass_squared = masses_squared[legs[0]]
    daughter_masses_squared = [masses_squared[n] for n in legs[1:]]
    gram = decay_gram(parent_mass_squared, daughter_masses_squared, invariants)
    directions = [vectors[:, n] for n in legs]
    unit = np.eye(4)
    contact = vertex(parameters, coordinates, directions, unit, gram)
    identity = np.eye(NUMBER_OF_FIELDS)
    exchange, radial_exchange, offdiagonal = 0.0, 0.0, 0.0
    for attached, (first, second) in PAIRINGS:
        momentum_squared = invariants[frozenset({first, second})]
        internal = -(unit[0] + unit[attached])
        left = np.array([vertex(parameters, coordinates, [directions[0], directions[attached], identity[n]],
                                [unit[0], unit[attached], internal], gram) for n in range(NUMBER_OF_FIELDS)])
        right = np.array([vertex(parameters, coordinates, [identity[n], directions[first], directions[second]],
                                 [-internal, unit[first], unit[second]], gram) for n in range(NUMBER_OF_FIELDS)])
        propagator = np.linalg.inv(mass_squared_matrix - momentum_squared * kinetic_matrix)
        exchange += left @ propagator @ right
        radial_exchange += left[0] * propagator[0, 0] * right[0]
        offdiagonal = max(offdiagonal, np.max(np.abs(propagator[0, 1:])))
    return {'total': contact + exchange, 'contact': contact, 'exchange': exchange,
            'radial_exchange': radial_exchange, 'radial_offdiagonal': offdiagonal}


# ---------------------------------------------------------------------------
# Analytic result: record amplitude without sigma + trace-form sigma exchange
# ---------------------------------------------------------------------------

def has_angular_component(vectors, index):
    return abs(vectors[1, index]) > ANGULAR_COMPONENT_THRESHOLD


def meson_matrix_of(vector):
    return np.einsum('k,kij->ij', vector[2:], LAMBDA)


def leg_column(parameters, mass_squared, vector):
    """Original polar column (Lambda_x, xi_x) of a pole: trace form if it has an ALP component.

    With ALP component: Lambda = (2 B0 f_pi/f_a) M_{m^2}^{-1}(P_01)/sqrt(N(m^2)), xi = 1/sqrt(N(m^2)).
    Without (exact meson pole): Lambda = its meson matrix, xi = 0.
    """
    if abs(vector[1]) <= ANGULAR_COMPONENT_THRESHOLD:
        return meson_matrix_of(vector), 0.0
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    root = np.sqrt(pole_normalization(parameters, mass_squared))
    inverse_image = inverse_mass_operator(parameters, mass_squared, parameters['sources']['P'][(0, 1)])
    return 2 * B0 * fpi / fa * inverse_image / root, 1 / root


def radial_vertex_trace_form(parameters, first, second, momentum_squared):
    """g_xy(q^2): the sigma -> PP coupling of the radial-decay record with m_sigma^2 -> q^2.

    first, second = (column, pole mass squared) with column = (Lambda, xi) from leg_column.
    f_a g = (B0/2)<S_10{Lambda_x,Lambda_y}> - (B0 f_pi/f_a)<P_11(xi_x Lambda_y + xi_y Lambda_x)>
            - (2 B0 f_pi^2/f_a^2) xi_x xi_y <S_12> + xi_x xi_y (q^2 - m_x^2 - m_y^2).
    """
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    S, P = parameters['sources']['S'], parameters['sources']['P']
    ((lambda_x, xi_x), mass_x), ((lambda_y, xi_y), mass_y) = first, second
    braces = (B0 / 2 * trace(S[(1, 0)] @ (lambda_x @ lambda_y + lambda_y @ lambda_x))
              - B0 * fpi / fa * trace(P[(1, 1)] @ (xi_x * lambda_y + xi_y * lambda_x))
              - 2 * B0 * fpi ** 2 / fa ** 2 * xi_x * xi_y * trace(S[(1, 2)])
              + xi_x * xi_y * (momentum_squared - mass_x - mass_y))
    return braces / fa


def radial_exchange_trace_form(parameters, masses_squared, vectors, legs, invariants):
    """sum over (k|ij) of g_ak(s_ij) g_ij(s_ij)/(m_sigma^2 - s_ij), trace-form columns."""
    radial_mass_squared = radial_mass_squared_trace_form(parameters)
    columns = [(leg_column(parameters, masses_squared[n], vectors[:, n]), masses_squared[n]) for n in legs]
    total = 0.0
    for attached, (first, second) in PAIRINGS:
        momentum_squared = invariants[frozenset({first, second})]
        total += (radial_vertex_trace_form(parameters, columns[0], columns[attached], momentum_squared)
                  * radial_vertex_trace_form(parameters, columns[first], columns[second], momentum_squared)
                  / (radial_mass_squared - momentum_squared))
    return total


def record_quadratic_blocks(parameters):
    """K = 1 and M^2 in the field order (lambda_0..lambda_8, a) from the trace forms.

    M^2_kl = (1/2)<lambda_k M_0(lambda_l)>, M^2_{k a} = -f_pi <C lambda_k>/f_a with C = B0 P_01,
    M^2_aa = -2 B0 f_pi^2 <S_02>/f_a^2 (mapping onto the finite-mixing record).
    """
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    coupling = B0 * parameters['sources']['P'][(0, 1)]
    mass_squared = np.zeros((10, 10))
    mass_squared[:9, :9] = mass_operator_matrix(parameters, 0.0)
    for k in range(9):
        mass_squared[k, 9] = mass_squared[9, k] = -fpi * trace(coupling @ LAMBDA[k]) / fa
    mass_squared[9, 9] = -2 * B0 * fpi ** 2 * trace(parameters['sources']['S'][(0, 2)]) / fa ** 2
    return np.eye(10), mass_squared


def _sympy_matrix(matrix):
    return sp.Matrix(3, 3, lambda i, j: sp.Float(matrix[i, j].real, 17) + sp.I * sp.Float(matrix[i, j].imag, 17))


def record_amplitude_without_radial(parameters, masses_squared, vectors, legs, invariants):
    """Boxed step-4 amplitude of the finite-mixing three-pseudoscalar record.

    Mapping (verified by the tests): widehat M_n = n! (S_0n + i P_0n) (x-derivatives), m_q -> S_00,
    C = B0 P_01, c_GA = 0, Theta_n = 0, m_0^2/N_f -> m_1^2/3, m_{a,0} = 0.  The pole-aligned
    coefficients S_n', P_n' follow from the record's Taylor law with c_A_hat = -M_{m_a^2}^{-1}(C).
    Uses compact_amplitude of scripts/derive_alp_finite_mixing_three_pseudoscalars.py.
    """
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    S, P = parameters['sources']['S'], parameters['sources']['P']
    zero = np.zeros((3, 3))
    parent_mass_squared = masses_squared[legs[0]]
    coupling = B0 * P[(0, 1)]
    hat_axial = -inverse_mass_operator(parameters, parent_mass_squared, coupling)
    parent_norm = pole_normalization(parameters, parent_mass_squared)
    dressed_scalar = [_sympy_matrix(math.factorial(n) * S.get((0, n), zero)) for n in range(5)]
    dressed_pseudoscalar = [_sympy_matrix(math.factorial(n) * P.get((0, n), zero)) for n in range(5)]
    aligned_scalar, aligned_pseudoscalar = aligned_source_coefficients(_sympy_matrix(hat_axial), dressed_scalar,
                                                                       dressed_pseudoscalar)
    daughters = {}
    for label, n in zip('ijk', legs[1:]):
        if has_angular_component(vectors, n):
            root = np.sqrt(pole_normalization(parameters, masses_squared[n]))
            column = (2 * fpi / fa * (parent_mass_squared - masses_squared[n])
                      * inverse_mass_operator(parameters, masses_squared[n], hat_axial) / root)
            daughters[label] = (_sympy_matrix(column), sp.Float(1 / root, 17))
        else:
            daughters[label] = (_sympy_matrix(meson_matrix_of(vectors[:, n])), 0)
    parent = (sp.zeros(3), sp.Float(1 / np.sqrt(parent_norm), 17))
    pair_of = {'i': frozenset({2, 3}), 'j': frozenset({1, 3}), 'k': frozenset({1, 2})}
    parent_dot = {label: sp.Float((parent_mass_squared + masses_squared[n] - invariants[pair_of[label]]) / 2, 17)
                  for label, n in zip('ijk', legs[1:])}
    record_invariants = {'jk': invariants[frozenset({2, 3})], 'ik': invariants[frozenset({1, 3})],
                         'ij': invariants[frozenset({1, 2})]}
    kinetic, mass_squared = record_quadratic_blocks(parameters)
    conversion = np.eye(10)
    shift = 2 * fpi / fa * (-hat_axial)
    conversion[:9, 9] = [0.5 * trace(LAMBDA[k] @ shift) for k in range(9)]
    aligned_kinetic = sp.Matrix(conversion.T @ kinetic @ conversion)
    aligned_mass = sp.Matrix(conversion.T @ mass_squared @ conversion)
    basis = [(_sympy_matrix(LAMBDA[k]), 0) for k in range(9)] + [(sp.zeros(3), 1)]
    expression = compact_amplitude(_sympy_matrix(hat_axial), aligned_pseudoscalar[1], sp.I * aligned_scalar[1],
                                   aligned_scalar, aligned_pseudoscalar, parent, daughters, parent_dot,
                                   record_invariants, basis, aligned_kinetic, aligned_mass)
    value = complex(sp.N(expression.subs({f_pi: fpi, f_a: fa, B_0: B0})))
    return value


def analytic_amplitude(parameters, masses_squared, vectors, legs, invariants):
    """(record amplitude without sigma, trace-form sigma exchange)."""
    return (record_amplitude_without_radial(parameters, masses_squared, vectors, legs, invariants).real,
            radial_exchange_trace_form(parameters, masses_squared, vectors, legs, invariants))


# ---------------------------------------------------------------------------
# Parameter sets
# ---------------------------------------------------------------------------

def random_parameters(seed=2026):
    """Random non-commuting Hermitian sources (GeV), r + s <= 4, obeying assumptions 1 and 2."""
    generator = np.random.default_rng(seed)
    unitary, _ = np.linalg.qr(generator.normal(size=(3, 3)) + 1j * generator.normal(size=(3, 3)))
    sources = {'S': {(0, 0): unitary @ np.diag([0.004, 0.008, 0.1]) @ unitary.conj().T}, 'P': {}}
    for (r, s) in itertools.product(range(5), range(5)):
        if r + s > MAXIMUM_SOURCE_ORDER or (r, s) == (0, 0):
            continue
        sources['S'][(r, s)] = random_hermitian(generator, 0.02)
        sources['P'][(r, s)] = random_hermitian(generator, 0.02)
    for key in [(0, 1), (1, 1)]:
        sources['S'][key] = sources['S'][key] - trace(sources['S'][key]) / 3 * np.eye(3)
    sources['P'][(1, 0)] = np.zeros((3, 3))
    sources['P'][(0, 0)] = np.zeros((3, 3))
    sources['S'][(0, 2)] = sources['S'][(0, 2)] - 0.3 * np.eye(3)  # positive angular mass
    parameters = {'B_0': 2.6, 'f_pi': 0.0924, 'f_a': 0.4, 'm_1_squared': 0.75, 'lambda': 0.8, 'sources': sources}
    return impose_stationarity(parameters)


def only_up_parameters(up_mass=0.0022, down_mass=0.0047, strange_mass=0.095, decay_constant=0.4, quartic_value=0.8):
    """M = diag[m_u (1 + sigma/f_a) e^{-i a/f_a}, m_d, m_s], coefficients with r + s <= 4."""
    projector = np.diag([1.0, 0.0, 0.0])
    sources = {'S': {(0, 0): np.diag([up_mass, down_mass, strange_mass])}, 'P': {}}
    for (r, s) in itertools.product(range(2), range(5)):
        if r + s > MAXIMUM_SOURCE_ORDER:
            continue
        if (r, s) != (0, 0):
            sources['S'][(r, s)] = up_mass * np.cos(np.pi * s / 2) / math.factorial(s) * projector
        sources['P'][(r, s)] = -up_mass * np.sin(np.pi * s / 2) / math.factorial(s) * projector
    parameters = {'B_0': 2.6, 'f_pi': 0.0924, 'f_a': decay_constant, 'm_1_squared': 0.75,
                  'lambda': quartic_value, 'sources': sources}
    return impose_stationarity(parameters)


def with_radial_mass(parameters, radial_mass_squared):
    """Copy with lambda chosen so that m_sigma^2 takes the given value (stationarity re-imposed)."""
    result = dict(parameters)
    S = parameters['sources']['S']
    shift = parameters['B_0'] * parameters['f_pi'] ** 2 / parameters['f_a'] ** 2 * (
        trace(S[(1, 0)]) - 2 * trace(S.get((2, 0), np.zeros((3, 3)))))
    result['lambda'] = (radial_mass_squared - shift) / (2 * parameters['f_a'] ** 2)
    return impose_stationarity(result)


# ---------------------------------------------------------------------------
# Comparison driver
# ---------------------------------------------------------------------------

def spectrum(parameters):
    """Polar quadratic forms, pole masses and vectors, ALP-like parent index."""
    gradient, kinetic, mass_squared = quadratic_forms(parameters, 'polar')
    masses_squared, vectors = pole_vectors(mass_squared)
    parent = int(np.argmax(np.abs(vectors[1])))
    return {'gradient': gradient, 'kinetic': kinetic, 'mass_squared': mass_squared,
            'masses_squared': masses_squared, 'vectors': vectors, 'parent': parent,
            'mixed': [n for n in range(len(masses_squared)) if has_angular_component(vectors, n)],
            'unmixed': [n for n in range(len(masses_squared)) if not has_angular_component(vectors, n)]}


def compare_point(parameters, data, daughters, s12, s13, include_cartesian=True):
    """Direct polar (and Cartesian) amplitude and the analytic result at one Dalitz point."""
    legs = (data['parent'],) + tuple(daughters)
    masses_squared = data['masses_squared']
    invariants = dalitz_invariants(masses_squared[legs[0]], [masses_squared[n] for n in legs[1:]], s12, s13)
    polar = direct_amplitude(parameters, 'polar', data['mass_squared'], data['kinetic'], masses_squared,
                             data['vectors'], legs, invariants)
    without_radial, radial = analytic_amplitude(parameters, masses_squared, data['vectors'], legs, invariants)
    result = {'polar': polar, 'record_without_radial': without_radial, 'radial_trace_form': radial,
              'analytic': without_radial + radial, 'legs': legs, 'invariants': invariants}
    if include_cartesian:
        _, cartesian_kinetic, cartesian_mass = quadratic_forms(parameters, 'cartesian')
        result['cartesian'] = direct_amplitude(parameters, 'cartesian', cartesian_mass, cartesian_kinetic,
                                               masses_squared, data['vectors'], legs, invariants)
        result['cartesian_hessian_difference'] = np.max(np.abs(cartesian_mass - data['mass_squared']))
    return result


def main():
    print('Units: GeV. Tree level, LO U(3), polar coordinates, no expansion in f_pi/f_a.\n')
    first_order, second_order = derive_symmetric_map_second_order()
    print('log(e^{tA} e^{tB} e^{tA}): first order - (2A + B) =', first_order, '; second order =', second_order)
    print('polar quadratic radial mixing: sigma a ->', derive_polar_radial_mixing()[0],
          '; sigma <P_rs Phi> ->', derive_polar_radial_mixing()[1])
    f_a_g, conjecture_difference, on_shell_difference, parent_pair = derive_off_shell_radial_vertex()
    print('f_a g_xy(q^2) =', f_a_g)
    print('  minus conjectured form:', conjecture_difference, '; minus on-shell g (m_sigma^2 -> q^2):',
          on_shell_difference, '; parent pair 2k_a.k_k - (s_ij - m_a^2 - m_k^2):', parent_pair)
    print('Cartesian minus polar quartic contacts at s = 0:')
    for key, value in derive_cartesian_quartic_contacts().items():
        print(f'  {key}: {value}')
    coupling, alp_entry = derive_record_mapping()
    print('mapping: C = (', coupling, ') P_01 ; record M^2_aa - polar M^2_aa =', alp_entry, '\n')

    parameters = random_parameters(2026)
    data = spectrum(parameters)
    mixed = [n for n in data['mixed'] if n != data['parent']]
    print('random sources (seed 2026): parent', data['parent'], 'mixed', data['mixed'], 'unmixed', data['unmixed'])
    print('  m_sigma^2 =', radial_mass_squared_trace_form(parameters))
    for daughters in ((mixed[0], mixed[1], mixed[2]), (mixed[0], mixed[0], data['unmixed'][0]),
                      (data['unmixed'][0], data['unmixed'][1], mixed[1])):
        point = compare_point(parameters, data, daughters, 0.9, 1.4)
        scale = abs(point['polar']['total'])
        print(f'  daughters {daughters}: polar {point["polar"]["total"]:.6e}, Cartesian '
              f'{point["cartesian"]["total"]:.6e}, analytic {point["analytic"]:.6e}; '
              f'relative |polar - analytic| {abs(point["polar"]["total"] - point["analytic"]) / scale:.1e}, '
              f'|polar - Cartesian| {abs(point["polar"]["total"] - point["cartesian"]["total"]) / scale:.1e}, '
              f'sigma part: direct {point["polar"]["radial_exchange"]:.4e} trace {point["radial_trace_form"]:.4e}')


if __name__ == '__main__':
    main()
