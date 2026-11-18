#!/usr/bin/env python3
"""On-shell equality of polar and Cartesian sigma -> P_i P_j vertices for general sources.

Conventions (those of the draft, section "Radial-mode decays", and Alda:2026nrd):

* Partonic source M(a, sigma) = S(a, sigma) + i P(a, sigma), with Hermitian S, P and
  S(a, sigma) = sum_{r,s} S_rs (sigma/f_a)^r (a/f_a)^s (same for P).  The coefficients
  S_rs, P_rs are general Hermitian 3x3 matrices of mass dimension one.
* Phi = sum_k P_k lambda_k, <lambda_k lambda_l> = 2 delta_kl, lambda_0 = sqrt(2/3) 1,
  u = exp(i Phi / (2 f_pi)), U = u^2, chi = 2 B0 M, meson Lagrangian
  (f_pi^2/4) <u_mu u^mu + chi_+> - (m_1^2/12) <Phi>^2.
* Polar fields sigma, a with kinetic term (1/2)(d sigma)^2 + (1/2)(1 + sigma/f_a)^2 (d a)^2.
  Cartesian fields s, p with f_a + s + i p = (f_a + sigma) exp(i a/f_a) (exact map) and
  kinetic term (1/2)(d s)^2 + (1/2)(d p)^2.
* Radial potential V_phi = (m_phi^2/2) rho^2 + (lambda/4) rho^4, rho = f_a + sigma.
* M_s(X) = B0 {S_00, X} + (m_1^2/3) <X> 1 - s X, acting on Hermitian 3x3 matrices.
* Coupling: L contains -(1/2) sigma sum_{ij} g_ij P_i P_j, outgoing on-shell momenta,
  angular component of every physical pole direction chosen positive.

Units: natural units, all masses, fields, B0, f_pi, f_a, S_rs, P_rs and g_ij in GeV.

Assumptions (all required for the trace form, see the memory record):
  P_00 = 0 and <S_01> = 0 (no tadpole in Phi or a); P_10 = 0 and <S_11> = 0 (no radial mixing
  with a or Phi at quadratic order); the anomalous singlet mass term does not depend on sigma or a;
  stationarity in sigma fixes m_phi^2; regular simple poles (M_{m_r^2} invertible).

The script contains
  (1) a sympy derivation of the coordinate map, of the general Cartesian cubic Lagrangian in
      terms of traces, and of the reduction of g^polar - g^Cart to the pole and radial-mass
      equations;
  (2) an independent numerical construction, with random non-commuting Hermitian sources, of
      both Lagrangians from the exact map, followed by the on-shell projection of the cubic
      terms on the numerical pole vectors.

Run:
    .venv/bin/python scripts/derive_radial_decay_polar_cartesian_general.py
"""
import itertools

import numpy as np
import sympy as sp

MAXIMUM_SOURCE_ORDER = 3  # r + s <= 3 is complete for cubic vertices

# ---------------------------------------------------------------------------
# (1) Symbolic derivation
# ---------------------------------------------------------------------------

f_a, f_pi, B_0, quartic = sp.symbols('f_a f_pi B_0 lambda', positive=True)
mass_phi_squared = sp.symbols('m_phi_squared', real=True)


def derive_coordinate_map():
    """Expand sigma(s, p), a(s, p) of the exact map to cubic order and check the kinetic metric.

    Returns (sigma_series, angular_series, metric_difference) where metric_difference is the
    simplified difference between the pull-back of the flat Cartesian metric and the polar metric.
    """
    radial, angular, scale = sp.symbols('s p epsilon', real=True)
    sigma_exact = sp.sqrt((f_a + radial) ** 2 + angular ** 2) - f_a
    angular_exact = f_a * sp.atan(angular / (f_a + radial))
    sigma_series = sp.expand(sp.series(sigma_exact.subs({radial: scale * radial, angular: scale * angular}),
                                       scale, 0, 4).removeO().subs(scale, 1))
    angular_series = sp.expand(sp.series(angular_exact.subs({radial: scale * radial, angular: scale * angular}),
                                         scale, 0, 4).removeO().subs(scale, 1))
    sigma, phase = sp.symbols('sigma a', real=True)
    cartesian_of_polar = sp.Matrix([(f_a + sigma) * sp.cos(phase / f_a) - f_a, (f_a + sigma) * sp.sin(phase / f_a)])
    jacobian = cartesian_of_polar.jacobian([sigma, phase])
    metric_difference = sp.simplify(jacobian.T * jacobian - sp.diag(1, (1 + sigma / f_a) ** 2))
    return sigma_series, angular_series, metric_difference


def _trace_symbols():
    """Scalar symbols standing for the traces <S_rs>, <P_rs Phi>, <S_rs Phi^2>."""
    orders = [(r, s) for r in range(MAXIMUM_SOURCE_ORDER + 1) for s in range(MAXIMUM_SOURCE_ORDER + 1)
              if r + s <= MAXIMUM_SOURCE_ORDER]
    trace_scalar = {o: sp.Symbol(f'tr_S_{o[0]}{o[1]}', real=True) for o in orders}
    trace_pseudoscalar_meson = {o: sp.Symbol(f'tr_P_{o[0]}{o[1]}_Phi', real=True) for o in orders}
    trace_scalar_meson_pair = {o: sp.Symbol(f'tr_S_{o[0]}{o[1]}_PhiPhi', real=True) for o in orders}
    return orders, trace_scalar, trace_pseudoscalar_meson, trace_scalar_meson_pair


def derive_general_cubic_lagrangians():
    """Derive the radial-pseudoscalar-pseudoscalar terms of both Lagrangians in trace form.

    The source part is B0 f_pi^2 <S> + B0 f_pi <P Phi> - (B0/2) <S Phi^2>, from
    (f_pi^2/4)<chi_+> with U expanded to second order in Phi.  A bookkeeping parameter counts
    the meson fields.  Returns a dict with the stationarity value of m_phi^2, m_sigma^2 and the
    cubic coefficients in both coordinates.
    """
    orders, trace_scalar, trace_pseudoscalar_meson, trace_scalar_meson_pair = _trace_symbols()
    radial, angular, scale = sp.symbols('s p epsilon', real=True)
    sigma, phase = sp.symbols('sigma a', real=True)

    def source_lagrangian(radial_field, angular_field):
        total = 0
        for (r, s) in orders:
            monomial = (radial_field / f_a) ** r * (angular_field / f_a) ** s
            total += (B_0 * f_pi ** 2 * trace_scalar[(r, s)]
                      + B_0 * f_pi * scale * trace_pseudoscalar_meson[(r, s)]
                      - B_0 / 2 * scale ** 2 * trace_scalar_meson_pair[(r, s)]) * monomial
        return total

    # Polar Lagrangian (potential part); the kinetic sigma (d a)^2 / f_a term is added separately.
    polar_potential_lagrangian = (source_lagrangian(sigma, phase)
                                  - mass_phi_squared / 2 * (f_a + sigma) ** 2 - quartic / 4 * (f_a + sigma) ** 4)
    vacuum_condition = sp.diff(polar_potential_lagrangian, sigma).subs({sigma: 0, phase: 0, scale: 0})
    mass_phi_squared_vacuum = sp.solve(vacuum_condition, mass_phi_squared)[0]
    radial_mass_squared = sp.expand(-sp.diff(polar_potential_lagrangian, sigma, 2).subs(
        {sigma: 0, phase: 0, scale: 0, mass_phi_squared: mass_phi_squared_vacuum}))

    # Cartesian Lagrangian: substitute the exact map, expanded to cubic order.
    sigma_series, angular_series, _ = derive_coordinate_map()
    count = sp.Symbol('field_count')
    substitution = {sigma: sigma_series.subs({radial: count * radial, angular: count * angular}),
                    phase: angular_series.subs({radial: count * radial, angular: count * angular})}
    cartesian_source = source_lagrangian(*[substitution[x] for x in (sigma, phase)])
    rho_squared = (f_a + count * radial) ** 2 + (count * angular) ** 2
    cartesian_lagrangian = cartesian_source - mass_phi_squared / 2 * rho_squared - quartic / 4 * rho_squared ** 2
    cartesian_lagrangian = sp.expand(cartesian_lagrangian.subs(mass_phi_squared, mass_phi_squared_vacuum))
    # Keep total field order three: the radial field, angular field and meson insertions (scale).
    cubic = sp.expand(sp.series(cartesian_lagrangian.subs(scale, count * scale), count, 0, 4).removeO()).coeff(count, 3)
    cartesian_s_p_p = sp.simplify(cubic.coeff(scale, 0).coeff(radial, 1).coeff(angular, 2))
    cartesian_s_p_meson = {o: sp.simplify(cubic.coeff(scale, 1).coeff(trace_pseudoscalar_meson[o]).coeff(radial, 1).coeff(angular, 1))
                           for o in orders}
    cartesian_s_meson_pair = {o: sp.simplify(cubic.coeff(scale, 2).coeff(trace_scalar_meson_pair[o]).coeff(radial, 1).coeff(angular, 0))
                              for o in orders}
    quadratic = sp.expand(sp.series(cartesian_lagrangian.subs(scale, count * scale), count, 0, 3).removeO()).coeff(count, 2)
    cartesian_radial_mass_squared = sp.simplify(-2 * quadratic.coeff(scale, 0).coeff(radial, 2))
    cartesian_s_p_mixing = sp.simplify(quadratic.coeff(scale, 0).coeff(radial, 1).coeff(angular, 1))
    cartesian_s_meson_mixing = {o: sp.simplify(quadratic.coeff(scale, 1).coeff(trace_pseudoscalar_meson[o]).coeff(radial, 1))
                                for o in orders}

    polar_cubic = sp.expand(sp.series(polar_potential_lagrangian.subs(
        {sigma: count * sigma, phase: count * phase, scale: count * scale}), count, 0, 4).removeO()).coeff(count, 3)
    polar_sigma_a_a = sp.simplify(polar_cubic.coeff(scale, 0).coeff(sigma, 1).coeff(phase, 2))
    polar_sigma_a_meson = {o: sp.simplify(polar_cubic.coeff(scale, 1).coeff(trace_pseudoscalar_meson[o]).coeff(sigma, 1).coeff(phase, 1))
                           for o in orders}
    polar_sigma_meson_pair = {o: sp.simplify(polar_cubic.coeff(scale, 2).coeff(trace_scalar_meson_pair[o]).coeff(sigma, 1).coeff(phase, 0))
                              for o in orders}
    return {
        'trace_scalar': trace_scalar,
        'mass_phi_squared': mass_phi_squared_vacuum,
        'radial_mass_squared': radial_mass_squared,
        'cartesian_radial_mass_squared': cartesian_radial_mass_squared,
        'cartesian_s_p_mixing': cartesian_s_p_mixing,
        'cartesian_s_meson_mixing': cartesian_s_meson_mixing,
        'cartesian_s_p_p': cartesian_s_p_p,
        'cartesian_s_p_meson': cartesian_s_p_meson,
        'cartesian_s_meson_pair': cartesian_s_meson_pair,
        'polar_sigma_a_a': polar_sigma_a_a,
        'polar_sigma_a_meson': polar_sigma_a_meson,
        'polar_sigma_meson_pair': polar_sigma_meson_pair,
        'polar_kinetic_sigma_da_da': 1 / f_a,
    }


def derive_vertex_difference():
    """Reduce (g^polar - g^Cart) sqrt(N_i N_j) to the pole and radial-mass relations.

    Symbols: Omega_r = (2 B0 f_pi/f_a) M_{m_r^2}^{-1}(P_01) enters only through the traces
    <S_10 {Omega_i, Omega_j}>, <P_11 Omega_r>, <P_01 Omega_r>.  Returns
    (polar_braces, cartesian_braces, difference, difference_after_relations).
    """
    data = derive_general_cubic_lagrangians()
    trace_scalar = data['trace_scalar']
    anticommutator_trace = sp.Symbol('tr_S10_anticomm_Omega_i_Omega_j', real=True)
    p11_omega = {r: sp.Symbol(f'tr_P11_Omega_{r}', real=True) for r in 'ij'}
    p01_omega = {r: sp.Symbol(f'tr_P01_Omega_{r}', real=True) for r in 'ij'}
    pole_mass_squared = {r: sp.Symbol(f'm_{r}_squared', real=True) for r in 'ij'}
    radial_mass_squared = sp.Symbol('m_sigma_squared', real=True)
    two_p_i_dot_p_j = radial_mass_squared - pole_mass_squared['i'] - pole_mass_squared['j']

    # Coefficient of sigma P_i P_j with Phi = Omega_i P_i + Omega_j P_j and a = P_i + P_j
    # (normalization 1/sqrt(N_i N_j) removed). g = minus that coefficient.
    def braces(meson_pair, angle_meson_01, angle_meson_11, angle_angle, kinetic):
        coefficient = (meson_pair[(1, 0)] * anticommutator_trace
                       + angle_meson_01 * (p01_omega['i'] + p01_omega['j'])
                       + angle_meson_11 * (p11_omega['i'] + p11_omega['j'])
                       + 2 * angle_angle
                       - kinetic * two_p_i_dot_p_j)  # (d P_i).(d P_j) -> -p_i.p_j, both orderings
        return sp.expand(-coefficient)

    polar_braces = braces(data['polar_sigma_meson_pair'], data['polar_sigma_a_meson'][(0, 1)],
                          data['polar_sigma_a_meson'][(1, 1)], data['polar_sigma_a_a'],
                          data['polar_kinetic_sigma_da_da'])
    cartesian_braces = braces(data['cartesian_s_meson_pair'], data['cartesian_s_p_meson'][(0, 1)],
                              data['cartesian_s_p_meson'][(1, 1)], data['cartesian_s_p_p'], 0)
    difference = sp.expand(polar_braces - cartesian_braces)
    # Angular row of the pole equation and radial mass from the stationary potential.
    angular_row = {pole_mass_squared[r]: -2 * B_0 * f_pi ** 2 / f_a ** 2 * trace_scalar[(0, 2)]
                   - B_0 * f_pi / f_a * p01_omega[r] for r in 'ij'}
    angular_row[radial_mass_squared] = data['radial_mass_squared']
    return polar_braces, cartesian_braces, difference, sp.simplify(difference.subs(angular_row))


# ---------------------------------------------------------------------------
# (2) Independent numerical construction
# ---------------------------------------------------------------------------

def gell_mann_basis():
    """lambda_0 = sqrt(2/3) 1 and the eight Gell-Mann matrices, <lambda_k lambda_l> = 2 delta_kl."""
    basis = np.zeros((9, 3, 3), dtype=complex)
    basis[0] = np.sqrt(2 / 3) * np.eye(3)
    basis[1][0, 1] = basis[1][1, 0] = 1
    basis[2][0, 1], basis[2][1, 0] = -1j, 1j
    basis[3][0, 0], basis[3][1, 1] = 1, -1
    basis[4][0, 2] = basis[4][2, 0] = 1
    basis[5][0, 2], basis[5][2, 0] = -1j, 1j
    basis[6][1, 2] = basis[6][2, 1] = 1
    basis[7][1, 2], basis[7][2, 1] = -1j, 1j
    basis[8] = np.diag([1, 1, -2]) / np.sqrt(3)
    return basis


LAMBDA = gell_mann_basis()
NUMBER_OF_FIELDS = 11  # radial, angular, P_0 ... P_8
ANGULAR_COMPONENT_THRESHOLD = 1e-8  # poles below this angular component have no ALP admixture
SUBSETS = range(8)  # bit masks of the three formal expansion variables t_1, t_2, t_3


class MultilinearExpansion:
    """Polynomial in t_1, t_2, t_3 with t_k^2 = 0 and scalar or matrix coefficients.

    The coefficient with bit mask 0b111 of a function evaluated at
    x = t_1 v_1 + t_2 v_2 + t_3 v_3 is the third derivative contracted with v_1, v_2, v_3.
    """

    def __init__(self, coefficients):
        self.coefficients = coefficients  # array of shape (8, ...) indexed by bit mask

    @classmethod
    def constant(cls, value):
        value = np.asarray(value, dtype=complex)
        coefficients = np.zeros((8,) + value.shape, dtype=complex)
        coefficients[0] = value
        return cls(coefficients)

    def __add__(self, other):
        other = other if isinstance(other, MultilinearExpansion) else MultilinearExpansion.constant(other)
        return MultilinearExpansion(self.coefficients + other.coefficients)

    __radd__ = __add__

    def __neg__(self):
        return MultilinearExpansion(-self.coefficients)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def _product(self, other, multiply):
        first = self.coefficients
        second = other.coefficients
        shape = np.shape(multiply(first[0], second[0]))
        result = np.zeros((8,) + shape, dtype=complex)
        for left, right in itertools.product(SUBSETS, SUBSETS):
            if left & right == 0:
                result[left | right] += multiply(first[left], second[right])
        return MultilinearExpansion(result)

    def __mul__(self, other):
        if not isinstance(other, MultilinearExpansion):
            return MultilinearExpansion(self.coefficients * other)
        return self._product(other, np.multiply)

    __rmul__ = __mul__

    def __truediv__(self, number):
        return MultilinearExpansion(self.coefficients / number)

    def matmul(self, other):
        return self._product(other, np.matmul)

    def times_matrix(self, matrix):
        """Scalar expansion times a constant matrix."""
        return MultilinearExpansion(np.einsum('m,ij->mij', self.coefficients, matrix))

    def trace(self):
        return MultilinearExpansion(np.trace(self.coefficients, axis1=1, axis2=2))

    def dagger(self):
        return MultilinearExpansion(np.conj(np.swapaxes(self.coefficients, 1, 2)))

    def apply(self, derivatives):
        """f(x0 + delta) from the values f(x0), f'(x0), f''(x0), f'''(x0) (scalar expansions)."""
        delta = self - self.coefficients[0]
        result = MultilinearExpansion.constant(derivatives[0])
        power = MultilinearExpansion.constant(1.0)
        factorial = 1.0
        for order in range(1, 4):
            power = power * delta
            factorial *= order
            result = result + power * (derivatives[order] / factorial)
        return result


def matrix_exponential(generator):
    """exp(X) for a matrix expansion with vanishing constant part (X^4 = 0)."""
    identity = MultilinearExpansion.constant(np.eye(3))
    square = generator.matmul(generator)
    return identity + generator + square / 2 + square.matmul(generator) / 6


def field_point(directions):
    """Field expansions x^A = sum_n t_n v_n^A for up to three direction vectors."""
    coefficients = np.zeros((8, NUMBER_OF_FIELDS))
    for position, direction in enumerate(directions):
        coefficients[1 << position] = direction
    return [MultilinearExpansion(coefficients[:, index].astype(complex)) for index in range(NUMBER_OF_FIELDS)]


def source_matrices(sources, radial_polar, angular_polar, decay_constant):
    """S(sigma, a) and P(sigma, a) as matrix expansions from the coefficient dictionaries."""
    scalar = MultilinearExpansion.constant(np.zeros((3, 3)))
    pseudoscalar = MultilinearExpansion.constant(np.zeros((3, 3)))
    radial_ratio = radial_polar / decay_constant
    angular_ratio = angular_polar / decay_constant
    for (r, s) in set(sources['S']) | set(sources['P']):
        monomial = MultilinearExpansion.constant(1.0)
        for _ in range(r):
            monomial = monomial * radial_ratio
        for _ in range(s):
            monomial = monomial * angular_ratio
        if (r, s) in sources['S']:
            scalar = scalar + monomial.times_matrix(sources['S'][(r, s)])
        if (r, s) in sources['P']:
            pseudoscalar = pseudoscalar + monomial.times_matrix(sources['P'][(r, s)])
    return scalar, pseudoscalar


def potential_lagrangian(fields, parameters, coordinates):
    """Non-derivative Lagrangian (source, anomalous singlet mass, radial potential) at the fields."""
    decay_constant = parameters['f_a']
    radial, angular = fields[0], fields[1]
    if coordinates == 'polar':
        radial_polar, angular_polar = radial, angular
        rho_squared = (radial + decay_constant) * (radial + decay_constant)
    else:
        rho_squared = (radial + decay_constant) * (radial + decay_constant) + angular * angular
        rho_value = np.sqrt(rho_squared.coefficients[0].real)
        radial_polar = rho_squared.apply([rho_value, 1 / (2 * rho_value), -1 / (4 * rho_value ** 3),
                                          3 / (8 * rho_value ** 5)]) - decay_constant
        inverse = (radial + decay_constant).apply([1 / decay_constant, -1 / decay_constant ** 2,
                                                   2 / decay_constant ** 3, -6 / decay_constant ** 4])
        angular_polar = decay_constant * (angular * inverse).apply([0.0, 1.0, 0.0, -2.0])
    meson_matrix = MultilinearExpansion.constant(np.zeros((3, 3)))
    for k in range(9):
        meson_matrix = meson_matrix + fields[2 + k].times_matrix(LAMBDA[k])
    chiral_field = matrix_exponential(meson_matrix * (1j / parameters['f_pi']))
    chiral_field_dagger = matrix_exponential(meson_matrix * (-1j / parameters['f_pi']))
    scalar, pseudoscalar = source_matrices(parameters['sources'], radial_polar, angular_polar, decay_constant)
    chi = (scalar + pseudoscalar * 1j) * (2 * parameters['B_0'])
    chi_plus_trace = (chi.matmul(chiral_field_dagger) + chi.dagger().matmul(chiral_field)).trace()
    singlet_trace = meson_matrix.trace()
    radial_potential = (rho_squared * (parameters['m_phi_squared'] / 2)
                        + rho_squared * rho_squared * (parameters['lambda'] / 4))
    return (chi_plus_trace * (parameters['f_pi'] ** 2 / 4) - singlet_trace * singlet_trace * (parameters['m_1_squared'] / 12)
            - radial_potential)


def kinetic_metric(fields, parameters, coordinates):
    """Field-dependent kinetic metric G_AB with L_kin = (1/2) G_AB d x^A d x^B."""
    metric = [[MultilinearExpansion.constant(float(row == column)) for column in range(NUMBER_OF_FIELDS)]
              for row in range(NUMBER_OF_FIELDS)]
    if coordinates == 'polar':
        factor = 1 + fields[0] / parameters['f_a']
        metric[1][1] = factor * factor
    return metric


def first_and_second_derivatives(parameters, coordinates):
    """Gradient and Hessian of the non-derivative Lagrangian at the vacuum (all fields)."""
    identity = np.eye(NUMBER_OF_FIELDS)
    gradient = np.array([potential_lagrangian(field_point([identity[a]]), parameters, coordinates).coefficients[1].real
                         for a in range(NUMBER_OF_FIELDS)])
    hessian = np.zeros((NUMBER_OF_FIELDS, NUMBER_OF_FIELDS))
    for a in range(NUMBER_OF_FIELDS):
        for b in range(a, NUMBER_OF_FIELDS):
            value = potential_lagrangian(field_point([identity[a], identity[b]]), parameters, coordinates)
            hessian[a, b] = hessian[b, a] = value.coefficients[3].real
    return gradient, hessian


def pole_vectors(hessian):
    """Masses and unit pole vectors of the angular and meson block (kinetic matrix is the identity).

    Vectors are embedded in the 11 field components; the angular component is made positive.
    """
    mass_matrix = -hessian[1:, 1:]
    masses_squared, vectors = np.linalg.eigh(mass_matrix)
    embedded = np.zeros((NUMBER_OF_FIELDS, len(masses_squared)))
    embedded[1:] = vectors * np.sign(vectors[0])
    return masses_squared, embedded


def direct_on_shell_vertex(parameters, coordinates, radial_mass_squared, masses_squared, vectors, i, j):
    """g_ij from the trilinear part of the Lagrangian contracted with the pole vectors."""
    radial_direction = np.eye(NUMBER_OF_FIELDS)[0]
    directions = [radial_direction, vectors[:, i], vectors[:, j]]
    fields = field_point(directions)
    potential_part = potential_lagrangian(fields, parameters, coordinates).coefficients[7].real
    # Outgoing momenta k_1 = -(p_i + p_j), k_2 = p_i, k_3 = p_j.
    p_i_dot_p_j = (radial_mass_squared - masses_squared[i] - masses_squared[j]) / 2
    momentum_products = {(1, 2): p_i_dot_p_j, (0, 1): -(masses_squared[i] + p_i_dot_p_j),
                         (0, 2): -(masses_squared[j] + p_i_dot_p_j)}
    metric = kinetic_metric(fields, parameters, coordinates)
    kinetic_part = 0.0
    for leg in range(3):
        others = [n for n in range(3) if n != leg]
        metric_derivative = np.array([[metric[a][b].coefficients[1 << leg].real for b in range(NUMBER_OF_FIELDS)]
                                      for a in range(NUMBER_OF_FIELDS)])
        contraction = directions[others[0]] @ metric_derivative @ directions[others[1]]
        # (1/2) * 2 orderings * (i k_m)(i k_n) = -k_m.k_n
        kinetic_part += -contraction * momentum_products[tuple(sorted(others))]
    return -(potential_part + kinetic_part)


def mass_operator_matrix(parameters, momentum_squared):
    """Matrix of M_s on the basis lambda_k: (M_s)_{kl} = (1/2) <lambda_k M_s(lambda_l)>."""
    operator = np.zeros((9, 9))
    S_00 = parameters['sources']['S'][(0, 0)]
    for column in range(9):
        image = (parameters['B_0'] * (S_00 @ LAMBDA[column] + LAMBDA[column] @ S_00)
                 + parameters['m_1_squared'] / 3 * np.trace(LAMBDA[column]) * np.eye(3)
                 - momentum_squared * LAMBDA[column])
        operator[:, column] = [0.5 * np.trace(LAMBDA[row] @ image).real for row in range(9)]
    return operator


def inverse_mass_operator(parameters, momentum_squared, matrix):
    """M_s^{-1}(X) for a Hermitian matrix X, returned as a matrix."""
    components = np.array([0.5 * np.trace(LAMBDA[k] @ matrix).real for k in range(9)])
    solution = np.linalg.solve(mass_operator_matrix(parameters, momentum_squared), components)
    return np.einsum('k,kij->ij', solution, LAMBDA)


def trace(matrix):
    return np.trace(matrix).real


def pole_shift_matrix(parameters, momentum_squared):
    """Omega_r = (2 B0 f_pi/f_a) M_{m_r^2}^{-1}(P_01)."""
    return (2 * parameters['B_0'] * parameters['f_pi'] / parameters['f_a']
            * inverse_mass_operator(parameters, momentum_squared, parameters['sources']['P'][(0, 1)]))


def pole_normalization(parameters, momentum_squared):
    """N(m^2) = 1 + (2 B0^2 f_pi^2/f_a^2) <M^{-1}(P_01)^2>."""
    inverse_image = inverse_mass_operator(parameters, momentum_squared, parameters['sources']['P'][(0, 1)])
    return 1 + 2 * parameters['B_0'] ** 2 * parameters['f_pi'] ** 2 / parameters['f_a'] ** 2 * trace(inverse_image @ inverse_image)


def pole_equation(parameters, momentum_squared):
    """s + (2 B0 f_pi^2/f_a^2)[<S_02> + B0 <P_01 M_s^{-1}(P_01)>] (vanishes at a pole)."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    P_01 = parameters['sources']['P'][(0, 1)]
    return momentum_squared + 2 * B0 * fpi ** 2 / fa ** 2 * (
        trace(parameters['sources']['S'][(0, 2)]) + B0 * trace(P_01 @ inverse_mass_operator(parameters, momentum_squared, P_01)))


def angular_row(parameters, momentum_squared):
    """m_r^2 + (2B0 f_pi^2/f_a^2)<S_02> + (B0 f_pi/f_a)<P_01 Omega_r> (vanishes at a pole)."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    return (momentum_squared + 2 * B0 * fpi ** 2 / fa ** 2 * trace(parameters['sources']['S'][(0, 2)])
            + B0 * fpi / fa * trace(parameters['sources']['P'][(0, 1)] @ pole_shift_matrix(parameters, momentum_squared)))


def radial_mass_squared_trace_form(parameters):
    """m_sigma^2 = 2 lambda f_a^2 + (B0 f_pi^2/f_a^2)(<S_10> - 2<S_20>)."""
    S = parameters['sources']['S']
    zero = np.zeros((3, 3))
    return (2 * parameters['lambda'] * parameters['f_a'] ** 2 + parameters['B_0'] * parameters['f_pi'] ** 2
            / parameters['f_a'] ** 2 * (trace(S.get((1, 0), zero)) - 2 * trace(S.get((2, 0), zero))))


def polar_vertex_trace_form(parameters, radial_mass_squared, mass_i_squared, mass_j_squared):
    """g^polar_ij of the trace formula (Omega_r expressed through M^{-1}_{m_r^2}(P_01))."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    S, P = parameters['sources']['S'], parameters['sources']['P']
    omega_i = pole_shift_matrix(parameters, mass_i_squared)
    omega_j = pole_shift_matrix(parameters, mass_j_squared)
    braces = (B0 / 2 * trace(S[(1, 0)] @ (omega_i @ omega_j + omega_j @ omega_i))
              - B0 * fpi / fa * trace(P[(1, 1)] @ (omega_i + omega_j))
              - 2 * B0 * fpi ** 2 / fa ** 2 * trace(S[(1, 2)])
              + radial_mass_squared - mass_i_squared - mass_j_squared)
    return braces / fa / np.sqrt(pole_normalization(parameters, mass_i_squared) * pole_normalization(parameters, mass_j_squared))


def cartesian_vertex_trace_form(parameters, mass_i_squared, mass_j_squared):
    """g^Cart_ij of the trace formula for general sources."""
    B0, fpi, fa = parameters['B_0'], parameters['f_pi'], parameters['f_a']
    S, P = parameters['sources']['S'], parameters['sources']['P']
    zero = np.zeros((3, 3))
    omega_i = pole_shift_matrix(parameters, mass_i_squared)
    omega_j = pole_shift_matrix(parameters, mass_j_squared)
    braces = (2 * parameters['lambda'] * fa
              + B0 / (2 * fa) * trace(S[(1, 0)] @ (omega_i @ omega_j + omega_j @ omega_i))
              - B0 * fpi / fa ** 2 * trace((P[(1, 1)] - P[(0, 1)]) @ (omega_i + omega_j))
              + B0 * fpi ** 2 / fa ** 3 * (trace(S[(1, 0)]) - 2 * trace(S.get((2, 0), zero))
                                          + 4 * trace(S[(0, 2)]) - 2 * trace(S[(1, 2)])))
    return braces / np.sqrt(pole_normalization(parameters, mass_i_squared) * pole_normalization(parameters, mass_j_squared))


def random_hermitian(generator, scale):
    matrix = generator.normal(size=(3, 3)) + 1j * generator.normal(size=(3, 3))
    return scale * (matrix + matrix.conj().T) / 2


def random_general_parameters(seed=2026):
    """Random non-commuting Hermitian sources (GeV) obeying the stated assumptions.

    P_00 = 0, <S_01> = 0, P_10 = 0, <S_11> = 0; every other coefficient with r + s <= 3 is a
    random Hermitian matrix.  m_phi^2 is fixed by stationarity in the radial direction.
    """
    generator = np.random.default_rng(seed)
    unitary, _ = np.linalg.qr(generator.normal(size=(3, 3)) + 1j * generator.normal(size=(3, 3)))
    sources = {'S': {}, 'P': {}}
    sources['S'][(0, 0)] = unitary @ np.diag([0.004, 0.008, 0.1]) @ unitary.conj().T
    for (r, s) in itertools.product(range(4), range(4)):
        if r + s > 3 or (r, s) == (0, 0):
            continue
        sources['S'][(r, s)] = random_hermitian(generator, 0.02)
        sources['P'][(r, s)] = random_hermitian(generator, 0.02)
    for key in [(0, 1), (1, 1)]:
        sources['S'][key] = sources['S'][key] - trace(sources['S'][key]) / 3 * np.eye(3)
    sources['P'][(1, 0)] = np.zeros((3, 3))
    sources['S'][(0, 2)] = sources['S'][(0, 2)] - 0.3 * np.eye(3)  # positive bare angular mass
    parameters = {'B_0': 2.6, 'f_pi': 0.0924, 'f_a': 0.4, 'm_1_squared': 0.75, 'lambda': 0.8, 'sources': sources}
    return impose_stationarity(parameters)


def only_up_parameters(up_mass=0.0022, down_mass=0.0047, strange_mass=0.095):
    """Sources of the only-up model, M = diag[m_u (1 + sigma/f_a) e^{-i a/f_a}, m_d, m_s] (GeV)."""
    projector = np.diag([1.0, 0.0, 0.0])
    sources = {'S': {(0, 0): np.diag([up_mass, down_mass, strange_mass])}, 'P': {}}
    for (r, s) in itertools.product(range(2), range(4)):
        if r + s > 3:
            continue
        if (r, s) != (0, 0):
            sources['S'][(r, s)] = up_mass * np.cos(np.pi * s / 2) / _factorial(s) * projector
        sources['P'][(r, s)] = -up_mass * np.sin(np.pi * s / 2) / _factorial(s) * projector
    parameters = {'B_0': 2.6, 'f_pi': 0.0924, 'f_a': 0.4, 'm_1_squared': 0.75, 'lambda': 0.8, 'sources': sources}
    return impose_stationarity(parameters)


def _factorial(number):
    result = 1
    for factor in range(2, number + 1):
        result *= factor
    return result


def impose_stationarity(parameters):
    """m_phi^2 = B0 f_pi^2 <S_10>/f_a^2 - lambda f_a^2 (radial stationarity)."""
    parameters['m_phi_squared'] = (parameters['B_0'] * parameters['f_pi'] ** 2 * trace(parameters['sources']['S'][(1, 0)])
                                   / parameters['f_a'] ** 2 - parameters['lambda'] * parameters['f_a'] ** 2)
    return parameters


def compare_vertices(parameters):
    """Return a dict of residuals and the vertex tables for one parameter set."""
    results = {}
    gradients, hessians = {}, {}
    for coordinates in ('polar', 'cartesian'):
        gradients[coordinates], hessians[coordinates] = first_and_second_derivatives(parameters, coordinates)
    results['tadpole'] = max(np.max(np.abs(g)) for g in gradients.values())
    results['hessian_difference'] = np.max(np.abs(hessians['polar'] - hessians['cartesian']))
    results['radial_mixing'] = np.max(np.abs(hessians['polar'][0, 1:]))
    radial_mass_squared = -hessians['polar'][0, 0]
    results['radial_mass_difference'] = abs(radial_mass_squared - radial_mass_squared_trace_form(parameters))
    masses_squared, vectors = pole_vectors(hessians['polar'])
    results['masses_squared'] = masses_squared
    results['radial_mass_squared'] = radial_mass_squared
    # Poles with no angular component (e.g. one state of each degenerate pair of B0{S_00, .})
    # are exact meson poles at which M_{m^2} is singular; the trace form does not apply there.
    mixed = [r for r in range(len(masses_squared)) if abs(vectors[1, r]) > ANGULAR_COMPONENT_THRESHOLD]
    results['mixed_poles'] = mixed
    results['pole_equation'] = max(abs(pole_equation(parameters, masses_squared[r])) for r in mixed)
    results['angular_row'] = max(abs(angular_row(parameters, masses_squared[r])) for r in mixed)
    shift_differences, normalization_differences = [], []
    for r in mixed:
        direction = vectors[2:, r] / vectors[1, r]
        shift = pole_shift_matrix(parameters, masses_squared[r])
        shift_differences.append(np.max(np.abs(np.einsum('k,kij->ij', direction, LAMBDA) - shift))
                                 / np.max(np.abs(shift)))
        normalization_differences.append(abs(1 / vectors[1, r] ** 2 - pole_normalization(parameters, masses_squared[r]))
                                          / pole_normalization(parameters, masses_squared[r]))
    results['pole_direction'] = max(shift_differences)
    results['normalization'] = max(normalization_differences)
    count = len(masses_squared)
    tables = {name: np.full((count, count), np.nan)
              for name in ('polar_direct', 'cartesian_direct', 'polar_trace', 'cartesian_trace')}
    for i in range(count):
        for j in range(i, count):
            for name, coordinates in (('polar_direct', 'polar'), ('cartesian_direct', 'cartesian')):
                tables[name][i, j] = tables[name][j, i] = direct_on_shell_vertex(
                    parameters, coordinates, radial_mass_squared, masses_squared, vectors, i, j)
            if i in mixed and j in mixed:
                tables['polar_trace'][i, j] = tables['polar_trace'][j, i] = polar_vertex_trace_form(
                    parameters, radial_mass_squared, masses_squared[i], masses_squared[j])
                tables['cartesian_trace'][i, j] = tables['cartesian_trace'][j, i] = cartesian_vertex_trace_form(
                    parameters, masses_squared[i], masses_squared[j])
    block = np.ix_(mixed, mixed)
    scale = np.max(np.abs(tables['polar_direct']))
    results['vertex_scale'] = scale
    results['polar_direct_vs_cartesian_direct'] = np.max(np.abs(tables['polar_direct'] - tables['cartesian_direct'])) / scale
    results['polar_direct_vs_trace'] = np.max(np.abs(tables['polar_direct'][block] - tables['polar_trace'][block])) / scale
    results['cartesian_direct_vs_trace'] = np.max(np.abs(tables['cartesian_direct'][block] - tables['cartesian_trace'][block])) / scale
    results['tables'] = tables
    return results


def main():
    print('Units: GeV. Tree level, LO U(3) chiral Lagrangian, no expansion in f_pi/f_a.\n')
    sigma_series, angular_series, metric_difference = derive_coordinate_map()
    print('sigma(s,p) =', sigma_series)
    print('a(s,p)     =', angular_series)
    print('metric difference =', metric_difference)
    data = derive_general_cubic_lagrangians()
    print('m_phi^2 (stationarity) =', data['mass_phi_squared'])
    print('m_sigma^2 =', data['radial_mass_squared'])
    print('Cartesian s p^2 coefficient =', data['cartesian_s_p_p'])
    print('Cartesian s p <P_rs Phi> coefficients =', {k: v for k, v in data['cartesian_s_p_meson'].items() if v != 0})
    print('Cartesian s <S_rs Phi^2> coefficients =', {k: v for k, v in data['cartesian_s_meson_pair'].items() if v != 0})
    _, _, difference, reduced = derive_vertex_difference()
    print('(g^polar - g^Cart) sqrt(N_i N_j) =', sp.factor(difference))
    print('after pole and radial-mass relations =', reduced, '\n')
    for label, parameters in (('random general sources', random_general_parameters()),
                              ('only-up model', only_up_parameters())):
        results = compare_vertices(parameters)
        print(label)
        print('  pole masses^2 [GeV^2] =', np.round(results['masses_squared'], 6), ' m_sigma^2 =', round(results['radial_mass_squared'], 6))
        print('  poles with an angular component:', results['mixed_poles'])
        for key in ('tadpole', 'hessian_difference', 'radial_mixing', 'radial_mass_difference', 'pole_equation',
                    'angular_row', 'pole_direction', 'normalization', 'polar_direct_vs_cartesian_direct',
                    'polar_direct_vs_trace', 'cartesian_direct_vs_trace'):
            print(f'  {key:34s} {np.max(results[key]):.2e}')


if __name__ == '__main__':
    main()
