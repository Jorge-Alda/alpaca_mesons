"""Field-redefinition invariant topological coefficient dressed by the axial connection.

Run (from the repository root):
    .venv/bin/python scripts/derive_alp_dressed_topological_invariant.py

Conventions and units
---------------------
Natural units; x = a/f_a is dimensionless and real. All flavor matrices are
noncommuting (3x3 for the exact rational checks). The quark mass function
M(x) = S(x) + i P(x) multiplies bar(q_R) M q_L and has mass dimension one.
C_R(x), C_L(x) are the dimensionless Hermitian coefficients of the
derivative currents (d_mu a / f_a) bar(q_eta) gamma^mu C_eta q_eta, with
C_V = (C_R + C_L)/2 and C_A = (C_R - C_L)/2. Theta(x) is the coefficient of
G tilde(G) in the Alda:2026nrd normalization, so Theta_1 = c_G. Taylor
coefficients X_n are derivatives d^n X/dx^n at x = 0, i.e.
X(x) = sum_n X_n x^n / n!.

Under q_eta' = g_eta(x) q_eta with g_eta(0) = 1:
    C_eta' = g_eta C_eta g_eta^dagger + i g_eta d_x g_eta^dagger,
    M' = g_R M g_L^dagger,
    Theta'(x) = Theta(x) - (1/(2i)) log(det g_R / det g_L),
with the logarithm continuous from x = 0. The transport matrices solve
dW_eta/dx = i C_eta W_eta, W_eta(0) = 1, and hat(M) = W_R^dagger M W_L.

Claim verified here: the dressed topological coefficient
    hat(Theta)(x) = Theta(x) + int_0^x tr C_A(x') dx'
is invariant under all such redefinitions, and its Taylor coefficients are
hat(Theta)_n = Theta_n + tr C_{A,n-1} (n >= 1), hat(Theta)_1 = c_GA.

Inputs
------
* memory/alp/shared/alp-general-chiral-redefinition-invariants.md:
  connection law for C_eta, homogeneous law for M, transport W_eta, hat(M),
  and hat(U) = W_R^dagger U W_L.
* memory/alp/shared/alp-qcd-topological-cp.md: anomaly law for Theta.
* memory/alp/shared/alp-dressed-mass-first-order-invariants.md:
  hat(P)_1 = -{c^AP, m}.
* memory/extracts/Alda:2026nrd/appendix-a-expansion.md: source
  parametrization M_1 = i[c^S, m] + i{c^P, m}, C_A(0) = c^A, c^AP = c^A - c^P,
  and the tadpole source C = -B_0 {m, c^AP} - (m_0^2/N_f) c_GA 1.
* memory/extracts/Alda:2026nrd/lo-topological-action.md:
  X_TC = 2i (tr log u + i c_G a/f_a), U = u^2.

Assumptions and approximations
------------------------------
* All series are truncated at O(x^4) (SERIES_ORDER); the exact rational
  checks are therefore checks through that order, not to all orders.
* Assumption (extension of Alda:2026nrd): in X_TC the linear term
  c_G a/f_a is replaced by the full function Theta(x), i.e.
  X_TC = 2i (tr log u + i Theta(x)), with tr log u := (1/2) log det U.
  Alda:2026nrd only defines the linear form.
* The physical vacuum angle is written theta_bar(x) = Theta(x) + (1/2) arg det M(x),
  the normalization fixed by the anomaly law above; arg det M(0) = 0 for the
  positive-definite m used here.

Each check prints PASS or FAIL; the process exits with status 1 on failure.
"""

import random
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from derive_alp_invariant_coefficient_functions import transport_series  # noqa: E402
from derive_alp_dressed_mass_first_order import (  # noqa: E402
    anticommutator,
    commutator,
    hermitian_matrix_symbol,
    hermitian_parts,
    random_rational_hermitian,
)

ALP_FIELD_RATIO = sp.Symbol("x", real=True)
INTEGRATION_VARIABLE = sp.Symbol("t", real=True)
SERIES_ORDER = 4
FLAVOR_NUMBER = 3
DEFAULT_SEED = 20261111

# Noncommuting Taylor coefficients C_0..C_{n-1} used to build W via transport_series.
TRANSPORT_SYMBOLS = sp.symbols(f"C_0:{SERIES_ORDER}", commutative=False)


# ---------------------------------------------------------------------------
# Truncated power series in x with scalar or explicit-matrix coefficients
# ---------------------------------------------------------------------------

def truncate(expression, order=SERIES_ORDER):
    """Drop powers of x above `order` from a polynomial in x."""
    expression = sp.expand(expression)
    return sp.Add(*(expression.coeff(ALP_FIELD_RATIO, degree) * ALP_FIELD_RATIO**degree
                    for degree in range(order + 1)))


def truncate_matrix(matrix, order=SERIES_ORDER):
    return matrix.applyfunc(lambda entry: truncate(entry, order))


def taylor_polynomial(coefficients):
    """sum_n X_n x^n / n! for scalar or matrix coefficients X_n."""
    return sum((coefficient * ALP_FIELD_RATIO**degree / sp.factorial(degree)
                for degree, coefficient in enumerate(coefficients)),
               start=0 * coefficients[0])


def taylor_coefficients(matrix_series, order=SERIES_ORDER):
    """Derivatives at x = 0 of a matrix polynomial, through `order`."""
    expanded = matrix_series.applyfunc(sp.expand)
    return [expanded.applyfunc(lambda entry: entry.coeff(ALP_FIELD_RATIO, degree))
            * sp.factorial(degree) for degree in range(order + 1)]


def integrate_from_zero(expression):
    """int_0^x expression(x') dx' for a polynomial in x."""
    return sp.integrate(expression.subs(ALP_FIELD_RATIO, INTEGRATION_VARIABLE),
                        (INTEGRATION_VARIABLE, 0, ALP_FIELD_RATIO))


def exponential_series(scalar_series, order=SERIES_ORDER):
    """exp(f) truncated, for a scalar series f with f(0) = 0."""
    result, power = sp.Integer(1), sp.Integer(1)
    for degree in range(1, order + 1):
        power = truncate(power * scalar_series, order)
        result += power / sp.factorial(degree)
    return truncate(result, order)


def logarithm_series(scalar_series, order=SERIES_ORDER):
    """log(f) truncated, for a scalar series f with f(0) = 1 (continuous branch)."""
    deviation = truncate(scalar_series - 1, order)
    result, power = sp.Integer(0), sp.Integer(1)
    for degree in range(1, order + 1):
        power = truncate(power * deviation, order)
        result += (-1)**(degree + 1) * power / degree
    return truncate(result, order)


def matrix_exponential_series(generator_series, order=SERIES_ORDER):
    """exp(i G(x)) truncated, for a matrix series G with G(0) = 0."""
    size = generator_series.shape[0]
    result, power = sp.eye(size), sp.eye(size)
    for degree in range(1, order + 1):
        power = truncate_matrix(power * sp.I * generator_series, order)
        result += power / sp.factorial(degree)
    return truncate_matrix(result, order)


def determinant_series(matrix_series, order=SERIES_ORDER):
    return truncate(matrix_series.det(method="berkowitz"), order)


def imaginary_part_series(scalar_series):
    """Imaginary part of a series in real x, coefficient by coefficient."""
    expanded = sp.expand(scalar_series)
    return sp.Add(*(sp.im(expanded.coeff(ALP_FIELD_RATIO, degree)) * ALP_FIELD_RATIO**degree
                    for degree in range(SERIES_ORDER + 1)))


def evaluate_noncommutative(expression, substitution, size):
    """Evaluate a polynomial in noncommuting symbols on explicit matrices."""
    expression = sp.expand(expression)
    if expression.is_Add:
        return sum((evaluate_noncommutative(term, substitution, size)
                    for term in expression.args), start=sp.zeros(size, size))
    if expression in substitution:
        return substitution[expression]
    if expression.is_Pow and expression.base in substitution:
        return substitution[expression.base] ** int(expression.exp)
    if expression.is_Mul:
        commutative_factors, ordered_factors = expression.args_cnc()
        result = sp.Mul(*commutative_factors) * sp.eye(size)
        for factor in ordered_factors:
            result = result * evaluate_noncommutative(factor, substitution, size)
        return result
    return expression * sp.eye(size)


def matrix_transport_series(current_coefficients, order=SERIES_ORDER):
    """W(x) with dW/dx = i C W, W(0) = 1, from transport_series evaluated on matrices."""
    size = current_coefficients[0].shape[0]
    symbolic_transport = transport_series(TRANSPORT_SYMBOLS[:order], ALP_FIELD_RATIO, order)
    substitution = dict(zip(TRANSPORT_SYMBOLS[:order], current_coefficients[:order]))
    return truncate_matrix(evaluate_noncommutative(symbolic_transport, substitution, size), order)


# ---------------------------------------------------------------------------
# Exact rational test data
# ---------------------------------------------------------------------------

def random_rational_complex(generator, size):
    return (random_rational_hermitian(generator, size)
            + sp.I * random_rational_hermitian(generator, size))


def positive_mass_matrix(generator, size):
    auxiliary = random_rational_hermitian(generator, size)
    return (auxiliary.H * auxiliary + sp.eye(size)).applyfunc(sp.expand)


def random_theory(seed=DEFAULT_SEED, size=FLAVOR_NUMBER):
    """Polynomial Hermitian C_R(x), C_L(x), mass function M(x), and Theta(x)."""
    generator = random.Random(seed)
    right_coefficients = [random_rational_hermitian(generator, size)
                          for _ in range(SERIES_ORDER + 1)]
    left_coefficients = [random_rational_hermitian(generator, size)
                         for _ in range(SERIES_ORDER + 1)]
    mass_coefficients = ([positive_mass_matrix(generator, size)]
                         + [random_rational_complex(generator, size)
                            for _ in range(SERIES_ORDER)])
    theta_coefficients = [sp.Rational(generator.randint(-9, 9), generator.randint(1, 5))
                          for _ in range(SERIES_ORDER + 1)]
    return {"right": right_coefficients, "left": left_coefficients,
            "mass": mass_coefficients, "theta": theta_coefficients}


def random_redefinition(seed, size=FLAVOR_NUMBER, traceless=False):
    """Hermitian G_R(x), G_L(x) = sum_{n>=1} G_n x^n/n!, noncommuting coefficients.

    For traceless=False the singlet part is nonlinear in x: extra identity
    pieces are added at orders x^2, x^3, x^4 with different weights for R and L.
    """
    generator = random.Random(seed)
    generators = {}
    for chirality in ("right", "left"):
        coefficients = [sp.zeros(size, size)]
        for degree in range(1, SERIES_ORDER + 1):
            coefficient = random_rational_hermitian(generator, size)
            if traceless:
                coefficient = coefficient - coefficient.trace() / size * sp.eye(size)
            elif degree >= 2:
                coefficient = coefficient + sp.Rational(generator.randint(-9, 9),
                                                        generator.randint(1, 5)) * sp.eye(size)
            coefficients.append(coefficient)
        generators[chirality] = taylor_polynomial(coefficients)
    return generators["right"], generators["left"]


def connection_transform(group_series, current_series):
    """g C g^dagger + i g d_x g^dagger, truncated."""
    return truncate_matrix(group_series * current_series * group_series.H
                           + sp.I * group_series * sp.diff(group_series.H, ALP_FIELD_RATIO))


def axial_trace_integral(right_series, left_series):
    """int_0^x tr C_A."""
    return sp.expand(integrate_from_zero((right_series.trace() - left_series.trace()) / 2))


def dressed_theta(theta_series, right_series, left_series):
    return truncate(theta_series + axial_trace_integral(right_series, left_series))


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def check_liouville(seed=DEFAULT_SEED):
    """Check 1: det W_eta = exp(i int_0^x tr C_eta) through O(x^4)."""
    theory = random_theory(seed)
    residuals = {}
    for chirality in ("right", "left"):
        transport = matrix_transport_series(theory[chirality])
        current_series = taylor_polynomial(theory[chirality])
        ode_residual = truncate_matrix(sp.diff(transport, ALP_FIELD_RATIO)
                                       - sp.I * current_series * transport,
                                       SERIES_ORDER - 1)
        determinant_residual = sp.expand(
            determinant_series(transport)
            - exponential_series(sp.I * integrate_from_zero(current_series.trace())))
        residuals[chirality] = {"ode_is_zero": ode_residual == sp.zeros(*ode_residual.shape),
                                "determinant_residual": determinant_residual}
    return residuals


def check_redefinition_invariance(seed=DEFAULT_SEED, redefinition_seed=7, traceless=False):
    """Check 2: hat(Theta) invariant under g_eta = exp(i G_eta(x)) through O(x^4)."""
    theory = random_theory(seed)
    right_series = taylor_polynomial(theory["right"])
    left_series = taylor_polynomial(theory["left"])
    theta_series = taylor_polynomial(theory["theta"])
    right_generator, left_generator = random_redefinition(redefinition_seed,
                                                          traceless=traceless)
    right_group = matrix_exponential_series(right_generator)
    left_group = matrix_exponential_series(left_generator)

    right_series_new = connection_transform(right_group, right_series)
    left_series_new = connection_transform(left_group, left_series)
    anomaly_shift = sp.expand(
        (logarithm_series(determinant_series(right_group))
         - logarithm_series(determinant_series(left_group))) / (2 * sp.I))
    theta_series_new = truncate(theta_series - anomaly_shift)

    # Covariance of the transport matrices, W' = g W.
    right_transport_new = matrix_transport_series(taylor_coefficients(right_series_new))
    right_transport = matrix_transport_series(theory["right"])
    covariance = truncate_matrix(right_transport_new - right_group * right_transport)

    singlet_axial_angle = sp.expand((right_generator.trace() - left_generator.trace()) / 2)
    return {
        "anomaly_shift_residual": sp.expand(anomaly_shift - truncate(singlet_axial_angle)),
        "theta_changes": sp.expand(theta_series_new - theta_series) != 0,
        "transport_covariance_is_zero": covariance == sp.zeros(*covariance.shape),
        "dressed_theta_residual": sp.expand(
            dressed_theta(theta_series_new, right_series_new, left_series_new)
            - dressed_theta(theta_series, right_series, left_series)),
    }


def symbolic_hermitian_matrix(prefix, size=FLAVOR_NUMBER):
    matrix = sp.zeros(size, size)
    for row in range(size):
        matrix[row, row] = sp.Symbol(f"{prefix}_{row}{row}", real=True)
        for column in range(row + 1, size):
            real_part = sp.Symbol(f"{prefix}r_{row}{column}", real=True)
            imaginary_part = sp.Symbol(f"{prefix}i_{row}{column}", real=True)
            matrix[row, column] = real_part + sp.I * imaginary_part
            matrix[column, row] = real_part - sp.I * imaginary_part
    return matrix


def check_taylor_form(order=6):
    """Check 3: hat(Theta)_n = Theta_n + tr C_{A,n-1}; hat(Theta)_1 = c_GA and its invariance."""
    theta_coefficients = sp.symbols(f"Theta_0:{order + 1}")
    axial_trace_coefficients = sp.symbols(f"trC_A0:{order + 1}")
    dressed = (taylor_polynomial(list(theta_coefficients))
               + integrate_from_zero(taylor_polynomial(list(axial_trace_coefficients))))
    taylor_residuals = [sp.expand(sp.diff(dressed, ALP_FIELD_RATIO, degree).subs(ALP_FIELD_RATIO, 0)
                                  - theta_coefficients[degree]
                                  - (axial_trace_coefficients[degree - 1] if degree else 0))
                        for degree in range(order + 1)]

    gluon_coupling = sp.Symbol("c_G", real=True)
    axial_coupling = symbolic_hermitian_matrix("cA")
    axial_generator = symbolic_hermitian_matrix("kappa")
    first_order = gluon_coupling + axial_coupling.trace()
    gluon_axial_invariant = sp.Symbol("c_G", real=True) + axial_coupling.trace()
    first_order_shifted = ((gluon_coupling - axial_generator.trace())
                           + (axial_coupling + axial_generator).trace())
    return {"taylor_residuals": taylor_residuals,
            "first_order_residual": sp.expand(first_order - gluon_axial_invariant),
            "shift_residual": sp.expand(first_order_shifted - first_order)}


def check_vacuum_angle(seed=DEFAULT_SEED, redefinition_seed=7):
    """Check 4 (exact rational): theta_bar invariant; hat(Theta) + arg det hat(M)/2 = theta_bar."""
    theory = random_theory(seed)
    right_series = taylor_polynomial(theory["right"])
    left_series = taylor_polynomial(theory["left"])
    theta_series = taylor_polynomial(theory["theta"])
    mass_series = taylor_polynomial(theory["mass"])
    mass_zero_determinant = theory["mass"][0].det()

    def phase_of_determinant(matrix_series):
        # arg det M(x) - arg det M(0), with arg det M(0) = 0 for positive m.
        return imaginary_part_series(logarithm_series(
            determinant_series(matrix_series) / mass_zero_determinant))

    vacuum_angle = truncate(theta_series + phase_of_determinant(mass_series) / 2)

    right_generator, left_generator = random_redefinition(redefinition_seed)
    right_group = matrix_exponential_series(right_generator)
    left_group = matrix_exponential_series(left_generator)
    mass_series_new = truncate_matrix(right_group * mass_series * left_group.H)
    theta_series_new = truncate(theta_series - sp.expand(
        (logarithm_series(determinant_series(right_group))
         - logarithm_series(determinant_series(left_group))) / (2 * sp.I)))
    vacuum_angle_new = truncate(theta_series_new + phase_of_determinant(mass_series_new) / 2)

    right_transport = matrix_transport_series(theory["right"])
    left_transport = matrix_transport_series(theory["left"])
    dressed_mass = truncate_matrix(right_transport.H * mass_series * left_transport)
    axial_integral = axial_trace_integral(right_series, left_series)
    dressed_phase_residual = sp.expand(phase_of_determinant(dressed_mass)
                                       - truncate(phase_of_determinant(mass_series)
                                                  - 2 * axial_integral))
    dressed_vacuum_angle = truncate(dressed_theta(theta_series, right_series, left_series)
                                    + phase_of_determinant(dressed_mass) / 2)

    # d/dx log det M at x = 0 equals tr(M_0^{-1} M_1).
    log_derivative_residual = sp.expand(
        sp.diff(logarithm_series(determinant_series(mass_series) / mass_zero_determinant),
                ALP_FIELD_RATIO).subs(ALP_FIELD_RATIO, 0)
        - (theory["mass"][0].inv() * theory["mass"][1]).trace())
    return {"vacuum_angle_residual": sp.expand(vacuum_angle_new - vacuum_angle),
            "dressed_phase_residual": dressed_phase_residual,
            "dressed_vacuum_angle_residual": sp.expand(dressed_vacuum_angle - vacuum_angle),
            "log_derivative_residual": log_derivative_residual}


def first_order_vacuum_angle(quark_mass, scalar, pseudoscalar, axial, gluon_coupling):
    """Alda parametrization at O(x): theta_bar_1 and the two invariant forms."""
    mass_one = (sp.I * (scalar * quark_mass - quark_mass * scalar)
                + sp.I * (pseudoscalar * quark_mass + quark_mass * pseudoscalar))
    inverse_mass = quark_mass.inv()
    log_det_derivative = sp.expand((inverse_mass * mass_one).trace())
    vacuum_angle_one = gluon_coupling + sp.im(log_det_derivative) / 2
    axial_pseudoscalar = axial - pseudoscalar
    gluon_axial_invariant = gluon_coupling + axial.trace()
    dressed_pseudoscalar_one = -(axial_pseudoscalar * quark_mass + quark_mass * axial_pseudoscalar)
    return {
        "pseudoscalar_trace_residual": sp.simplify(
            vacuum_angle_one - gluon_coupling - pseudoscalar.trace()),
        "invariant_form_residual": sp.simplify(
            vacuum_angle_one - gluon_axial_invariant + axial_pseudoscalar.trace()),
        "trace_from_dressed_residual": sp.simplify(
            axial_pseudoscalar.trace()
            + (inverse_mass * dressed_pseudoscalar_one).trace() / 2),
    }


def check_first_order_vacuum_angle(seed=DEFAULT_SEED):
    """Check 4 (first order): theta_bar_1 = c_G + tr c^P = c_GA - tr c^AP."""
    gluon_coupling = sp.Symbol("c_G", real=True)
    diagonal_mass = sp.diag(*sp.symbols("m_u m_d m_s", positive=True))
    symbolic = first_order_vacuum_angle(
        diagonal_mass, symbolic_hermitian_matrix("cS"), symbolic_hermitian_matrix("cP"),
        symbolic_hermitian_matrix("cA"), gluon_coupling)
    generator = random.Random(seed)
    rational = first_order_vacuum_angle(
        positive_mass_matrix(generator, FLAVOR_NUMBER),
        random_rational_hermitian(generator, FLAVOR_NUMBER),
        random_rational_hermitian(generator, FLAVOR_NUMBER),
        random_rational_hermitian(generator, FLAVOR_NUMBER), gluon_coupling)
    return {"symbolic_diagonal_mass": symbolic, "rational_generic_mass": rational}


def check_meson_topological_combination(seed=DEFAULT_SEED, goldstone_seed=11):
    """Check 5: det hat(U) = det U det(W_R)^{-1} det W_L and X_TC equality, O(x^4)."""
    theory = random_theory(seed)
    right_series = taylor_polynomial(theory["right"])
    left_series = taylor_polynomial(theory["left"])
    theta_series = taylor_polynomial(theory["theta"])
    hermitian = random_rational_hermitian(random.Random(goldstone_seed), FLAVOR_NUMBER)
    identity = sp.eye(FLAVOR_NUMBER)
    goldstone_matrix = ((identity + sp.I * hermitian)
                        * (identity - sp.I * hermitian).inv()).applyfunc(sp.simplify)
    unitarity = (goldstone_matrix.H * goldstone_matrix - identity).applyfunc(sp.simplify)

    right_transport = matrix_transport_series(theory["right"])
    left_transport = matrix_transport_series(theory["left"])
    dressed_goldstone = truncate_matrix(right_transport.H * goldstone_matrix * left_transport)
    goldstone_determinant = goldstone_matrix.det()
    right_determinant = determinant_series(right_transport)
    unitary_determinant_residual = truncate(determinant_series(right_transport.H)
                                            * right_determinant - 1)
    factorization_residual = truncate(
        determinant_series(dressed_goldstone) * right_determinant
        - goldstone_determinant * determinant_series(left_transport))

    # tr log hat(u) - tr log u = (1/2) log(det hat(U) / det U), continuous from x = 0.
    trace_log_difference = logarithm_series(
        sp.expand(determinant_series(dressed_goldstone) / goldstone_determinant)) / 2
    dressed = dressed_theta(theta_series, right_series, left_series)
    topological_difference = truncate(
        2 * sp.I * (trace_log_difference + sp.I * dressed)
        - 2 * sp.I * sp.I * theta_series)
    return {"goldstone_is_unitary": unitarity == sp.zeros(FLAVOR_NUMBER, FLAVOR_NUMBER),
            "unitary_determinant_residual": sp.expand(unitary_determinant_residual),
            "factorization_residual": sp.expand(factorization_residual),
            "trace_log_shift_residual": truncate(
                trace_log_difference + sp.I * axial_trace_integral(right_series, left_series)),
            "topological_combination_residual": sp.expand(topological_difference)}


def check_tadpole_source():
    """Check 6: B_0 hat(P)_1 - (m_0^2/N_f) hat(Theta)_1 1 equals Alda's C."""
    quark_mass = hermitian_matrix_symbol("m_q")
    vector, axial, scalar, pseudoscalar = (hermitian_matrix_symbol(name)
                                           for name in ("c^V", "c^A", "c^S", "c^P"))
    axial_pseudoscalar = hermitian_matrix_symbol("c^{AP}")
    gluon_coupling = sp.Symbol("c_G", real=True)
    axial_trace = sp.Symbol(r"\mathrm{tr}\,c^A", real=True)
    chiral_constant, topological_mass_squared, flavor_number = sp.symbols(
        "B_0 m_0^2 N_f", positive=True)

    right_transport = transport_series([vector + axial], ALP_FIELD_RATIO, 1)
    left_transport = transport_series([vector - axial], ALP_FIELD_RATIO, 1)
    mass_function = quark_mass + ALP_FIELD_RATIO * (
        sp.I * commutator(scalar, quark_mass) + sp.I * anticommutator(pseudoscalar, quark_mass))
    dressed_mass = sp.expand(sp.adjoint(right_transport) * mass_function * left_transport)
    _, dressed_pseudoscalar_one = hermitian_parts(dressed_mass.coeff(ALP_FIELD_RATIO, 1))
    dressed_pseudoscalar_one = sp.expand(dressed_pseudoscalar_one.subs(
        axial, axial_pseudoscalar + pseudoscalar))
    dressed_theta_one = gluon_coupling + axial_trace  # Theta_1 + tr C_{A,0}
    source = (chiral_constant * dressed_pseudoscalar_one
              - topological_mass_squared / flavor_number * dressed_theta_one)
    alda_source = (-chiral_constant * anticommutator(quark_mass, axial_pseudoscalar)
                   - topological_mass_squared / flavor_number
                   * (gluon_coupling + axial_trace))
    return {"dressed_pseudoscalar_residual": sp.expand(
                dressed_pseudoscalar_one + anticommutator(axial_pseudoscalar, quark_mass)),
            "source_residual": sp.expand(source - alda_source)}


def report(label, passed):
    print(f"[{'PASS' if passed else 'FAIL'}] {label}")
    return passed


def main():
    results = []
    liouville = check_liouville()
    for chirality in ("right", "left"):
        results.append(report(f"1 {chirality}: dW/dx = iCW and det W = exp(i int tr C) to O(x^4)",
                              liouville[chirality]["ode_is_zero"]
                              and liouville[chirality]["determinant_residual"] == 0))

    for label, traceless in (("nonlinear singlet part", False), ("traceless generators", True)):
        invariance = check_redefinition_invariance(traceless=traceless)
        results.append(report(f"2 {label}: (1/2i) log(det g_R/det g_L) = tr(G_R-G_L)/2",
                              invariance["anomaly_shift_residual"] == 0))
        results.append(report(f"2 {label}: Theta changes = {invariance['theta_changes']}"
                              f" (expected {not traceless})",
                              invariance["theta_changes"] == (not traceless)))
        results.append(report(f"2 {label}: W_R' = g_R W_R to O(x^4)",
                              invariance["transport_covariance_is_zero"]))
        results.append(report(f"2 {label}: hat Theta' = hat Theta to O(x^4)",
                              invariance["dressed_theta_residual"] == 0))

    taylor = check_taylor_form()
    results.append(report("3 hat Theta_n = Theta_n + tr C_{A,n-1} (n = 0..6)",
                          all(residual == 0 for residual in taylor["taylor_residuals"])))
    results.append(report("3 hat Theta_1 = c_G + tr c^A = c_GA",
                          taylor["first_order_residual"] == 0))
    results.append(report("3 hat Theta_1 invariant under c_G -> c_G - tr kappa, c^A -> c^A + kappa",
                          taylor["shift_residual"] == 0))

    vacuum = check_vacuum_angle()
    results.append(report("4 d/dx log det M |_0 = tr(M_0^{-1} M_1)",
                          vacuum["log_derivative_residual"] == 0))
    results.append(report("4 theta_bar = Theta + arg det M / 2 invariant to O(x^4)",
                          vacuum["vacuum_angle_residual"] == 0))
    results.append(report("4 arg det hat M = arg det M - 2 int tr C_A to O(x^4)",
                          vacuum["dressed_phase_residual"] == 0))
    results.append(report("4 hat Theta + arg det hat M / 2 = theta_bar to O(x^4)",
                          vacuum["dressed_vacuum_angle_residual"] == 0))
    first_order = check_first_order_vacuum_angle()
    for case, values in first_order.items():
        results.append(report(f"4 first order ({case}): theta_bar_1 = c_G + tr c^P",
                              values["pseudoscalar_trace_residual"] == 0))
        results.append(report(f"4 first order ({case}): theta_bar_1 = c_GA - tr c^AP",
                              values["invariant_form_residual"] == 0))
        results.append(report(f"4 first order ({case}): tr c^AP = -(1/2) tr(m^-1 hat P_1)",
                              values["trace_from_dressed_residual"] == 0))

    meson = check_meson_topological_combination()
    results.append(report("5 Cayley U is unitary", meson["goldstone_is_unitary"]))
    results.append(report("5 det(W_R^dagger) det(W_R) = 1 to O(x^4)",
                          meson["unitary_determinant_residual"] == 0))
    results.append(report("5 det hat U det W_R = det U det W_L to O(x^4)",
                          meson["factorization_residual"] == 0))
    results.append(report("5 tr log hat u = tr log u - i int tr C_A to O(x^4)",
                          meson["trace_log_shift_residual"] == 0))
    results.append(report("5 2i(tr log u + i Theta) = 2i(tr log hat u + i hat Theta) to O(x^4)",
                          meson["topological_combination_residual"] == 0))

    tadpole = check_tadpole_source()
    results.append(report("6 hat P_1 = -{c^AP, m} from the dressing",
                          tadpole["dressed_pseudoscalar_residual"] == 0))
    results.append(report("6 B_0 hat P_1 - (m_0^2/N_f) hat Theta_1 = -B_0{m,c^AP} - (m_0^2/N_f) c_GA",
                          tadpole["source_residual"] == 0))

    if not all(results):
        sys.exit(1)


if __name__ == "__main__":
    main()
